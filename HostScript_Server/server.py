"""
Server for test data streaming and receiving depth images.


Protocol overview (server-side expectations only):
- Client sends 1 byte: Request type
- 0: End/Close
- 1: Request: Image + Calibration
- 2: Request: Image (without calibration)
- 3: Client sends calculated depth image (float32 per pixel, line by line)


The server delivers the requested data (header + image data) and evaluates the received depth images against ground truth.

Parameters: 
HOST: IP address of the server
PORT: Port on which the server should listen
WIDTH: Width to which the images should be scaled and sent.
HEIGHT: Height to which the images should be scaled and sent.
"""
HOST = '0.0.0.0'
PORT = 5001
WIDTH = 1280
HEIGHT = 720

import socket
import threading
import struct
import numpy as np
import time
import os
from pathlib import Path
from dataclasses import dataclass
import math
import logging
import json

import stereo_image
import stereo_calibration
import depth_image

# configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(client_id)s] %(message)s')
class GlobalClientIDFilter(logging.Filter):
    def filter(self, record):
        if not hasattr(record, "client_id"):
            record.client_id = "global"
        return True

global_logger = logging.getLogger('stereo_matching_server')
global_logger.addFilter(GlobalClientIDFilter())

DATA_BASE = Path(__file__).parent / "data" / "all"
RESULT_BASE = Path('result')
RESULT_STATS = RESULT_BASE / 'stats'
RESULT_STATS.mkdir(parents=True, exist_ok=True)

@dataclass
class Measurement:
    seq: int
    start_time: float = float('nan')
    end_time: float = float('nan')
    rmse: float = float('nan')
    bpr: float = float('nan')
    n_valid: int = 0

    def duration(self) -> float:
        """
        Returns the duration of the measurement in seconds (end - start).
        If start or end are NaN, NaN is returned.
        """
        if math.isnan(self.start_time) or math.isnan(self.end_time):
            return float('nan')
        return self.end_time - self.start_time

# List of paths to all test data folders
test_data = []

def get_path_for_index(index: int) -> Path:
    return test_data[index]

def load_all_test_data(base_folder: Path) -> None:
    """
    Recursively search `base_folder` and add every directory that contains a `calib.txt` file.
    Args:
        base_folder: Base folder (path) containing the test data (e.g., data/all).
    """
    base = Path(base_folder)
    if not base.exists():
        global_logger.warning('Base folder %s does not exist', base)
        return
    for dirpath, dirnames, filenames in os.walk(base):
        if 'calib.txt' in filenames:
            test_data.append(dirpath)
    global_logger.info('Loaded test datasets: %d', len(test_data))

# ------------------------------
# Network helper functions
# ------------------------------

def send_close_status(conn: socket.socket) -> None:
    """Send a close status (Type 0) to the client."""
    conn.send(struct.pack('<B', 0))

def send_image(conn : socket.socket, type_id: int, seq: int, img_left: np.ndarray, img_right: np.ndarray, calib=None, logger: logging.LoggerAdapter = None) -> None:
    """
    Send header + (optional) calibration data + image data (B,G,R line by line) to the client.

    Warning: The current protocol sends the calibration bytes directly behind the header —
    the client must know how many bytes to read (or know the format).
    """
    height, width = img_left.shape[:2]
    # Header: type (1B), seq (4B), width (2B), height (2B)
    conn.send(struct.pack('<BiHH', type_id, seq, width, height))
    
    # Optional send StereoCalib data
    if type_id == 1 and calib is not None:
        logger.debug('Sending calibration for seq=%d (w=%d, h=%d)', seq, width, height)
        conn.send(calib.pack())

    # For testing purposes: optionally fill the images with white (the original code filled the images)
    img_left.fill(255)
    img_right.fill(255)

    # send Imagedata: left image B,G,R line by line, then right image
    for ch in range(3):
        for y in range(height):
            conn.send(img_left[y, :, ch].tobytes())
    for ch in range(3):
        for y in range(height):
            conn.send(img_right[y, :, ch].tobytes())

def recv_exact(conn: socket.socket, size: int) -> bytes:
    """
    Receives exactly `size` bytes or throws ConnectionError.

    This helper function ensures that `recv` is called repeatedly until
    all bytes have been received.
    """
    buf = b''
    while len(buf) < size:
        chunk = conn.recv(size - len(buf))
        if not chunk:
            raise ConnectionError('Connection interrupted during recv_exact')
        buf += chunk
    return buf

def receive_depth_image(conn: socket.socket, logger: logging.LoggerAdapter) -> tuple:
    """
    Receives a depth image (float32) from the client.

    Already expects the type byte; then reads seq (int32), width (uint16), height (uint16), and
    then width*height float32 values (line by line).

    Returns:
        seq, img, width, height
    """
    # receive Header:  seq (4B) + width (2B) + height (2B) (Type (1B) already received)
    header_bytes = recv_exact(conn, 8)
    seq, width, height = struct.unpack('<iHH', header_bytes)
    logger.info('Receive depth image: seq=%d, size=%dx%d', seq, width, height)

    # Create image buffer
    img = np.zeros((height, width), dtype=np.float32)

    # line by line image reception
    # per line width * 4 bytes
    row_size = width * 4
    for y in range(height):
        row = recv_exact(conn, row_size)
        img[y, :] = np.frombuffer(row, dtype=np.float32)
        
        # if(np.any(img[y, :] < 254) and len(np.where(img[y, :] < 254) ) >1):
        #    print(f"Row {y} received with invalid depth values.")
        #    print(f"Bytes: {row[(np.where(img[y, :] < 254)[0]-1)*4 : (np.where(img[y, :] < 254)[1]+2)*4]}")

    return seq, img, width, height

# ------------------------------
# Client-Handler
# ------------------------------

def handle_client(conn, addr):
    """
    Handle a connected client.
    
     Args:
        conn: Socket connection to the client.
        addr: Address of the client.
    """
    client_id = f"{addr[0]}:{addr[1]}"
    
    logger = logging.LoggerAdapter(
        logging.getLogger(__name__),
        extra={"client_id": client_id}
    )

    logger.info('[+] New Connection from %s', client_id)

    seq = 0
    measurements = []

    try:
        while True:
            # Warte auf Typ-Anfrage vom Client
            try:
                request_byte = conn.recv(1)
                if not request_byte:
                    logger.info('Client %s disconnected', client_id)
                    break
                request = struct.unpack('<B', request_byte)[0]
            except ConnectionResetError:
                logger.warning('Connection lost with client %s', client_id)
                break

            if (request == 1 or request == 2) and seq >= len(test_data):
                logger.info("Maximum number (%d) of test data reached, send Ende-Status", seq)
                send_close_status(conn)
                break
            
            if request == 0:
                logger.info("Client stopped connection")
                break
            elif request == 1:
                logger.info("Client ask for img and calib (seq=%d)", seq)
                base = get_path_for_index(seq)
                img_left = stereo_image.read_image(base+"/im0.png", WIDTH, HEIGHT, logger=logger)
                img_right = stereo_image.read_image(base+"/im1.png", WIDTH, HEIGHT, logger=logger)
                if img_left is None or img_right is None:
                    logger.error("Error reading images for seq=%d, send close status", seq)
                    send_close_status(conn)
                    break
                calib = stereo_calibration.StereoCalib(base+"/calib.txt")
                calib.scale_calib(WIDTH, HEIGHT)

                send_image(conn, type_id=1, seq=seq, img_left=img_left, img_right=img_right, calib=calib, logger=logger)
                meas = Measurement(seq=seq, start_time=time.time())
                measurements.append(meas)
                seq += 1
            elif request == 2:
                logger.info("Client ask for img (seq=%d)", seq)
                base = get_path_for_index(seq)
                img_left = stereo_image.read_image(base+"/im0.png", WIDTH, HEIGHT, logger=logger)
                img_right = stereo_image.read_image(base+"/im1.png", WIDTH, HEIGHT, logger=logger)
                if img_left is None or img_right is None:
                    logger.error("Error reading images for seq=%d, send close status", seq)
                    send_close_status(conn)
                    break
                send_image(conn, type_id=2, seq=seq, img_left=img_left, img_right=img_right, logger=logger)
                meas = Measurement(seq=seq, start_time=time.time())
                measurements.append(meas)
                seq += 1
            elif request == 3:
                logger.info("Client send image")
                end_time = time.time()
                seq_comp, img, width, height = receive_depth_image(conn, logger)
                if (width, height) != (WIDTH, HEIGHT):
                    logger.warning(f"Receive Image size does not match expected size (received: {width}x{height}, expected: {WIDTH}x{HEIGHT})")
                    send_close_status(conn)
                    # ToDO: notify client about size mismatch
                    break

                base = get_path_for_index(seq_comp)
                calib = stereo_calibration.StereoCalib(base+"/calib.txt")
                ground_truth = depth_image.get_depth_image(base+"/disp0.pfm", base+"/disp1.pfm", calib, WIDTH, HEIGHT)
                rmse, bpr, n_valid = depth_image.compare_img(ground_truth, img, logger=logger)

                if 0 <= seq_comp < len(measurements) and measurements[seq_comp] is not None:
                    if measurements[seq_comp].end_time is not None and not math.isnan(measurements[seq_comp].end_time):
                        logger.warning("Measurement for seq %d already completed", seq_comp)
                    else:
                        measurements[seq_comp].end_time = end_time
                        measurements[seq_comp].rmse = rmse
                        measurements[seq_comp].bpr = bpr
                        measurements[seq_comp].n_valid = n_valid
                else:
                    logger.warning(f"Measurement for seq {seq_comp} not found")
                logger.info(f"Compare result: SEQ={seq_comp} RMSE={rmse}, BPR={bpr}, N_valid={n_valid}")
            else:
                logger.warning(f"Unknown Request: {request}")
    except ConnectionError as e:
        logger.error('Connection error with %s: %s', client_id, e)
    finally:
        conn.close()
    
    # save stats
    frame_count = len([m for m in measurements if m.end_time is not None and not math.isnan(m.end_time)])
    duration = sum([m.duration() for m in measurements if m.end_time is not None and not math.isnan(m.end_time)])
    fps = frame_count / duration if duration > 0 else 0
    logger.info('Result Client %s: Frames: %d, Time: %.2fs, FPS: %.2f', client_id, frame_count, duration, fps)

    with open(f"{RESULT_STATS}/{addr[0]}_{addr[1]}.txt", 'w', encoding='utf-8') as f:
        stats = {
            "client": client_id,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            "total_frames": frame_count,
            "total_time_s": duration,
            "fps": fps,
            "measurements": [
                {
                    "seq": m.seq,
                    "duration_s": m.duration(),
                    "rmse_mm": m.rmse,
                    "bpr": m.bpr,
                    "n_valid": m.n_valid
                } for m in measurements if m.end_time is not None and not math.isnan(m.end_time)
            ]
        }

        json.dump(stats, f, ensure_ascii=False, indent=2)
    

# =========================
# start Server 
# =========================
def server_main():
    """
    Main server function to accept client connections.
    """ 
    load_all_test_data(DATA_BASE)
    
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind((HOST, PORT))
    sock.listen(1)
    global_logger.info('Server listening on %s:%d', HOST, PORT)

    try:
        while True:
            conn, addr = sock.accept()
            t = threading.Thread(target=handle_client, args=(conn, addr), daemon=True)
            t.start()
    except KeyboardInterrupt:
        global_logger.info("Server stopping due to KeyboardInterrupt")
    finally:
        sock.close()

if __name__ == "__main__":
    server_main()