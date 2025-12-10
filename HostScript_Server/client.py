import socket
import struct
import numpy as np
import cv2
import time
from matplotlib import pyplot as plt

HOST = '127.0.0.1'   # Server IP
PORT = 5001          # Server Port

def recv_all(sock, n):
    buf = b''
    while len(buf) < n:
        chunk = sock.recv(n - len(buf))
        if not chunk:
            raise ConnectionError("Verbindung beendet (recv_all)")
        buf += chunk
    return buf

def recv_image_rows(sock, width, height, channels=3):
    # Server sendet zeilenweise: jede Zeile als raw bytes of width*channels
    row_bytes = width * channels
    data = bytearray()
    for _ in range(height):
        data += recv_all(sock, row_bytes)
    arr = np.frombuffer(data, dtype=np.uint8)
    arr = arr.reshape((height, width, channels))
    return arr

def parse_calib_bytes(b):
    # erwartet 20 floats: cam0(9), cam1(9), doffs, baseline
    vals = struct.unpack('<20f', b)
    cam0 = np.array(vals[0:9], dtype=np.float32).reshape((3,3))
    cam1 = np.array(vals[9:18], dtype=np.float32).reshape((3,3))
    doffs = float(vals[18])
    baseline = float(vals[19])
    return cam0, cam1, doffs, baseline

def compute_depth_from_disp(disp, fx, baseline, doffs):
    # disp: float32 in pixel units, np.nan for invalid
    denom = disp + doffs
    depth = np.full(disp.shape, np.nan, dtype=np.float32)
    valid = np.isfinite(denom) & (denom != 0.0)
    depth[valid] = (fx * baseline) / denom[valid]
    return depth

def run_client(host=HOST, port=PORT, max_seq=10):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((host, port))
    try:
        seq = 0
        while seq < max_seq:
            # 1 = ask for img + calib (server implementation)
            s.send(struct.pack('<B', 1))
            # Receive header: type(1B) + seq(int32) + width(int32) + height(int32) => 13 bytes
            hdr = recv_all(s, 9)
            type_id, seq_recv, width, height = struct.unpack('<BiHH', hdr)
            if type_id != 1: # type_id = 2 for img only
                print("Unerwarteter Type-ID vom Server:", type_id)
                break
            print(f"Header empfangen: seq={seq_recv}, width={width}, height={height}")
            # receive calib (20 floats = 80 bytes)
            calib_bytes = recv_all(s, 80)
            cam0, cam1, doffs, baseline = parse_calib_bytes(calib_bytes)
            fx = float(cam0[0,0])

            # receive left and right images (server sends color rows)
            left = recv_image_rows(s, width, height, channels=3)
            right = recv_image_rows(s, width, height, channels=3)

            # optional: show/save
            # cv2.imwrite(f"client/left_{seq_recv}.png", left)
            # cv2.imwrite(f"right_{seq_recv}.png", right)
            print(f"Received seq={seq_recv} (w={width},h={height})")
            # convert to gray for SGM
            left_gray = cv2.cvtColor(left, cv2.COLOR_BGR2GRAY)
            right_gray = cv2.cvtColor(right, cv2.COLOR_BGR2GRAY)
            # cv2.imwrite(f"client/left_{seq_recv}_gray.png", left_gray)
            # cv2.imwrite(f"client/right_{seq_recv}_gray.png", right_gray)
            
            stereo = cv2.StereoSGBM_create(numDisparities=176, blockSize=1)
            disparity = stereo.compute(left_gray, right_gray)

            disparity = disparity.astype(np.float32) / 16.0

            # visualize disparity
            # disparity_img = disparity.copy()
            # disparity_img[disparity_img < 0] = 0  # ungültige Disparitäten
            # img_norm = (disparity_img - np.min(disparity_img)) / (np.max(disparity_img) - np.min(disparity_img))
            # img_uint16 = np.clip(img_norm * 65535, 0, 65535).astype(np.uint16)
            # cv2.imwrite(f"client/disparity_{seq_recv}.png", img_uint16)
            
            # compute depth using calib parameters
            disparity[disparity < 0] = np.nan  # ungültige Disparitäten
            depth = compute_depth_from_disp(disparity, fx=fx, baseline=baseline, doffs=doffs).astype(np.float32)
            # visualize depth
            # depth_sanitized = np.nan_to_num(depth, nan=0.0, posinf=0.0, neginf=0.0)
            # depth_uint16 = np.round(depth_sanitized).astype(np.uint16)
            # cv2.imwrite(f"client/depth_{seq_recv}.png", depth_uint16)

            # send back depth: request type 3 then header and rows (float32 rows)
            # first send single request byte 3 (server already expects request byte separately)
            # s.send(struct.pack('<B', 3))
            # now payload header: type(1B)=3, seq(int32), width(int32), height(int32)
            payload_hdr = struct.pack('<BiHH', 3, seq_recv, width, height)
            s.send(payload_hdr)
            # send rows as float32
            for y in range(height):
                row = depth[y, :].astype(np.float32).tobytes()
                s.send(row)

            print(f"Sent depth for seq={seq_recv} (w={width},h={height}), fx={fx:.2f}, baseline={baseline:.2f}")
            seq += 1
            # small delay to avoid busy loop
            time.sleep(0.1)
    except Exception as e:
        print("Fehler aufgetreten:", e)
    finally:
        print("Schließe Verbindung zum Server.")
        s.close()

if __name__ == "__main__":
    print(cv2.__version__)
    run_client()