import socket
import threading
import struct
import numpy as np
import time
import os
import re
from pathlib import Path
import cv2
from dataclasses import dataclass
import math

# =========================
# TCP-Server
# =========================
HOST = '0.0.0.0'
PORT = 5001

WIDTH = 1280
HEIGHT = 720

# Statistik (Thread-sicher)
stats_lock = threading.Lock()
client_stats = {}

# all test data
test_data = []

@dataclass
class Measurement:
    seq: int
    start_time: float = float('nan')
    end_time: float = float('nan')
    rmse: float = float('nan')
    bpr: float = float('nan')
    n_valid: int = 0

    def duration(self) -> float:
        if math.isnan(self.start_time) or math.isnan(self.end_time):
            return float('nan')
        return self.end_time - self.start_time


# =========================
# Stereo Calibration Struktur
# =========================
class StereoCalib:
    def __init__(self, cam0, cam1, doffs, baseline):
        self.cam0 = cam0  # 3x3 np.array float32
        self.cam1 = cam1  # 3x3 np.array float32
        self.doffs = doffs
        self.baseline = baseline

    def pack(self):
        """Return bytes to send via TCP."""
        data = struct.pack('<18f', *(self.cam0.flatten().tolist() + self.cam1.flatten().tolist()))
        data += struct.pack('<2f', self.doffs, self.baseline)
        return data


def parse_matrix(text):
    # erwartet etwas wie "[a b c; d e f; g h i]"
    inner = text.strip()
    inner = inner.lstrip('[').rstrip(']')
    rows = [r.strip() for r in inner.split(';') if r.strip()]
    mat = []
    for r in rows:
        # mehrere Trennzeichen handhaben (Spaces / commas)
        parts = re.split(r'[,\s]+', r.strip())
        mat.append([float(x) for x in parts if x != ''])
    return np.array(mat, dtype=np.float32)

def get_path_for_index(index):
    return test_data[index]

def parse_calib_file(index):
    data = {}
    with open(Path(get_path_for_index(index)) / "calib.txt", 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                continue
            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip()
            if val.startswith('[') and val.endswith(']'):
                # Matrix
                try:
                    data[key] = parse_matrix(val)
                except Exception:
                    data[key] = val
            else:
                # numerische Werte oder Strings
                # versuche int, dann float, sonst raw string
                if re.fullmatch(r'-?\d+', val):
                    data[key] = int(val)
                else:
                    try:
                        data[key] = float(val)
                    except Exception:
                        data[key] = val

    if not all(k in data for k in ("cam0", "cam1", "doffs", "baseline")):
        raise ValueError("Kalibrierungsdatei unvollständig")
    
    calib = StereoCalib(data["cam0"], data["cam1"], data["doffs"], data["baseline"])
    return calib

def scale_calib(calib: StereoCalib, sx: float, sy: float) -> StereoCalib:
    cam0 = calib.cam0.copy()
    cam1 = calib.cam1.copy()

    cam0[0, 0] *= sx  # fx
    cam0[1, 1] *= sy  # fy
    cam0[0, 2] *= sx  # cx
    cam0[1, 2] *= sy  # cy

    cam1[0, 0] *= sx
    cam1[1, 1] *= sy
    cam1[0, 2] *= sx
    cam1[1, 2] *= sy

    doffs = calib.doffs * sx
    return StereoCalib(cam0, cam1, doffs, calib.baseline)


def resize_image(img, new_size):
    h, w = img.shape[:2]

    if h == new_size[1] and w == new_size[0]:
        return img

    new_w, new_h = int(new_size[0]), int(new_size[1])
    scale_x = new_w / w
    scale_y = new_h / h

    # Bild skalieren
    interp_img = cv2.INTER_AREA if (scale_x < 1 or scale_y < 1) else cv2.INTER_LINEAR
    img_rs = cv2.resize(img, (new_w, new_h), interpolation=interp_img)

    return img_rs

def resize_depth_image(img, new_size):
    h, w = img.shape[:2]

    if h == new_size[1] and w == new_size[0]:
        return img

    new_w, new_h = int(new_size[0]), int(new_size[1])
    scale_x = new_w / w
    scale_y = new_h / h

    # Bild skalieren
    img_rs = cv2.resize(img, (new_w, new_h), interpolation=cv2.INTER_NEAREST)

    return img_rs

def resize_disparity(disp, new_size):
    h, w = disp.shape[:2]

    if h == new_size[1] and w == new_size[0]:
        return disp

    new_w, new_h = int(new_size[0]), int(new_size[1])
    scale_x = new_w / w

    # Disparität: Maske für gültige Werte
    valid_mask = np.isfinite(disp).astype(np.uint8)

    # nearest neighbor für Disparität (keine Glättung)
    disp_nn = cv2.resize(disp.astype(np.float32), (new_w, new_h), interpolation=cv2.INTER_NEAREST)

    # Skalieren der Disparitätswerte in neue Pixel-Einheit
    disp_rs = disp_nn * scale_x

    # Maske ebenfalls nearest, Invalids zurücksetzen
    mask_rs = cv2.resize(valid_mask, (new_w, new_h), interpolation=cv2.INTER_NEAREST).astype(bool)
    disp_rs[~mask_rs] = np.nan

    return disp_rs


def read_test_image(index, num=0):
    img_path = Path(get_path_for_index(index)) / f"im{num}.png"
    img = cv2.imread(str(img_path)).astype(np.uint8)
    img = resize_image(img, (WIDTH, HEIGHT))
    return img

def read_disp_image(index, num=0):
    img_path = Path(get_path_for_index(index)) / f"disp{num}.pfm"
    """ Read Middlebury disparity map and calculate depth map.
        http://davis.lbl.gov/Manuals/NETPBM/doc/pfm.html
        https://github.com/singer-yang/Middlebury-Depth-Map
    """
    # ==> Read disparity map from *.pfm file
    disp = cv2.imread(img_path, cv2.IMREAD_UNCHANGED)
    with open(img_path, 'rb') as pfm_file:
        header = pfm_file.readline().decode().rstrip()
        channels = 3 if header == 'PF' else 1

        dim_match = re.match(r'^(\d+)\s(\d+)\s$', pfm_file.readline().decode('utf-8'))
        if dim_match:
            width, height = map(int, dim_match.groups())
        else:
            raise Exception("Malformed PFM header.")

        scale = float(pfm_file.readline().decode().rstrip())    # read disparity scale factor
        if scale < 0:
            endian = '<' # littel endian
            scale = -scale
        else:
            endian = '>' # big endian

    disp = disp * scale
    return disp

def load_all_test_data(base_folder):
    base = Path(base_folder)
    for dirpath, dirnames, filenames in os.walk(base):
        if 'calib.txt' in filenames:
            test_data.append(dirpath)            

def send_close_status(conn):
    # Sende einen Abschlussstatus an den Client
    conn.send(struct.pack('<B', 0))  # Type 0 = Ende

def send_image(conn, type_id, seq, img_left, img_right, calib=None):
    height, width = img_left.shape[:2]
    # Header senden
    conn.send(struct.pack('<BiHH', type_id, seq, width, height))
    
    # Optional StereoCalib senden
    if type_id == 1 and calib is not None:
        print(f"Sending calib for seq={seq} (w={width},h={height})")
        conn.send(calib.pack())
    """
    # testcase
    for y in range(100):
        # change some pixels to white for testing 
        img_left[0, y, 0] = 255
        img_left[0, y, 1] = 255
        img_left[0, y, 2] = 255
        img_left[1, y, 0] = 0
        img_left[1, y, 1] = 0
        img_left[1, y, 2] = 0
        img_left[2, y, 0] = 255
        img_left[2, y, 1] = 255
        img_left[2, y, 2] = 255

        img_left[0, 1279-y, 0] = 255
        img_left[0, 1279-y, 1] = 255
        img_left[0, 1279-y, 2] = 255
        img_left[1, 1279-y, 0] = 0
        img_left[1, 1279-y, 1] = 0
        img_left[1, 1279-y, 2] = 0
        img_left[2, 1279-y, 0] = 255
        img_left[2, 1279-y, 1] = 255
        img_left[2, 1279-y, 2] = 255

        img_left[717, 1279-y, 0] = 255
        img_left[717, 1279-y, 1] = 255
        img_left[717, 1279-y, 2] = 255
        img_left[718, 1279-y, 0] = 0
        img_left[718, 1279-y, 1] = 0
        img_left[718, 1279-y, 2] = 0
        img_left[719, 1279-y, 0] = 255
        img_left[719, 1279-y, 1] = 255
        img_left[719, 1279-y, 2] = 255
    """
    # fülle das ganze bild mit weiß für test
    img_left.fill(255)
    img_right.fill(255)
    # Zeilenweise linkes Bild senden
    # blue, green, red
    for y in range(height):
        conn.send(img_left[y, :, 0].tobytes())
    for y in range(height):
        conn.send(img_left[y, :, 1].tobytes())
    for y in range(height):
        conn.send(img_left[y, :, 2].tobytes())
    
    # Zeilenweise rechtes Bild senden
    # blue, green, red
    for y in range(height):
        conn.send(img_right[y, :, 0].tobytes())
    for y in range(height):
        conn.send(img_right[y, :, 1].tobytes())
    for y in range(height):
        conn.send(img_right[y, :, 2].tobytes())

def receive_depth_image(conn):
    """
    Empfängt eine Anfrage vom Client:
    Type (1 Byte)
    seq  (int32)
    width (int16)
    height(int16)
    dann Zeilenweise ein Graustufenbild (uint8)
    
    Rückgabe:
        type_id, seq, img (np.array), width, height
    """
    # Header empfangen:  seq (4B) + width (2B) + height (2B) (Type (1B) bereits empfangen)
    header_bytes = conn.recv(8)
    if len(header_bytes) < 8:
        raise ConnectionError("Verbindung unterbrochen oder unvollständiger Header")

    seq, width, height = struct.unpack('<iHH', header_bytes)
    print(f"Empfangen: seq={seq}, size={width}x{height}")

    # Bildpuffer erstellen
    # img = np.full((height, width), 255, dtype=np.float32)
    img = np.zeros((height, width), dtype=np.float32)

    # Zeilenweise Bild empfangen
    for y in range(height):
        row = b''
        while len(row) < width * 4:  # float32 = 4 Bytes
            chunk = conn.recv(width * 4 - len(row))
            if not chunk:
                raise ConnectionError("Verbindung unterbrochen beim Bildempfang")
            row += chunk
        if len(row) > width * 4:
            print("zu viele Daten für eine Zeile empfangen")
        # if y < 3:
        #    print(row.hex())
        img[y, :] = np.frombuffer(row, dtype=np.float32)
        row = b''
        if(np.any(img[y, :] < 254) and len(np.where(img[y, :] < 254) ) >1):
            print(f"Row {y} received with invalid depth values.")
            print(f"Bytes: {row[(np.where(img[y, :] < 254)[0]-1)*4 : (np.where(img[y, :] < 254)[1]+2)*4]}")

    return seq, img, width, height

def disparity_to_depth(disp: np.ndarray, calib: StereoCalib, num = 0) -> np.ndarray:
    # ==> Calculate depth map
    depth = calib.baseline * (calib.cam0[0, 0] if num == 0 else calib.cam1[0, 0]) / (disp + calib.doffs) # depth in [mm]

    return depth


def depth_from_left_right(disp_left: np.ndarray, disp_right: np.ndarray, calib: StereoCalib) -> np.ndarray:
    """
    Erzeugt ein Tiefenbild aus beiden Disparitätskarten.
    - Zuerst Depth aus left berechnen.
    - Fehlende Pixel (NaN) aus right verwenden (berechnet mit cam1[0,0]).
    """
    depth = disparity_to_depth(disp_left, calib, 0)

    # depth aus right (verwende cam1 fx falls vorhanden)
    # temporär ein StereoCalib mit cam0=cam1 um die gleiche Funktion zu nutzen:
    depth_r = disparity_to_depth(disp_right, calib, 1)

    # Merge: links bevorzugen, rechts wo links NaN
    mask_fill = ~np.isfinite(depth) & np.isfinite(depth_r)
    depth[mask_fill] = depth_r[mask_fill]
    return depth

def convert_and_save_depth_image(seq, pofix, depth, min_depth=None, max_depth=None):
    """Konvertiert das Tiefenbild in ein 8-bit Schwarz-Weiß-Bild und speichert es als PNG."""
    # ==> Write depth map into 16-bit .png
    depth_path = Path(get_path_for_index(seq)) / f"depth_{seq}_{pofix}.png"
    
    # normalize to 16-bit
    if min_depth is not None and max_depth is not None:
        depth_sanitized = np.clip(depth, min_depth, max_depth)
        depth_sanitized = (depth_sanitized - min_depth) / (max_depth - min_depth) * 65535
    else:
        if np.any(np.isfinite(depth)):
            depth_sanitized = (depth - np.nanmin(depth)) / (np.nanmax(depth) - np.nanmin(depth)) * 65535
        else:
            print(f"Warnung: Kein gültiger Tiefenwert für Seq {seq}, speichere leeres Bild")
    

    depth_sanitized = np.nan_to_num(depth_sanitized, nan=0, posinf=0, neginf=0)
    depth_uint16 = np.round(depth_sanitized).astype(np.uint16)
    #with open(depth_path, 'wb') as f:
        # writer = png.Writer(width=1920, height=1080, bitdepth=16, greyscale=True)
    cv2.imwrite(str(depth_path), np.reshape(depth_uint16, (-1, depth.shape[1])))
    return np.nanmin(depth), np.nanmax(depth)

def compare_img(seq, test_img, abs_thresh: float = 10):
    """Vergleicht das empfangene Bild mit dem Testbild und gibt den mittleren absoluten Fehler zurück."""
    u = test_img.astype(np.uint8, copy=True)
    cv2.imwrite(f"test2.png", u)

    ref_disp_left = read_disp_image(seq, 0)
    ref_disp_right = read_disp_image(seq, 1)
    # ref_disp_left = resize_disparity(ref_disp_left, (1280, 720))
    # ref_disp_right = resize_disparity(ref_disp_right, (1280, 720))
    calib = parse_calib_file(seq)
    ref_depth = depth_from_left_right(ref_disp_left, ref_disp_right, calib)
    ref_depth = resize_depth_image(ref_depth, (WIDTH, HEIGHT))
    min_value, max_value = convert_and_save_depth_image(seq, "ref", ref_depth)
    # test_img[test_img > 10000] = 10000  # max 10m
    convert_and_save_depth_image(seq, "test", test_img, min_value, max_value)
    # Root Mean Square Error (RMSE)
    # Bad Pixel Rate (BPR)

    if test_img.shape != ref_depth.shape:
        print(f"Fehler: Bildgrößen stimmen nicht überein (empfangen: {test_img.shape}, referenz: {ref_depth.shape})")

    valid = np.isfinite(test_img) & np.isfinite(ref_depth)
    n_valid = int(np.count_nonzero(valid))
    if n_valid == 0:
        return float('nan'), float('nan'), 0

    diff = test_img[valid] - ref_depth[valid]
    rmse = float(np.sqrt(np.mean(np.square(diff))))
    me = np.mean(np.abs(diff))
    b1m = float(np.count_nonzero(np.abs(diff) > 1000) / n_valid)
    b1dm = float(np.count_nonzero(np.abs(diff) > 100) / n_valid)
    bpr = float(np.count_nonzero(np.abs(diff) > abs_thresh) / n_valid)
    print(f"SEQ={seq}: RMSE={rmse:.2f}mm, ME={me:.2f}mm, B1m={(1-b1m)*100:.2f}%, B1dm={(1-b1dm)*100:.2f}%, BPR={(1-bpr)*100:.2f}%, N_valid={n_valid}")
    return rmse, bpr, n_valid

def handle_client(conn, addr):
    print(f"Verbindung von {addr}")

    """Wird pro Client gestartet."""
    client_id = f"{addr[0]}:{addr[1]}"
    print(f"[+] Neue Verbindung von {client_id}")

    seq = 0
    measurements = []

    while True:
        # Warte auf Typ-Anfrage vom Client
        try:
            request_byte = conn.recv(1)
            if not request_byte:
                print(f"Client getrennt")
                break
            request = struct.unpack('<B', request_byte)[0]
        except ConnectionResetError:
            print(f"Verbindung verloren")
            break
        if (request == 1 or request == 2) and seq >= len(test_data):
            print(f"Maximum number ({seq}) of test data reached, send Ende-Status")
            send_close_status(conn)
            break
        if request == 0:
            print(f"Client beendet Übertragung")
            break
        elif request == 1:
            print(f"Client ask for img and calib")
            send_image(conn, type_id=1, seq=seq, img_left=read_test_image(seq, 0), img_right=read_test_image(seq, 1), calib=scale_calib(parse_calib_file(seq), WIDTH/1920, HEIGHT/1080))
            if seq >= len(measurements) or measurements[seq] is None:
                meas = Measurement(seq=seq, start_time=time.time())
                measurements.append(meas)
            else:
                print(f"Warnung: Seq {seq} wurde bereits gesendet")
            seq += 1
        elif request == 2:
            print(f"Client ask for img")
            send_image(conn, type_id=2, seq=seq, img_left=read_test_image(seq, 0), img_right=read_test_image(seq, 1))
            if seq >= len(measurements) or measurements[seq] is None:
                meas = Measurement(seq=seq, start_time=time.time())
                measurements.append(meas)
            else:
                print(f"Warnung: Seq {seq} wurde bereits gesendet")
            seq += 1
        elif request == 3:
            print(f"Client send image")
            end_time = time.time()
            seq_comp, img, width, height = receive_depth_image(conn)
            rmse, bpr, n_valid = compare_img(seq_comp, img)
            if measurements[seq_comp] is not None:
                if measurements[seq_comp].end_time is not None and not math.isnan(measurements[seq_comp].end_time):
                    print(f"Warnung: Messung für seq {seq_comp} bereits abgeschlossen")
                else:
                    measurements[seq_comp].end_time = end_time
                    measurements[seq_comp].rmse = rmse
                    measurements[seq_comp].bpr = bpr
                    measurements[seq_comp].n_valid = n_valid
            else:
                print(f"Warnung: Messung für seq {seq_comp} nicht gefunden")
            print(f"Compare result: SEQ={seq_comp} RMSE={rmse}, BPR={bpr}, N_valid={n_valid}")
        else:
            print(f"Unknown Request: {request}")

    conn.close()
    with stats_lock:
        #client_stats[client_id] = {
        #    "frames": frame_count,
        #    "bytes": bytes_received,
        #    "duration": duration,
        #    "fps": frame_count / duration if duration > 0 else 0
        #}
        pass
    # print(f"[-] {client_id} getrennt (Frames: {frame_count}, Zeit: {duration:.2f}s)")



# =========================
# Server starten u = f.astype(np.uint8, copy=True) cv2.imwrite(f"client/disparity_{seq_recv}.png", img_uint16)
# =========================
def server_main():
    base = Path(__file__).parent / "data" / "all"
    load_all_test_data(base)
    
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind((HOST, PORT))
    sock.listen(1)
    print(f"Server wartet auf Port {PORT}")

    try:
        while True:
            conn, addr = sock.accept()
            t = threading.Thread(target=handle_client, args=(conn, addr), daemon=True)
            t.start()
    except KeyboardInterrupt:
        print("Server stoppt...")
    finally:
        sock.close()

if __name__ == "__main__":
    server_main()