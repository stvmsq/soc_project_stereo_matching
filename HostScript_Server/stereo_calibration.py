import struct
import re
import numpy as np
from pathlib import Path

# ------------------------------
# Helpher functions
# ------------------------------

def parse_3x3_float_matrix(text):
    """
    Parse a 3x3 float matrix from a calibration file.
    The function expects a string representation formatted like:
        "[a b c; d e f; g h i]"
    It supports flexible separators such as spaces or commas.

    Parameters
    ----------
        text : str
            The matrix string to parse.

    Returns
    -------
        np.ndarray : A 3x3 NumPy array of dtype float32.

    Raises
    ------
        ValueError : If the matrix cannot be parsed into 3 rows of floats.
    """
    inner = text.strip()
    inner = inner.lstrip('[').rstrip(']')
    rows = [r.strip() for r in inner.split(';') if r.strip()]
    mat = []
    for r in rows:
        # split by spaces or commas
        parts = re.split(r'[,\s]+', r.strip())
        mat.append([float(x) for x in parts if x != ''])
    return np.array(mat, dtype=np.float32)

# ------------------------------
# Stereo Calibration Structure
# ------------------------------

class StereoCalib:
    """
    Represents stereo camera calibration data including both intrinsic
    camera matrices, image size, disparity offset, and baseline.
    The class can either be created manually by passing all parameters
    or loaded from a calibration file.

    Parameters
    ----------
    path : str or Path, optional
        Path to the calibration file to load.
    width : int, optional
        Image width in pixels.
    height : int, optional
        Image height in pixels.
    cam0 : np.ndarray, optional
        3×3 intrinsic matrix of the left camera.
    cam1 : np.ndarray, optional
        3×3 intrinsic matrix of the right camera.
    doffs : float, optional
        Disparity offset.
    baseline : float, optional
        Baseline distance between the two cameras.

    Raises
    ------
    ValueError
        If neither a file path nor the complete calibration parameters
        are provided.
    """

    def __init__(self, path = None, width=None, height=None, cam0=None, cam1=None, doffs=None, baseline=None):
        if path is not None:
            self.__init_from_file(path)
            return
        elif None in (width, height, cam0, cam1, doffs, baseline):
            raise ValueError("Either provide a calibration file path or all calibration parameters.")
        self.cam0 = cam0  # 3x3 np.array float32
        self.cam1 = cam1  # 3x3 np.array float32
        self.doffs = doffs
        self.baseline = baseline
        self.width = width
        self.height = height

    def __init_from_file(self, path):
        """
        Load and parse stereo calibration parameters from a text file.

        The file is expected to contain key-value pairs such as:
            width = 1240
            height = 376
            cam0 = [fx 0 cx; 0 fy cy; 0 0 1]
            cam1 = [fx 0 cx; 0 fy cy; 0 0 1]
            doffs = 124.2
            baseline = 0.54

        Parameters
        ----------
        path : str or Path
            The path to the calibration file.

        Raises
        ------
        ValueError
            If required fields are missing.
        """
        data = {}
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if '=' not in line:
                    continue
                key, val = line.split('=', 1)
                key = key.strip()
                val = val.strip()
                # Detect matrix fields
                if val.startswith('[') and val.endswith(']'):
                    try:
                        data[key] = parse_3x3_float_matrix(val)
                    except Exception:
                        data[key] = val
                else:
                    # Try parsing ints, then floats, else raw strings
                    if re.fullmatch(r'-?\d+', val):
                        data[key] = int(val)
                    else:
                        try:
                            data[key] = float(val)
                        except Exception:
                            data[key] = val

        if not all(k in data for k in ("width", "height", "cam0", "cam1", "doffs", "baseline")):
            raise ValueError("Calibration file is missing required parameters.")
        self.width = data["width"]
        self.height = data["height"]
        self.cam0 = data["cam0"]
        self.cam1 = data["cam1"]
        self.doffs = data["doffs"]
        self.baseline = data["baseline"]
        
    def scale_calib(self, width, height):
        """
        Scale the calibration parameters to a new image resolution.

        The intrinsic matrices are scaled proportionally, and the
        disparity offset is scaled horizontally.

        Parameters
        ----------
        width : int
            New target image width.
        height : int
            New target image height.
        """
        sx = width / self.width
        sy = height / self.height
        self.cam0[0, 0] *= sx  # fx
        self.cam0[1, 1] *= sy  # fy
        self.cam0[0, 2] *= sx  # cx
        self.cam0[1, 2] *= sy  # cy

        self.cam1[0, 0] *= sx
        self.cam1[1, 1] *= sy
        self.cam1[0, 2] *= sx
        self.cam1[1, 2] *= sy

        self.width = int(self.width * sx)
        self.height = int(self.height * sy)
        self.doffs = self.doffs * sx


    def pack(self):
        """
        Serialize the stereo calibration parameters for transmission
        over a TCP connection.

        The data layout is as follows (little-endian floats):
            - 9 floats: cam0 matrix flattened row-major
            - 9 floats: cam1 matrix flattened row-major
            - 1 float: disparity offset (doffs)
            - 1 float: baseline

        Returns
        -------
        bytes
            Raw byte buffer containing 20 floats (18 + 2).
        """
        data = struct.pack('<18f', *(self.cam0.flatten().tolist() + self.cam1.flatten().tolist()))
        data += struct.pack('<2f', self.doffs, self.baseline)
        return data