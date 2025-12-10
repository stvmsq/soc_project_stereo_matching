import cv2
import numpy as np
from pathlib import Path
import re
import logging

from stereo_calibration import StereoCalib

# ------------------------------
# Disparity Image Loading (Middlebury PFM Files)
# ------------------------------

def read_disp_image(disp_path : str | Path) -> np.ndarray:
    """
    Read a Middlebury-style PFM disparity file and return the disparity map
    scaled to the correct floating-point disparity values.
    
    Inspired by:
        http://davis.lbl.gov/Manuals/NETPBM/doc/pfm.html
        https://github.com/singer-yang/Middlebury-Depth-Map

    Parameters
    ----------
    disp_path : str or Path
        Path to the disparity PFM file.

    Returns
    -------
    np.ndarray
        The disparity image as a float32 array.
    """
    disp = cv2.imread(disp_path, cv2.IMREAD_UNCHANGED)

    # Read the PFM header to extract scale and dimension info
    with open(disp_path, 'rb') as pfm_file:
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

# ------------------------------
# Image Resizing
# ------------------------------

def resize_depth_image(img : np.ndarray, new_size: tuple[int, int]) -> np.ndarray:
    """
    Resize a depth map using nearest-neighbor interpolation.

    Parameters
    ----------
    img : np.ndarray
        Input depth map.
    new_size : tuple(int, int)
        Target size (width, height).

    Returns
    -------
    np.ndarray
        Resized depth map.
    """
    h, w = img.shape[:2]

    if h == new_size[1] and w == new_size[0]:
        return img

    new_w, new_h = int(new_size[0]), int(new_size[1])
    scale_x = new_w / w
    scale_y = new_h / h

    # scale image using nearest neighbor (no smoothing)
    img_rs = cv2.resize(img, (new_w, new_h), interpolation=cv2.INTER_NEAREST)

    return img_rs

def resize_disparity(disp : np.ndarray, new_size: tuple[int, int]) -> np.ndarray:
    """
    Resize a disparity image while preserving disparity correctness.

    Notes
    -----
    - Disparity values must be scaled horizontally when resizing.
    - Invalid pixels (NaN) remain invalid after resizing.
    - Nearest-neighbor interpolation is used to avoid smoothing errors.

    Parameters
    ----------
    disp : np.ndarray
        Disparity map.
    new_size : tuple(int, int)
        Target size (width, height).

    Returns
    -------
    np.ndarray
        Resized disparity map with correct value scaling.
    """
    h, w = disp.shape[:2]

    if h == new_size[1] and w == new_size[0]:
        return disp

    new_w, new_h = int(new_size[0]), int(new_size[1])
    scale_x = new_w / w

    # Valid mask
    valid_mask = np.isfinite(disp).astype(np.uint8)

    # Resize raw disparity (nearest-neighbor)
    disp_nn = cv2.resize(disp.astype(np.float32), (new_w, new_h), interpolation=cv2.INTER_NEAREST)

    # Scale disparity values to new pixel units
    disp_rs = disp_nn * scale_x

    # Resize mask and reapply NaNs
    mask_rs = cv2.resize(valid_mask, (new_w, new_h), interpolation=cv2.INTER_NEAREST).astype(bool)
    disp_rs[~mask_rs] = np.nan

    return disp_rs

# ------------------------------
# Depth Calculation
# ------------------------------

def disparity_to_depth(disp: np.ndarray, calib: StereoCalib, num = 0) -> np.ndarray:
    """
    Convert a disparity map into a depth map using stereo calibration.

    Depth is computed as:

        depth = baseline * fx / (disparity + doffs)

    Parameters
    ----------
    disp : np.ndarray
        Disparity map.
    calib : StereoCalib
        Stereo calibration structure.
    num : int, optional
        Camera index:
            0 = left camera (cam0)
            1 = right camera (cam1)

    Returns
    -------
    np.ndarray
        Depth map with the same dimensions as `disp`.
    """
    fx = calib.cam0[0, 0] if num == 0 else calib.cam1[0, 0]
    depth = calib.baseline * fx / (disp + calib.doffs) # depth in [mm]

    return depth

def depth_from_left_and_right_disp(disp_left: np.ndarray, disp_right: np.ndarray, calib: StereoCalib) -> np.ndarray:
    """
    Compute a depth map by combining left and right disparity maps.

    Strategy
    --------
    - Compute depth from left disparity.
    - Compute depth from right disparity (using cam1 fx).
    - Fill missing pixels (NaNs) from the right depth map.

    Parameters
    ----------
    disp_left : np.ndarray
        Left disparity map.
    disp_right : np.ndarray
        Right disparity map.
    calib : StereoCalib
        Calibration data.

    Returns
    -------
    np.ndarray
        Combined (filled) depth map.
    """
    depth = disparity_to_depth(disp_left, calib, 0)
    depth_r = disparity_to_depth(disp_right, calib, 1)

    # Fill missing values from the right camera
    mask_fill = ~np.isfinite(depth) & np.isfinite(depth_r)
    depth[mask_fill] = depth_r[mask_fill]
    return depth

def get_depth_image(disp_path_left, disp_path_right, calib: StereoCalib, width, height) -> np.ndarray:
    """
    Load disparity maps, compute depth, and resize to target resolution.

    Parameters
    ----------
    disp_path_left : str or Path
        Path to left disparity PFM.
    disp_path_right : str or Path
        Path to right disparity PFM.
    calib : StereoCalib
        Stereo calibration.
    width : int
        Target width.
    height : int
        Target height.

    Returns
    -------
    np.ndarray
        The resized depth map.
    """
    disp_left = read_disp_image(disp_path_left)
    disp_right = read_disp_image(disp_path_right)
    depth = depth_from_left_and_right_disp(disp_left, disp_right, calib)
    depth = resize_depth_image(depth, (width, height))
    return depth

# ------------------------------
# Depth Image Saving
# ------------------------------

def convert_and_save_depth_image(depth_path: Path, depth : np.ndarray, logger : logging.Logger, min_depth : float = None, max_depth: float = None) -> tuple[float, float]:
    """Konvertiert das Tiefenbild in ein 8-bit Schwarz-Weiß-Bild und speichert es als PNG."""
    """
    Convert a depth map into a 16-bit grayscale PNG image and save it.

    Depth values can be normalized explicitly using min/max values,
    or automatically normalized based on the depths range.

    Parameters
    ----------
    depth_path : str or Path
        Output image path (.png recommended).
    depth : np.ndarray
        Depth map.
    min_depth : float, optional
        Minimum depth for normalization.
    max_depth : float, optional
        Maximum depth for normalization.

    Returns
    -------
    tuple(float, float)
        (min_depth, max_depth) actually used after normalization.
    """
    if min_depth is not None and max_depth is not None:
        depth_sanitized = np.clip(depth, min_depth, max_depth)
        depth_sanitized = (depth_sanitized - min_depth) / (max_depth - min_depth) * 65535
    else:
        if np.any(np.isfinite(depth)):
            depth_sanitized = (depth - np.nanmin(depth)) / (np.nanmax(depth) - np.nanmin(depth)) * 65535
        else:
            logger.warning(f"No valid depth values, saving empty image")
    

    depth_sanitized = np.nan_to_num(depth_sanitized, nan=0, posinf=0, neginf=0)
    depth_uint16 = np.round(depth_sanitized).astype(np.uint16)

    cv2.imwrite(str(depth_path), np.reshape(depth_uint16, (-1, depth.shape[1])))
    return np.nanmin(depth), np.nanmax(depth)

# ------------------------------
# Depth Comparison
# ------------------------------


def compare_img(ground_truth : np.ndarray, test_img: np.ndarray, logger : logging.Logger, abs_thresh: float = 10) -> tuple[float, float, int]:
    """
    Compare a test depth map against ground truth using several metrics.

    Metrics
    -------
    - RMSE (Root Mean Square Error)
    - BPR (Bad Pixel Rate): |error| > abs_thresh
    - Count of valid pixels

    Parameters
    ----------
    ground_truth : np.ndarray
        Reference depth map.
    test_img : np.ndarray
        Depth map to evaluate.
    logger : logging.Logger
        Logger for error messages.
    abs_thresh : float, optional
        Threshold in millimeters for BPR.

    Returns
    -------
    tuple(float, float, int)
        (rmse, bpr, valid_pixel_count)
    """
    if test_img.shape != ground_truth.shape:
        logger.error(f"Image size doesn't match (received: {test_img.shape}, reference: {ground_truth.shape})")

    valid = np.isfinite(test_img) & np.isfinite(ground_truth)
    n_valid = int(np.count_nonzero(valid))
    if n_valid == 0:
        return float('nan'), float('nan'), 0

    diff = test_img[valid] - ground_truth[valid]
    rmse = float(np.sqrt(np.mean(np.square(diff))))
    me = np.mean(np.abs(diff))
    b1m = float(np.count_nonzero(np.abs(diff) > 1000) / n_valid)
    b1dm = float(np.count_nonzero(np.abs(diff) > 100) / n_valid)
    bpr = float(np.count_nonzero(np.abs(diff) > abs_thresh) / n_valid)
    
    # logger.info(f"SEQ={seq}: RMSE={rmse:.2f}mm, ME={me:.2f}mm, B1m={(1-b1m)*100:.2f}%, B1dm={(1-b1dm)*100:.2f}%, BPR={(1-bpr)*100:.2f}%, N_valid={n_valid}")
    
    return rmse, bpr, n_valid

