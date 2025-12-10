import cv2
import numpy as np
import logging
def resize_image(img: np.ndarray, new_size: tuple[int, int]) -> np.ndarray:
    """
    Resize an image to a new resolution.

    Parameters:
        img : np.ndarray
            Input image as a NumPy array (H x W x C).
        new_size : tuple[int, int]
            Target size in the format (width, height).
    Returns:
        np.ndarray
            The resized image. If the size already matches, the input
            image is returned unchanged.
    """
    h, w = img.shape[:2]

    # No resizing needed
    if h == new_size[1] and w == new_size[0]:
        return img

    new_w, new_h = int(new_size[0]), int(new_size[1])
    scale_x = new_w / w
    scale_y = new_h / h

    # Choose interpolation:
    # - INTER_AREA for downscaling
    # - INTER_LINEAR for upscaling
    interp_img = cv2.INTER_AREA if (scale_x < 1 or scale_y < 1) else cv2.INTER_LINEAR
    img_rs = cv2.resize(img, (new_w, new_h), interpolation=interp_img)

    return img_rs


def read_image(img_path: str, width: int, height: int, logger : logging.Logger) -> np.ndarray | None:
    """
    Load an image from disk and resize it to the target resolution.

    Parameters:
        img_path : Path to the image file.
        width : Target width after resizing.
        height : Target height after resizing.

    Returns:
        np.ndarray | None : The loaded and resized image in uint8 format, or None if
        loading fails.
    """
    try:
        img = cv2.imread(str(img_path)).astype(np.uint8)
        img = resize_image(img, (width, height))
        return img
    except Exception as e: 
        logger.error(f"Error reading image {img_path}: {e}")
        return None