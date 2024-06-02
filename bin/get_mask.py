#!/usr/bin/env python3

import fire
import numpy as np
import tifffile as tf
from scipy import ndimage
from skimage import filters
import workflow_utils as ut
from PIL import Image


def img_arr_to_grayscale(img_arr: np.array) -> np.array:
    return np.dot(img_arr[..., :3], np.array([0.299, 0.587, 0.114])).astype(np.uint8)

def blur_img(img_arr: np.array, k: int=54) -> np.array:    
    blurred_img = ndimage.uniform_filter(img_arr, size=k)
    return blurred_img

def get_otsu_mask(grayscale_img_arr: np.array) -> np.array:
    """
    Apply Otsu thresholding and return the mask separating background from foreground.
    True values in the mask correspond to background (high intensity pixels).

    Parameters
    ----------
    grayscale_img_arr : np.ndarray
        Single channel input image.

    Returns
    -------
    otsu_mask : ndarray
        Background mask of the input image.
    """

    otsu_threshold = filters.threshold_otsu(grayscale_img_arr)
    otsu_mask = grayscale_img_arr > otsu_threshold
    return otsu_mask.astype(np.bool_)

# TODO: combine these functions in one?
# TODO: add downscaling method option
def downscale(
    img_path: str, 
    height: int, 
    width: int, 
    savefile:str, 
    **kwargs) -> None:

    img = ut.read_img(img_path)
    img = Image.fromarray(img)
    img = img.resize((width, height), **kwargs) # dimensions are swapped as compared to numpy
    img = np.asarray(img)

    np.save(savefile, img)

# TODO: add downscaling method option
def downscale_using_patchsize(
    img_path: str, 
    patchsize: int,
    savefile:str, 
    **kwargs) -> None:

    img = ut.read_img(img_path)

    # PIL dimensions are swapped as compared to numpy arr
    height, width = img.shape[0] // patchsize, img.shape[1] // patchsize

    img = Image.fromarray(img)
    img = img.resize((width, height), **kwargs)
    img = np.asarray(img)

    np.save(savefile, img)

def get_mask(
    img_path: str,
    blur_k: int=54,
    savefile: str='pixelmask.npy') -> None:

    img = ut.read_img(img_path)
    img_grayscale = img_arr_to_grayscale(img)
    img_blurred = blur_img(img_grayscale, blur_k)
    background_mask = get_otsu_mask(img_blurred)
    tissue_mask = ~background_mask
    np.save(savefile, tissue_mask)

if __name__ == '__main__':
    fire.Fire()