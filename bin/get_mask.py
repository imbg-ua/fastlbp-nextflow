#!/usr/bin/env python3

import fire
import numpy as np
import tifffile as tf
from scipy import ndimage
from skimage import filters
import workflow_utils as ut
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

def annot_to_tissue_binmask(annot_img: np.array, background_val: list | None = None) -> np.array:
    annot_img = annot_img.reshape((annot_img.shape[0], annot_img.shape[1], -1))
    channels_num = annot_img.shape[2]

    # default background value is 0 for each channel
    if background_val is None:
        background_val = np.repeat(0, channels_num)
    else:
        background_val = np.array(background_val)

    background_mask = np.ones((annot_img.shape[0], annot_img.shape[1]), dtype=np.bool_)

    for channel in range(channels_num):
        channel_mask = annot_img[..., channel] == background_val[channel]
        background_mask &= channel_mask

    return ~background_mask

# TODO: make this step optional in the pipeline
# (if annotations are already 2D, don't duplicate them in binmask.npy)
def check_annotations(annotations: str, background_val_str: str = "", savefile: str = 'binmask.npy') -> None:
    annot = ut.read_img(annotations)

    # TODO: refactor this
    if background_val_str == "":
        background_val = None
    else:
        background_val = ut.convert_to_proper_type(background_val_str)

    if len(annot.shape) > 2:
        tissue_mask = annot_to_tissue_binmask(annot, background_val)
    elif len(annot.shape) == 2:
        tissue_mask = annot
    else:
        raise ValueError('Annnotations shape is not valid.')

    np.save(savefile, tissue_mask)
    

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
    background: str='light',
    savefile: str='pixelmask.npy') -> None:

    img = ut.read_img(img_path)
    img_grayscale = img_arr_to_grayscale(img)
    img_blurred = blur_img(img_grayscale, blur_k)
    background_mask = get_otsu_mask(img_blurred)

    if background == 'light':
        tissue_mask = ~background_mask
    elif background == 'dark':
        tissue_mask = background_mask

    np.save(savefile, tissue_mask)

if __name__ == '__main__':
    fire.Fire()