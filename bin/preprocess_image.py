#!/usr/bin/env python3

import fire
import numpy as np
import tifffile as tf


def crop_to_superpixel(img_arr: np.array, pixelsize: int) -> np.array:
    """
    Crop image to superpixels given a pixel size.
    
    Parameters
    ----------
    imr_arr : np.ndarray
        Input image.

    pixelsize : int
        Super pixel size. 

    Returns
    -------
    cropped_image : ndarray
        Cropped image.

    top, bottom, left, right : int
        Number of pixels that were cropped from the top, bottom, left and right side of the image respectively.
    """
    remainder_0 = img_arr.shape[0] % int(pixelsize)
    remainder_1 = img_arr.shape[1] % int(pixelsize)

    x0 = img_arr.shape[0] - remainder_0
    x1 = img_arr.shape[1] - remainder_1
    cropped_image = img_arr[0:x0, 0:x1]
    
    return cropped_image, remainder_0, remainder_1


def main(img_path, patchsize, savepath):
    img = tf.imread(img_path)
    img, rem_1, rem_2 = crop_to_superpixel(img, patchsize)

    tf.imwrite(img, savepath)


if __name__ == '__main__':
    fire.Fire(main)