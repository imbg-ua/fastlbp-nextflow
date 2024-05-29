#!/usr/bin/env python3

import fire
import numpy as np
import tifffile as tf
import workflow_utils as ut

def get_superpixel_mask(img_path: str, pixelsize: int) -> np.array:
    """
    Crop image to superpixels.
    
    Parameters
    ----------
    imr_path : str
        Input image path.

    pixelsize : int
        Super pixel size. 

    Returns
    -------
    cropped_image : ndarray
        Cropped image.

    top, bottom, left, right : int
        Number of pixels that were cropped from the top, bottom, left and right side of the image respectively.
    """
    
    img = ut.read_img(img_path)

    h = img.shape[0] # rows
    w = img.shape[1] # columns
    mask = np.zeros((h, w), dtype=np.uint8) # change to np._bool?

    mask[0:(h // pixelsize), 0:(w // pixelsize)] = 1

    return mask




def main(img_path, patchsize, savepath):
    img = tf.imread(img_path)
    img, rem_1, rem_2 = crop_to_superpixel(img, patchsize)
    tf.imwrite(savepath, img)



if __name__ == '__main__':
    fire.Fire(main)