#!/usr/bin/env python3

import fire
import numpy as np
import workflow_utils as ut


def masked_labels_to_patch_img(
    np_labels_path: str, 
    np_mask_path: str,
    increment: int = 2, 
    savefile: str='patch_labels.npy') -> None:

    labels = np.load(np_labels_path) # dtype = np.int_
    mask = np.load(np_mask_path).astype(np.bool_) # patchmask

    labels_img = np.zeros(shape=mask.shape, dtype=np.uint8)
    labels += increment
    labels_img[mask] = labels

    np.save(savefile, labels_img)

def labels_to_patch_img( 
    np_labels_path: str,
    lbp_output_path: str,
    increment: int = 2, 
    savefile: str='patch_labels.npy') -> None:

    labels = np.load(np_labels_path) # dtype = np.int_
    lbp_res_shape = np.load(lbp_output_path).shape

    labels_img = labels.reshape((lbp_res_shape[0], lbp_res_shape[1]))
    labels_img += increment

    np.save(savefile, labels_img)    

def patch_to_pixel_img(
    np_patch_img_path: str,
    patchsize: int,
    remainder_x0: int,
    remainder_x1: int,
    dtype = np.uint8, 
    savefile: str='img_labels.npy') -> None:

    patch_img = ut.read_img(np_patch_img_path)

    result = np.empty((patch_img.shape[0] * patchsize + remainder_x0, 
                       patch_img.shape[1] * patchsize + remainder_x1), dtype=dtype)

    for row in range(patch_img.shape[0]):
        for col in range(patch_img.shape[1]):
            patch = np.full((patchsize, patchsize), fill_value=patch_img[row, col])
            result[row*patchsize:(row + 1)*patchsize, col*patchsize:(col + 1)*patchsize] = patch

    np.save(savefile, result)

if __name__ == '__main__':
    fire.Fire()