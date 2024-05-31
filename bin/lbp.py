#!/usr/bin/env python3

import fastlbp_imbg as fastlbp
import numpy as np
import fire
import workflow_utils as ut

import os
from PIL import Image

def get_radii_list(patchsize: int, a: float=1.499, b: float=1.327) -> list:
    """
    Progressively increase radii up to the size of patch.
    """
    radii = []
    i = 0
    while True:
        radius = round(a * b**(float(i)))
        if radius <= patchsize:
            radii.append(radius)
            i += 1
        else:
            break
    return radii

# TODO: move to separate module + separate process
def downscale_mask(pixelmask: np.array, np_shape_0, np_shape_1) -> np.array:
    # TODO: change downscaling approach
    pil = Image.fromarray(pixelmask)
    pil = pil.resize((np_shape_1, np_shape_0)) # PIL dimensions are swapped as compared to np array
    return np.array(pil, dtype=np.bool_)

def run_lbp(
    img_path: str,
    patchsize: int = 100, 
    ncpus: int = 1, 
    img_mask: str = None, 
    outfile_name: str = 'lbp_result.py', 
    img_name: str = 'lbp_result', 
    save_intermediate_results: bool = True) -> None:

    img = ut.read_img(img_path)

    radii = get_radii_list(patchsize)
    npoints_list = [fastlbp.get_p_for_r(r) for r in radii]

    if img_mask:
        img_mask = ut.read_img(img_mask)
        fastlbp.run_fastlbp(img, radii, npoints_list, patchsize, ncpus, 
        outfile_name=outfile_name, img_name=img_name, save_intermediate_results=save_intermediate_results,
        img_mask=img_mask.astype(np.uint8))

        # TODO: implement separately (not here)
        lbp_result = np.load(f'data/out/{outfile_name}')
        patch_mask = downscale_mask(img_mask, lbp_result.shape[0], lbp_result.shape[1])

        np.save('patchmask.npy', patch_mask)
        np.save('pixelmask.npy', img_mask)

        lbp_result_flattened = lbp_result[patch_mask]
        
        np.save(f'data/out/{os.path.splitext(outfile_name)[0]}_flattened.npy', lbp_result_flattened)

    else:
        fastlbp.run_fastlbp(img, radii, npoints_list, patchsize, ncpus,
        outfile_name=outfile_name, img_name=img_name, save_intermediate_results=save_intermediate_results)

        # TODO: refactor this
        lbp_result = np.load(f'data/out/{outfile_name}')
        lbp_result_flattened = lbp_result.reshape((-1, lbp_result.shape[-1])) # reshape to 2d <pixel, lbp codes>
        np.save(f'data/out/{os.path.splitext(outfile_name)[0]}_flattened.npy', lbp_result_flattened)


if __name__ == '__main__':
    fire.Fire(run_lbp)
    