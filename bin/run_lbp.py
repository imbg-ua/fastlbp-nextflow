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

# # TODO: move to get_mask.py, maybe even as a pipeline step
# def convert_annot_to_binmask(annot: np.array) -> np.array:
#     pass

def run_lbp(
    img_path: str,
    patchsize: int = 100, 
    radii: list[int] | None = None,
    npoints: list[int] | None = None,
    ncpus: int = 1, 
    img_mask: str = None,
    patch_mask: str = None, 
    outfile_name: str = 'lbp_result.py', 
    img_name: str = 'lbp_result', 
    save_intermediate_results: bool = True) -> None:

    if radii is None:
        radii = get_radii_list(patchsize)

    if npoints is None:
        npoints = [fastlbp.get_p_for_r(r) for r in radii]

    img = ut.read_img(img_path)

    if img_mask:
        img_mask = ut.read_img(img_mask).astype(np.bool_)
        patch_mask = ut.read_img(patch_mask).astype(np.bool_)

        fastlbp.run_fastlbp(img, radii, npoints, patchsize, ncpus, 
        outfile_name=outfile_name, img_name=img_name, save_intermediate_results=save_intermediate_results,
        img_mask=img_mask.astype(np.uint8))

        lbp_result = np.load(f'data/out/{outfile_name}')
        lbp_result_flattened = lbp_result[patch_mask]
        
        np.save(f'data/out/{os.path.splitext(outfile_name)[0]}_flattened.npy', lbp_result_flattened)

    else:
        fastlbp.run_fastlbp(img, radii, npoints, patchsize, ncpus,
        outfile_name=outfile_name, img_name=img_name, save_intermediate_results=save_intermediate_results)

        # TODO: refactor this
        lbp_result = np.load(f'data/out/{outfile_name}')
        lbp_result_flattened = lbp_result.reshape((-1, lbp_result.shape[-1])) # reshape to 2d <pixel, lbp codes>
        np.save(f'data/out/{os.path.splitext(outfile_name)[0]}_flattened.npy', lbp_result_flattened)


# TODO: make these 2 functions a single one
def main_grid_search(
    img_path: str, 
    params_str: str, 
    ncpus: int = 1, 
    img_mask: str = None,
    patch_mask: str = None, 
    outfile_name: str = 'lbp_result.py', 
    img_name: str = 'lbp_result', 
    save_intermediate_results: bool = True) -> None:

    params_dict = ut.parse_params_str(params_str)
    run_lbp(img_path=img_path, 
            ncpus=ncpus, 
            img_mask=img_mask, 
            patch_mask=patch_mask, 
            outfile_name=outfile_name, 
            img_name=img_name, 
            save_intermediate_results=save_intermediate_results, 
            **params_dict)

def main(
    img_path: str,
    params_str: str,
    img_mask: str = None,
    patch_mask: str = None,
    save_intermediate_results: bool = True) -> None:

    params_dict = ut.parse_params_str(params_str)
    run_lbp(img_path=img_path,
            img_mask=img_mask,
            patch_mask=patch_mask,
            save_intermediate_results=save_intermediate_results,
            **params_dict)

if __name__ == '__main__':
    fire.Fire()
    