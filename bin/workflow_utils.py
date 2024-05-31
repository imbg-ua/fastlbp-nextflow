#!/usr/bin/env python3

from PIL import Image
import os
import fire
import tifffile as tf
import numpy as np
import ast

PIL_EXTENSIONS = ('.jpeg', '.jpg', '.png')
TIFFFILE_EXTENSIONS = ('.tif', '.tiff', '.ndpi')
NUMPY_EXTENSIONS = ('.npy')

def read_img(path: str) -> np.array:
    """
    Generic image reading function.
    """
    path_lower = path.lower()

    if path_lower.endswith(PIL_EXTENSIONS):
        return np.asarray(Image.open(path))

    elif path_lower.endswith(TIFFFILE_EXTENSIONS):
        return tf.imread(path) # returns numpy array

    elif path_lower.endswith(NUMPY_EXTENSIONS):
        return np.load(path)

    else:
        raise ValueError(f'Image format {os.path.splitext(path)[1]} is not supported')

def convert_to_proper_type(val: str):
    try:
        val = ast.literal_eval(val)
    except ValueError:
        pass
    return val

def parse_params_str(params_str: str) -> dict:
    """
    Input example: 
    params_str = [param_1, 10, param_2, 20]
    """

    params_list = list(map(lambda x: x.strip("[ ]\"'"), str(params_str).split(',')))
    params_dict = {params_list[i]: convert_to_proper_type(params_list[i + 1]) for i in range(0, len(params_list), 2)}

    return params_dict


if __name__ == '__main__':
    fire.Fire()