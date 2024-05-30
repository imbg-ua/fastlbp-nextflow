#!/usr/bin/env python3

from PIL import Image
import os
import fire
import tifffile as tf
import numpy as np

PIL_EXTENSIONS = ('.jpeg', '.jpg', '.png')
TIFFFILE_EXTENSIONS = ('.tif', '.tiff', '.ndpi')
NUMPY_EXTENSIONS = ('.npy')

DIGITS = set('0123456789') # cringe

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
    # This is the stupidest unimaginable piece of crap 
    if DIGITS.intersection(set(val)):
        if '.' in val:
            return float(val)
        else:
            return int(val)
    elif val.lower() == 'true' or val.lower() == 'false':
        return bool(val)
    return val

def parse_params_str(params_str: str) -> dict:
    """
    Input example: 
    params_str = [param_1, 10, param_2, 20]
    """
    # params_list = list(params_str.strip('[]').split(','))
    params_list = list(map(lambda x: x.strip("[ ]\"'"), str(params_str).split(',')))
    params_dict = {params_list[i]: convert_to_proper_type(params_list[i + 1]) for i in range(0, len(params_list), 2)}

    print(f'{params_str = } ; {params_list = } ; {params_dict = }')
    return params_dict


if __name__ == '__main__':
    fire.Fire()