#!/usr/bin/env python3

from PIL import Image
import os
import fire
import tifffile as tf
import numpy as np
import ast
import re

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

# TODO: make this better
def parse_params_str(params_str: str) -> dict:
    """
    Input example: 
    params_str = [param_1, 10, param_2, 20]
    """
    params_str = str(params_str)

    # TODO: this is not versatile, need to cover all possible options for str params
    words_without_quotes = r"(?<!')\b[a-zA-Z_](\w*-*\.?)*\w*\b(?!')"
    params_str = re.sub(words_without_quotes, r"'\g<0>'", params_str) # single quote every string

    bool_regex = r"'\b(true|false)\b'"

    def bool_to_capital(match):
        word = match.group(1)
        if word.lower() == 'true':
            return 'True'
        elif word.lower() == 'false':
            return 'False'

    params_str = re.sub(bool_regex, bool_to_capital, params_str, flags=re.IGNORECASE)

    params_list = convert_to_proper_type(params_str)
    params_dict = {param_name: param_val for param_name, param_val in params_list}

    return params_dict


if __name__ == '__main__':
    fire.Fire()