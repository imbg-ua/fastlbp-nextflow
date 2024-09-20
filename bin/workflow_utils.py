#!/usr/bin/env python3

from PIL import Image
import os
import fire
import tifffile as tf
import numpy as np
import ast
import re

Image.MAX_IMAGE_PIXELS = None

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


# === Jaccard score calculation utils === #

import pandas as pd
from sklearn.metrics import jaccard_score

def load_and_resize(gt_path: str, patch_labels_path: str) -> tuple[np.array, np.array]:
    gt = read_img(gt_path)
    pred = read_img(patch_labels_path)
    
    target_shape = (pred.shape[1], pred.shape[0])
    gt_pil = Image.fromarray(gt)
    gt_pil = gt_pil.resize(target_shape)
    gt_np = np.asarray(gt_pil)
    return pred, gt_np

def make_pair_jacc(gt_res, pred, savefile: str = None):
    
    # get gt and pred labels
    gt_lbls = np.unique(gt_res).astype(np.uint8)  
    pred_lbls = np.unique(pred).astype(np.uint8)

    gt_lbls = gt_lbls[1:]
    pred_lbls = pred_lbls[1:]

    # dataframe where jacc scores for all possible gt-pd labels are stored
    df = pd.DataFrame(columns=tuple([f'pd_{pred_lbl}' for pred_lbl in pred_lbls]))
    count = 0
    jaccs = []

    for gt_lbl in gt_lbls:
        gr_tr = np.where(gt_res==gt_lbl, 1, 0).astype(np.uint8)   

        for pd_lbl in pred_lbls:
            predict = np.where(pred==pd_lbl, 1, 0).astype(np.uint8)       
            
            jacc = jaccard_score(gr_tr, predict, average='micro') # calc jacc score labels=[1], average=None average=None
            jaccs.append(round(jacc, 3))   

        df.loc[count] = jaccs
        jaccs = []
        count += 1

    df['index'] = [f'gt_{gt_lab}' for gt_lab in gt_lbls]
    df = df.set_index('index', drop=True)
    df.index.name = None
    # save df
    if savefile:
        df.to_csv(savefile)
    return df


# the number of times we need to iterate over the df
def calc_times(df):
    if df.shape[0] == df.shape[1]:
        times = df.shape[0]
    else:
        times = min(df.shape)
    return times

# dataframe to store gt-pd pairs as well as max Jaccard values
def make_df_max(times):
    df_max = pd.DataFrame(columns=['gt', 'pd', 'max_Jacc'])
    df_max['index'] = range(times)
    df_max = df_max.set_index('index', drop=True)
    df_max.index.name = None
    return df_max

# find out WHAT the the max val is
def find_max_elem(df):
    max_elem = df.max().max()
    return max_elem

# find out WHERE the max val is
def find_max_elem_ind_col(df, max_elem):
    for col in df.columns:
        for ind in df.index:
            if df.at[ind, col] == max_elem:
                return ind, col
            
# add the pair and the max value to df_max
def add_to_max(ind, col, df, max_elem, df_max, time):
    ind, col = find_max_elem_ind_col(df, max_elem)
    df_max.loc[time] = [ind, col, max_elem]

# drop row and column where the max_elem is located
def drop_max(df, ind ,col):
    df.drop(col, axis=1, inplace=True)
    df.drop(str(ind), axis=0, inplace=True)
            
def run_paiwise_max_jacc(gt_path: str, 
                         patch_labels_path: str, 
                         savefile_all_pairs: str = 'all_pairs_jacc.csv',
                         savefile_max_pairs: str = 'pairs_max_jacc.csv'):
    """
    Calculate pairwise greedy cluster matching based off Jaccard scores between classes/clusters. 

    Parameters
    ----------
    gt_path : str
        Ground truth (annotation) path.

    patch_labels_path : str
        Path to the clustering result.

    save_dir : str
        Directory to save resulting csv with Jaccard scores.

    Returns
    -------
    None
    """
    pred, gt_res = load_and_resize(gt_path, patch_labels_path) # load
    
    df = make_pair_jacc(gt_res, pred, savefile=savefile_all_pairs) # original df is saved here
    times = calc_times(df) # the number of times we need to iterate over the df
    df_max = make_df_max(times) # dataframe to store gt-pd pairs as well as max Jaccard values
            
    for time in range(times):
        max_elem = find_max_elem(df)
        ind, col = find_max_elem_ind_col(df, max_elem)
        add_to_max(ind, col, df, max_elem, df_max, time)
        drop_max(df, ind ,col)
    
    df_max.to_csv(savefile_max_pairs) # max parwise Jacc df is saved here


if __name__ == '__main__':
    fire.Fire()