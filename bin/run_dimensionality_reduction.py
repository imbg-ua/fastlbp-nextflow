#!/usr/bin/env python3

import fire
import umap
import numpy as np
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

import workflow_utils as ut


def run_umap(
    np_data_path: str, 
    n_components=2, 
    savefile: str='embeddings.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    embeddings = umap.UMAP(n_components=n_components, **kwargs).fit_transform(data)
    np.save(savefile, embeddings)

def run_tsne(
    np_data_path: str, 
    n_components=2, 
    savefile: str='embeddings.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    embeddings = TSNE(**kwargs).fit_transform(data)
    np.save(savefile, embeddings)

def run_pca(
    np_data_path: str,  
    savefile: str='embeddings.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    embeddings = PCA(**kwargs).fit_transform(data)
    np.save(savefile, embeddings)

def main(np_data_path: str, params_str: str) -> None:

    params_dict = ut.parse_params_str(params_str)
    try:
        method = params_dict.pop('method')
    except KeyError:
        # default dimred method
        # TODO: make explicit only?
        method = 'umap'    

    if method == 'umap':
        run_umap(np_data_path=np_data_path, 
                **params_dict)
    elif method == 'pca':
        run_pca(np_data_path=np_data_path,
                    **params_dict)
    elif method == 'tsne':
        run_tsne(np_data_path=np_data_path,
                    **params_dict)
    else:
        raise ValueError(f'Dimensionality reduction method \"{method}\" is not a valid option.')

    


if __name__ == '__main__':
    fire.Fire(main)
    