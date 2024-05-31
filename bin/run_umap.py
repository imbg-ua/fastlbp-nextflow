#!/usr/bin/env python3

import umap
import fire
import numpy as np
import workflow_utils as ut

def run_umap(
    np_data_path: str, 
    n_components=2, 
    savefile: str='umap_embeddings.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    embeddings = umap.UMAP(n_components=n_components, **kwargs).fit_transform(data)
    np.save(savefile, embeddings)


def main(np_data_path: str, params_str: str) -> None:
    params_dict = ut.parse_params_str(params_str)
    run_umap(np_data_path=np_data_path, 
             **params_dict)

if __name__ == '__main__':
    fire.Fire(main)
    