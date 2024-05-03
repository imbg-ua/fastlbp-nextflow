#!/usr/bin/env python3

import umap
import fire
import numpy as np

def run_umap(
    np_data_path: str, 
    n_components=2, 
    savefile: str='umap_embeddings.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    embeddings = umap.UMAP(n_components=n_components, **kwargs).fit_transform(data)
    np.save(savefile, embeddings)


if __name__ == '__main__':
    fire.Fire(run_umap)
    