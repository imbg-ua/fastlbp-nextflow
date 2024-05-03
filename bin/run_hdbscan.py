#!/usr/bin/env python3

import hdbscan
import fire
import numpy as np

def run_hdbscan(
    np_data_path: str, 
    min_samples=10, 
    min_cluster_size=200, 
    cluster_selection_epsilon=0.0,
    gen_min_span_tree=True,
    savefile: str='hdbscan_labels.npy', 
    **kwargs) -> None:

    data = np.load(np_data_path)
    hdbscan_labels = hdbscan.HDBSCAN(min_samples=min_samples, min_cluster_size=min_cluster_size, **kwargs).fit_predict(data)
    np.save(savefile, hdbscan_labels)


if __name__ == '__main__':
    fire.Fire(run_hdbscan)
