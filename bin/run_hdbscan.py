#!/usr/bin/env python3

import hdbscan
from sklearn.cluster import KMeans, AffinityPropagation
from sklearn.mixture import GaussianMixture
import anndata as adata
import scanpy as sc
import fire
import numpy as np
import workflow_utils as ut


def run_hdbscan(
    np_data_path: str, 
    min_samples=10, 
    min_cluster_size=200, 
    cluster_selection_epsilon=0.0,
    gen_min_span_tree=True,
    savefile: str='clustering_labels.npy', 
    **kwargs) -> None:

    data = np.load(np_data_path)
    hdbscan_labels = hdbscan.HDBSCAN(min_samples=min_samples, min_cluster_size=min_cluster_size, **kwargs).fit_predict(data)
    np.save(savefile, hdbscan_labels)

def run_kmeans(
    np_data_path: str, 
    savefile: str='clustering_labels.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    kmeans_labels = KMeans(**kwargs).fit_predict(data)
    np.save(savefile, kmeans_labels)
    
def run_affinity_propagation(
    np_data_path: str, 
    savefile: str='clustering_labels.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    affprop_labels = AffinityPropagation(**kwargs).fit_predict(data)
    np.save(savefile, affprop_labels)

def run_gaussian_mixture(
    np_data_path: str, 
    savefile: str='clustering_labels.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    gauss_mix_labels = GaussianMixture(**kwargs).fit_predict(data)
    np.save(savefile, gauss_mix_labels)
    
def run_leiden(
    np_data_path: str, 
    savefile: str='clustering_labels.npy', 
    **kwargs) -> None:
    data = np.load(np_data_path)
    ad = adata.AnnData(X=embeddings)
    sc.pp.neighbors(ad, n_neighbors=10, use_rep='X')
    sc.tl.leiden(ad, resolution=1.0)
    lied_clust = ad.obs['leiden']
    leid_clust = np.array([int(x) for x in lied_clust])
    np.save(savefile, leid_clust)

def main(np_data_path: str, params_str: str) -> None:
    params_dict = ut.parse_params_str(params_str)
    try:
        method = params_dict.pop('method')
    except KeyError:
        # default clustering method
        # TODO: remove default value?
        method = 'hdbscan'    

    if method == 'k_means':
        run_kmeans(np_data_path=np_data_path,
                    **params_dict)
    elif method == 'hdbscan':
        run_hdbscan(np_data_path=np_data_path,
                    **params_dict)
    elif method == 'leiden':
        run_leiden(np_data_path=np_data_path,
                    **params_dict)
    elif method == 'affinity_propagation':
        run_affinity_propagation(np_data_path=np_data_path, 
                                 **params_dict)
    elif method == 'gaussian_mixture':
        run_gaussian_mixture(np_data_path=np_data_path, 
                             **params_dict)
    elif method == 'leiden':
        run_leiden(np_data_path=np_data_path, 
                   **params_dict)
    else:
        raise ValueError(f'Clustering method \"{method}\" is not a valid option.')

if __name__ == '__main__':
    fire.Fire(main)
