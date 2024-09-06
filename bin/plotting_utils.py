#!/usr/bin/env python3

import numpy as np
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from PIL import Image
import os
import workflow_utils as ut
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as mcolors
from copy import deepcopy
import base64
from io import BytesIO


def dimred_2D_html(dimred_path: str, labels_path: str, component_1: int=0, component_2: int=1,
                   include_plotlyjs: str | bool = 'cdn') -> str:
    dimred = np.load(dimred_path)
    labels = np.load(labels_path)

    data = dimred_and_labels_to_dataframe(dimred, labels, component_1, component_2)

    fig = plot_dimred_2D_labelled(data)

    plot_html = fig.to_html(full_html = False, include_plotlyjs=include_plotlyjs)

    return plot_html

def img_html(img_path: str, include_plotlyjs: str | bool = 'cdn') -> str:
    img = ut.read_img(img_path)

    fig = plot_img(img)

    plot_html = fig.to_html(full_html = False, include_plotlyjs=include_plotlyjs)

    return plot_html

def img_downscaled_html(img_path: str, size=(400, 400), include_plotlyjs: str | bool = 'cdn') -> str:
    img = ut.read_img(img_path)

    fig = plot_img_downscaled(img, size)

    plot_html = fig.to_html(full_html = False, include_plotlyjs=include_plotlyjs)

    return plot_html

# TODO: move to workflow utils?
def dimred_and_labels_to_dataframe(dimred_embeddings: np.array, labels: np.array, component_1: int=0, component_2: int=1):
    comp_1 = dimred_embeddings[:, component_1]
    comp_2 = dimred_embeddings[:, component_2]
    
    dat = np.array([comp_1, comp_2, labels])
    dat = dat.T

    # TODO: refactor hardcoded names
    dat = pd.DataFrame(dat, columns = ["component_1", "component_2", "label"])
    return dat

# TODO: probably don't need this one in reports
def plot_dimred_2D(dimred_embeddings: np.array, component_1: int=0, component_2: int=1):
    fig = px.scatter(x=dimred_embeddings[:, component_1], y=dimred_embeddings[:, component_2])
    return fig

def plot_dimred_2D_labelled(dat: pd.DataFrame):
    fig = px.scatter(dat, x="component_1", y="component_2", color="label",
                    labels={
                        "component_1": "Component 1",
                        "component_2": "Component 2",
                        "label": "Cluster Value"
                    },
                    title="Projections")
    return fig

# TODO: downscale without PIL
def plot_img_downscaled(img_arr: np.array, size=(400, 400)):
    img_pil = Image.fromarray(img_arr)
    img_pil.thumbnail(size)
    img_pil = np.asarray(img_pil)
    fig = px.imshow(img_pil)
    return fig

# Can be used for patches, individual LBP features, etc
def plot_img(img: np.array):
    fig = px.imshow(img)
    return fig



def get_preconstructed_cmap(id: str = 'viridis') -> mpl.colors.Colormap:
    """Get colormap from a set if options.

    Get a colormap by the name. If there is no existing matplotlib colormap with a given name,
    return a custom one with the color sets specified below.

    Parameters
    ----------
    id : str, default: 'viridis'
        Colormap name. Available customs colormaps include:
        - 'omer_colours'
        - 'ben_colorblind'
        - 'ben_colorblind_large'
    
    Returns
    -------
    cmap : mpl.colors.Colormap
        Maplotlib colormap that can be used further.
    """
    try:
        cmap = mpl.colormaps[id]
        return cmap
    except KeyError:
        if id == 'omer_colours':
            # adopted from https://github.com/spectralnanodiamond/ImageTextureFinder/blob/2907bc2782dffb95ed4e5930c5c0ed0b02e15134/useful_functions.py#L1185
            out = ['blue', 'darkorange', 'deeppink', 'yellow', 'lightskyblue', 'peachpuff']
            return LinearSegmentedColormap.from_list('omer_cmap', out)
        elif id == 'ben_colorblind':
            # adopted from https://github.com/spectralnanodiamond/ImageTextureFinder/blob/2907bc2782dffb95ed4e5930c5c0ed0b02e15134/useful_functions.py#L1205
            out = ['lightgray', 'blue', 
                   'darkorange', 'deeppink', 
                    'yellow', 'lightskyblue', 
                    'peachpuff', 'purple', 
                    'darkgreen', 'black', 
                    'dimgray']
            return LinearSegmentedColormap.from_list('ben_colorblind_cmap', out)
        elif id == 'ben_colorblind_large':
            # adopted from https://github.com/spectralnanodiamond/ImageTextureFinder/blob/2907bc2782dffb95ed4e5930c5c0ed0b02e15134/useful_functions.py#L1216
            colors3 = list(mcolors.CSS4_COLORS.keys())
            colors3.remove('rebeccapurple')
            colors3.remove('teal')
            colors3.remove('white')
            colors3.remove('whitesmoke')
            colors3.remove('ghostwhite')
            colors3.remove('black')
            colors3.remove('darkgrey')
            colors3.remove('dimgrey')
            colors3.remove('lightgray')
            colors3.remove('lightgrey')
            colors3.remove('grey')
            colors3.remove('lightslategrey')
            colors3.remove('slategrey')
            colors3.remove('darkslategrey')
            colors_to_repeat = deepcopy(colors3)
            out = ['lightgray', 'blue', 
                   'darkorange', 'deeppink', 
                   'yellow', 'lightskyblue', 
                   'peachpuff', 'purple', 
                   'darkgreen', 'black', 
                   'dimgray'] + colors_to_repeat
            return LinearSegmentedColormap.from_list('ben_colorblind_large_cmap', out)


def get_allen_brain_int_annot_legend(allen_brain_csv_metadata: str) -> dict:
    dat = pd.read_csv(allen_brain_csv_metadata)
    return dict(zip(dat.color_indices, dat.region_names))

def get_allen_brain_rgb_annot_legend(allen_brain_csv_metadata: str) -> dict:
    dat = pd.read_csv(allen_brain_csv_metadata)
    return dict(zip(dat.color_indices, dat.region_names))


# TODO: add subfigures
def plot_jaccard_grid(img_path: str, annot_path: str, annot_legend_dict: dict, df_to_plot: pd.DataFrame,
                      fig_width: int=15, cmap: str | mpl.colors.Colormap = 'viridis'):
    """
    Return a summary grid of unsupervised clustering runs ranked by the jaccard score given the ground truth.
    """
    
    img = ut.read_img(img_path)
    annot = ut.read_img(annot_path)

    # found experimentally, can't help but use a magic number
    PRETTY_HEIGH_TO_WIDTH_RATIO = 45 / len(df_to_plot) 
    
    # fig, axs = plt.subplots(nrows=len(df_to_plot), ncols=4, figsize=(fig_width, fig_width * PRETTY_HEIGH_TO_WIDTH_RATIO)) # img, annot, clustering result, individual classes when their amount is fixed
    fig = plt.figure(figsize=(fig_width, fig_width * PRETTY_HEIGH_TO_WIDTH_RATIO), constrained_layout=True)
    subfigs = fig.subfigures(nrows=len(df_to_plot), ncols=1, hspace=.11)

    clustering_result = df_to_plot.patch_img
    run_hash = df_to_plot.hash
    
    annot_values = np.unique(annot.ravel())
    num_annot_values = len(annot_values)
    annot_min = np.min(annot_values)
    annot_max = np.max(annot_values)

    if isinstance(cmap, str):
        cmap_by_key = get_preconstructed_cmap(cmap)
        cmap_annot = cmap_by_key.resampled(annot_max - annot_min + 1)
    elif isinstance(cmap, mpl.colors.Colormap):
        cmap_annot = cmap.resampled(annot_max - annot_min + 1)
    
    # TODO: make modular and flexible

    for idx, subfig in enumerate(subfigs):

        subfig.suptitle(f'Run {run_hash[idx]}', fontsize='x-large')
        axs_row = subfig.subplots(nrows=1, ncols=4)

        clustering_result_img = ut.read_img(clustering_result[idx])
        total_num_found_clusters = len(np.unique(clustering_result_img.ravel()))

        # remap cluster labels to match ground truth based on Jaccard scores
        run_outdir = os.path.dirname(clustering_result[idx])
        max_pairwise_jacc_csv = os.path.join(run_outdir, 'pairs_max_jacc.csv') # TODO: parameterise
        labels_remapping_dict = ut.get_class_remapping(max_pairwise_jacc_csv)
        clustering_result_img = ut.remap_patchimg(clustering_result_img, labels_remapping_dict)

        # remapping done

        # highlight class with the largest matching score
        labels_remapping_list_sorted = sorted(list(labels_remapping_dict.values()), 
                                              key=lambda x: x[1], reverse=True)
        largest_class, largest_jacc_score = labels_remapping_list_sorted[0]

        largest_class_img = np.zeros(shape=(clustering_result_img.shape[0], clustering_result_img.shape[1]), 
                                     dtype=np.bool_)

        largest_class_mask = clustering_result_img == int(largest_class)
        largest_class_img[largest_class_mask] = True

        # hightlighting done
       

        axs_row[0].imshow(img)
        axs_row[0].axis("off")
        axs_row[0].set_title(f'Image')
    
        ax1 = axs_row[1].imshow(annot, cmap=cmap_annot, vmin=annot_min - 0.5, 
                                vmax=annot_max + 0.5)
        axs_row[1].axis("off")
        axs_row[1].set_title(f'Annotation')
        
        # add metadata legend for annotations
        annot_colors = [ax1.cmap(ax1.norm(annot_value)) for annot_value in annot_values]
        patches = [mpatches.Patch(color=annot_colors[i], label=f'{annot_values[i]} - {annot_legend_dict[i]}') for i in range(len(annot_values))]
        axs_row[1].legend(handles=patches, bbox_to_anchor=(1, 1), loc='best', borderaxespad=0, framealpha=0.6)

        clustering_result_unique = np.unique(clustering_result_img.ravel())
        clustering_result_max = np.max(clustering_result_unique)
        clustering_result_min = np.min(clustering_result_unique)
        
        if isinstance(cmap, str):
            # cmap_clustering = plt.get_cmap(cmap, clustering_result_max - clustering_result_min + 1)
            cmap_by_key = get_preconstructed_cmap(cmap)
            cmap_clustering = cmap_by_key.resampled(clustering_result_max - clustering_result_min + 1)
        elif isinstance(cmap, mpl.colors.Colormap):
            cmap_clustering = cmap.resampled(clustering_result_max - clustering_result_min + 1)

        ax2 = axs_row[2].imshow(clustering_result_img, cmap=cmap_annot, vmin=annot_min - 0.5, 
                                vmax=annot_max + 0.5)
        axs_row[2].axis("off")

        # indicate how many clusters there were before remapping to match annotations
        num_found_clusters = len(np.unique(clustering_result_img.ravel()))
        found_less_clusters_than_classes = f'A total of {total_num_found_clusters - 1} clusters were found' # without background
        found_more_clusters_than_classes = f'Showing top {num_found_clusters - 1} clusters out of {total_num_found_clusters - 1}'
        unsup_clust_title = found_less_clusters_than_classes if total_num_found_clusters < num_annot_values else found_more_clusters_than_classes
        axs_row[2].set_title(f'Unsupervised Clustering\n{unsup_clust_title}')

        divider = make_axes_locatable(axs_row[2])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(ax2, cax=cax, ticks=np.arange(clustering_result_min, clustering_result_max + 1),
                       orientation='vertical', pad=.1, fraction=0.05)

        axs_row[3].imshow(largest_class_img, cmap=cmap_clustering)
        axs_row[3].axis("off")
        axs_row[3].set_title(f'Best matching cluster {largest_class}\nJaccard: {largest_jacc_score}') # this is a placeholder

    return fig

def mpl_figure_to_base64(fig: mpl.figure.Figure) -> bytes:
    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png')
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
    return encoded

def base64img_to_html(encoded_img: bytes) -> str:
    # TODO: remove responsiveness
    return f'<img src=\'data:image/png;base64,{encoded_img}\' class=\"img-fluid\" alt=\"Centered Image\" style=\"max-width: 100%; height: auto;\">'

def mpl_figure_to_html(fig: mpl.figure.Figure) -> str:
    return base64img_to_html(mpl_figure_to_base64(fig))

def main():
    pass

if __name__ == '__main__':
    main()