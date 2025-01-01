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

from typing import List


PLOTTING_BACKENDS = ('matplotlib', 'plotly')

# TODO: change `discrete image` to `label image`
def shrink_discrete_img_to_uniform_range(img_arr: np.array) -> np.array:
    """Remap values in an image with integer labels to a uniform range `[0, n]` where
    `n` is the total number of unique values.
    """

    all_values = np.unique(img_arr)

    img_arr_uniform = np.empty(shape=(img_arr.shape[0], img_arr.shape[1]), dtype=np.uint16)
    for idx, val in enumerate(all_values):
        val_mask = img_arr == val
        img_arr_uniform[val_mask] = idx

    return img_arr_uniform, all_values

# TODO: improve function naming
def plot_img_discrete_uniform_and_create_colorbar_axis(img_arr: np.array, ax: mpl.axes.Axes, cmap: mpl.colors.Colormap = 'viridis'):
    img_arr_uniform, all_values = shrink_discrete_img_to_uniform_range(img_arr)
    total_values_num = len(all_values)
    cmap = plt.get_cmap('viridis', total_values_num) # cmap restricted to the total number of clusters

    # set limits 0.5 outside the true range
    img = ax.imshow(img_arr_uniform, cmap=cmap, vmin=-.5, vmax=total_values_num - .5)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05) # TODO: make tuneable

    return img, cax, all_values

# TODO: replace and/or combine with `plot_img_discrete_uniform_and_create_colorbar_axis()`
# generic function to plot images with discrete colorbars
def plot_img_discrete(img_arr: np.array, figsize: tuple=(8, 8)):
    """
    Display image with discrete values.

    Parameters
    ----------
    img_arr : np.array
        Image as numpy array.

    figsize : tuple, default: (8, 8)
        Size of the matplotlib figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure with the displayed image.
    """
    # get discrete colormap
    img_arr = img_arr.copy().astype('int') # TODO: don't copy
    all_values = np.unique(img_arr)
    total_clust_num = len(all_values)
    cmap = plt.get_cmap('viridis', total_clust_num) # cmap restricted to the total number of clusters
    fig, ax = plt.subplots(figsize=figsize)

    # FIXME: remap cluster values to a uniform range to simplify visualization and not 
    # implement a custom colormap normalization method
    # TODO: implement a custom colormap normalization method
    # TODO: use `shrink_discrete_img_to_uniform_range()` function
    img_arr_uniform = np.empty(shape=(img_arr.shape[0], img_arr.shape[1]), dtype=np.uint16)
    for idx, val in enumerate(all_values):
        val_mask = img_arr == val
        img_arr_uniform[val_mask] = idx        

    # set limits 0.5 outside the true range
    img = ax.imshow(img_arr_uniform, cmap=cmap, vmin=-.5, vmax=total_clust_num - .5)

    # tell the colorbar to tick at integers
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(img, cax=cax, ticks=np.arange(total_clust_num),
                    orientation='vertical', pad=.1, fraction=0.05)

    cbar.ax.set_yticklabels(all_values) # set real cluster labels as ticks

    return fig




def dimred_2D_html(dimred_path: str, labels_path: str, component_1: int=0, component_2: int=1, backend: str='plotly', *args, **kwargs):
    assert backend in PLOTTING_BACKENDS, f'Plotting backend {backend} is not supported. Please \
    choose from the available options: {PLOTTING_BACKENDS}'

    match backend:
        case 'plotly':
            return dimred_2D_html_plotly(dimred_path, labels_path,
                                             component_1, component_2,
                                             *args, **kwargs)
        
        case 'matplotlib':
            return dimred_2D_html_matplotlib(dimred_path, labels_path,
                                             component_1, component_2,
                                             *args, **kwargs)
        
def dimred_2D_html_matplotlib(dimred_path: str, labels_path: str, component_1: int=0, component_2: int=1) -> str:
    dimred = np.load(dimred_path)
    labels = np.load(labels_path)

    data = dimred_and_labels_to_dataframe(dimred, labels, component_1, component_2)
    fig = plot_dimred_2D_labelled_matplotlib(data)

    return mpl_figure_to_html(fig)

def dimred_2D_html_plotly(dimred_path: str, labels_path: str, component_1: int=0, component_2: int=1,
                   include_plotlyjs: str | bool = 'cdn') -> str:
    dimred = np.load(dimred_path)
    labels = np.load(labels_path)

    data = dimred_and_labels_to_dataframe(dimred, labels, component_1, component_2)
    fig = plot_dimred_2D_labelled(data)
    plot_html = fig.to_html(full_html = False, include_plotlyjs=include_plotlyjs)

    return plot_html

# --- image plotting function --- #
def img_html(img_path: str, backend: str = 'plotly', *args, **kwargs) -> str:
    """
    Available backends:
    """

    assert backend in PLOTTING_BACKENDS, f'Plotting backend {backend} is not supported. Please \
    choose from the available options: {PLOTTING_BACKENDS}'

    match backend:
        case 'plotly':
            return img_html_plotly(img_path, *args, **kwargs)
        case 'matplotlib':
            return img_html_matplotlib(img_path, *args, **kwargs)


def img_html_matplotlib(img_path: str) -> str:
    img = ut.read_img(img_path)
    fig = plot_img_discrete(img) # TODO: adjust function and use in other places too

    fig_html = mpl_figure_to_html(fig)

    plt.close(fig)

    return fig_html

def img_html_plotly(img_path: str, include_plotlyjs: str | bool = 'cdn') -> str:
    img = ut.read_img(img_path)
    fig = plot_img(img)
    plot_html = fig.to_html(full_html = False, include_plotlyjs=include_plotlyjs)

    return plot_html


def img_downscaled_html(img_path: str, backend: str='plotly', *args, **kwargs) -> str:
    assert backend in PLOTTING_BACKENDS, f'Plotting backend {backend} is not supported. Please \
    choose from the available options: {PLOTTING_BACKENDS}'

    match backend:
        case 'plotly':
            return img_downscaled_html_plotly(img_path, *args, **kwargs)
        case 'matplotlib':
            return img_downscaled_html_matplotlib(img_path, *args, **kwargs)


def img_downscaled_html_plotly(img_path: str, size=(400, 400), include_plotlyjs: str | bool = 'cdn') -> str:
    img = ut.read_img(img_path)
    fig = plot_img_downscaled(img, size) # TODO: rename plot_img_downscaled -> plot_img_downscaled_plotly
    plot_html = fig.to_html(full_html = False, include_plotlyjs=include_plotlyjs)

    return plot_html

def img_downscaled_html_matplotlib(img_path: str, size=(400, 400)) -> str:
    img = ut.read_img(img_path)
    img_downscaled = downscale_image(img)
    fig, ax = plt.subplots()
    ax.imshow(img_downscaled)
    return mpl_figure_to_html(fig)


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


def plot_dimred_2D_labelled_matplotlib(dat: pd.DataFrame):
    fig, ax = plt.subplots()
    ax.scatter(x=dat['component_1'], y=dat['component_2'], c=dat['label'])
    ax.set_xlabel('Component 1')
    ax.set_ylabel('Component 2')
    fig.suptitle('Projections')
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


def downscale_image(img_arr: np.array, size: tuple=(400, 400)) -> np.array:
    img_pil = Image.fromarray(img_arr)
    img_pil.thumbnail(size)
    img_pil = np.asarray(img_pil)

    return img_pil

# TODO: downscale without PIL
def plot_img_downscaled(img_arr: np.array, size=(400, 400)):
    img_downscaled = downscale_image(img_arr)
    fig = px.imshow(img_downscaled)
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


def read_annotation_legend(annotation_legend_csv: str, class_value_column: str="id", 
                          class_name_column: str="name") -> dict:
    """Read annotation legend from a csv file.

    The csv file **must** contain 2 mandatory columns:
        - `id` with the integer value of the class
        - `name` name of the corresponding class
    
    Alternatively, you can explicitly specify other names for the `id` and `name` columns using function parameters.

    Parameters
    ----------
    annotation_legend_csv : str
        Path to a csv file with the annotation legend.

    class_value_column : str, default: 'id'
        Name of the column storing integer values representing 
        different classes in the annotations. 
    
    class_name_column : str, default: 'name'
        Name of the column storing class names.

    Returns
    -------
    result : dict
        Python dictionary with the annotation legend. 
        Keys are class integer values, values are names of the classes.

    """

    dat = pd.read_csv(annotation_legend_csv)
    return dict(zip(dat[class_value_column], dat[class_name_column]))
    
# TODO: replace with a generic function `read_annotation_legend`
def get_allen_brain_int_annot_legend(allen_brain_csv_metadata: str) -> dict:
    dat = pd.read_csv(allen_brain_csv_metadata)
    return dict(zip(dat.color_indices, dat.region_names))

# TODO: delete as not used
def get_allen_brain_rgb_annot_legend(allen_brain_csv_metadata: str) -> dict:
    dat = pd.read_csv(allen_brain_csv_metadata)
    return dict(zip(dat.color_indices, dat.region_names))

    

# TODO: add subfigures for each run
def plot_jaccard_grid(img_path: str, annot_path: str, clustering_imgs: List[str], run_hashes: List[str], annot_legend_dict: dict | None = None,
                      class_mapping_csv: str='pairs_max_jacc.csv',
                      fig_width: int=15, cmap: str | mpl.colors.Colormap = 'viridis'):
    """
    Return a summary grid of the unsupervised clustering results mapped to annotation with multiple classes. 
    The runs are ranked by the Jaccard score based on the greedy cluster matching algorithm described in `workflow_utils.py`
    """

    assert len(clustering_imgs) == len(run_hashes), "Length of the array of clustering results doesn't match length of the array of run hashes."

    num_runs = len(run_hashes)
    
    img = ut.read_img(img_path)
    annot = ut.read_img(annot_path)

    # found experimentally, the rationale is that each row subfigure contains 4 image plots
    # TODO: find a better ratio?
    PRETTY_HEIGHT_TO_WIDTH_RATIO = num_runs * 0.25
    
    # fig, axs = plt.subplots(nrows=num_runs, ncols=4, figsize=(fig_width, fig_width * PRETTY_HEIGHT_TO_WIDTH_RATIO)) # img, annot, clustering result, individual classes when their amount is fixed
    fig = plt.figure(figsize=(fig_width, fig_width * PRETTY_HEIGHT_TO_WIDTH_RATIO), constrained_layout=True)
    subfigs = fig.subfigures(nrows=num_runs, ncols=1, hspace=.11)
    
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

        subfig.suptitle(f'Run {run_hashes[idx]}', fontsize='x-large')
        axs_row = subfig.subplots(nrows=1, ncols=4)

        clustering_result_img = ut.read_img(clustering_imgs[idx])
        total_num_found_clusters = len(np.unique(clustering_result_img.ravel()))

        # remap cluster labels to match ground truth based on Jaccard scores
        run_outdir = os.path.dirname(clustering_imgs[idx])

        max_pairwise_jacc_csv = os.path.join(run_outdir, class_mapping_csv)
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
        if annot_legend_dict:
            # TODO: add dict key type check or choose a dtype to cast to in read_annotation_legend()
            patches = [mpatches.Patch(color=annot_colors[i], label=f'{annot_values[i]} - {annot_legend_dict[annot_values[i]]}') for i in range(len(annot_values))]
        else:
            patches = [mpatches.Patch(color=annot_colors[i], label=f'Class {annot_values[i]}') for i in range(len(annot_values))]

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
        num_found_clusters = len(clustering_result_unique)
        found_less_clusters_than_classes = f'A total of {total_num_found_clusters - 1} clusters were found' # without background
        found_same_num_of_clusters_as_classes = f'Found {total_num_found_clusters - 1} clusters with the following matching'
        found_more_clusters_than_classes = f'Showing top {num_found_clusters - 1} clusters out of {total_num_found_clusters - 1}'

        if total_num_found_clusters < num_annot_values:
            unsup_clust_title = found_less_clusters_than_classes
        elif total_num_found_clusters > num_annot_values:
            unsup_clust_title = found_more_clusters_than_classes
        else:
            unsup_clust_title = found_same_num_of_clusters_as_classes

        axs_row[2].set_title(f'Unsupervised Clustering\n{unsup_clust_title}')
        # TODO: make discrete colorbar show actual number of clusters
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