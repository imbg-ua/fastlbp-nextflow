#!/usr/bin/env python3

import numpy as np
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from PIL import Image
import os
import workflow_utils as ut


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

if __name__ == '__main__':
    main()