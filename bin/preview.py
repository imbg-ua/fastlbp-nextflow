#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import fire

def generate_img_thumbnail(img_path: str, max_size: tuple[int] = (500, 500)):
    print(f"creating thumbnail for {img_path}")
    input_img = Image.open(img_path)
    input_size = input_img.size
    input_img.thumbnail(max_size)
    input_img.save('thumbnail.jpg')
    print(f"Done. Original size {input_size}, thumbnail size {max_size}.")

def generate_labels_preview(patchsize: int, img_path: str, np_labels_path: str, max_size: tuple[int] = (500, 500)):
    print(f"creating labels preview for {np_labels_path}")
    with Image.open(img_path) as input_img:
        input_size = input_img.size
    patch_shape = input_size[0] // patchsize, input_size[1] // patchsize
    labels = np.load(np_labels_path).reshape(patch_shape)
    fig = plt.figure()
    plt.imshow(labels)
    fig.savefig('labels.jpg', dpi=fig.dpi)
    plt.close()
    print(f"Done. Original size {input_size}, thumbnail size {max_size}.")


if __name__ == "__main__":
    fire.Fire()
