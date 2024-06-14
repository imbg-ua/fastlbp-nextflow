#!/usr/bin/env python3

import pandas as pd
import os
import workflow_utils as ut
import plotting_utils as pt
import fire
from jinja2 import Environment, FileSystemLoader


def params_str_to_list_of_dicts(params_str: str) -> str:
    # in metadata tsv semicolons are used instead of commas
    params_str = params_str.replace(';', ',')
    steps = params_str.split('/') # separates pipeline steps
    steps = list(map(lambda x: ut.parse_params_str(x), steps))
    return steps

def parse_metadata(tsv: str) -> pd.DataFrame:
    dat = pd.read_csv(tsv, sep='\t')
    for column in dat.columns:
        dat[f'{column}'] = dat[f'{column}'].astype('string')

    dat['parameters_list_of_dicts'] = dat['parameters'].apply(params_str_to_list_of_dicts)
    return dat


# TODO: use dataclass
def get_runs_data(outdir: str, dimred_outname: str="embeddings.npy",
                  clustering_outname: str="clustering_labels.npy",
                  patchimg_outname: str="patch_labels.npy"):
    # TODO: remove hardcoded values
    hash_runs = []
    umap_runs = []
    hdbscan_runs = []
    patchimg_runs = []

    for run_dir in os.listdir(outdir):
        if os.path.isdir(abs_run_dir := os.path.join(outdir, run_dir)):
            hash_runs.append(run_dir)
            umap_runs.append(os.path.join(abs_run_dir, dimred_outname))
            hdbscan_runs.append(os.path.join(abs_run_dir, clustering_outname))
            patchimg_runs.append(os.path.join(abs_run_dir, patchimg_outname))

    res_dict = { 
        'hash': hash_runs,
        'dimred': umap_runs,
        'clustering': hdbscan_runs,
        'patch_img': patchimg_runs
    }

    return pd.DataFrame(res_dict, dtype='string')

def process_outdir(outdir: str, metadata_tsv: str='hash_to_combination.tsv') -> pd.DataFrame:
    runs_data = get_runs_data(outdir)

    metadata_abs_path = os.path.join(outdir, metadata_tsv)
    metadata = parse_metadata(metadata_abs_path)

    # contains hash, paths to the output files and params dict
    runs_info = runs_data.merge(metadata, on='hash')

    runs_info['dimred_html'] = runs_info.apply(lambda dat: pt.dimred_2D_html(dat.dimred, dat.clustering), axis=1)
    runs_info['patchimg_html'] = runs_info.apply(lambda dat: pt.img_html(dat.patch_img), axis=1)

    # TODO: hardcoded values
    essential_info = runs_info[['hash', 'dimred_html', 'patchimg_html', 'parameters_list_of_dicts']]

    return essential_info


def render_template(img_path: str, runs_info: pd.DataFrame,
                    annot_path: str='', template: str='template_report.html.jinja') -> str:
    pipeline_bin_dir = os.path.dirname(os.path.abspath(__file__))
    template_abs_path = os.path.join(pipeline_bin_dir, template)

    # jinja2 setup
    env = Environment(loader=FileSystemLoader(pipeline_bin_dir))
    template = env.get_template(template) # FIXME: or template_abs_path?

    # generate downscaled image plot
    img = pt.img_downscaled_html(img_path)

    # generate downscaled annotations if possible
    annot = None
    if annot_path:
        annot = pt.img_downscaled_html(annot_path)
    
    runs_info_list = list(runs_info.T.to_dict().values())

    html_out = template.render(img=img, annot=annot, runs_info_list=runs_info_list)

    return html_out


# TODO: add other default parameters
def create_report(outdir: str, img_path: str, 
                  metadata_tsv: str='hash_to_combination.tsv',
                  annot_path: str='', 
                  template: str='template_report.html.jinja',
                  savefile: str='report.html') -> None:
    
    runs_info = process_outdir(outdir, metadata_tsv)
    report_html = render_template(img_path, runs_info, annot_path, template)

    with open(savefile, 'w') as f:
        f.write(report_html)

if __name__ == '__main__':
    fire.Fire(create_report)