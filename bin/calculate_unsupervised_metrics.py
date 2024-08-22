#!/usr/bin/env python3

import workflow_utils as ut
import fire

def calculate_jaccard_score_for_run(patch_img_path: str, 
                                    annotation_path: str, 
                                    savefile_all_pairs: str = 'all_pairs_jacc.csv',
                                    savefile_max_pairs: str = 'pairs_max_jacc.csv'):

    ut.run_paiwise_max_jacc(gt_path=annotation_path, 
                            patch_labels_path=patch_img_path, 
                            savefile_all_pairs=savefile_all_pairs,
                            savefile_max_pairs=savefile_max_pairs)

if __name__ == '__main__':
    fire.Fire(calculate_jaccard_score_for_run)
