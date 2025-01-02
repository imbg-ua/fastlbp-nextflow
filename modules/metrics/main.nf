#!/usr/bin/env/ nextflow

process calculate_unsupervised_clustering_score {
    debug params.debug_flag
    tag "analysis"
    // publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(annot), path(patch_img)

    output:
    tuple val(step_outdir), path(all_pairs_csv), path(pairs_max_csv)

    script:
    // TODO: fix inconsistent namings
    all_pairs_csv = params.all_pairs_csv
    pairs_max_csv = params.pairs_max_csv
    """
    calculate_unsupervised_metrics.py \
        --patch_img_path ${patch_img} \
        --annotation_path ${annot} \
        --savefile_all_pairs ${all_pairs_csv} \
        --savefile_max_pairs ${pairs_max_csv}
    """
}