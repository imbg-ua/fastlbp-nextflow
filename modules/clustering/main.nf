#!/usr/bin/env/ nextflow

// TODO: add workflows for the case when clustering is not the final step
// (RunClusteringAndPrepareForNextStep)

params.args.clustering.method = 'k_means' // default clustering method
params.debug_flag = true

// uses image basename as run id
process clustering_publish {
    tag "${img_id}"
    debug params.debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple val(img_id), path(data), val(params_str)

    output:
    tuple val(img_id), path("clustering_labels.npy")

    script:
    """
    run_clustering.py \
        --np_data_path ${data} \
        --params_str "${params_str}"
    """
}

process clustering {
    debug debug_flag
    tag "${params.args.clustering.method}"

    input:
    tuple val(run_id), path(data), val(params_str)

    output:
    tuple val(run_id), path("clustering_labels.npy")

    script:
    """
    run_clustering.py \
        --np_data_path ${data} \
        --params_str "${params_str}"
    """
}


workflow RunClustering {
    take:
    clustering_inputs_channel

    main:
    clustering(clustering_inputs_channel)

    emit:
    all_clustering_outputs = clustering.out
}


workflow {
    // TODO: read needed parameters from tsv to run as a separate module
    RunClustering()
}