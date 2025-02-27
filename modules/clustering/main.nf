#!/usr/bin/env/ nextflow

// TODO: add workflows for the case when clustering is not the final step
// (RunClusteringAndPrepareForNextStep)

// params.args.clustering.method = 'k_means' // default clustering method
// params.debug_flag = true

// it is only used in the tag in the clustering process, so maybe the default value is not needed at all here
// params.args = [:]
// params.args.clustering = [:]
// params.args.clustering.method = 'k_means'

if ( !params.args )
    params.args = [:]
if ( !params.args.clustering )
    params.args.clustering = [:]
if ( !params.args.clustering.method )
    params.args.clustering.method = 'k_means'

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
    debug params.debug_flag
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

// TODO: merge with the workflow above to have a single entry point?
workflow RunClusteringPublish {
    take:
    clustering_publish_inputs_channel

    main:
    clustering_publish(clustering_publish_inputs_channel)

    emit:
    all_clustering_outputs = clustering_publish.out
}


workflow {
    // TODO: move binary scripts the corresponding module bin folders

    // TODO: describe required input tsv structure in the module README
    def inputs_tsv_ch = Channel.fromPath(params.inputs_tsv)

    // required tsv file structure: run_id, input_np_array, params_str
    inputs_tsv_ch
        .splitCsv(header:true, sep:'\t')
        .map { row -> 
        tuple( row.run_id, row.input_np_array, row.params_str )
        }
        .set{ module_inputs_ch }

    RunClusteringPublish(module_inputs_ch)
}
