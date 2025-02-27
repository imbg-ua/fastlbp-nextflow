#!/usr/bin/env/ nextflow


include { createCombinations; 
          checkNestedParameterCombinations; 
          getFinalOutdirFromParamsCombinations; 
          getAllOutputFilesFromParamsCombinations } from '../../lib/nf/utils'


// generic utils module

// == deconvolve labels to patch image == //

process labels_to_patch_img {
    debug params.debug_flag
    tag "generating segmented image"
    // publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(labels_data), path(patchmask), path(lbp_output)

    output:
    tuple val(step_outdir), path(savefile)

    script:
    savefile = "patch_labels.npy"
    command = patchmask ? "masked_labels_to_patch_img --np_mask_path ${patchmask}" : \
    "labels_to_patch_img --lbp_output_path ${lbp_output}"
    """
    deconvolve_masked_labels.py ${command} \
        --np_labels_path ${labels_data} \
        --savefile ${savefile}
    """
}


workflow LabelsToPatchImage {
    take:
    labels_to_patch_img_inputs_channel

    main:

    labels_to_patch_img(labels_to_patch_img_inputs_channel)

    emit:
    all_deconvolved_images_outputs = labels_to_patch_img.out
}

// == copy files process == //

process copy_files_to_target_dir {
    tag "wrapping up"
    debug params.debug_flag
    publishDir "${outdir}", mode: "copy" // TODO: switch to `move`?

    input:
    tuple val(outdir), path(data)

    output:
    tuple val(outdir), path(data)
    
    script:
    """
    """
}

//  == grid search combinations utils == //

process combinations_metadata_to_tsv {
    tag "wrapping up"
    debug params.debug_flag
    publishDir "${params.outdir}", mode: "copy"

    input:
    val(params_and_hash_str)

    output:
    path(metadata_csv)

    script:
    metadata_csv = "hash_to_combination.tsv"
    """
    params_to_hash.py \
        --params_list_str "${params_and_hash_str}" \
        --savefile ${metadata_csv}
    """
}

// TODO: implement repeating parts as subworkflows
workflow GetParameterCombinations {
    main:
    
    // TODO: don't resort to config params inside the workflow

    // TODO: find a more concise way to do this
    // check if params for each step are given explicitly or as a path to a tsv file
    if ( params.args_tsv ) {
        if ( params.args_tsv.lbp ) {
            def tsv_lbp_args = Channel.fromPath(params.args_tsv.lbp)
            tsv_lbp_args
                .splitCsv(header:true, sep:'\t')
                .map { row -> row.collect { param_k, param_v -> tuple(param_k, param_v) } }
                .set { lbp_combinations }
        } else {
            def lbp_step = params.args.lbp
            lbp_combinations = createCombinations(lbp_step)
        }

        if ( params.args_tsv.dimred ) {
            def tsv_umap_args = Channel.fromPath(params.args_tsv.dimred)
            tsv_umap_args
                .splitCsv(header:true, sep:'\t')
                .map { row -> row.collect { param_k, param_v -> tuple(param_k, param_v) } }
                .set { umap_combinations }
        } else {
            def umap_step = params.args.dimred
            umap_combinations = createCombinations(umap_step)
        }
        if ( params.args_tsv.clustering ) {
            def tsv_hdbscan_args = Channel.fromPath(params.args_tsv.clustering)
            tsv_hdbscan_args
                .splitCsv(header:true, sep:'\t')
                .map { row -> row.collect { param_k, param_v -> tuple(param_k, param_v) } }
                .set { hdbscan_combinations }
        } else {
            def hdbscan_step = params.args.clustering
            hdbscan_combinations = createCombinations(hdbscan_step)
        }

    } else {
        // no tsv arguments provided
        lbp_combinations = createCombinations(params.args.lbp)
        umap_combinations = createCombinations(params.args.dimred)
        hdbscan_combinations = createCombinations(params.args.clustering)
    }

    emit:
    lbp_combinations
    umap_combinations
    hdbscan_combinations

}

workflow ParameterCombinationsToMap {
    take:
    parameter_combinations_list

    main:

    parameter_combinations_list
        .map { it -> 
        tuple("${it.toString()}", it) }
        .set { parameter_combinations_map }

    emit:
    parameter_combinations_map
}

workflow AllCombinationsToHashes {
    take:
    all_input_channels // each channel must have the following structure: [step_param_list_str, step_output]

    main:

    // combine all channels to get all possible parameter combinations used in different steps
    def first_channel = all_input_channels[0]
    def residual_channels = all_input_channels[1..-1]

    def combined_channel = residual_channels
                                .inject(first_channel, { accumulated_channels, curr_channel ->
                                accumulated_channels.combine(curr_channel) 
                                }
                            )

    // filter the resulting combinations to retain the final parameter combination strings
    // i.e. the channel currently contains combinations for lbp, for lbp+umap, for lbp+umap+hdbscan
    // we want to leave only the entries describing the full pipeline, that is only lbp+umap+hdbscan
    combined_channel
        .filter { outputs -> checkNestedParameterCombinations(outputs) }
        .map { outputs -> 
        tuple(getFinalOutdirFromParamsCombinations(outputs).replaceAll(',', ';'), 
        getAllOutputFilesFromParamsCombinations(outputs), 
        getFinalOutdirFromParamsCombinations(outputs).replaceAll(',', ';').md5()) }
        .multiMap { final_combination_for_outdir, all_outputs_tuple, final_combination_for_outdir_hash ->
        // replace commas with semicolons to distinguish between inner parameters list and
        // outer combination-to-hash list 
        data_to_save: tuple(final_combination_for_outdir_hash, all_outputs_tuple)
        dir_and_hash: tuple(final_combination_for_outdir, final_combination_for_outdir_hash) }
        .set { final_data_and_folders }

    emit:
    data_to_save = final_data_and_folders.data_to_save
    dir_and_hash = final_data_and_folders.dir_and_hash
}
