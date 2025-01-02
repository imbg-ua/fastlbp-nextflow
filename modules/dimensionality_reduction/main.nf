#!/usr/bin/env/ nextflow



process dimred {
    debug debug_flag
    tag "${params.args.dimred.method}"
    // publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(lbp_result_flattened), val(params_str)

    output:
    tuple val(step_outdir), path("embeddings.npy")

    script:
    """
    run_dimensionality_reduction.py \
        --np_data_path ${lbp_result_flattened} \
        --params_str "${params_str}"
    """
}


// shorthand for Dimensinality Reduction
workflow RunDimRed {
    take:
    dimred_inputs_channel

    main:
    dimred(dimred_inputs_channel)

    emit:
    all_dimred_outputs = dimred.out
}

workflow RunDimRedAndPrepareForNextStep {
    take:
    dimred_inputs_channel
    next_step_combinations_map

    main:
    RunDimRed(dimred_inputs_channel)

    all_dimred_outputs = RunDimRed.out.all_dimred_outputs

    all_dimred_outputs
        .combine(next_step_combinations_map)
        .map { dimred_outdir, embeddings, next_step_params_str, next_step_params ->
        tuple("${dimred_outdir}/${next_step_params_str}", embeddings, next_step_params) }
        .set { next_step_inputs_channel }

    emit:
    all_dimred_outputs
    next_step_inputs_channel
}