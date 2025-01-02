#!/usr/bin/env/ nextflow

// TODO: add default parameter values here


process fastlbp {
    debug params.debug_flag
    tag "fastlbp"
    // publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(img), path(mask), path(patchmask), val(params_str)

    output:
    tuple val(step_outdir), path("data"), emit: lbp_data_folder
    tuple val(step_outdir), path("data/out/${file(params.constargs.lbp.outfile_name).getBaseName()}_flattened.npy"), emit: lbp_result_file_flattened
    tuple val(step_outdir), path("data/out/${file(params.constargs.lbp.outfile_name).getBaseName()}.npy"), emit: lbp_result_file_img

    script:
    mask_path = mask ? "--img_mask ${mask}" : ""
    patchmask_path = patchmask ? "--patch_mask ${patchmask}" : ""
    """
    run_lbp.py main_grid_search \
        --img_path ${img} \
        --params_str "${params_str}" \
        --ncpus ${params.constargs.lbp.ncpus} \
        --outfile_name ${params.constargs.lbp.outfile_name} \
        ${mask_path} \
        ${patchmask_path} \
        --img_name ${params.constargs.lbp.img_name}
    """
}


workflow RunFastLBP {
    take:
    fastlbp_inputs_channel

    main:
    fastlbp(fastlbp_inputs_channel)

    emit:
    all_lbp_outputs = fastlbp.out.lbp_data_folder
    lbp_result_flattened = fastlbp.out.lbp_result_file_flattened
    lbp_result_img = fastlbp.out.lbp_result_file_img
}

workflow RunFastLBPAndPrepareForNextStep {
    take:
    fastlbp_inputs_channel
    next_step_combinations_map

    main:
    RunFastLBP(fastlbp_inputs_channel)

    all_lbp_outputs = RunFastLBP.out.all_lbp_outputs
    lbp_result_flattened = RunFastLBP.out.lbp_result_flattened
    lbp_result_img = RunFastLBP.out.lbp_result_img

    lbp_result_flattened
        .combine(next_step_combinations_map)
        .map { lbp_outdir, lbp_result_flattened_cur, next_step_params_str, next_step_params ->
        tuple("${lbp_outdir}/${next_step_params_str}", lbp_result_flattened_cur, next_step_params) }
        .set {next_step_inputs_channel}

    emit:
    all_lbp_outputs
    next_step_inputs_channel
    lbp_result_img
}


