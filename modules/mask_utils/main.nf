#!/usr/bin/env/ nextflow

include { getValueFromParamList } from '../../lib/nf/utils'

process get_tissue_mask {
    tag "preprocessing"
    debug params.debug_flag
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple path(img), val(background_type)

    output:
    tuple path(img), path('pixelmask.npy')

    script:
    background_type_param = background_type ? "--background ${background_type}" : ""
    """
    get_mask.py get_mask \
        --img_path ${img} \
        ${background_type_param}
    """
}


process downscale_mask {
    debug params.debug_flag
    tag "preprocessing"
    // publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(img), path(pixelmask), val(patchsize)

    output:
    tuple val(step_outdir), path(patchmask_savefile)

    script:
    patchmask_savefile = "patchmask.npy"
    """
    get_mask.py downscale_using_patchsize \
        --img_path ${pixelmask} \
        --patchsize ${patchsize} \
        --savefile ${patchmask_savefile}
    """
}

process convert_annotations_to_binmask {
    tag "mask preprocessing"
    debug params.debug_flag
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple path(img), path(mask), val(background_val_str)

    output:
    tuple path(img), path(binmask)

    script:
    binmask = "binmask.npy"
    background_val = background_val_str ? "--background_val_str \"${background_val_str}\"" : ""
    """
    get_mask.py check_annotations \
        --annotations ${mask} \
        ${background_val} \
        --savefile ${binmask}
    """

}


workflow GetTissueMask {
    take:
    get_tissue_mask_inputs_channel

    main:
    get_tissue_mask(get_tissue_mask_inputs_channel)

    get_tissue_mask.out
        .map { imgg, pixelmask -> pixelmask }
        .set { pixel_mask_only }

    emit:
    img_and_pixel_mask = get_tissue_mask.out
    pixel_mask_only
}



workflow GetTissuePatchMask {
    take:
    get_tissue_mask_inputs_channel
    lbp_combinations_hash_outdir

    main:
    GetTissueMask(get_tissue_mask_inputs_channel)
    img_and_pixel_mask = GetTissueMask.out.img_and_pixel_mask
    pixel_mask_only = GetTissueMask.out.pixel_mask_only

    pixel_mask_only
        .set { mask_for_report }

    // combine the pixel mask with the lbp parameters channel
    // to get the patchsize value which is needed to properly downscale 
    // pixel mask into patch mask
    img_and_pixel_mask
        .combine(lbp_combinations_hash_outdir)
        .map { imgg, pixelmask, lbp_outdir, lbp_params ->
        tuple(lbp_outdir, imgg, pixelmask, getValueFromParamList(lbp_params, param_name = 'patchsize')) }
        .set { downscale_me }

    downscale_mask(downscale_me)

    downscale_mask.out
        .set { all_downscale_outputs }

    downscale_mask.out
        .combine(lbp_combinations_hash_outdir, by:0)
        .combine(img_and_pixel_mask)
        .map { downscale_outdir, patchmask_path, lbp_params, imgg, pixelmask_path ->
        tuple(downscale_outdir, imgg, pixelmask_path, patchmask_path, lbp_params) }
        .set { feed_me_into_lbp }

    emit:
    feed_me_into_lbp
    all_downscale_outputs
    mask_for_report

}

workflow GetTissuePatchMaskFromAnnnotations {
    take:
    convert_to_binmask_ch
    lbp_combinations_hash_outdir

    main:
    convert_annotations_to_binmask(convert_to_binmask_ch)

    convert_annotations_to_binmask.out
        .set { converted_mask_ch }

    convert_annotations_to_binmask.out
        .combine(lbp_combinations_hash_outdir)
        .map { imgg, maskk, lbp_outdirr, lbp_paramss ->
        tuple(lbp_outdirr, imgg, maskk, getValueFromParamList(lbp_paramss, param_name = 'patchsize')) }
        .set { downscale_me }

    downscale_mask(downscale_me)

    downscale_mask.out
        .set { all_downscale_outputs }

    downscale_mask.out
        .combine(lbp_combinations_hash_outdir, by:0)
        .combine(converted_mask_ch)
        .map { downscale_outdir, patchmask_path, lbp_params, imgg, binmaskk ->
        tuple(downscale_outdir, imgg, binmaskk, patchmask_path, lbp_params) }
        .set { feed_me_into_lbp }

    emit:
    feed_me_into_lbp
    all_downscale_outputs

}