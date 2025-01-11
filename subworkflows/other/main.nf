#!/usr/bin/env/ nextflow

// this workflow works correctly only with the fair directive is set to true
params.fair = true
params.lbp_runs_tsv = ''

include { checkNextflowVersion;
          get_param_value_from_param_str } from '../../lib/nf/utils'

// TODO: make sure there is no circular dependency
include { convert_annotations_to_binmask;
          fastlbp;
          get_tissue_mask;
          downscale_mask } from '../../subworkflows/normal'

nextflow.enable.dsl = 2
pipeline_version = "0.0.2"

checkNextflowVersion()



// TODO: this doesn't allow processing same images with different lbp parameters
// well, it allows, but the outdirs are the same and outputs for different runs interfere 
// and get overwritten by each other. To fix this, need to create run hashes and save outputs 
// similarly to the grid search mode
workflow MultiImageLBPTSV {
    if ( !params.lbp_runs_tsv ) {
        log.error("""
        Provide lbp_runs_tsv parameter to run the pipeline in the LBP-only tsv mode. 
        """)
    }

    def lbp_runs_tsv = Channel.fromPath(params.lbp_runs_tsv)
    lbp_runs_tsv
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple( row.image, row.mask, row.background_color, row.lbp_params_str)
        }
        .set{ lbp_runs_params }

    lbp_runs_params
        .branch { 
            without_mask: (it[1] == "" || it[1] == "no")
                return tuple(it[0], it[3]) // (image, lbp_params_str)
            otsu_mask: it[1] == "auto"
                return it // the whole tuple
            with_mask: (it[1] != "" && it[1] != "no")
                return it // the whole tuple
        }
        .set { lbp_runs_params_split_by_mask }


    // ------------------------- //
    // process images with masks //
    // ------------------------- //
    lbp_runs_params_split_by_mask.with_mask
        .multiMap { image, mask, background_color, lbp_params_str ->
            convert_to_binmask_ch: tuple(image, mask, background_color)
            lbp_params_ch: lbp_params_str
        }
        .set { convert_me_to_binmask }
    
    convert_annotations_to_binmask(convert_me_to_binmask.convert_to_binmask_ch)
    
    // patchsize parameter value from parameter string extraction pattern
    def patchsize_pattern = /patchsize, (\d+)/

    convert_annotations_to_binmask.out
        .merge(convert_me_to_binmask.lbp_params_ch)
        .map { imgg, binmaskk, lbp_params_strr ->
        tuple(imgg, binmaskk, get_param_value_from_param_str(lbp_params_strr, patchsize_pattern)) }
        .set { downscale_me_with_mask }
    
    // downscale_mask(downscale_me_with_mask)

    // downscale_mask.out
    //     .merge(convert_me_to_binmask.lbp_params_ch)
    //     .set { feed_me_into_lbp }

    // fastlbp(feed_me_into_lbp)

    // ---------------------------- //
    // process images without masks //
    // ---------------------------- //
    lbp_runs_params_split_by_mask.without_mask
        .map { image, lbp_params_str ->
        tuple(image, [], [], lbp_params_str) }
        .set { feed_me_into_lbp_no_mask }

    // fastlbp(feed_me_into_lbp_no_mask)

    // ---------------------------------- //
    // process images using Otsu's method //
    // ---------------------------------- //

    lbp_runs_params_split_by_mask.otsu_mask
        .multiMap { image, mask, background_color, lbp_params_str ->
            get_otsu_mask_ch: tuple(image, background_color)
            lbp_params_otsu_ch: lbp_params_str
        }
        .set { img_to_get_masks }

    get_tissue_mask(img_to_get_masks.get_otsu_mask_ch)

    get_tissue_mask.out
        .merge(img_to_get_masks.lbp_params_otsu_ch)
        .map { imgg, binmaskk, lbp_params_strr ->
        tuple(imgg, binmaskk, get_param_value_from_param_str(lbp_params_strr, patchsize_pattern)) }
        .set { downscale_me_otsu }
    
    // downscale_mask(downscale_me_otsu)

    // downscale_mask.out
    //     .merge(img_to_get_masks.lbp_params_otsu_ch)
    //     .set { feed_me_into_lbp_otsu }


    // downscale provided masks and otsu masks and create a combined channel with fastlbp inputs

    downscale_me_otsu.concat(downscale_me_with_mask)
        .set { downscale_me_otsu_and_with_mask }

    img_to_get_masks.lbp_params_otsu_ch.concat(convert_me_to_binmask.lbp_params_ch)
        .set { lbp_params_str_otsu_and_with_mask }

    downscale_mask(downscale_me_otsu_and_with_mask)
    downscale_mask.out
        .merge(lbp_params_str_otsu_and_with_mask)
        .set { feed_me_into_lbp_otsu_and_with_mask } 


    // fastlbp
    feed_me_into_lbp_no_mask.concat(feed_me_into_lbp_otsu_and_with_mask)
        .set { all_feed_me_into_lbp }
    fastlbp(all_feed_me_into_lbp)

}