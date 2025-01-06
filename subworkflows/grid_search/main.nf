#!/usr/bin/env/ nextflow


// TODO: add script-specific default params

pipeline_version = "0.0.2"
nextflow.enable.dsl = 2

params.debug_flag = true

include { infoLog; 
          checkNextflowVersion; 
          getValueFromParamList; 
          checkNestedParameterCombinations;
          getFinalOutdirFromParamsCombinations;
          getAllOutputFilesFromParamsCombinations;
          createCombinations } from '../../lib/nf/utils'


checkNextflowVersion()


include { RunFastLBP;
          RunFastLBPAndPrepareForNextStep } from '../../modules/feature_extraction'

include { RunDimRed;
          RunDimRedAndPrepareForNextStep } from '../../modules/dimensionality_reduction'

include { RunClustering } from '../../modules/clustering'

include { GetTissueMask;
          GetTissuePatchMask;
          GetTissuePatchMaskFromAnnnotations } from '../../modules/mask_utils'

include { calculate_unsupervised_clustering_score } from '../../modules/metrics'
include { generate_report } from '../../modules/reporting'

include { LabelsToPatchImage;
          copy_files_to_target_dir;
          combinations_metadata_to_tsv;
          GetParameterCombinations;
          ParameterCombinationsToMap } from '../../modules/utils'

workflow NoMaskWorkflow {

    take:
    lbp_combinations_hash_outdir
    umap_combinations_hash_outdir
    hdbscan_combinations_hash_outdir

    main:
    infoLog("No mask mode")

    lbp_combinations_hash_outdir
        .map { lbp_params_str, lbp_params ->
        tuple(lbp_params_str, params.img_path, [], [], lbp_params) }
        .set { feed_me_into_lbp }
    
    RunFastLBPAndPrepareForNextStep(feed_me_into_lbp, umap_combinations_hash_outdir)

    all_lbp_outputs = RunFastLBPAndPrepareForNextStep.out.all_lbp_outputs
    feed_me_into_umap = RunFastLBPAndPrepareForNextStep.out.next_step_inputs_channel
    lbp_img = RunFastLBPAndPrepareForNextStep.out.lbp_result_img

    RunDimRedAndPrepareForNextStep(feed_me_into_umap, hdbscan_combinations_hash_outdir)
    all_umap_outputs = RunDimRedAndPrepareForNextStep.out.all_dimred_outputs
    feed_me_into_hdbscan = RunDimRedAndPrepareForNextStep.out.next_step_inputs_channel


    RunClustering(feed_me_into_hdbscan)
    all_hdbscan_outputs = RunClustering.out.all_clustering_outputs

    // connect LBP outputs to corresponding masks used

    // TODO: check if this is really necessary as the outputs
    // should retain the order of the respective inputs
    lbp_img
        .join(feed_me_into_lbp)
        .set { lbp_and_masks }

    // connect clustering result to the LBP result and used mask
    // TODO: switch from combine to merge? is the order the same? 
    // prolly not, clustering combinations were calculated separately from lbp combinations 
    // (but still on top of the previous steps). Anyway, TODO remains valid.
    all_hdbscan_outputs
        .combine(lbp_and_masks)
        .filter { it[0].toString().contains(it[2].toString()) } // clustering_outdir, clustering_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
        .map { clustering_outdir, clustering_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
        tuple(clustering_outdir, clustering_labels_path, patchmask_path, lbp_img_path) }
        .set { convert_me_back_to_image }

    LabelsToPatchImage(convert_me_back_to_image)
    all_deconvolved_images_outputs = LabelsToPatchImage.out.all_deconvolved_images_outputs

    AllCombinationsToHashes([all_lbp_outputs, all_umap_outputs, all_hdbscan_outputs, all_deconvolved_images_outputs])
    data_to_save = AllCombinationsToHashes.out.data_to_save
    dir_and_hash = AllCombinationsToHashes.out.dir_and_hash

    data_to_save
        .transpose()
        .map { params_hash, step_dataa ->
        tuple("${params.outdir}/${params_hash}", step_dataa) }
        .set { dirs_and_data_to_save }

    copy_files_to_target_dir(dirs_and_data_to_save)

    dir_and_hash
        .toList()
        .set { combinations_metadata }
    
    combinations_metadata_to_tsv(combinations_metadata)

    // wait until all processes finish before genering the report
    combinations_metadata_to_tsv.out.collect()
        .map { tuple(params.img_path, [], params.outdir) }
        .set { data_for_report }

    generate_report(data_for_report)

}


workflow OtsuMaskWorkflow {
    take:
    lbp_combinations_hash_outdir
    umap_combinations_hash_outdir
    hdbscan_combinations_hash_outdir

    main:
    infoLog("Otsu mask mode")

    def img_and_background_color_channel = Channel.of(tuple(params.img_path, params.background_color))

    GetTissuePatchMask(img_and_background_color_channel, lbp_combinations_hash_outdir)
    feed_me_into_lbp = GetTissuePatchMask.out.feed_me_into_lbp
    all_downscale_outputs = GetTissuePatchMask.out.all_downscale_outputs
    mask_for_report = GetTissuePatchMask.out.mask_for_report

    RunFastLBPAndPrepareForNextStep(feed_me_into_lbp, umap_combinations_hash_outdir)
    all_lbp_outputs = RunFastLBPAndPrepareForNextStep.out.all_lbp_outputs
    feed_me_into_umap = RunFastLBPAndPrepareForNextStep.out.next_step_inputs_channel
    lbp_img = RunFastLBPAndPrepareForNextStep.out.lbp_result_img


    RunDimRedAndPrepareForNextStep(feed_me_into_umap, hdbscan_combinations_hash_outdir)
    all_umap_outputs = RunDimRedAndPrepareForNextStep.out.all_dimred_outputs
    feed_me_into_hdbscan = RunDimRedAndPrepareForNextStep.out.next_step_inputs_channel

    RunClustering(feed_me_into_hdbscan)
    all_hdbscan_outputs = RunClustering.out.all_clustering_outputs

    lbp_img
        .join(feed_me_into_lbp)
        .set { lbp_and_masks }

    all_hdbscan_outputs
        .combine(lbp_and_masks)
        .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
        .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
        tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
        .set { convert_me_back_to_image }

    LabelsToPatchImage(convert_me_back_to_image)
    all_deconvolved_images_outputs = LabelsToPatchImage.out.all_deconvolved_images_outputs

    AllCombinationsToHashes([all_downscale_outputs, all_lbp_outputs, 
                             all_umap_outputs, all_hdbscan_outputs, 
                             all_deconvolved_images_outputs])
    data_to_save = AllCombinationsToHashes.out.data_to_save
    dir_and_hash = AllCombinationsToHashes.out.dir_and_hash

    data_to_save
        .transpose()
        .map { params_hash, step_dataa ->
        tuple("${params.outdir}/${params_hash}", step_dataa) }
        .set { dirs_and_data_to_save }

    copy_files_to_target_dir(dirs_and_data_to_save)

    dir_and_hash
        .toList()
        .set { combinations_metadata }

    combinations_metadata_to_tsv(combinations_metadata)

    combinations_metadata_to_tsv.out.collect()
        .map { tuple(params.img_path, params.outdir) }
        .set { data_for_report }

    data_for_report
        .combine(mask_for_report)
        .map { imgg, outdirr, maskk ->
        tuple(imgg, maskk, outdirr) }
        .set { data_for_report_auto_mask }

    generate_report(data_for_report_auto_mask)
}

workflow MaskWorkflow {
    take:
    lbp_combinations_hash_outdir
    umap_combinations_hash_outdir
    hdbscan_combinations_hash_outdir

    main:
    infoLog("Mask/annotations mode")

    convert_to_binmask_ch = Channel.of([params.img_path, params.mask, 
    params.background_color])

    GetTissuePatchMaskFromAnnnotations(convert_to_binmask_ch, lbp_combinations_hash_outdir)
    feed_me_into_lbp = GetTissuePatchMaskFromAnnnotations.out.feed_me_into_lbp
    all_downscale_outputs = GetTissuePatchMaskFromAnnnotations.out.all_downscale_outputs

    RunFastLBPAndPrepareForNextStep(feed_me_into_lbp, umap_combinations_hash_outdir)
    all_lbp_outputs = RunFastLBPAndPrepareForNextStep.out.all_lbp_outputs
    feed_me_into_umap = RunFastLBPAndPrepareForNextStep.out.next_step_inputs_channel
    lbp_img = RunFastLBPAndPrepareForNextStep.out.lbp_result_img

    RunDimRedAndPrepareForNextStep(feed_me_into_umap, hdbscan_combinations_hash_outdir)
    all_umap_outputs = RunDimRedAndPrepareForNextStep.out.all_dimred_outputs
    feed_me_into_hdbscan = RunDimRedAndPrepareForNextStep.out.next_step_inputs_channel

    RunClustering(feed_me_into_hdbscan)
    all_hdbscan_outputs = RunClustering.out.all_clustering_outputs

    lbp_img
        .join(feed_me_into_lbp)
        .set { lbp_and_masks }

    all_hdbscan_outputs
        .combine(lbp_and_masks)
        .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
        .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
        tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
        .set { convert_me_back_to_image }

    LabelsToPatchImage(convert_me_back_to_image)
    all_deconvolved_images_outputs = LabelsToPatchImage.out.all_deconvolved_images_outputs


    // if integer annotations are provided
    if ( params.integer_annot ) {

        // Calculate Jaccard scores for runs with respect to annotations
        // and do cluster matching in a greedy fashion
        all_deconvolved_images_outputs
            .map { patch_img_outdir, patch_img ->
            tuple(patch_img_outdir, params.integer_annot, patch_img) }
            .set { runs_to_calculate_eval_metric }

        calculate_unsupervised_clustering_score(runs_to_calculate_eval_metric)
        calculate_unsupervised_clustering_score.out
            .set { calculated_metrics_outputs }

        // TODO: generalise AllCombinationsToHashes to use here
        all_downscale_outputs
                .combine(all_lbp_outputs)
                .combine(all_umap_outputs)
                .combine(all_hdbscan_outputs)
                .combine(all_deconvolved_images_outputs)
                .combine(calculated_metrics_outputs)
                .filter { it[2].contains(it[0]) && it[4].contains(it[2]) &&
                it[6].contains(it[4]) && it[8].contains(it[6]) && it[10].contains(it[8])} // TOOD: as there are additional steps as compared to no mask / Otsu
                .multiMap { downscale_outdirr, downscale_dataa, lbp_outdirr, lbp_dataa, umap_outdirr, umap_dataa, hdbscan_outdirr, hdbscan_dataa, 
                deconvolve_outdirr, deconvolve_dataa, metrics_outdirr, metrics_all_pairs_dataa, metrics_max_pairs_data -> // the difference from the AllCombinationsToHashes is in the metric calculcation outputs
                data_to_save: tuple(deconvolve_outdirr.replaceAll(",", ";").md5(), tuple(downscale_dataa, lbp_dataa, umap_dataa, hdbscan_dataa, deconvolve_dataa, metrics_all_pairs_dataa, metrics_max_pairs_data))
                dir_and_hash: tuple(deconvolve_outdirr.replaceAll(",", ";"), deconvolve_outdirr.replaceAll(",", ";").md5()) }
                .set { final_data_and_folders }

        // TODO: I don't like it
        final_data_and_folders.data_to_save
            .set { data_to_save }
        final_data_and_folders.dir_and_hash
            .set { dir_and_hash }

        
    } else {
        // TODO: create a process that converts any type of mask into integer annotations
        // for now it behaves similarly to the Otsu mask mode if no integer annotations were provided

        AllCombinationsToHashes([all_downscale_outputs, all_lbp_outputs, 
                                all_umap_outputs, all_hdbscan_outputs, 
                                all_deconvolved_images_outputs])
        data_to_save = AllCombinationsToHashes.out.data_to_save
        dir_and_hash = AllCombinationsToHashes.out.dir_and_hash
    }

    data_to_save
        .transpose()
        .map { params_hash, step_dataa ->
        tuple("${params.outdir}/${params_hash}", step_dataa) }
        .set { dirs_and_data_to_save }

    copy_files_to_target_dir(dirs_and_data_to_save)

    dir_and_hash
        .toList()
        .set { combinations_metadata }

    combinations_metadata_to_tsv(combinations_metadata)

    combinations_metadata_to_tsv.out.collect()
        .map { tuple(params.img_path, params.mask, params.outdir) }
        .set { data_for_report }

    generate_report(data_for_report)
}


workflow SingleImage {

    GetParameterCombinations()

    lbp_combinations = GetParameterCombinations.out.lbp_combinations
    umap_combinations = GetParameterCombinations.out.umap_combinations
    hdbscan_combinations = GetParameterCombinations.out.hdbscan_combinations

    // generate output dir names for each combination of parameters
    
    lbp_combinations_hash_outdir = ParameterCombinationsToMap(lbp_combinations)
    umap_combinations_hash_outdir = ParameterCombinationsToMap(umap_combinations)
    hdbscan_combinations_hash_outdir = ParameterCombinationsToMap(hdbscan_combinations)

    if ( !params.mask ) {

        NoMaskWorkflow(lbp_combinations_hash_outdir, 
                       umap_combinations_hash_outdir, 
                       hdbscan_combinations_hash_outdir)

    } else if ( params.mask == "auto" ) {

        OtsuMaskWorkflow(lbp_combinations_hash_outdir, 
                       umap_combinations_hash_outdir, 
                       hdbscan_combinations_hash_outdir)

    } else {

        MaskWorkflow(lbp_combinations_hash_outdir, 
                     umap_combinations_hash_outdir, 
                     hdbscan_combinations_hash_outdir)

    }
}

// workflow entry point
workflow Pipeline {
    if ( params.img_path && !params.imgs_dir )    
        SingleImage()
    else {
        log.error("""
        Parameter search is available only in SingleImage mode.
        Make sure you provided correct path to the image.
        Processing multiple images with multiple sets of parameters
        will be soon available through TSV file-based interface.
        """)
    }
}
