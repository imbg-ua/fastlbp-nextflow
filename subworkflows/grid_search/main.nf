#!/usr/bin/env/ nextflow

include { infoLog; 
          checkNextflowVersion; 
          getValueFromParamList; 
          checkNestedParameterCombinations;
          getFinalOutdirFromParamsCombinations;
          getAllOutputFilesFromParamsCombinations;
          createCombinations } from '../../lib/nf/utils'

include { RunFastLBP;
          RunFastLBPAndPrepareForNextStep } from '../../modules/feature_extraction'
include { RunDimRed;
          RunDimRedAndPrepareForNextStep } from '../../modules/dimensionality_reduction'
include { RunClustering } from '../../modules/clustering'


nextflow.enable.dsl = 2

pipeline_version = "0.0.2"

// Default params
// TODO: move to config
// TODO: move before include statements
params.background_color = ""
params.pairs_max_csv = 'pairs_max_jacc.csv'
// params.pairs_max_csv = '' // TODO: BUG: Make dependent oт the availability of annotaitons
params.annot_legend_path = ""
params.plots_backend = 'matplotlib'

debug_flag = true

checkNextflowVersion()

process copy_files_to_target_dir {
    tag "wrapping up"
    debug debug_flag
    publishDir "${outdir}", mode: "copy"

    input:
    tuple val(outdir), path(data)

    output:
    tuple val(outdir), path(data)
    
    script:
    """
    """
}

process combinations_metadata_to_tsv {
    tag "wrapping up"
    debug debug_flag
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

process convert_annotations_to_binmask {
    tag "mask preprocessing"
    debug debug_flag
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


process get_tissue_mask {
    tag "preprocessing"
    debug debug_flag
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
    debug debug_flag
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

process labels_to_patch_img {
    debug debug_flag
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

process calculate_unsupervised_clustering_score {
    debug debug_flag
    tag "analysis"
    // publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(annot), path(patch_img)

    output:
    tuple val(step_outdir), path(all_pairs_csv), path(pairs_max_csv)

    script:
    // TODO: fix inconsistent namings
    all_pairs_csv = 'all_pairs_jacc.csv'
    // pairs_max_csv = 'pairs_max_jacc.csv'
    pairs_max_csv = params.pairs_max_csv
    """
    calculate_unsupervised_metrics.py \
        --patch_img_path ${patch_img} \
        --annotation_path ${annot} \
        --savefile_all_pairs ${all_pairs_csv} \
        --savefile_max_pairs ${pairs_max_csv}
    """
}

process generate_report {
    debug debug_flag
    tag "wrapping up"
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple path(img), path(annot), path(outdir)

    output:
    path(html_report)

    script:
    html_report = "report.html"
    annot_path = annot ? "--annot_path ${annot}" : ""
    integer_annot_path = params.integer_annot ? "--integer_annot_path ${params.integer_annot}" : ""
    annot_legend_path = params.annot_legend_path ? "--annot_legend_path ${params.annot_legend_path}" : "" // TODO: use as input?
    class_mapping_csv = params.pairs_max_csv ? "--class_mapping_csv ${params.pairs_max_csv}" : "" // TODO: fix inconsistent names
    plots_backend = params.plots_backend ? "--plots_backend ${params.plots_backend}" : "" // TODO: if I leave empty str as default will it work correctly?
    """
    generate_report.py \
        --outdir ${outdir} \
        --img_path ${img} \
        ${annot_path} \
        ${integer_annot_path} \
        ${annot_legend_path} \
        ${class_mapping_csv} \
        ${plots_backend} \
        --savefile ${html_report}
    """
}

// TODO: implement repeating parts as subworkflows

workflow GetParameterCombinations {
    main:
    
    // TODO: don't resort to config params inside the workflow

    println "getting all paramter combinations" // DEBUG
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


workflow LabelsToPatchImage {
    take:
    labels_to_patch_img_inputs_channel

    main:

    labels_to_patch_img(labels_to_patch_img_inputs_channel)

    emit:
    all_deconvolved_images_outputs = labels_to_patch_img.out
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

    println "finished getting all param combs" // DEBUG

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
        """)
    }
}
