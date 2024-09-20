#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

pipeline_version = "0.0.2"

// Default params
params.background_color = ""

debug_flag = true

if ( !nextflow.version.matches(">=24.04") ) {
    error "The workflow requires Nextflow version 24.04 or greater, your current version is ${nextflow.version}"
}

def info_log(msg) {
    // let's make it yellow on a black bg
    log.info("\u001B[93;40;1m" + msg + "\u001B[0m")
}


def extract_parameter_val_from_params_list(param_list, parameter_to_extract = 'patchsize') {
    // params list of the following structure: [[param_1, value_1], [param_2, value_2], ...]
    def res = param_list.find { it[0] == parameter_to_extract }
    res[1]
}

// TODO: refactor
// BUG: works incorrectly with string parameters
def createCombinations(method_params) {
    // example: 
    /*
    [method:umap, n_components:[10, 20, 30], n_neighbors:[20, 30], min_dist:[values:[0.0, 0.1, 0.3], bind_to:n_components], random_state:42]
    */

    // total number of parameters for the method
    def method_params_num = method_params.size()
    // example: 5

    println "method_params_num = ${method_params_num}"

    def collected_params_list = method_params.collect { key, values -> values instanceof Map ? \
    tuple(key, values.values, values.bind_to) : tuple(key, values, key) }
    // collected_params_list example: 
    /*
    [[method, umap, method], [n_components, [10, 20, 30], n_components], [n_neighbors, [20, 30], n_neighbors], [min_dist, [0.0, 0.1, 0.3], n_components], [random_state, 42, random_state]]
    */


    def collected_params = Channel.fromList(collected_params_list)

    def combinations_entities = collected_params
        .groupTuple(by:2)
        .map { params_list, values_list, bind_to_key ->
        [params_list, values_list].transpose() }

    
    // combinations_entities example:

    /* [1] groupTuple(by:2)
    [[n_components, min_dist], [[10, 20, 30], [0.0, 0.1, 0.3]], n_components]
    [[n_neighbors], [[20, 30]], n_neighbors]
    */

    /* [2] map {params_list, values_list, bind_to_key -> [params_list, values_list].transpose()}
    [[n_components, [10, 20, 30]], [min_dist, [0.0, 0.1, 0.3]]]
    [[n_neighbors, [20, 30]]]
    */

    def entities_grouped = combinations_entities
        .map { bound_parameter_values ->
        
        def transformParams = { data ->
            // parentheses needed for groovy key evaluation
            def output = data.collectEntries { iit -> [(iit[0]): iit[1]] }
            // example
            /*
            output = [n_components:[10, 20, 30], min_dist:[0.0, 0.1, 0.3]]
            */

            // values in the map above are guaranteed to be of the same length
            // due to the config parsing strategy 
            // (different params will appear in the same map only if bind_to parameter was used)

            // now combine parameters pairwise
            def res = []
            output.values()[0].eachWithIndex {_, idx ->
                // iterating over the indices of the values in the list (in the example: 3)
                def combined_res = []
                output.each { kk, vv ->
                    // for each [parameter, values_list] pair

                    // append to the combined result parameter name
                    combined_res << kk

                    // and parameter value from the list (if it's a list)
                    if (vv instanceof List)
                        combined_res << vv[idx]
                    else
                        combined_res << vv

                    // example: 
                    // combined_res = [n_components, 10, min_dist, 0.0]
                }

                // append to the total combinations combined_res
                res << combined_res
            }
            return res
        }
        // example:
        /*
        bound_parameter_values = [[n_components, [10, 20, 30]], [min_dist, [0.0, 0.1, 0.3]]]

        output = [n_components:[10, 20, 30], min_dist:[0.0, 0.1, 0.3]]

        res = [[n_components, 10, min_dist, 0.0], [n_components, 20, min_dist, 0.1], [n_components, 30, min_dist, 0.3]]
        */
        return transformParams(bound_parameter_values) 
        }

    // now we need to flatten obtained channel of parameter combinations
    // where some elements consist of multiple parameters bound to each other

    // so we define a flattening function
    // function
    flatten_in_my_way = { my_list ->
        // example of my_list:
        /*
        [[[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]], 
        [[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]], ...]
        */

        def flatList = []
        my_list.each { itemm ->
            // for each element in the list:
            // example of such itemm:
            /*
            [[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]]
            */
            itemm.each { bound_params_list -> 
            bound_params_list.each {alternating_param_name_value -> 
            flatList.add(alternating_param_name_value) } 
            }
        }
        return flatList
    }


    // now we need to find all combinations of parameters
    // properly taking care of bound parameters in entitied_grouped
    def all_combs_flat = entities_grouped // example: [[n_components, 10, min_dist, 0.0], [n_components, 20, min_dist, 0.1], [n_components, 30, min_dist, 0.3]]
        .map { kk -> [kk] } // enclose in additional list
        .collect() // collect in a single list 
        .combinations() // groovy method to find all combinations of elements in list
        .map { ll -> flatten_in_my_way(ll) } // flatten elements of type: [[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]]

    // example
    /*
    all_combs_flat = [method, umap, n_components, 10, min_dist, 0.0, n_neighbors, 20, random_state, 42, method, umap, n_components, 10, min_dist, 0.0, n_neighbors, 20, random_state, 42, method, umap, ...]
    */
    
    all_combs_flat.view()
    // now as we flattened the combinations list we need to group back parameter combinations for different runs
    def all_combs_final = all_combs_flat.flatMap{ kkk -> kkk } // flatten all parameter lists 
        .collate(method_params_num * 2) // group individual run combinations (if there are n parameters -> there are n values, so we group every 2*n elements)
        .map { run_parameters -> run_parameters.collate(2) } // group [parameter, value] pairs in a flattened list of all parameters for a single run

    // TODO: the method above produces duplicates, urgently need optimisation
    return all_combs_final.distinct() // TODO: optimise this
}


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

process fastlbp {
    debug debug_flag
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
    run_umap.py \
        --np_data_path ${lbp_result_flattened} \
        --params_str "${params_str}"
    """
}

process clustering {
    debug debug_flag
    tag "${params.args.clustering.method}"
    // publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(data), val(params_str)

    output:
    tuple val(step_outdir), path("clustering_labels.npy")

    script:
    """
    run_hdbscan.py \
        --np_data_path ${data} \
        --params_str "${params_str}"
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
    pairs_max_csv = 'pairs_max_jacc.csv'
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
    integer_annot_path = params.integer_annot_path ? "--integer_annot_path ${params.integer_annot_path}" : ""
    annot_legend_path = params.annot_legend_path ? "--annot_legend_path ${params.annot_legend_path}" : "" // TODO: use as input?
    """
    generate_report.py \
        --outdir ${outdir} \
        --img_path ${img} \
        ${annot_path} \
        ${integer_annot_path} \
        ${annot_legend_path} \
        --savefile ${html_report}
    """
}

// TODO: implement repeating parts as subworkflows
workflow SingleImage {

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

    // generate output dir names for each combination of parameters
    lbp_combinations
        .map { it -> 
        tuple("${it.toString()}", it) }
        .set { lbp_combinations_hash_outdir }

    umap_combinations
        .map { it ->
        tuple("${it.toString()}", it) }
        .set { umap_combinations_hash_outdir }

    hdbscan_combinations
        .map { it ->
        tuple("${it.toString()}", it) }
        .set { hdbscan_combinations_hash_outdir }

    if ( !params.mask ) {

        info_log("No mask mode")

        lbp_combinations_hash_outdir
            .map { lbp_params_str, lbp_params ->
            tuple(lbp_params_str, params.img_path, [], [], lbp_params) }
            .set { feed_me_into_lbp }
        
        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_data_folder
            .set { all_lbp_outputs }

        fastlbp.out.lbp_result_file_flattened
            .combine(umap_combinations_hash_outdir)
            .map { lbp_outdir, lbp_result_flattened_cur, umap_params_str, umap_params ->
            tuple("${lbp_outdir}/${umap_params_str}", lbp_result_flattened_cur, umap_params) }
            .set {feed_me_into_umap}

        fastlbp.out.lbp_result_file_img
            .set { lbp_img }

        dimred(feed_me_into_umap)

        dimred.out
            .set { all_umap_outputs }

        dimred.out
            .combine(hdbscan_combinations_hash_outdir)
            .map { umap_outdir, umap_embeddings, hdbscan_params_str, hdbscan_params ->
            tuple("${umap_outdir}/${hdbscan_params_str}", umap_embeddings, hdbscan_params) }
            .set { feed_me_into_hdbscan }

        clustering(feed_me_into_hdbscan)

        clustering.out
            .set { all_hdbscan_outputs }

        lbp_img
            .join(feed_me_into_lbp)
            .set { lbp_and_masks }

        clustering.out
            .combine(lbp_and_masks)
            .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
            .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
            tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
            .set { convert_me_back_to_image }

        labels_to_patch_img(convert_me_back_to_image)

        labels_to_patch_img.out
            .set { all_deconvolves_images_outputs }

        all_lbp_outputs
            .combine(all_umap_outputs)
            .combine(all_hdbscan_outputs)
            .combine(all_deconvolves_images_outputs)
            .filter { it[2].contains(it[0]) && it[4].contains(it[2]) && it[6].contains(it[4]) }
            .multiMap { lbp_outdirr, lbp_dataa, umap_outdirr, umap_dataa, hdbscan_outdirr, hdbscan_dataa, 
            deconvolve_outdirr, deconvolve_dataa ->
            // replace commas with semicolons to distinguish between inner parameters list and
            // outer combination-to-hash list 
            data_to_save: tuple(deconvolve_outdirr.replaceAll(",", ";").md5(), tuple(lbp_dataa, umap_dataa, hdbscan_dataa, deconvolve_dataa))
            dir_and_hash: tuple(deconvolve_outdirr.replaceAll(",", ";"), deconvolve_outdirr.replaceAll(",", ";").md5()) }
            .set { final_data_and_folders }

        final_data_and_folders.data_to_save
            .transpose()
            .map { params_hash, step_dataa ->
            tuple("${params.outdir}/${params_hash}", step_dataa) }
            .set { dirs_and_data_to_save }

        copy_files_to_target_dir(dirs_and_data_to_save)

        final_data_and_folders.dir_and_hash
            .toList()
            .set { combinations_metadata }
        
        combinations_metadata_to_tsv(combinations_metadata)

        combinations_metadata_to_tsv.out.collect()
            .map { tuple(params.img_path, [], params.outdir) }
            .set { data_for_report }

        generate_report(data_for_report)

    } else if ( params.mask == "auto" ) {

        info_log("Otsu mask mode")

        get_tissue_mask(tuple(params.img_path, params.background_color))

        get_tissue_mask.out
            .set { pixel_mask_ch }

        get_tissue_mask.out
            .map { imgg, pixelmask ->
            pixelmask }
            .toList()
            .set { mask_for_report }

        get_tissue_mask.out
            .combine(lbp_combinations_hash_outdir)
            .map { imgg, pixelmask, lbp_outdir, lbp_params ->
            tuple(lbp_outdir, imgg, pixelmask, extract_parameter_val_from_params_list(lbp_params, parameter_to_extract = 'patchsize')) }
            .set {downscale_me}

        downscale_mask(downscale_me)

        downscale_mask.out
            .set { all_downscale_outputs }

        downscale_mask.out
            .combine(lbp_combinations_hash_outdir, by:0)
            .combine(pixel_mask_ch)
            .map { downscale_outdir, patchmask_path, lbp_params, imgg, pixelmask_path ->
            tuple(downscale_outdir, imgg, pixelmask_path, patchmask_path, lbp_params) }
            .set { feed_me_into_lbp }

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_data_folder
            .set { all_lbp_outputs }

        fastlbp.out.lbp_result_file_flattened
            .combine(umap_combinations_hash_outdir)
            .map { lbp_outdir, lbp_result_flattened_cur, umap_params_str, umap_params ->
            tuple("${lbp_outdir}/${umap_params_str}", lbp_result_flattened_cur, umap_params) }
            .set {feed_me_into_umap}

        fastlbp.out.lbp_result_file_img
            .set { lbp_img }

        dimred(feed_me_into_umap)

        dimred.out
            .set { all_umap_outputs }

        dimred.out
            .combine(hdbscan_combinations_hash_outdir)
            .map { umap_outdir, umap_embeddings, hdbscan_params_str, hdbscan_params ->
            tuple("${umap_outdir}/${hdbscan_params_str}", umap_embeddings, hdbscan_params) }
            .set { feed_me_into_hdbscan }

        clustering(feed_me_into_hdbscan)

        clustering.out
            .set { all_hdbscan_outputs }

        lbp_img
            .join(feed_me_into_lbp)
            .set { lbp_and_masks }

        clustering.out
            .combine(lbp_and_masks)
            .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
            .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
            tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
            .set { convert_me_back_to_image }

        labels_to_patch_img(convert_me_back_to_image)

        labels_to_patch_img.out
            .set { all_deconvolves_images_outputs }

        all_downscale_outputs
            .combine(all_lbp_outputs)
            .combine(all_umap_outputs)
            .combine(all_hdbscan_outputs)
            .combine(all_deconvolves_images_outputs)
            .filter { it[2].contains(it[0]) && it[4].contains(it[2]) && 
            it[6].contains(it[4]) && it[8].contains(it[6]) } // as there is additional step as compared to no mask 
            .multiMap { downscale_outdirr, downscale_dataa, lbp_outdirr, lbp_dataa, umap_outdirr, umap_dataa, hdbscan_outdirr, hdbscan_dataa, 
            deconvolve_outdirr, deconvolve_dataa ->
            // replace commas with semicolons to distinguish between inner parameters list and
            // outer combination-to-hash list 
            data_to_save: tuple(deconvolve_outdirr.replaceAll(",", ";").md5(), tuple(downscale_dataa, lbp_dataa, umap_dataa, hdbscan_dataa, deconvolve_dataa))
            dir_and_hash: tuple(deconvolve_outdirr.replaceAll(",", ";"), deconvolve_outdirr.replaceAll(",", ";").md5()) }
            .set { final_data_and_folders }

        final_data_and_folders.data_to_save
            .transpose()
            .map { params_hash, step_dataa ->
            tuple("${params.outdir}/${params_hash}", step_dataa) }
            .set { dirs_and_data_to_save }

        copy_files_to_target_dir(dirs_and_data_to_save)

        final_data_and_folders.dir_and_hash
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

    } else {

        info_log("Mask/annotations mode")

        convert_to_binmask_ch = Channel.of([params.img_path, params.mask, 
        params.background_color])

        convert_annotations_to_binmask(convert_to_binmask_ch)

        convert_annotations_to_binmask.out
            .set { converted_mask_ch }

        convert_annotations_to_binmask.out
            .combine(lbp_combinations_hash_outdir)
            .map { imgg, maskk, lbp_outdirr, lbp_paramss ->
            tuple(lbp_outdirr, imgg, maskk, extract_parameter_val_from_params_list(lbp_paramss, parameter_to_extract = 'patchsize')) }
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

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_data_folder
            .set { all_lbp_outputs }

        fastlbp.out.lbp_result_file_flattened
            .combine(umap_combinations_hash_outdir)
            .map { lbp_outdir, lbp_result_flattened_cur, umap_params_str, umap_params ->
            tuple("${lbp_outdir}/${umap_params_str}", lbp_result_flattened_cur, umap_params) }
            .set {feed_me_into_umap}

        fastlbp.out.lbp_result_file_img
            .set { lbp_img }

        dimred(feed_me_into_umap)

        dimred.out
            .set { all_umap_outputs }

        dimred.out
            .combine(hdbscan_combinations_hash_outdir)
            .map { umap_outdir, umap_embeddings, hdbscan_params_str, hdbscan_params ->
            tuple("${umap_outdir}/${hdbscan_params_str}", umap_embeddings, hdbscan_params) }
            .unique() // FIXME: not optimal hotfix
            .set { feed_me_into_hdbscan }

        clustering(feed_me_into_hdbscan)

        clustering.out
            .set { all_hdbscan_outputs }

        lbp_img
            .join(feed_me_into_lbp)
            .set { lbp_and_masks }

        clustering.out
            .combine(lbp_and_masks)
            .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
            .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
            tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
            .set { convert_me_back_to_image }

        labels_to_patch_img(convert_me_back_to_image)

        labels_to_patch_img.out
            .set { all_deconvolves_images_outputs }
        
        // Calculate Jaccard scores for runs with respect to annotations
        // and do cluster matching in a greedy fashion
        all_deconvolves_images_outputs
            .map { patch_img_outdir, patch_img ->
            tuple(patch_img_outdir, params.integer_annot, patch_img) }
            .set { runs_to_calculate_eval_metric }

        calculate_unsupervised_clustering_score(runs_to_calculate_eval_metric)
        calculate_unsupervised_clustering_score.out
            .set { calculated_metrics_outputs }

        all_downscale_outputs
            .combine(all_lbp_outputs)
            .combine(all_umap_outputs)
            .combine(all_hdbscan_outputs)
            .combine(all_deconvolves_images_outputs)
            .combine(calculated_metrics_outputs)
            .filter { it[2].contains(it[0]) && it[4].contains(it[2]) && // TODO: don't embarass yourself with this
            it[6].contains(it[4]) && it[8].contains(it[6]) && it[10].contains(it[8])} // as there are additional steps as compared to no mask 
            // .view()
            .multiMap { downscale_outdirr, downscale_dataa, lbp_outdirr, lbp_dataa, umap_outdirr, umap_dataa, hdbscan_outdirr, hdbscan_dataa, 
            deconvolve_outdirr, deconvolve_dataa, metrics_outdirr, metrics_all_pairs_dataa, metrics_max_pairs_data ->
            // replace commas with semicolons to distinguish between inner parameters list and
            // outer combination-to-hash list 
            data_to_save: tuple(deconvolve_outdirr.replaceAll(",", ";").md5(), tuple(downscale_dataa, lbp_dataa, umap_dataa, hdbscan_dataa, deconvolve_dataa, metrics_all_pairs_dataa, metrics_max_pairs_data))
            dir_and_hash: tuple(deconvolve_outdirr.replaceAll(",", ";"), deconvolve_outdirr.replaceAll(",", ";").md5()) }
            .set { final_data_and_folders }

        final_data_and_folders.data_to_save
            .transpose()
            .map { params_hash, step_dataa ->
            tuple("${params.outdir}/${params_hash}", step_dataa) }
            .set { dirs_and_data_to_save }


        copy_files_to_target_dir(dirs_and_data_to_save)

        final_data_and_folders.dir_and_hash
            .toList()
            .set { combinations_metadata }
        
        combinations_metadata_to_tsv(combinations_metadata)

        combinations_metadata_to_tsv.out.collect()
            .map { tuple(params.img_path, params.mask, params.outdir) }
            .set { data_for_report }

        generate_report(data_for_report)

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
