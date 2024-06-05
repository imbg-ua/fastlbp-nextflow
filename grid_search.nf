#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

debug_flag = true
params.annot_suffix = "annotation"

def info_log(msg) {
    // let's make it yellow on a black bg
    log.info("\u001B[93;40;1m" + msg + "\u001B[0m")
}

def extract_patchsize_from_lbp_params_list(lbp_params_list) {
    // params list of the following structure: [[param_1, value_1], [param_2, value_2], ...]
    def res = lbp_params_list.find { it[0] == "patchsize" }
    res[1]
}

def createCombinations(method_params) {
    def method_params_num = method_params.size()
    def channels = method_params.collect { key, values -> values instanceof List ? \
    Channel.of([key, values]).transpose() : Channel.of([key, values])} // Channel.from()

    // println "$method_params DEBUG"
    // def channels_test = method_params.findAll { it.value.containsKey("values") } // Channel.from()
    // def no_channels_test = method_params.findAll { !it.value.containsKey("values") } // Channel.from()
    def channels_test = method_params.findAll { it.value instanceof Map } // Channel.from()
    def no_channels_test = method_params.findAll { !(it.value instanceof Map) } // Channel.from()

    // def testest_supertest = method_params.collect { key, values -> values instanceof Map ? \
    // tuple("parameter": key, "values": values.values, "bind_to": values.bind_to) : tuple("parameter": key, "values": values, "bind_to": key) }


    def testest_supertest = method_params.collect { key, values -> values instanceof Map ? \
    tuple(key, values.values, values.bind_to) : tuple(key, values, key) }

    def testest_superttest_ch = Channel.fromList(testest_supertest)

    // testest_superttest_ch
    //     .map { param, vals, bind_to ->
    //     bind_to }
    //     .unique()
    //     .set { unique_set_of_params }

    // def unique_params = unique_set_of_params.toList()

    testest_superttest_ch
        .groupTuple(by:2)
        .map { params_l, values_l, join_key ->
        [params_l, values_l].transpose() }
        .set { combinations_entities }

    combinations_entities
        .map { it ->
        
        def transformParams = { data ->
            def output = data.collectEntries { iit -> [(iit[0]): iit[1]] }
            def res = []
            output.values()[0].eachWithIndex {_, idx ->
                def combined_res = []
                output.each { kk, vv ->
                    combined_res << kk
                    if (vv instanceof List)
                        combined_res << vv[idx]
                    else
                        combined_res << vv
                }
                res << combined_res
            }

            return res
        }
        return transformParams(it) 
        }
        .set { entities_grouped }

    flatten_in_my_way = { my_list ->
        def flatList = []
        my_list.each { itemm ->
            itemm.each { otem -> otem.each {totem -> flatList.add(totem) } }
        }
        return flatList
    }
    // entities_grouped.view()
    entities_grouped
        .map { kk -> [kk] }
        .collect()
        .combinations()
        .map { ll -> flatten_in_my_way(ll) }
        .set { all_combs_flat }

    all_combs_flat.flatMap{ kkk -> kkk }.collate(method_params_num * 2)
        .map { aboba -> aboba.collate(2) }
        .set { all_combs_final }

    // all_combs_final.view()

    return all_combs_final



    // entities_grouped.flatMap { lists -> lists }.view()
    // entities_grouped.collect().combinations().view()
    // entities_grouped.flatten().view()
    // println(mmmeega_list.toString())
    // def entities_grouped_ch_ch = entities_grouped_ch.collect { ww -> Channel.of(ww) }


    // entities_grouped_ch_ch[0].view()

    
    // def first_entity = entities_grouped_ch.head()
    // def other_entities = entities_grouped_ch.tail()

    // other_entities.view()

    // def entities_grouped_ch_ch = entities_grouped_ch.collect { ww -> Channel.of([ww]) }

    // entities_grouped_ch_ch.view()

    // other_entities
    //     .map { its -> first_entity = first_entity.combine(its) }
    
    // first_entity.view()

        // .view()
    // def combinations_entities_list = combinations_entities.toList()

    // def new_nomap = no_channels_test

    // def new_map = [:]

    // channels_test.each { key, value -> 
    //     new_map["${value.bind_to}_${key}"] = tuple(["${value.bind_to}", method_params["${value.bind_to}"], "${key}", value.values])
    //     new_nomap = new_nomap.findAll { it.key != value.bind_to }
    // }

    // def new_map_combined_entries = new_map.collect { key, values -> tuple(values) }

    // def combined_entries_ch = Channel.fromList(new_map_combined_entries)

    // combined_entries_ch.view()

    // println "${new_map_combined_entries} DEBUG ALL PARAMS"
    // println "${new_nomap} DEBUG ALL PARAMS NOMAP"
    // println "$channels_test DEBUG MAPS"
    // println "$no_channels_test DEBUG NO MAPS"

    // channels_test[0].view()
    //     // .map { k, v ->
    //     // [k: v] }
    //     // .view()


    def combinedChannel = channels[0]
    combinedChannel = combinedChannel.map { [it] }

    // https://github.com/nextflow-io/nextflow/issues/1403
    channels.tail().each {
        combinedChannel = combinedChannel.combine(it.map { itit -> [itit] })
    }
    // combinedChannel.view()
    combinedChannel
}

def params_args_list = params.args.collect { k, v -> [k, v] }

def lbp_step = params.args.lbp
def lbp_step_param_num = params.args.lbp.size()
lbp_params_names = params.args.lbp.collect { k, v -> k }

def dimred_step = params.args.umap
def dimred_step_param_num = params.args.umap.size()
dimred_params_names = params.args.umap.collect { k, v -> k }

def clust_step = params.args.hdbscan
def clust_step_param_num = params.args.hdbscan.size()
clust_params_names = params.args.hdbscan.collect { k, v -> k }

def lbp_combinations = createCombinations(lbp_step)
def umap_combinations = createCombinations(dimred_step)
def hdbscan_combinations = createCombinations(clust_step)


process get_tissue_mask {
    tag "preprocessing"
    debug debug_flag
    publishDir "${params.outdir}", mode: "copy"

    input:
    path(img)

    output:
    path('pixelmask.npy')

    script:
    """
    get_mask.py get_mask \
        --img_path ${img} \
    """
}

// process convert_to_binmask {
//     tag "${run_id}"
//     debug debug_flag
//     publishDir "${params.outdir}/${run_id}", mode: "copy"

//     input:
//     tuple val(run_id), path(pixelmask)

//     output:
//     tuple path(pixelmask), path('pixelmask.npy')

//     script:
//     """
//     get_mask.py convert_to_binmask \
//         --annot_path ${pixelmask} \
//     """
// }

process downscale_mask {
    debug debug_flag
    publishDir "${step_outdir}", mode: "copy"

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
    publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(img), path(mask), path(patchmask), val(params_str)

    output:
    path("data")
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

process umap {
    debug debug_flag
    publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(lbp_result_flattened), val(params_str)

    output:
    tuple val(step_outdir), path("umap_embeddings.npy")

    script:
    """
    run_umap.py \
        --np_data_path ${lbp_result_flattened} \
        --params_str "${params_str}"
    """
}

process hdbscan {
    debug debug_flag
    publishDir "${step_outdir}", mode: "copy"

    input:
    tuple val(step_outdir), path(data), val(params_str)

    output:
    tuple val(step_outdir), path("hdbscan_labels.npy")

    script:
    """
    run_hdbscan.py \
        --np_data_path ${data} \
        --params_str "${params_str}"
    """
}

process labels_to_patch_img {
    debug debug_flag
    publishDir "${step_outdir}", mode: "copy"

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

workflow SingleImage {
    // [string, hash] is the correct order

    // lbp_combinations.view()
    // umap_combinations.view()
    // hdbscan_combinations.view()

    lbp_combinations
        .map { it -> 
        tuple("${params.outdir}/${it.flatten().collect().join('_').toString()}", it) }
        .set { lbp_combinations_hash_outdir }

    umap_combinations
        .map { it ->
        tuple("${it.flatten().collect().join('_').toString()}", it) }
        .set { umap_combinations_hash_outdir }

    hdbscan_combinations
        .map { it ->
        tuple("${it.flatten().collect().join('_').toString()}", it) }
        .set { hdbscan_combinations_hash_outdir }

    if ( !params.mask ) {

        info_log("No mask mode")

        lbp_combinations_hash_outdir
            .map { lbp_params_str, lbp_params ->
            tuple(lbp_params_str, params.img_path, [], [], lbp_params) }
            .set { feed_me_into_lbp }
        
        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_result_file_flattened
            .combine(umap_combinations_hash_outdir)
            .map { lbp_outdir, lbp_result_flattened_cur, umap_params_str, umap_params ->
            tuple("${lbp_outdir}/${umap_params_str}", lbp_result_flattened_cur, umap_params) }
            .set {feed_me_into_umap}

        fastlbp.out.lbp_result_file_img
            .set { lbp_img }

        umap(feed_me_into_umap)

        umap.out
            .combine(hdbscan_combinations_hash_outdir)
            .map { umap_outdir, umap_embeddings, hdbscan_params_str, hdbscan_params ->
            tuple("${umap_outdir}/${hdbscan_params_str}", umap_embeddings, hdbscan_params) }
            .set { feed_me_into_hdbscan }

        hdbscan(feed_me_into_hdbscan)

        lbp_img
            .join(feed_me_into_lbp)
            .set { lbp_and_masks }

        hdbscan.out
            .combine(lbp_and_masks)
            .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
            .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
            tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
            .set { convert_me_back_to_image }

        labels_to_patch_img(convert_me_back_to_image)

    } else if ( params.mask == "auto" ) {

        info_log("Otsu mask mode")

        get_tissue_mask(params.img_path)

        get_tissue_mask.out
            .set { pixel_mask_ch }

        get_tissue_mask.out
            .combine(lbp_combinations_hash_outdir)
            .map { pixelmask, lbp_outdir, lbp_params ->
            tuple(lbp_outdir, params.img_path, pixelmask, extract_patchsize_from_lbp_params_list(lbp_params)) }
            .set {downscale_me}

        downscale_mask(downscale_me)
        downscale_mask.out
            .combine(lbp_combinations_hash_outdir, by:0)
            .combine(pixel_mask_ch)
            .map { downscale_outdir, patchmask_path, lbp_params, pixelmask_path ->
            tuple(downscale_outdir, params.img_path, pixelmask_path, patchmask_path, lbp_params) }
            .set { feed_me_into_lbp }

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_result_file_flattened
            .combine(umap_combinations_hash_outdir)
            .map { lbp_outdir, lbp_result_flattened_cur, umap_params_str, umap_params ->
            tuple("${lbp_outdir}/${umap_params_str}", lbp_result_flattened_cur, umap_params) }
            .set {feed_me_into_umap}

        fastlbp.out.lbp_result_file_img
            .set { lbp_img }

        umap(feed_me_into_umap)

        umap.out
            .combine(hdbscan_combinations_hash_outdir)
            .map { umap_outdir, umap_embeddings, hdbscan_params_str, hdbscan_params ->
            tuple("${umap_outdir}/${hdbscan_params_str}", umap_embeddings, hdbscan_params) }
            .set { feed_me_into_hdbscan }

        hdbscan(feed_me_into_hdbscan)

        lbp_img
            .join(feed_me_into_lbp)
            .set { lbp_and_masks }

        hdbscan.out
            .combine(lbp_and_masks)
            .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
            .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
            tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
            .set { convert_me_back_to_image }

        labels_to_patch_img(convert_me_back_to_image)
    } else {

        info_log("Provided mask mode")

        lbp_combinations_hash_outdir
            .map { lbp_outdir, lbp_params ->
            tuple(lbp_outdir, params.img_path, params.mask, extract_patchsize_from_lbp_params_list(lbp_params)) }
            .set { downscale_me }

        downscale_mask(downscale_me)

        downscale_mask.out
            .combine(lbp_combinations_hash_outdir, by:0)
            .map { downscale_outdir, patchmask_path, lbp_params ->
            tuple(downscale_outdir, params.img_path, params.mask, patchmask_path, lbp_params) }
            .set { feed_me_into_lbp }

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_result_file_flattened
            .combine(umap_combinations_hash_outdir)
            .map { lbp_outdir, lbp_result_flattened_cur, umap_params_str, umap_params ->
            tuple("${lbp_outdir}/${umap_params_str}", lbp_result_flattened_cur, umap_params) }
            .set {feed_me_into_umap}

        fastlbp.out.lbp_result_file_img
            .set { lbp_img }

        umap(feed_me_into_umap)

        umap.out
            .combine(hdbscan_combinations_hash_outdir)
            .map { umap_outdir, umap_embeddings, hdbscan_params_str, hdbscan_params ->
            tuple("${umap_outdir}/${hdbscan_params_str}", umap_embeddings, hdbscan_params) }
            .set { feed_me_into_hdbscan }

        hdbscan(feed_me_into_hdbscan)

        lbp_img
            .join(feed_me_into_lbp)
            .set { lbp_and_masks }

        hdbscan.out
            .combine(lbp_and_masks)
            .filter { it[0].toString().contains(it[2].toString()) } // hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params
            .map { hdbscan_outdir, hdbscan_labels_path, lbp_outdir, lbp_img_path, img_path, pixelmask_path, patchmask_path, lbp_params ->
            tuple(hdbscan_outdir, hdbscan_labels_path, patchmask_path, lbp_img_path) }
            .set { convert_me_back_to_image }

        labels_to_patch_img(convert_me_back_to_image)
    }
}

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
