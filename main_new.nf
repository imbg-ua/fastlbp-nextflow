#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

debug_flag = true
params.annot_suffix = "annotation"

def info_log(msg) {
    // let's make it yellow on a black bg
    log.info("\u001B[93;40;1m" + msg + "\u001B[0m")
}

def extract_patchsize_from_lbp_params_list(lbp_params_list) {
    def res = lbp_params_list.find { it[0] == "patchsize" }
    res[1]
}

def createCombinations(method_params) {
    def channels = method_params.collect { key, values -> values instanceof List ? \
    Channel.of([key, values]).transpose() : Channel.of([key, values])} // Channel.from()
    def combinedChannel = channels[0]
    combinedChannel = combinedChannel.map { [it] }

    // https://github.com/nextflow-io/nextflow/issues/1403
    channels.tail().each {
        combinedChannel = combinedChannel.combine(it.map { itit -> [itit] })
    }
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

lbp_combinations.combine(umap_combinations.combine(hdbscan_combinations)).set {all_parameters_combinations}

process generate_params_combinations_dir {
    tag "setup"
    debug debug_flag
    publishDir "${params.outdir}", mode: "copy"

    input:
    val(params_combination_string)

    output:
    path(outfile_name)

    script:
    outfile_name = "hash_to_combination.tsv"
    """
    params_to_hash.py \
        --params_list_str "${params_combination_string}" \
        --savefile ${outfile_name}
    """
}

process get_tissue_mask {
    tag "preprocessing"
    debug debug_flag
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple path(img)

    output:
    tuple path(img), path('pixelmask.npy')

    script:
    """
    get_mask.py get_mask \
        --img_path ${img} \
    """
}

process downscale_mask {
    tag "${img_id}"
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple val(img_id), path(img), path(pixelmask), val(patchsize)

    output:
    tuple val(img_id), path(img), path(pixelmask), path(savefile)

    script:
    savefile = "patchmask.npy"
    """
    get_mask.py downscale_using_patchsize \
        --img_path ${pixelmask} \
        --patchsize ${patchsize} \
        --savefile ${savefile}
    """
}

process fastlbp {
    tag "${img_id}"
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple val(img_id), path(img), path(mask), path(patchmask), val(params_str)

    output:
    path("data")
    tuple val(img_id), path("data/out/${file(params.constargs.lbp.outfile_name).getBaseName()}_flattened.npy"), emit: lbp_result_file_flattened
    tuple val(img_id), path("data/out/${file(params.constargs.lbp.outfile_name).getBaseName()}.npy"), emit: lbp_result_file_img

    script:
    mask_path = mask ? "--img_mask ${mask}" : ""
    patchmask_path = patchmask ? "--patch_mask ${patchmask}" : ""
    """
    run_lbp.py \
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
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple val(img_id), path(lbp_result_flattened), val(params_str)

    output:
    tuple val(img_id), path("umap_embeddings.npy")

    script:
    """
    run_umap.py \
        --np_data_path ${lbp_result_flattened} \
        --params_str "${params_str}"
    """
}

process hdbscan {
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple val(img_id), path(data), val(params_str)

    output:
    tuple val(img_id), path("hdbscan_labels.npy")

    script:
    """
    run_hdbscan.py \
        --np_data_path ${data} \
        --params_str "${params_str}"
    """
}

process labels_to_patch_img {
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple val(img_id), path(labels_data), path(patchmask), path(lbp_output)

    output:
    tuple val(img_id), path(savefile)

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
    take:
    parameters_set

    main:
    // [string, hash] is the correct order
    parameters_set
        .map { [it, it.flatten().collect().join('_').toString(),  it.flatten().collect().join('_').toString().md5()] }
        .tap {hash_and_str}
        .map { [it[2], it[0]] }
        .set { hash_and_comb }

    hash_and_str
        .map { [it[1], it[2]] }
        .toList().set { combinations_to_save }

    // generate tsv file mappign hash to parameter combination
    generate_params_combinations_dir(combinations_to_save)

    hash_and_comb.transpose().transpose().collate(2)
        .map { tuple(it[0][0], tuple(it[0][1], it[1][1])) }
        .branch {
            lbp: it[1][0] in lbp_params_names
            dimred: it[1][0] in dimred_params_names
            clust: it[1][0] in clust_params_names
        }
        .set { hash_and_combs_split } 
    
    hash_and_combs_split.dimred.groupTuple().set { dimred_runs }
    hash_and_combs_split.lbp.groupTuple().set { lbp_runs }
    hash_and_combs_split.clust.groupTuple().set { clust_runs }

    if ( !params.mask ) {

        info_log("No mask mode")

        lbp_runs
            .map { hash, params_list -> [hash, params.img_path, [], [], params_list.toString()] }
            .set { no_mask }
        fastlbp(no_mask)

        fastlbp.out.lbp_result_file_flattened
            .join(dimred_runs)
            .set { feed_me_into_umap }

        fastlbp.out.lbp_result_file_img
            .set { lbp_runs }

        umap(feed_me_into_umap)

        umap.out
            .join(clust_runs)
            .set { feed_me_into_hdbscan }
        hdbscan(feed_me_into_hdbscan)

        hdbscan.out
            .join(no_mask)
            .join(lbp_runs)
            .map { img_id, clust_labels, img_path, pixelmask, patchmask, lbp_run_params, lbp_result_img ->
            tuple(img_id, clust_labels, patchmask, lbp_result_img)}
            .set { convert_my_labels_to_img }

        labels_to_patch_img(convert_my_labels_to_img)
    } else if ( params.mask == "auto" ) {

        info_log("Otsu mask mode")

        img = Channel.of([params.img_path])
        get_tissue_mask(img)

        get_tissue_mask.out
            .combine(lbp_runs)
            .map { img, pixelmask, img_id, lbp_params ->
            tuple(img_id, img, pixelmask, extract_patchsize_from_lbp_params_list(lbp_params))
            }
            .set {downscale_me}

        downscale_mask(downscale_me)

        downscale_mask.out
            .join(lbp_runs)
            .set { feed_me_into_lbp }

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_result_file_flattened
            .join(dimred_runs)
            .set { feed_me_into_umap }

        fastlbp.out.lbp_result_file_img
            .set { lbp_runs }

        umap(feed_me_into_umap)

        umap.out
            .join(clust_runs)
            .set { feed_me_into_hdbscan }
        hdbscan(feed_me_into_hdbscan)

        hdbscan.out
            .join(feed_me_into_lbp)
            .map { img_id, clust_labels, img, pixelmask, patchmask, lbp_params ->
            tuple(img_id, clust_labels, patchmask) }
            .join(lbp_runs)
            .set { convert_my_labels_to_img }

        labels_to_patch_img(convert_my_labels_to_img)
    } else {

        info_log("Provided mask mode")
        img_and_mask = Channel.of([params.img_path, file(params.mask)])
        
        img_and_mask
            .combine(lbp_runs)
            .map { img, annot, img_id, lbp_params_str ->
            tuple(img_id, img, annot, extract_patchsize_from_lbp_params_list(lbp_params_str)) }
            .set { downscale_me }
        
        downscale_mask(downscale_me)

        downscale_mask.out
            .join(lbp_runs)
            .set { feed_me_into_lbp }

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_result_file_flattened
            .join(dimred_runs)
            .set { feed_me_into_umap }

        fastlbp.out.lbp_result_file_img
            .set { lbp_runs }

        umap(feed_me_into_umap)

        umap.out
            .join(clust_runs)
            .set { feed_me_into_hdbscan }
        hdbscan(feed_me_into_hdbscan)

        hdbscan.out
            .join(feed_me_into_lbp)
            .map { img_id, clust_labels, img, pixelmask, patchmask, lbp_params ->
            tuple(img_id, clust_labels, patchmask) }
            .join(lbp_runs)
            .set { convert_my_labels_to_img }

        labels_to_patch_img(convert_my_labels_to_img)
    }
}

workflow Pipeline {
    if ( params.img_path && !params.imgs_dir )    
        SingleImage(all_parameters_combinations)
    else {
        log.error("""
        Parameter search is available only in SingleImage mode.
        Make sure you specified correct image path.
        """)
    }
}
