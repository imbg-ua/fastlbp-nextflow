#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

debug_flag = true
params.annot_suffix = "annotation"


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
println "${lbp_params_names} HELLLOOO"

println "${lbp_step_param_num} lbp params num"

def dimred_step = params.args.umap
def dimred_step_param_num = params.args.umap.size()
dimred_params_names = params.args.umap.collect { k, v -> k }
println "${dimred_params_names} HELLLOvccccccOO"

println "${dimred_step_param_num} umap params num"

def clust_step = params.args.hdbscan
def clust_step_param_num = params.args.hdbscan.size()
clust_params_names = params.args.hdbscan.collect { k, v -> k }
println "${clust_params_names} HELLLOOsdasdasdO"

println "${clust_step_param_num} scan params num"

def lbp_combinations = createCombinations(lbp_step)
def umap_combinations = createCombinations(dimred_step)
def hdbscan_combinations = createCombinations(clust_step)

hdbscan_combinations.view()

lbp_combinations.combine(umap_combinations.combine(hdbscan_combinations)).set {all_parameters_combinations}

process generate_params_combinations_dir {
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
    tag "${img_id}"
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple val(img_id), path(img)

    output:
    tuple val(img_id), path(img), path('pixelmask.npy')

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
    tuple val(img_id), path(lbp_output), path(patchmask), path(labels_data)

    output:
    tuple val(img_id), path(savefile)

    script:
    savefile = "patch_labels.npy"
    command = "masked_labels_to_patch_img --np_mask_path ${patchmask}" ? patchmask : \
    "labels_to_patch_img --lbp_output_path ${lbp_output}"
    """
    deconvolve_masked_labels.py ${command} \
        --np_labels_path ${labels_data} \
        --savefile ${savefile}
    """
}

process split_params_combinations_str_in_separate_channels {
    fair true

    input:
    val(param_str)

    output:
    val(lbp_params), emit: lbp_params_combs
    val(dimred_params), emit: dimred_params_combs
    val(clust_params), emit: clust_params_combs

    script:
    // TODO: this is the ugliest thing ever
    lbp_params = param_str.subList(0, 2 * lbp_step_param_num)

    dimred_params = param_str.subList(2 * lbp_step_param_num,
    2 * (dimred_step_param_num + lbp_step_param_num))

    clust_params = param_str.subList(2 * (dimred_step_param_num + lbp_step_param_num),
    2 * (dimred_step_param_num + lbp_step_param_num + clust_step_param_num))
    """
    """
}

workflow SingleImage {
    take:
    parameters_set

    main:

    // parameters_set.view()

    // [string, hash] is the correct order
    parameters_set
        .map { [it, it.flatten().collect().join('_').toString(),  it.flatten().collect().join('_').toString().md5()] }
        .tap {testik}
        .map { [it[2], it[0]] }
        .set { hash_and_comb }

    testik
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

        log.info("No mask mode")

        lbp_runs
            .map { hash, params_list -> [hash, params.img_path, [], [], params_list.toString()] }
            .set { no_mask }
        fastlbp(no_mask)

        fastlbp.out.lbp_result_file_flattened
            .join(dimred_runs)
            .set { feed_me_into_umap }
        umap(feed_me_into_umap)

        umap.out
            .join(clust_runs)
            .set { feed_me_into_hdbscan }
        hdbscan(feed_me_into_hdbscan)

        // hdbscan.out
        //     .join(no_mask)
        //     .view()

        // labels_to_patch_img()

    }

    // } else if ( params.mask == "auto" ) {

    //     log.info("Otsu mask mode")

    //     img = Channel.of([params.img_path])
    //     get_tissue_mask(img) | downscale_mask

    //     downscale_mask.out
    //         .map { img_path, pixelmask_path, patchmask_path ->
    //             [
    //                 patchmask_path, 
    //                 img_path.getBaseName()
    //             ]
    //          }
    //         .set { img_id_and_mask_ch }

    //     fastlbp(downscale_mask.out)
    //     umap(fastlbp.out.lbp_result_file_flattened) | hdbscan
        
    //     img_id_and_mask_ch
    //         .join(hdbscan.out, by:1)
    //         .map {
    //             imd_id, patchmask_path, clustering_labels_path -> 
    //                 [
    //                     patchmask_path, 
    //                     clustering_labels_path, 
    //                     imd_id
    //                 ]
    //         } 
    //         .set { deconvolve_ch }

    //     labels_to_patch_img(deconvolve_ch)
        
    // } else {
        
    //     log.info("Provided mask mode")

    //     img_and_mask = Channel.of([params.img_path, file(params.mask)])
    //     downscale_mask(img_and_mask) | fastlbp
    //     umap(fastlbp.out.lbp_result_file_flattened) | hdbscan

    // }
}

workflow MultiImage {
    // imgs = files("${params.imgs_dir}/*")

    // if ( !params.masks ) {

    //     log.info("No mask mode")

    //     imgs_and_masks = Channel.fromList(imgs)
    //         .map {it -> [
    //             it, [], []
    //             ]}
    //     fastlbp(imgs_and_masks)
    //     umap(fastlbp.out.lbp_result_file_flattened) | hdbscan

    // } else if ( params.masks == "auto" ) {

    //     log.info("Otsu mask mode")

    //     imgs_and_masks = Channel.fromList(imgs)
    //     get_tissue_mask(imgs_and_masks) | downscale_mask | fastlbp
    //     umap(fastlbp.out.lbp_result_file_flattened) | hdbscan
        
    // } else {
        
    //     log.info("Provided mask mode")

    //     imgs_and_masks = Channel.fromList(imgs)
    //         .map {it -> [
    //             it,
    //             file(params.masks) \
    //             / it.getBaseName() + "_${params.annot_suffix}." + it.getExtension()
    //             ]}
        
    //     downscale_mask(imgs_and_masks) | fastlbp
    //     umap(fastlbp.out.lbp_result_file_flattened) | hdbscan

    // }
}

workflow Pipeline {
    if ( params.img_path && !params.imgs_dir )    
        SingleImage(all_parameters_combinations)
        
    else
        MultiImage()
}
