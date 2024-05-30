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

println "${lbp_step_param_num} lbp params num"

def dimred_step = params.args.umap
def dimred_step_param_num = params.args.umap.size()

println "${dimred_step_param_num} umap params num"

def clust_step = params.args.hdbscan
def clust_step_param_num = params.args.hdbscan.size()

println "${clust_step_param_num} scan params num"

def lbp_combinations = createCombinations(lbp_step)
def umap_combinations = createCombinations(dimred_step)
def hdbscan_combinations = createCombinations(clust_step)

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
    path(img), val(img_id)

    output:
    tuple path(img), path('pixelmask.npy')

    script:
    img_id = img_id ? img_id : img.getBaseName()
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
    tuple path(img), path(pixelmask)
    val(patchsize)

    output:
    tuple path(img), path(pixelmask), path(savefile)

    script:
    img_id = img_id ? img_id : img.getBaseName()
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
    tuple path(img), path(mask), path(patchmask), val(img_id)
    val(patchsize)

    output:
    path("data")
    tuple path("data/out/${file(params.constargs.lbp.outfile_name).getBaseName()}_flattened.npy"), val(img_id), emit: lbp_result_file_flattened

    script:
    img_id = img_id ? img_id : img.getBaseName()
    mask_path = mask ? "--img_mask ${mask}" : ""
    patchmask_path = patchmask ? "--patch_mask ${patchmask}" : ""
    """
    run_lbp.py \
        --img_path ${img} \
        --patchsize ${patchsize} \
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
    tuple path(lbp_result_flattened), val(img_id)
    val(n_components)

    output:
    tuple path("umap_embeddings.npy"), val(img_id)

    script:
    """
    run_umap.py \
        --np_data_path ${lbp_result_flattened} \
        --n_components ${n_components}
    """
}

process hdbscan {
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple path(data), val(img_id)
    tuple val(min_samples), val(min_cluster_size), val(cluster_selection_epsilon), val(gen_min_span_tree)

    output:
    tuple path("hdbscan_labels.npy"), val(img_id)

    script:
    """
    run_hdbscan.py \
        --np_data_path ${data} \
        --min_samples ${min_samples} \
        --min_cluster_size ${min_cluster_size} \
        --cluster_selection_epsilon ${cluster_selection_epsilon} \
        --gen_min_span_tree ${gen_min_span_tree}
    """
}

process labels_to_patch_img {
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple path(patchmask), path(labels_data), val(img_id)

    output:
    tuple path(savefile), val(img_id)

    script:
    savefile = "patch_labels.npy"
    """
    deconvolve_masked_labels.py labels_to_patch_img \
        --np_labels_path ${labels_data} \
        --np_mask_path ${patchmask} \
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
    lbp_params = param_str.subList(0, 2 * lbp_step_param_num)
    dimred_params = param_str.subList(2 * lbp_step_param_num, 2 * (dimred_step_param_num + lbp_step_param_num))
    clust_params = param_str.subList(2 * (dimred_step_param_num + lbp_step_param_num), 2 * (dimred_step_param_num + lbp_step_param_num + clust_step_param_num))
    """
    """
}

workflow SingleImage {
    take:
    parameters_set

    main:

    parameters_set
        .map { it -> it.flatten().toList() }
        .set { testik }

    split_params_combinations_str_in_separate_channels(testik)

    split_params_combinations_str_in_separate_channels.out.lbp_params_combs
        .set { lbp_params_combss }
    split_params_combinations_str_in_separate_channels.out.dimred_params_combs
        .set { dimred_params_combss }
    split_params_combinations_str_in_separate_channels.out.clust_params_combs
        .set { clust_params_combss }

    lbp_params_combss.view()
    dimred_params_combss.view()
    clust_params_combss.view()

    parameters_set
        .flatten()
        .collate(2) // pairs <parameter_name, parameters_value>
        .set { collated_parameters_set_flat }

    // collated_parameters_set_flat.view()
    // collated_parameters_set_flat.transpose().groupTuple()
    //     .map { nn -> [nn[0], nn[1].unique()] }
    //     .set {  }

    collated_parameters_set_flat
        .branch { param_def ->
            patchsize: param_def[0] == "patchsize"
                return param_def[1]

            n_components: param_def[0] == "n_components"
                return param_def[1]
            n_neighbors: param_def[0] == "n_neighbors"
                return param_def[1]
            min_dist: param_def[0] == "min_dist"
                return param_def[1]
            random_state: param_def[0] == "random_state"
                return param_def[1]

            min_samples: param_def[0] == "min_samples"
                return param_def[1]
            min_cluster_size: param_def[0] == "min_cluster_size"
                return param_def[1]
            cluster_selection_epsilon: param_def[0] == "cluster_selection_epsilon"
                return param_def[1]
            gen_min_span_tree: param_def[0] == "gen_min_span_tree"
                return param_def[1]
        }
        .set {parameters_split}


    parameters_split.min_samples
            .combine(parameters_split.min_cluster_size)
            .combine(parameters_split.cluster_selection_epsilon)
            .combine(parameters_split.gen_min_span_tree)
            .distinct()
            .set { all_params_hdbscan_ch }

    
    parameters_set
        .map { list -> list.flatten().collect().join('_') }
        .map { it.toString() }
        .tap { params_combinations_strings }
        .map { it.md5() }
        .set { run_id_ch }
    
    params_combinations_strings
        .map { [1, it] }
        .set { params_combinations_strings_with_key }

    run_id_ch
        .map { [1, it] }
        .set { run_id_ch_with_key }
    
    params_combinations_strings_with_key.cross(run_id_ch_with_key)
        .map { [it[0][1], it[1][1]] }
        .collect()
        .set { combinations_to_save }

    generate_params_combinations_dir(combinations_to_save)

    if ( !params.mask ) {

        log.info("No mask mode")

        run_id_ch
            .map { [params.img_path, [], [], it] }
            .set { no_mask }

        fastlbp(no_mask, parameters_split.patchsize)

        umap(fastlbp.out.lbp_result_file_flattened, parameters_split.n_components)

        hdbscan(umap.out, all_params_hdbscan_ch)

    } else if ( params.mask == "auto" ) {

        log.info("Otsu mask mode")

        img = Channel.of([params.img_path])
        get_tissue_mask(img) | downscale_mask

        downscale_mask.out
            .map { img_path, pixelmask_path, patchmask_path ->
                [
                    patchmask_path, 
                    img_path.getBaseName()
                ]
             }
            .set { img_id_and_mask_ch }

        fastlbp(downscale_mask.out)
        umap(fastlbp.out.lbp_result_file_flattened) | hdbscan
        
        img_id_and_mask_ch
            .join(hdbscan.out, by:1)
            .map {
                imd_id, patchmask_path, clustering_labels_path -> 
                    [
                        patchmask_path, 
                        clustering_labels_path, 
                        imd_id
                    ]
            } 
            .set { deconvolve_ch }

        labels_to_patch_img(deconvolve_ch)
        
    } else {
        
        log.info("Provided mask mode")

        img_and_mask = Channel.of([params.img_path, file(params.mask)])
        downscale_mask(img_and_mask) | fastlbp
        umap(fastlbp.out.lbp_result_file_flattened) | hdbscan

    }
}

workflow MultiImage {
    imgs = files("${params.imgs_dir}/*")

    if ( !params.masks ) {

        log.info("No mask mode")

        imgs_and_masks = Channel.fromList(imgs)
            .map {it -> [
                it, [], []
                ]}
        fastlbp(imgs_and_masks)
        umap(fastlbp.out.lbp_result_file_flattened) | hdbscan

    } else if ( params.masks == "auto" ) {

        log.info("Otsu mask mode")

        imgs_and_masks = Channel.fromList(imgs)
        get_tissue_mask(imgs_and_masks) | downscale_mask | fastlbp
        umap(fastlbp.out.lbp_result_file_flattened) | hdbscan
        
    } else {
        
        log.info("Provided mask mode")

        imgs_and_masks = Channel.fromList(imgs)
            .map {it -> [
                it,
                file(params.masks) \
                / it.getBaseName() + "_${params.annot_suffix}." + it.getExtension()
                ]}
        
        downscale_mask(imgs_and_masks) | fastlbp
        umap(fastlbp.out.lbp_result_file_flattened) | hdbscan

    }
}

workflow Pipeline {
    if ( params.img_path && !params.imgs_dir )    
        SingleImage(all_parameters_combinations)
        
    else
        MultiImage()
}
