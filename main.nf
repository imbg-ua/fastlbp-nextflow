#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

debug_flag = true
params.annot_suffix = "annotation"

// Default params
params.masks_dir = ""
params.imgs_dir = ""
params.img_path = ""

process fastlbp {
    tag "${img.getBaseName()}"
    debug debug_flag
    publishDir "${params.outdir}/${img.getBaseName()}", mode: "copy"

    input:
    tuple path(img), path(mask)

    output:
    path("data")
    path("data/out/${params.lbp.outfile_name}"), emit: lbp_result_file
    tuple path("data/out/${file(params.lbp.outfile_name).getBaseName()}_flattened.npy"), val("${img.getBaseName()}"), emit: lbp_result_file_flattened
    path("patchmask.npy"), emit: lbp_patchmask, optional: true
    path("pixelmask.npy"), emit: lbp_pixelmask, optional: true

    script:
    mask_path = mask ? "--img_mask ${mask}" : ""
    """
    lbp.py \
        --img_path ${img} \
        --patchsize ${params.lbp.patchsize} \
        --ncpus ${params.lbp.ncpus} \
        --outfile_name ${params.lbp.outfile_name} \
        ${mask_path} \
        --img_name ${params.lbp.img_name}
    """
}

process umap {
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple path(lbp_result_flattened), val(img_id)

    output:
    tuple path("umap_embeddings.npy"), val(img_id)

    script:
    """
    run_umap.py \
        --np_data_path ${lbp_result_flattened} \
        --n_components ${params.umap.n_components}
    """
}

process hdbscan {
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple path(data), val(img_id)

    output:
    path("hdbscan_labels.npy")

    script:
    """
    run_hdbscan.py \
        --np_data_path ${data} \
        --min_samples ${params.hdbscan.min_samples} \
        --min_cluster_size ${params.hdbscan.min_cluster_size} \
        --cluster_selection_epsilon ${params.hdbscan.cluster_selection_epsilon} \
        --gen_min_span_tree ${params.hdbscan.gen_min_span_tree}
    """
}

workflow SingleImage {
    fastlbp([params.img_path, params.mask.binmask_path ? file(params.mask.binmask_path) : []])
    umap(fastlbp.out.lbp_result_file_flattened) | hdbscan
}

workflow MultiImage {
    imgs = files("${params.imgs_dir}/*")
    
    imgs_and_masks = Channel.fromList(imgs)
        .map {it -> [
            it,
            params.masks_dir ? file(params.masks_dir) \
            / it.getBaseName() + "_${params.annot_suffix}." + it.getExtension() : []
            ]}
    fastlbp(imgs_and_masks)
    umap(fastlbp.out.lbp_result_file_flattened) | hdbscan
}

workflow Pipeline {
    if ( params.img_path && !params.imgs_dir )
        SingleImage()
    else
        MultiImage()
}
