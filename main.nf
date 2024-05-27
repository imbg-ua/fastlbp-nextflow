#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

debug_flag = true
params.annot_suffix = "annotation"

// Default params
params.masks = ""
params.imgs_dir = ""
params.img_path = ""
params.mask = ""


process get_tissue_mask {
    tag "${img_id}"
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    path(img)

    output:
    tuple path(img), path('pixelmask.npy')

    script:
    img_id = img.getBaseName()
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

    output:
    tuple path(img), path(pixelmask), path(savefile)

    script:
    img_id = img.getBaseName()
    savefile = "patchmask.npy"
    """
    get_mask.py downscale_using_patchsize \
        --img_path ${pixelmask} \
        --patchsize ${params.lbp.patchsize} \
        --savefile ${savefile}
    """

}

process fastlbp {
    tag "${img_id}"
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple path(img), path(mask), path(patchmask)

    output:
    path("data")
    tuple path("data/out/${file(params.lbp.outfile_name).getBaseName()}_flattened.npy"), val(img_id), emit: lbp_result_file_flattened

    script:
    img_id = img.getBaseName()
    mask_path = mask ? "--img_mask ${mask}" : ""
    patchmask_path = patchmask ? "--patch_mask ${patchmask}" : ""
    """
    run_lbp.py \
        --img_path ${img} \
        --patchsize ${params.lbp.patchsize} \
        --ncpus ${params.lbp.ncpus} \
        --outfile_name ${params.lbp.outfile_name} \
        ${mask_path} \
        ${patchmask_path} \
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
    path("hdbscan_labels.npy"), emit: labels_path
    val(img_id), emit: img_id

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
    if ( !params.mask ) {

        log.info("No mask mode")

        no_mask = Channel.of([params.img_path, [], []])
        fastlbp(no_mask)
        umap(fastlbp.out.lbp_result_file_flattened) | hdbscan

    } else if ( params.mask == "auto" ) {

        log.info("Otsu mask mode")

        img = Channel.of([params.img_path])
        get_tissue_mask(img) | downscale_mask | fastlbp
        umap(fastlbp.out.lbp_result_file_flattened) | hdbscan
        
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
        SingleImage()
    else
        MultiImage()
}
