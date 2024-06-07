#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

debug_flag = true

if ( !nextflow.version.matches(">=24.04") ) {
    error "The workflow requires Nextflow version 24.04 or greater, your current version is ${nextflow.version}"
}

params.annot_suffix = "annotation"

def info_log(msg) {
    // let's make it yellow on a black bg
    log.info("\u001B[93;40;1m" + msg + "\u001B[0m")
}

def extract_patchsize_from_lbp_params_list(lbp_params_list) {
    def res = lbp_params_list.find { it[0] == "patchsize" }
    res[1]
}

def params_args_list = params.args.collect { k, v -> [k, v] }

def lbp_params = params.args.lbp.collect { k, v -> [k, v] }
def umap_params = params.args.umap.collect { k, v -> [k, v] }
def hdbscan_params = params.args.hdbscan.collect { k, v -> [k, v] }


process convert_annotations_to_binmask {
    tag "mask preprocessing"
    debug debug_flag
    publishDir "${params.outdir}/${img_id}", mode: "copy"

    input:
    tuple path(img), path(mask), val(background_val_str)

    output:
    tuple path(img), path(binmask)

    script:
    img_id = img.getBaseName()
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
    tuple path(img), path(pixelmask), val(patchsize)

    output:
    tuple path(img), path(pixelmask), path(savefile)

    script:
    img_id = img.getBaseName()
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
    tuple path(img), path(mask), path(patchmask), val(params_str)

    output:
    path("data")
    tuple val(img_id), path("data/out/${file(params.args.lbp.outfile_name).getBaseName()}_flattened.npy"), emit: lbp_result_file_flattened
    tuple val(img_id), path("data/out/${file(params.args.lbp.outfile_name).getBaseName()}.npy"), emit: lbp_result_file_img

    script:
    img_id = img.getBaseName()
    mask_path = mask ? "--img_mask ${mask}" : ""
    patchmask_path = patchmask ? "--patch_mask ${patchmask}" : ""
    """
    run_lbp.py main \
        --img_path ${img} \
        --params_str "${params_str}" \
        ${mask_path} \
        ${patchmask_path}
    """
}

process umap {
    tag "${img_id}"
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
    tag "${img_id}"
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
    tag "${img_id}"
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

// TODO: refactor redundant code
workflow SingleImage {

    if ( !params.mask ) {

        info_log("No mask mode")

        lbp_input_ch = Channel.of([params.img_path, [], [], lbp_params])

        fastlbp(lbp_input_ch)

        fastlbp.out.lbp_result_file_flattened
            .combine([[umap_params]])
            .set { feed_me_into_umap }

        fastlbp.out.lbp_result_file_img
            .set { lbp_runs }
        
        umap(feed_me_into_umap)

        umap.out
            .combine([[hdbscan_params]])
            .set { feed_me_into_hdbscan }
        
        hdbscan(feed_me_into_hdbscan)

        hdbscan.out
            .join(lbp_runs)
            .map { img_id, clust_labels, lbp_result ->
            tuple(img_id, clust_labels, [], lbp_result) }
            .set { convert_my_labels_to_img }

        labels_to_patch_img(convert_my_labels_to_img)

    } else if ( params.mask == "auto" ) {
        
        info_log("Otsu mask mode")

        get_tissue_mask(params.img_path)

        get_tissue_mask.out
            .combine([[lbp_params]])
            .map { img, pixelmask, lbp_params_list_cur ->
            tuple(img, pixelmask, extract_patchsize_from_lbp_params_list(lbp_params_list_cur))}
            .set { downscale_me }
        
        downscale_mask(downscale_me)

        downscale_mask.out
            .combine([[lbp_params]])
            .set { feed_me_into_lbp }

        feed_me_into_lbp
            .map { img, pixelmask, patchmask, lbp_params_list_cur ->
            tuple(img.getBaseName(), patchmask) }
            .set { image_and_patchmask_ch }

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_result_file_flattened
            .combine([[umap_params]])
            .set { feed_me_into_umap }

        fastlbp.out.lbp_result_file_img
            .set { lbp_runs }
        
        umap(feed_me_into_umap)

        umap.out
            .combine([[hdbscan_params]])
            .set { feed_me_into_hdbscan }
        
        hdbscan(feed_me_into_hdbscan)

        hdbscan.out
            .join(image_and_patchmask_ch)
            .join(lbp_runs)
            .set { convert_my_labels_to_img }

        labels_to_patch_img(convert_my_labels_to_img)
    } else {

        // Mask mode
        
        info_log("Provided mask mode")

        annot_to_convert = Channel.of([params.img_path, file(params.mask), 
        params.background_color ? params.background_color : ""])

        convert_annotations_to_binmask(annot_to_convert)

        convert_annotations_to_binmask.out
            .set { img_and_mask }
        
        img_and_mask
            .combine([[lbp_params]])
            .map { img, annot, lbp_params_list_cur ->
            tuple(img, annot, extract_patchsize_from_lbp_params_list(lbp_params_list_cur)) }
            .set { downscale_me }

        downscale_mask(downscale_me)

        downscale_mask.out
            .combine([[lbp_params]])
            .set { feed_me_into_lbp }

        feed_me_into_lbp
            .map { img, pixelmask, patchmask, lbp_params_list_cur ->
            tuple(img.getBaseName(), patchmask) }
            .set { image_and_patchmask_ch }

        fastlbp(feed_me_into_lbp)

        fastlbp.out.lbp_result_file_flattened
            .combine([[umap_params]])
            .set { feed_me_into_umap }

        fastlbp.out.lbp_result_file_img
            .set { lbp_runs }
        
        umap(feed_me_into_umap)

        umap.out
            .combine([[hdbscan_params]])
            .set { feed_me_into_hdbscan }
        
        hdbscan(feed_me_into_hdbscan)

        hdbscan.out
            .join(image_and_patchmask_ch)
            .join(lbp_runs)
            .set { convert_my_labels_to_img }

        labels_to_patch_img(convert_my_labels_to_img)
    }
}

workflow OtsuWorkflow {
    take:
    imgs_ch

    main:
    get_tissue_mask(imgs_ch)

    get_tissue_mask.out
        .combine([[lbp_params]])
        .map { img, pixelmask, lbp_params_list_cur ->
        tuple(img, pixelmask, extract_patchsize_from_lbp_params_list(lbp_params_list_cur))}
        .set { downscale_me }
    
    downscale_mask(downscale_me)

    downscale_mask.out
        .combine([[lbp_params]])
        .set { feed_me_into_lbp }

    feed_me_into_lbp
        .map { img, pixelmask, patchmask, lbp_params_list_cur ->
        tuple(img.getBaseName(), patchmask) }
        .set { image_and_patchmask_ch }

    fastlbp(feed_me_into_lbp)

    fastlbp.out.lbp_result_file_flattened
        .combine([[umap_params]])
        .set { feed_me_into_umap }

    fastlbp.out.lbp_result_file_img
        .set { lbp_runs }
    
    umap(feed_me_into_umap)

    umap.out
        .combine([[hdbscan_params]])
        .set { feed_me_into_hdbscan }
    
    hdbscan(feed_me_into_hdbscan)

    hdbscan.out
        .join(image_and_patchmask_ch)
        .join(lbp_runs)
        .set { convert_my_labels_to_img }

    labels_to_patch_img(convert_my_labels_to_img)
}

workflow NoMaskWorkFlow {
    take:
    imgs_and_lbp_params

    main:

    fastlbp(imgs_and_lbp_params)

    fastlbp.out.lbp_result_file_flattened
        .combine([[umap_params]])
        .set { feed_me_into_umap }

    fastlbp.out.lbp_result_file_img
        .set { lbp_runs }

    umap(feed_me_into_umap)

    umap.out
        .combine([[hdbscan_params]])
        .set { feed_me_into_hdbscan }

    hdbscan(feed_me_into_hdbscan)

    hdbscan.out
        .join(lbp_runs)
        .map { img_id, clust_labels, lbp_result ->
        tuple(img_id, clust_labels, [], lbp_result) }
        .set { convert_my_labels_to_img }

    labels_to_patch_img(convert_my_labels_to_img)
}

workflow ProvidedMaskWorkflow {
    take:
    imgs_and_masks

    main:

    downscale_mask(imgs_and_masks)

    downscale_mask.out
        .combine([[lbp_params]])
        .set { feed_me_into_lbp }

    feed_me_into_lbp
        .map { img, pixelmask, patchmask, lbp_params_list_cur ->
        tuple(img.getBaseName(), patchmask) }
        .set { image_and_patchmask_ch }

    fastlbp(feed_me_into_lbp)

    fastlbp.out.lbp_result_file_flattened
        .combine([[umap_params]])
        .set { feed_me_into_umap }

    fastlbp.out.lbp_result_file_img
        .set { lbp_runs }

    umap(feed_me_into_umap)

    umap.out
        .combine([[hdbscan_params]])
        .set { feed_me_into_hdbscan }

    hdbscan(feed_me_into_hdbscan)

    hdbscan.out
        .join(image_and_patchmask_ch)
        .join(lbp_runs)
        .set { convert_my_labels_to_img }

    labels_to_patch_img(convert_my_labels_to_img)

}

// TODO: organize repetitive code into workflows
workflow MultiImage {

    if ( params.imgs_and_masks_tsv ) {
        def imgs_and_masks = Channel.fromPath(params.imgs_and_masks_tsv)
        // tsv structure: img_path - mask_path
        imgs_and_masks
            .splitCsv(header:true, sep:'\t')
            .map { row -> 
            tuple(row.image, row.mask, row.background_color ? row.background_color : "" ) }
            .set { imgs_and_masks_ch }


        imgs_and_masks_ch
            .branch { 
                without_mask: it[1] == ""
                    return it[0]
                otsu_mask: it[1] == "auto"
                    return it[0]
                with_mask: it[1] != ""
                    return it
             }
             .set { imgs_and_masks_ch_split }

        // ------------------------- //
        // process images with masks //
        // ------------------------- //
        imgs_and_masks_ch_split.with_mask
            .set { convert_me_to_binmask }

        convert_annotations_to_binmask(convert_me_to_binmask)
        
        convert_annotations_to_binmask.out
            .map { imgg, binmaskk ->
            tuple(imgg, binmaskk,  extract_patchsize_from_lbp_params_list(lbp_params)) }
            .set { downscale_me_with_mask }
        
        ProvidedMaskWorkflow(downscale_me_with_mask)

        // ---------------------------- //
        // process images without masks //
        // ---------------------------- //
        imgs_and_masks_ch_split.without_mask
            .map { imgg -> 
            tuple(imgg, [], [], lbp_params) }
            .set { imgs_to_feed_into_lbp }
        NoMaskWorkFlow(imgs_to_feed_into_lbp)

        // ---------------------------------- //
        // process images using Otsu's method //
        // ---------------------------------- //
        imgs_and_masks_ch_split.otsu_mask
            .set { imgs_to_get_masks }
        OtsuWorkflow(imgs_to_get_masks)
           
    } else {

        imgs = files("${params.imgs_dir}/*")

        if ( !params.masks ) {
            info_log("No mask mode")

            imgs_and_masks = Channel.fromList(imgs)
                .map { it -> 
                [it, [], [], lbp_params] }
                .set { lbp_input_ch }
            
            fastlbp(lbp_input_ch)

            fastlbp.out.lbp_result_file_flattened
                .combine([[umap_params]])
                .set { feed_me_into_umap }

            fastlbp.out.lbp_result_file_img
                .set { lbp_runs }
            
            umap(feed_me_into_umap)

            umap.out
                .combine([[hdbscan_params]])
                .set { feed_me_into_hdbscan }
            
            hdbscan(feed_me_into_hdbscan)

            hdbscan.out
                .join(lbp_runs)
                .map { img_id, clust_labels, lbp_result ->
                tuple(img_id, clust_labels, [], lbp_result) }
                .set { convert_my_labels_to_img }

            labels_to_patch_img(convert_my_labels_to_img)
        } else if ( params.masks == "auto" ) {
            info_log("Otsu mask mode")
            
            imgs_ch = Channel.fromList(imgs)
            get_tissue_mask(imgs_ch)

            get_tissue_mask.out
                .combine([[lbp_params]])
                .map { img, pixelmask, lbp_params_list_cur ->
                tuple(img, pixelmask, extract_patchsize_from_lbp_params_list(lbp_params_list_cur))}
                .set { downscale_me }

            downscale_mask(downscale_me)

            downscale_mask.out
                .combine([[lbp_params]])
                .set { feed_me_into_lbp }

            feed_me_into_lbp
                .map { img, pixelmask, patchmask, lbp_params_list_cur ->
                tuple(img.getBaseName(), patchmask) }
                .set { image_and_patchmask_ch }

            fastlbp(feed_me_into_lbp)

            fastlbp.out.lbp_result_file_flattened
                .combine([[umap_params]])
                .set { feed_me_into_umap }

            fastlbp.out.lbp_result_file_img
                .set { lbp_runs }
            
            umap(feed_me_into_umap)

            umap.out
                .combine([[hdbscan_params]])
                .set { feed_me_into_hdbscan }
            
            hdbscan(feed_me_into_hdbscan)

            hdbscan.out
                .join(image_and_patchmask_ch)
                .join(lbp_runs)
                .set { convert_my_labels_to_img }

            labels_to_patch_img(convert_my_labels_to_img)
        } else {
            info_log("Masks directory mode")

            imgs_and_masks = Channel.fromList(imgs)
                .map {it -> 
                    tuple(it, file(params.masks) \
                    / it.getBaseName() + "_${params.annot_suffix}." + it.getExtension())
                }
                .combine([[lbp_params]])
                .map { img, pixelmask, lbp_params_list_cur ->
                tuple(img, pixelmask, extract_patchsize_from_lbp_params_list(lbp_params_list_cur))}
                .set { downscale_me }

            downscale_mask(downscale_me)

            downscale_mask.out
                .combine([[lbp_params]])
                .set { feed_me_into_lbp }

            feed_me_into_lbp
                .map { img, pixelmask, patchmask, lbp_params_list_cur ->
                tuple(img.getBaseName(), patchmask) }
                .set { image_and_patchmask_ch }

            fastlbp(feed_me_into_lbp)

            fastlbp.out.lbp_result_file_flattened
                .combine([[umap_params]])
                .set { feed_me_into_umap }

            fastlbp.out.lbp_result_file_img
                .set { lbp_runs }
            
            umap(feed_me_into_umap)

            umap.out
                .combine([[hdbscan_params]])
                .set { feed_me_into_hdbscan }
            
            hdbscan(feed_me_into_hdbscan)

            hdbscan.out
                .join(image_and_patchmask_ch)
                .join(lbp_runs)
                .set { convert_my_labels_to_img }

            labels_to_patch_img(convert_my_labels_to_img)
        }

    }    
}

workflow Pipeline {
    if ( params.img_path && !params.imgs_dir )
        SingleImage()
    else
        MultiImage()
}
