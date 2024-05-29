#!/usr/bin/env/ nextflow

params.output_dir = './preprocessing'

process cropTiffToSuperPixel {
    publishDir params.output_dir, mode: 'move'

    input:
    each path(image)
    val(patchsize)

    output:
    path(out)

    script:
    out = "${image.baseName}_cropped_to_superpixel_patchsize_${patchsize}.tif"
    """
    preprocess_image.py \
    --img_path ${image} \
    --patchsize ${patchsize} \
    --savepath ${out}
    """


}

workflow {
    data = channel.fromPath('/nfs/lbp_project/GBM/AT10-BRA-5-FO-1-A1/spatial/AT10-BRA-5-FO-1-A1.ndpi')
    patchsize = channel.of(20, 100)
    cropTiffToSuperPixel(data, patchsize)
}
