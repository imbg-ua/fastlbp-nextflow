#!/usr/bin/env/ nextflow

process generate_report {
    debug params.debug_flag
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