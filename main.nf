#!/usr/bin/env/ nextflow

// TODO: nested arguments are not correctly overriden by the parameters file?
// params.args.clustering.method = 'k_means' // default clustering method
params.debug_flag = true
params.mode = 'grid' // [normal, lbp_only, lbp_tsv]
params.fair = false
if ( params.mode == 'lbp_tsv' ) {
    params.fair = true // for lbp_tsv fair directive is needed. TODO: optimise this
}

include { checkNextflowVersion } from './lib/nf/utils'

// TODO: change module file names
include { Pipeline as GridPipeline } from './subworkflows/grid_search'
include { Pipeline as NormalPipeline} from './subworkflows/normal'
include { MultiImageLBPTSV } from './subworkflows/other'

nextflow.enable.dsl = 2

pipeline_version = "0.0.2"

// // Default params
// // TODO: move to config
// // TODO: move before imports
// params.background_color = ""
// params.pairs_max_csv = 'pairs_max_jacc.csv'
// // params.pairs_max_csv = '' // TODO: BUG: Make dependent oт the availability of annotaitons
// params.annot_legend_path = ""
// params.plots_backend = 'matplotlib'

// debug_flag = true

// params.constargs = ''

println "fair in main : ${params.fair} DEBUG "
checkNextflowVersion()

// entry point for the whole tool
workflow {
    // TODO: choose better mode selection condition
    if ( params.mode == 'grid' ) {
        GridPipeline()
    } else if ( params.mode == 'lbp_tsv' ) {
        // TODO: this is more like MultiImage Grid mode
        MultiImageLBPTSV()
    } else if ( params.mode == 'normal' || params.mode == 'lbp_only' ) {
        NormalPipeline()
    } else {
        log.error("""
        Used unknown or unsupported mode: `${params.mode}`. 
        Currently available execution modes are: lbp_only, lbp_tsv, normal, grid 
        """)
    }
}