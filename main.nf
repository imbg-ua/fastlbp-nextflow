#!/usr/bin/env/ nextflow

// TODO: nested arguments are not correctly overriden by the parameters file?
// params.args.clustering.method = 'k_means' // default clustering method
params.debug_flag = true
params.mode = 'normal' // [normal, lbp_only, lbp_tsv, grid]
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

// params.constargs = ''

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