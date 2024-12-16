def info_log(msg) {
    // let's make it yellow on a black bg
    log.info("\u001B[93;40;1m" + msg + "\u001B[0m")
}

def check_nextflow_version() {
    if ( !nextflow.version.matches(">=24.04") ) {
        error "The workflow requires Nextflow version 24.04 or greater, your current version is ${nextflow.version}"
    }
}

def get_value_from_param_list(params_list, param_name = 'patchsize') {
    // params list of the following structure: [[param_1, value_1], [param_2, value_2], ...]
    def res = params_list.find { it[0] == param_name }
    res[1]
}

def check_nested_parameter_combinations(param_combinations_and_outputs) {
    // Mandatory structure of the input param_combinations_and_outputs: [step1_parameter_combination_hash, step1_output, 
    // step2_parameter_combination_hash, step2_output, ...]

    for (int i = 0; i < param_combinations_and_outputs.size(); i++)
        if (i >= 2 && i % 2 == 0 && 
            !param_combinations_and_outputs[i].contains(param_combinations_and_outputs[i - 2]))
            return false
    return true

    // it[2].contains(it[0]) && it[4].contains(it[2]) && it[6].contains(it[4])
}

def get_final_outdir_from_params_combinations(param_combinations_and_outputs) {
    // assuming that all the parameters and combinations go like this:
    /*
        .multiMap { lbp_outdirr, lbp_dataa, 
                    umap_outdirr, umap_dataa, 
                    hdbscan_outdirr, hdbscan_dataa, 
                    deconvolve_outdirr, deconvolve_dataa -> ... }
    */
    // i.e. param_combinations_and_outputs has the following structure: [step1_parameter_combination_hash, step1_output, 
    // step2_parameter_combination_hash, step2_output, ...]

    // where `step1_parameter_combination_hash` also represents the temporary output folder for the step, before 
    // the outputs get moved to the combination hash-based directory in the final outdir

    // in this setting, the the second to last elemen of the list represents the final output folder of the whole pipeline

    return param_combinations_and_outputs[-2]
}


def get_all_output_files_from_params_combinations(param_combinations_and_outputs) {
    def even_elements = param_combinations_and_outputs.findAll { param_combinations_and_outputs.indexOf(it) % 2 != 0 }
    // println "$even_elements qweqweqwesss"
    return even_elements as Tuple
}