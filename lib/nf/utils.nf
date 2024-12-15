def info_log(msg) {
    // let's make it yellow on a black bg
    log.info("\u001B[93;40;1m" + msg + "\u001B[0m")
}

def check_nextflow_version() {
    if ( !nextflow.version.matches(">=24.04") ) {
        error "The workflow requires Nextflow version 24.04 or greater, your current version is ${nextflow.version}"
    }
}

def get_value_from_param_list(params_list, param_name) {
    // params list with strucute [[key_1, val_1], [key_2, val_2], ...]
    def res = params_list.find { it[0] == param_name }
    res[1]
}
