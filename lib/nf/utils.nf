def infoLog(msg) {
    // let's make it yellow on a black bg
    log.info("\u001B[93;40;1m" + msg + "\u001B[0m")
}

def checkNextflowVersion() {
    if ( !nextflow.version.matches(">=24.04") ) {
        error "The workflow requires Nextflow version 24.04 or greater, your current version is ${nextflow.version}"
    }
}

def getValueFromParamList(params_list, param_name = 'patchsize') {
    // params list of the following structure: [[param_1, value_1], [param_2, value_2], ...]
    def res = params_list.find { it[0] == param_name }
    res[1]
}

def checkNestedParameterCombinations(param_combinations_and_outputs) {
    // Mandatory structure of the input param_combinations_and_outputs: [step1_parameter_combination_hash, step1_output, 
    // step2_parameter_combination_hash, step2_output, ...]

    for (int i = 0; i < param_combinations_and_outputs.size(); i++)
        if (i >= 2 && i % 2 == 0 && 
            !param_combinations_and_outputs[i].contains(param_combinations_and_outputs[i - 2]))
            return false
    return true

    // it[2].contains(it[0]) && it[4].contains(it[2]) && it[6].contains(it[4])
}

def getFinalOutdirFromParamsCombinations(param_combinations_and_outputs) {
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

// TODO: switch to camel case to enhance consistency with Groovy conventions
def getAllOutputFilesFromParamsCombinations(param_combinations_and_outputs) {
    def even_elements = param_combinations_and_outputs.findAll { param_combinations_and_outputs.indexOf(it) % 2 != 0 }
    // println "$even_elements qweqweqwesss"
    return even_elements as Tuple
}

// TODO: refactor and optimise
def createCombinations(method_params) {
    // example: 
    /*
    [method:umap, n_components:[10, 20, 30], n_neighbors:[20, 30], min_dist:[values:[0.0, 0.1, 0.3], bind_to:n_components], random_state:42]
    */

    // total number of parameters for the method
    def method_params_num = method_params.size()
    // example: 5

    def collected_params_list = method_params.collect { key, values -> values instanceof Map ? \
    tuple(key, values.values, values.bind_to) : tuple(key, values, key) }
    // collected_params_list example: 
    /*
    [[method, umap, method], [n_components, [10, 20, 30], n_components], [n_neighbors, [20, 30], n_neighbors], [min_dist, [0.0, 0.1, 0.3], n_components], [random_state, 42, random_state]]
    */


    def collected_params = Channel.fromList(collected_params_list)

    def combinations_entities = collected_params
        .groupTuple(by:2)
        .map { params_list, values_list, bind_to_key ->
        [params_list, values_list].transpose() }

    
    // combinations_entities example:

    /* [1] groupTuple(by:2)
    [[n_components, min_dist], [[10, 20, 30], [0.0, 0.1, 0.3]], n_components]
    [[n_neighbors], [[20, 30]], n_neighbors]
    */

    /* [2] map {params_list, values_list, bind_to_key -> [params_list, values_list].transpose()}
    [[n_components, [10, 20, 30]], [min_dist, [0.0, 0.1, 0.3]]]
    [[n_neighbors, [20, 30]]]
    */

    def entities_grouped = combinations_entities
        .map { bound_parameter_values ->
        
        def transformParams = { data ->
            // parentheses needed for groovy key evaluation
            def output = data.collectEntries { iit -> [(iit[0]): iit[1]] }
            // example
            /*
            output = [n_components:[10, 20, 30], min_dist:[0.0, 0.1, 0.3]]
            */

            // values in the map above are guaranteed to be of the same length
            // due to the config parsing strategy 
            // (different params will appear in the same map only if bind_to parameter was used)

            // now combine parameters pairwise
            def res = []

            if (output.values()[0] instanceof String) {
                res << [output.keySet()[0], output.values()[0]]
            } else {
                output.values()[0].eachWithIndex {_, idx ->
                    // iterating over the indices of the values in the list (in the example: 3)
                    def combined_res = []
                    output.each { kk, vv ->
                        // for each [parameter, values_list] pair

                        // append to the combined result parameter name
                        combined_res << kk

                        // and parameter value from the list (if it's a list)
                        if (vv instanceof List)
                            combined_res << vv[idx]
                        else
                            combined_res << vv

                        // example: 
                        // combined_res = [n_components, 10, min_dist, 0.0]
                    }

                    // append to the total combinations combined_res
                    res << combined_res
                }
                
            }
            
            return res
        }
        // example:
        /*
        bound_parameter_values = [[n_components, [10, 20, 30]], [min_dist, [0.0, 0.1, 0.3]]]

        output = [n_components:[10, 20, 30], min_dist:[0.0, 0.1, 0.3]]

        res = [[n_components, 10, min_dist, 0.0], [n_components, 20, min_dist, 0.1], [n_components, 30, min_dist, 0.3]]
        */
        return transformParams(bound_parameter_values) 
        }

    // now we need to flatten obtained channel of parameter combinations
    // where some elements consist of multiple parameters bound to each other

    // so we define a flattening function
    // function
    flatten_in_my_way = { my_list ->
        // example of my_list:
        /*
        [[[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]], 
        [[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]], ...]
        */

        def flatList = []
        my_list.each { itemm ->
            // for each element in the list:
            // example of such itemm:
            /*
            [[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]]
            */
            itemm.each { bound_params_list -> 
            bound_params_list.each {alternating_param_name_value -> 
            flatList.add(alternating_param_name_value) } 
            }
        }
        return flatList
    }


    // now we need to find all combinations of parameters
    // properly taking care of bound parameters in entitied_grouped
    def all_combs_flat = entities_grouped // example: [[n_components, 10, min_dist, 0.0], [n_components, 20, min_dist, 0.1], [n_components, 30, min_dist, 0.3]]
        .map { kk -> [kk] } // enclose in additional list
        .collect() // collect in a single list 
        .combinations() // groovy method to find all combinations of elements in list
        .map { ll -> flatten_in_my_way(ll) } // flatten elements of type: [[method, umap], [n_components, 10, min_dist, 0.0], [n_neighbors, 20], [random_state, 42]]

    // example
    /*
    all_combs_flat = [method, umap, n_components, 10, min_dist, 0.0, n_neighbors, 20, random_state, 42, method, umap, n_components, 10, min_dist, 0.0, n_neighbors, 20, random_state, 42, method, umap, ...]
    */
    
    // now as we flattened the combinations list we need to group back parameter combinations for different runs
    def all_combs_final = all_combs_flat.flatMap{ kkk -> kkk } // flatten all parameter lists 
        .collate(method_params_num * 2) // group individual run combinations (if there are n parameters -> there are n values, so we group every 2*n elements)
        .map { run_parameters -> run_parameters.collate(2) } // group [parameter, value] pairs in a flattened list of all parameters for a single run

    // TODO: the method above produces duplicates, urgently need optimisation
    return all_combs_final.distinct() // TODO: optimise this
}
