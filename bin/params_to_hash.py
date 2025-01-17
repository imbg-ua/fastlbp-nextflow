#!/usr/bin/env python3

import fire

import workflow_utils as ut


def main(params_list_str: str, 
    savefile: str='hash_to_combination.tsv') -> None:
    # params_list = [ [[[param_1, value_1], [param_2, value_2]]/[[method_2_param_1, value_1], ...], hash_1], ... ]

    # TODO: make this versatile and not hardcoded
    params_list_str = params_list_str.strip('[ ]') # remove one leading and one trailing bracket
    params_list = list(str(params_list_str).split(','))

    assert len(params_list) % 2 == 0

    tsv_content = 'hash\tparameters\n'

    # TODO: make this versatile and not hardcoded
    for i in range(0, len(params_list), 2):
        # add exactly 2 opening and closing brackets for complete list syntax after stripping
        param_string = '[[' + params_list[i].strip('[ ]') + ']]'
        # remove bracket leftovers enclosing hash values
        param_hash = params_list[i + 1].strip('[ ]')
        tsv_content += f'{param_hash}\t{param_string}\n'

    with open(savefile, 'w') as f:
        f.write(tsv_content)

if __name__ == '__main__':
    fire.Fire(main)