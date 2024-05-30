#!/usr/bin/env python3

import fire

def main(params_list_str: str, 
    savefile: str='hash_to_combination.tsv') -> None:
    # params_list = list(str(params_list_str).strip('[]').split(',')) # obsolete version
    params_list = list(map(lambda x: x.strip('[ ]'), str(params_list_str).split(',')))

    assert len(params_list) % 2 == 0

    tsv_content = 'hash\tparameters\n'

    for i in range(0, len(params_list), 2):
        param_string = params_list[i]
        param_hash = params_list[i + 1]
        tsv_content += f'{param_hash}\t{param_string}\n'

    with open(savefile, 'w') as f:
        f.write(tsv_content)

if __name__ == '__main__':
    fire.Fire(main)