import re
import sys
import os
import json
import argparse
from os import listdir
from os.path import isfile, join
from typing import Pattern, Dict, Optional
import anndata
import scanpy as sc
from anndata import AnnData

from src.utils.parsing_utils import *

def get_files(dir_input: str) -> list:
    '''
    A function to get a list of all files stored
    in the certain directory.

    Parameters:
    -----------
    dir_input : str
        A path to the directory which contains files.
    '''
    def extract_num(s: str, p: Pattern, ret: int | float = 0) -> int | float:
        '''
        A function to extract the pattern (ex. numbers) 
        from a string.

        Parameters:
        -----------
        s : str
            An input string.
        p : Pattern
            A pattern for extraction.
        ret : int | float
            A dummy value to return 
            in case we have not found a number.
        '''
        search = p.search(s)
        if search:
            return int(search.groups()[0])
        return ret

    p = re.compile(r'(\d+)(\.h5ad)')

    files = [file for file in listdir(dir_input) if isfile(join(dir_input, file))]
    return sorted(files, key=lambda s: extract_num(s, p, float('inf')))

def unite_adatas(dir_input: str) -> AnnData:
    '''
    A function to read and unite the pseudobulk
    data stored in separate files (ex. one file - one plate).

    Parameters:
    -----------
    dir_input : str
        A path to the directory which contains files.
    '''
    logger.info('Unite datasets from several plates')
    files = get_files(dir_input)
    adatas = []
    for i, file in enumerate(files):
        adata = sc.read_h5ad(join(dir_input, file))
        adatas.append(adata)
    return anndata.concat(adatas, merge='same')



def main():
    '''
    A function to process the arguments entered in the console/
    read from the config file and to run the downstream pipeline
    with the set parameters.
    '''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str)
    parser.add_argument('--input', type=str)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']

    if not args.config is None:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)

    for key in ['input', 'output']:
        if not d_args.get(key):
            raise Exception(f'The argument {key} is not set')


    adata = unite_adatas(d_args['input'])
    logger.info(f'Save the united dataset to {d_args["output"]}')

    dir_output = os.path.dirname(d_args['output'])
    if not os.path.exists(dir_output):
        try:
            os.makedirs(dir_output)
        except FileExistsError as e:
            logger.warning('%s', str(e))

    adata.write_h5ad(d_args['output'], compression='gzip')


if __name__ == '__main__':
    main()
