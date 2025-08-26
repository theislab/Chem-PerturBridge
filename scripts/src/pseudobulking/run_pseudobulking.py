import json
import argparse
import psutil
import resource
import math
from typing import List, Dict, Optional

from src.utils.parsing_utils import *
from .standardization import *
from .pseudobulk import Pseudobulk

MEMORY_RATIO = 0.9

def setMemoryLimit(n_bytes: int):
    '''
    Force Python to raise an exception when it uses more than
    n_bytes bytes of memory.
    source: https://medium.com/metabob/chasing-memory-spikes-and-leaks-in-python-172ae99290d3
    
    Parameters:
    -----------
    n_bytes : int
        Limitations for the memory
    '''
    if n_bytes <= 0: return
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (n_bytes, hard))
    soft, hard = resource.getrlimit(resource.RLIMIT_DATA)
    if n_bytes < soft*1024:
        resource.setrlimit(resource.RLIMIT_DATA, (n_bytes, hard))


setMemoryLimit(math.ceil(psutil.virtual_memory().total * MEMORY_RATIO))


def main():
    '''
    A function to process the arguments entered in the console/
    read from the config file and to run the downstream pipeline
    with the set parameters.
    '''
    def check_sub_arg_values(d_args: Dict[str, Optional[str]], 
                             arg_names: List[str]) -> Dict[str, Optional[str]]:
        '''
        A function to set int type to the dictionary values
        associated with the filtering parameters (e.g. min_counts, 
        min_genes, etc.)

        Parameters:
        -----------
        d_args : Dict[str, Optional[str]]
            A dictionary containing the information about arguments
        arg_names : List[str]
            A list of filtering parameters to check
        '''
        for name in arg_names:
            if d_args.get(name):
                for par in d_args[name].keys():
                    try:
                        d_args[name][par] = int(d_args[name][par])
                    except:
                        raise Exception("{par}'s value should be integer")
        return d_args.copy()

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str)
    parser.add_argument('--dataset_name', type=str)
    parser.add_argument('--input', nargs='+', action=ParseKW)
    parser.add_argument('--output', type=str)
    parser.add_argument('--filter_malat1', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_low_counts', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_nans', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_cells_params', nargs='+', action=ParseKW)
    parser.add_argument('--filter_genes_params', nargs='+', action=ParseKW)
    parser.add_argument('--groupby_fields', nargs='+', default=['plate', 'well', 'perturbagen', 'cell_type', 'guide'])
    parser.add_argument('--drug_col', type=str, default='perturbagen')
    parser.add_argument('--ignore_cell_lines', nargs='+', default=[])
    required_args = ['dataset_name', 'input', 'output', 'groupby_fields', 'drug_col']
    required_sub_args = {'input': ['path2adata', 'path2obs', 'path2var'],
                         'filter_cells_params': ['min_counts', 'min_genes', 'max_counts', 'max_genes'],
                         'filter_genes_params': ['min_counts', 'min_cells', 'max_counts', 'max_cells']}

    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']

    if not args.config is None:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)

    for key in required_args:
        if not d_args.get(key):
            if key == 'dataset_name':
                raise Exception(f"The argument {key} is not set. Chose the dataset name from {dataset_standardization.keys()}")
            else:
                raise Exception(f"The argument {key} is not set")
    
    check_sub_args(d_args, required_sub_args)
    check_sub_arg_values(d_args, ['filter_cells_params', 'filter_genes_params'])

    pseudo = Pseudobulk(
                        files_input=d_args['input'],
                        file_output=d_args['output'],
                        groupby_fields=d_args['groupby_fields'],
                        dataset_standardization_=dataset_standardization[d_args['dataset_name']],
                        filter_malat1=d_args['filter_malat1'],
                        filter_low_counts=d_args['filter_low_counts'],
                        filter_nans=d_args['filter_nans'],
                        filter_cells_params=d_args['filter_cells_params'],
                        filter_genes_params=d_args['filter_genes_params'],
                        ignore_cell_lines=d_args['ignore_cell_lines'],
                        )
    pseudo.run_pseudobulking()
    logger.info("Mapping drugs to PubChem idx")
    pseudo.add_pubchem_cids_to_padata({}, drug_col=d_args['drug_col'])
    logger.info("Process pseudobulk")
    pseudo.process_pseudobulk()
    logger.info(f"Save dataset to {pseudo.file_output}")
    pseudo.save()



if __name__ == "__main__":
    main()
