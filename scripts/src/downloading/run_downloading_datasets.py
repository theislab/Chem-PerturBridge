import os
import time
import argparse
import json
from typing import Dict, Optional
import numpy as np
import lamindb as ln
from django.db.utils import OperationalError
from src.utils.parsing_utils import *


def download_save_from_lamindb(key_adata: str = None,
                           key_obs: str = None, 
                           key_var: str = None,
                           path2adata: str = None,
                           path2obs: str = None,
                           path2var: str = None,
                           subsampling: bool = False,
                           subsample_size: int = 5000,
                           instance: str = 'laminlabs/pertdata') -> None:
    
    '''
    Download datasets from lamindb and save them.

    Parameters:
    -----------
    key_adata : str
        A key of AnnData object stored in lamindb.
    key_obs : str
        A key of obs DataFrame stored in lamindb.
    key_var : str
        A key of var DataFrame stored in lamindb.
    path2adata : str = None
        Path to save AnnData object
    path2obs : str = None
        Path to save obs DataFrame
    path2var : str = None
        Path to save var DataFrame
    subsampling : bool = False
        A flag, True value means to make a subsample
        from AnnData object. It has no effect on obs 
        and var datasets.
    subsample_size : int = 5000
        The number of observations in the subsample.
    instance : str
        Instance in the lamindb. 
    '''
    def connect2lamin_with_retry(instance,
                                 max_retries=10,
                                 delay=10,
                                 ):
        retries = 0
        while retries < max_retries:
            try:
                ln.connect(instance)
                return None
            except OperationalError as e:
                logger.error('Connection failed: %s, retried %d times', str(e), retries)
                retries += 1
                time.sleep(delay)
        raise Exception('Max retries reached. Could not connect to lamindb.')

    def append_subsample_suffix(filename: str):
        '''
        Change a file name in case subsampling.

        Parameters:
        -----------
        filename : str
            A path to save AnnData object
        '''
        name, ext = os.path.splitext(filename)
        return f'{name}_subsample{ext}'
    
    if (key_adata) and (path2adata):
        if subsampling:
            path2adata = append_subsample_suffix(path2adata)
        if path2adata.startswith('~'):
            path2adata = os.path.expanduser(path2adata)
        if not os.path.isfile(path2adata):
            dir_path = os.path.dirname(path2adata)
            if not os.path.exists(dir_path):
                try:
                    os.makedirs(dir_path)
                except FileExistsError as e:
                    logger.warning('%s', str(e))
            
            logger.info('Download and save an anndata file to %s', path2adata)
            connect2lamin_with_retry(instance)
            adata = ln.Artifact.get(key=key_adata).load()
            ln.setup.disconnect()
            if subsampling:
                idx = np.random.choice(adata.obs.index, replace=False, size=subsample_size)
                adata = adata[idx].copy()
            adata.write_h5ad(path2adata, compression='gzip')
        else:
            logger.info('An anndata file already exists')
    else:
        logger.info('Skipped downloading an anndata file as some arguments are empty')
    
    if (key_obs) and (path2obs):
        if path2obs.startswith('~'):
            path2obs = os.path.expanduser(path2obs)
        if not os.path.isfile(path2obs):
            dir_path = os.path.dirname(path2obs)
            if not os.path.exists(dir_path):
                try:
                    os.makedirs(dir_path)
                except FileExistsError as e:
                    logger.warning('%s', str(e))
            
            logger.info('Download and save an .obs file to %s', path2obs)
            connect2lamin_with_retry(instance)
            ln.Artifact.get(key=key_obs).load().to_parquet(path2obs)
            ln.setup.disconnect()
        else:
            logger.info('An .obs file already exists')
    else:
        logger.info('Skipped downloading an .obs file as some arguments are empty')

    if (key_var) and (path2var):
        if path2var.startswith('~'):
            path2var = os.path.expanduser(path2var)
        if not os.path.isfile(path2var):
            dir_path = os.path.dirname(path2var)
            if not os.path.exists(dir_path):
                try:
                    os.makedirs(dir_path)
                except FileExistsError as e:
                    logger.warning('%s', str(e))
            
            logger.info('Download and save a .var file to %s', path2var)
            connect2lamin_with_retry(instance)
            ln.Artifact.get(key=key_var).load().to_parquet(path2var)
            ln.setup.disconnect()
        else:
            logger.info('A .var file already exists')
    else:
        logger.info('Skipped downloading an .var file as some arguments are empty')

    logger.info('All data is saved')



def main():
    '''
    A function to process the arguments entered in the console/
    read from the config file and to run the downstream pipeline
    with the set parameters.
    '''
    def check_input_output_sub_args(d_args: Dict[str, Optional[str]],
                                    required_sub_args: Dict[str, str]) -> None:
        '''
        A function to check if the user provided both input path to the dataset
        in the lamin database and the path to save the dataset.

        Parameters:
        -----------
        d_args : Dict[str, Optional[str]]
            Input arguments represented as a dictionary.
        required_sub_args : Dict[str, str]
            A dictionary of predetermined parameter names which
            should be used with arguments.
        '''
        for i, isub_args in enumerate(required_sub_args['input']):
            osub_args = required_sub_args['output'][i]
            if d_args['input'].get(isub_args) and (not d_args['output'].get(osub_args)):
                raise Exception(f'Both {isub_args} and {osub_args} sub-arguments should be specified')

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str)
    parser.add_argument('--input', nargs='+', action=ParseKW)
    parser.add_argument('--output', nargs='+', action=ParseKW)
    parser.add_argument('--subsampling', type=bool, default=False, action=argparse.BooleanOptionalAction)
    required_args = ['input', 'output']
    required_sub_args = {'input': ['key_adata', 'key_obs', 'key_var'],
                         'output': ['path2adata', 'path2obs', 'path2var']}

    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']

    if not args.config is None:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)

    for key in required_args:
        if not d_args.get(key):
            raise Exception(f'The argument {key} is not specified')

    check_sub_args(d_args, required_sub_args)
    check_input_output_sub_args(d_args, required_sub_args)
    download_save_from_lamindb(**{**d_args['input'], **d_args['output']}, subsampling=d_args['subsampling'])

if __name__ == '__main__':
    main()
