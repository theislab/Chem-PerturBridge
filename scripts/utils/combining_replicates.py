import re
import sys
import json
import logging
import argparse
import anndata
import pandas as pd
import scanpy as sc
from os import listdir
from os.path import isfile, join
from typing import Pattern, Dict, Optional
from anndata import AnnData

#Init logger
fmt = '%(asctime)s | [%(levelname)s] %(message)s'
datefmt = '%Y-%m-%d %H:%M:%S'
formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)

h1 = logging.StreamHandler(sys.stdout)
h1.setLevel(logging.INFO)
h1.addFilter(lambda log: log.levelno == logging.INFO)
h1.setFormatter(formatter)

h2 = logging.StreamHandler(sys.stderr)
h2.setLevel(logging.WARNING)
h2.setFormatter(formatter)

logger = logging.getLogger(__name__)
logger.propagate = False
logger.setLevel(logging.DEBUG)
logger.handlers = [h1, h2]

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
        ret : int
            A dummy value to return 
            in case we have not found a number.
        '''
        search = p.search(s)
        if search:
            return int(search.groups()[0])
        else:
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
    files = get_files(dir_input)
    adatas = []
    for i, file in enumerate(files):
        adata = sc.read_h5ad(join(dir_input, file))
        adata.obs['plate'] = i + 1
        adatas.append(adata)
    return anndata.concat(adatas)

def filter_cells(ads: AnnData,
                 groupby: str,
                 min_cells: int) -> AnnData:
    '''
    A function to remove aggregated pseudobulk samples
    which have the number of cells lower than threshold.

    Parameters:
    -----------
    ads : AnnData
        An AnnData object.
    groupby : groupby
        A column to aggregate pseudobulk samples by.
    min_cells : int
        A threshold to filter samples (the number of cells).
    '''
    vc = ads.obs[[groupby, 'psbulk_cells']].groupby(groupby, observed=False).sum()
    c_kept = vc[vc.psbulk_cells>min_cells].index.values
    adss = ads[ads.obs.drugname_drugconc.isin(c_kept)].copy()
    return adss

def save2csv(adss: AnnData, 
             cell_line: str, 
             dir_output: str) -> None:
    '''
    A function to save AnnData object containing the
    information about different bulk experiments within
    one cell line as separate .csv files (_X.csv, _var.csv, _obs.csv).

    Parameters:
    -----------
    adss : AnnData
        An AnnData object sliced by a cell line before.
    cell_line : str
        A cell line which the dataset is sliced by.
    dir_output : str
        A path to the directory to save files to.
    '''
    if adss.n_obs > 0:
        pd.DataFrame(adss.X).to_csv(join(dir_output, f"{cell_line}_X.csv"), index=False)
        adss.var.to_csv(join(dir_output, f"{cell_line}_var.csv"))
        adss.obs.to_csv(join(dir_output, f"{cell_line}_obs.csv"))
    else:
        logger.warning(f'Cell line {cell_line} does not contain observations after filtering')

def save_by_cell_lines(adata: AnnData, 
                       dir_output: str, 
                       groupby: str = 'drugname_drugconc', 
                       min_cells: int = 100) -> None:
    '''
    A function to slice the united pseudobulk data 
    from different plates by cell lines, filter it and save the
    sliced parts in the output directory.

    Parameters:
    -----------
    adata : 
        A united dataset.
    dir_output : 
        A path to the output directory.
    groupby : 
        A column to group by during a filtering step.
    min_cells : 
        A threshold for the number of cells in the aggregated
        pseudobulk to filter.
    '''
    cell_lines = sorted(set(adata.obs.cell_name))
    for cell_line in cell_lines:
        ads = adata[adata.obs['cell_name']==cell_line].copy()
        cell_line = re.sub(r'[^a-zA-Z0-9]', '', cell_line)
        adss = filter_cells(ads, groupby, min_cells)
        save2csv(adss, cell_line, dir_output)

def main():
    '''
    A function to process the arguments entered in the console/
    read from the config file and to run the downstream pipeline
    with the set parameters.
    '''
    def merge_args(d_args: Dict[str, Optional[int]], 
                   config: Dict[str, Optional[int]]) -> Dict[str, Optional[int]]:
        '''
        A function to unite the parameters entered as the arguments
        from the console and the parameters loaded from a config file.

        Parameters:
        -----------
        d_args : Dict
            input arguments represented as a dictionary.
        config : Dict
            parameters loaded from a config file.
        '''
        for key in config.keys():
            if (not key in d_args.keys()) or (not d_args[key]):
                d_args[key] = config[key]
        return d_args.copy()

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str)
    parser.add_argument('--input', type=str)
    parser.add_argument('--output', type=str)
    parser.add_argument('--groupby', type=str, default='drugname_drugconc')
    parser.add_argument('--min_cells', type=int, default=100)
    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']

    if not args.config is None:
        with open(args.config) as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)

    for key in ['input', 'output']:
        if not d_args[key]:
            raise Exception(f"The argument {key} is not set")
    
    logger.info("Unite datasets from several plates")
    
    adata = unite_adatas(d_args['input'])
    logger.info("Save by cell lines")
    save_by_cell_lines(adata, 
                       d_args['output'], 
                       min_cells=d_args['min_cells'], 
                       groupby=d_args['groupby'])

if __name__ == "__main__":
    main()
