import os
import sys
import json
import logging
import argparse
from typing import List, Dict, Optional, Sequence, Any
import pandas as pd
import scanpy as sc
import decoupler as dc
import pubchempy as pcp
from pandas.api.types import is_float_dtype
from anndata import AnnData

#Init logger
FMT = '%(asctime)s | [%(levelname)s] %(message)s'
DATEFMT = '%Y-%m-%d %H:%M:%S'
formatter = logging.Formatter(fmt=FMT, datefmt=DATEFMT)

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


class ParseKW(argparse.Action):
    '''
    A class to parse the dictionary-like input args.
    From https://sumit-ghosh.com/posts/parsing-dictionary-key-value-pairs-kwargs-argparse-python/

    Parameters:
    -----------
    parser: argparse.ArgumentParser
        The parser object. 
    namespace: argparse.Namespace
        An object, which holds attributes and returns it.
    values: List[str]
        A list of values which follow the argument.
    option_string: Optional[List[str]] = None
        The option string that is used to invoke the action.     
    '''
    def __call__(self,
                 parser: argparse.ArgumentParser,
                 namespace: argparse.Namespace,
                 values: str | Sequence[Any] | None,
                 option_string: str | None = None
                 ):

        setattr(namespace, self.dest, {})
        if values:
            for value in values:
                key, value = value.split('=')
                getattr(namespace, self.dest)[key] = int(value)

class Pseudobulk:
    '''
    A class to process raw scRNA-seq AnnData and
    aggregate counts to pseudobulk.

    Parameters:
    -----------
    file_input : str
        A path to the input file containing raw data.
    file_output : str
        A path to the output file containing pseudobulk data.
    bulk_fields : List[str]
        A list of columns to construct sample_id based on them.
    filter_cells_params : Optional[Dict[str, int]]
        A dictionary containing the params for cell filtering.
    filter_genes_params : Optional[Dict[str, int]]
        A dictionary containing the params for gene filtering.
    dropna : Optional[List[str]]
        A list of columns containing NaN rows to filter.
    drop_cols : Optional[List[str]]
        A list of columns to drop.
    select_cols : Optional[List[str]]
        A list of columns to select.

    Attributes:
    -----------
    file_input : str
        A path to the input file containing raw data.
    file_output : str
        A path to the output file containing pseudobulk data.
    bulk_fields : List[str]
        A list of columns to construct sample_id based on them.
    filter_cells_params : Optional[Dict[str, int]]
        A dictionary containing the params for cell filtering.
    filter_genes_params : Optional[Dict[str, int]]
        A dictionary containing the params for gene filtering.
    dropna : Optional[List[str]]
        A list of columns containing NaN rows to filter.
    drop_cols : Optional[List[str]]
        A list of columns to drop.
    select_cols : Optional[List[str]]
        A list of columns to select.
    padata : AnnData
        A pseudobulk dataset.
    '''
    def __init__(
                 self,
                 file_input: str,
                 file_output: str,
                 bulk_fields: List[str],
                 filter_cells_params: Optional[Dict[str, int]] = None,
                 filter_genes_params: Optional[Dict[str, int]] = None,
                 dropna: Optional[List[str]] = None,
                 drop_cols: Optional[List[str]] = None,
                 select_cols: Optional[List[str]] = None
                 ):

        self.file_input = file_input
        dir_output = os.path.dirname(file_output)
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)
        self.file_output = file_output

        if filter_cells_params is None:
            self.filter_cells_params = {}
        else:
            self.filter_cells_params = filter_cells_params

        if filter_genes_params is None:
            self.filter_genes_params = {}
        else:
            self.filter_genes_params = filter_genes_params

        self.bulk_fields = bulk_fields

        if dropna is None:
            self.dropna = []
        else:
            self.dropna = dropna

        if drop_cols is None:
            self.drop_cols = []
        else:
            self.drop_cols = drop_cols

        if select_cols is None:
            self.select_cols = []
        else:
            self.select_cols = select_cols

        self.padata = None

    def run_pseudobulking(self, ) -> None:
        '''
        Run pseudobulking pipeline.
        '''
        logger.info("Read anndata from %s", self.file_input)
        adata = sc.read_h5ad(self.file_input)
        logger.info("Process anndata")
        adata = self.preprocess(adata)
        logger.info("Build pseudobulk")
        self.padata = self.build_pseudoulk(adata)
        logger.info("Filter columns")
        self.filter_cols()

    def preprocess(self, adata: AnnData) -> AnnData:
        '''
        A function to preprocess a single cell dataset:
        - remove observations containing nan values in the
          specified columns
        - filter cell/gene outliers

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        adata = self.dropna_obs(adata)
        self.filter_cells_(adata)
        self.filter_genes_(adata)
        if adata.obs.empty:
            logger.warning("The AnnData object is empty after processing!")
        return adata

    def build_pseudoulk(self, adata: AnnData) -> AnnData:
        '''
        A function to construct a pseudobulk dataset from
        scRNA-seq data by applying the Decoupler package; 
        and to filter samples with low number of cells 
        from the obtained dataset.

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        self.define_sample_group_cols(adata)
        padata = dc.pp.pseudobulk(adata,
                                  sample_col='sample_id',
                                  groups_col='pseudo_group',
                                  layer='counts')
        self.set_sample_idx(padata)
        del padata.layers["psbulk_props"]
        return padata


    def filter_cols(self) -> None:
        '''
        Drop or select columns in adata.obs.
        '''
        if not self.padata is None:
            if len(self.drop_cols) > 0:
                self.padata.obs.drop(columns=self.drop_cols, inplace=True)
            if len(self.select_cols) > 0:
                self.padata.obs = self.padata.obs[
                                self.select_cols + ['psbulk_cells', 'psbulk_counts']
                                ]
        else:
            raise Exception("The pseudobulk dataset is empty")

    def dropna_obs(self, adata: AnnData) -> AnnData:
        '''
        Drop rows (observations) which have NaN values in the selected columns
        stored in self.dropna (before pseudobulking).

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        # TODO: look at counts field
        adata.layers['counts'] = adata.X.copy()
        if len(self.dropna) != 0:
            adata = adata[~adata.obs[self.dropna].isna().all(axis=1)].copy()
        return adata

    def filter_cells_(self, adata: AnnData) -> None:
        '''
        A function to filter cell outliers by specified params
        (min_counts, min_genes, max_counts, max_genes).

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        for key in self.filter_cells_params.keys():
            kwarg = dict({key: self.filter_cells_params[key]})
            sc.pp.filter_cells(adata, **kwarg, inplace=True)

    def filter_genes_(self, adata: AnnData) -> None:
        '''
        A function to filter gene outliers by specified params
        (min_counts, min_cells, max_counts, max_cells).

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        for key in self.filter_genes_params.keys():
            kwarg = dict({key: self.filter_genes_params[key]})
            sc.pp.filter_genes(adata, **kwarg, inplace=True)

    def define_sample_group_cols(self, adata: AnnData) -> None:
        '''
        A function to construct 'sample_id' column for pseudobulking 
        procedure based on the list of chosen columns (self.bulk_fields).

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        for f in self.bulk_fields:
            if f not in adata.obs:
                raise KeyError(f"adata.obs lacks '{f}'")
            if is_float_dtype(adata.obs[f]) and not (adata.obs[f].isnull().values.any()):
                adata.obs[f] = adata.obs[f].astype(int)

        adata.obs['sample_id'] = adata.obs[self.bulk_fields]\
            .astype(str)\
            .apply(lambda row: '_'.join(row.values), axis=1, result_type='reduce')
        adata.obs['pseudo_group'] = 'all'

    def set_sample_idx(self, padata: AnnData) -> None:
        '''
        A function to set index to 'sample_id'.

        Parameters:
        -----------
        padata : AnnData
            A pseudobulk dataset.
        '''
        padata.obs['sample_id'] = padata.obs['sample_id'].astype(str)
        padata.obs.set_index('sample_id', inplace=True)
        padata.obs.drop(columns=['pseudo_group'], inplace=True)

    def save(self, ) -> None:
        '''
        Save the pseudobulk dataset.
        '''
        if not self.padata is None:
            self.padata.write_h5ad(self.file_output, compression="gzip")


    def add_pubchem_cids_to_padata(
                                   self,
                                   cache: Dict[str, Optional[int]],
                                   drug_col: str = 'perturbation'
                                   ) -> None:
        '''
        Add a 'pubchem_cid' column to adata.obs based on the drug_col.
        
        Parameters:
        -----------
        cache : Dict[str, Optional[int]]
            A dictionary storing the mapping of drugs to PubChem CIDs.
        drug_col : str
            A name of column containing drug names.
        '''
        def get_pubchem_cid(
                            drug_name: str,
                            cache: Dict[str, Optional[int]]
                            ) -> Optional[int]:
            '''
            Fetch PubChem CID for a given drug name, using cache to skip repeat lookups.
            
            Parameters:
            -----------
            drug_name : str
                A drug name for mapping to PubChem CID.
            cache : Dict[str, Optional[int]]
                A dictionary storing the mapping of drugs to PubChem CIDs.      
            '''
            if pd.isna(drug_name) or not drug_name:
                return None

            if drug_name in cache:
                return cache[drug_name]
            try:
                compounds = pcp.get_compounds(drug_name, 'name')
                cid = compounds[0].cid if compounds else None
                logger.debug("CID for '%s': %d", drug_name, cid)
            except Exception as e:
                logger.warning("PubChem lookup failed for '%s': %d", drug_name, e)
                cid = None
            cache[drug_name] = cid
            return cid

        if not self.padata is None:
            unique_drugs = self.padata.obs[drug_col].dropna().unique().tolist()
            to_fetch = [d for d in unique_drugs if d not in cache]
            if to_fetch:
                logger.info("Looking up %d new perturbations on PubChem...", len(to_fetch))
                for drug in to_fetch:
                    get_pubchem_cid(drug, cache)

            self.padata.obs['pubchem_cid'] = self.padata.obs[drug_col].map(cache)
        else:
            raise Exception("The pseudobulk dataset is empty")



def main():
    '''
    A function to process the arguments entered in the console/
    read from the config file and to run the downstream pipeline
    with the set parameters.
    '''
    def merge_args(
                   d_args: Dict[str, Optional[str | Dict[str, int] | List[str]]],
                   config: Dict[str, Optional[str  |Dict[str, int] | List[str]]]
                   ) -> Dict[str, Optional[str | Dict[str, int] | List[str]]]:
        '''
        A function to unite the parameters entered as the arguments
        from the console and the parameters loaded from a config file.

        Parameters:
        -----------
        d_args : Dict[str, Optional[str | Dict[str, int] | List[str]]]
            input arguments represented as a dictionary.
        config : Dict[str, Optional[str  |Dict[str, int] | List[str]]]
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
    parser.add_argument('--filter_cells_params', nargs='+', action=ParseKW)
    parser.add_argument('--filter_genes_params', nargs='+', action=ParseKW)
    parser.add_argument('--bulk_fields', nargs='+')
    parser.add_argument('--dropna', nargs='+', default=[])
    parser.add_argument('--drop_cols', nargs='+', default=[])
    parser.add_argument('--select_cols', nargs='+', default=[])
    parser.add_argument('--drug_col', type=str)
    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']

    if not args.config is None:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)

    if 'filter_cells_params' in d_args.keys():
        if d_args['filter_cells_params']:
            pars_appr = ['min_counts', 'min_genes', 'max_counts', 'max_genes']
            for par in d_args['filter_cells_params'].keys():
                if not par in pars_appr:
                    raise KeyError(
                        f"The parameter's name {par} isn't appropriate, "\
                        f"should be one of {pars_appr}"
                        )

    if 'filter_genes_params' in d_args.keys():
        if d_args['filter_genes_params']:
            pars_appr = ['min_counts', 'min_cells', 'max_counts', 'max_cells']
            for par in d_args['filter_genes_params'].keys():
                if not par in pars_appr:
                    raise KeyError(
                        f"The parameter's name {par} isn't appropriate, "\
                        f"should be one of {pars_appr}"
                        )


    for key in ['input', 'output', 'bulk_fields', 'drug_col']:
        if not d_args[key]:
            raise KeyError(f"The argument {key} is not set")


    pseudo = Pseudobulk(
                        file_input = d_args['input'],
                        file_output = d_args['output'],
                        bulk_fields = d_args['bulk_fields'],
                        filter_cells_params = d_args['filter_cells_params'],
                        filter_genes_params = d_args['filter_genes_params'],
                        dropna = d_args['dropna'],
                        drop_cols = d_args['drop_cols'],
                        select_cols = d_args['select_cols'],
                        )
    pseudo.run_pseudobulking()
    logger.info("Mapping drugs to PubChem idx")
    pseudo.add_pubchem_cids_to_padata({}, drug_col=d_args['drug_col'])
    logger.info("Save dataset")
    pseudo.save()




if __name__ == "__main__":
    main()
