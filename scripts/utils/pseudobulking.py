import os
#import sys
import json
#import logging
import argparse
from IPython import get_ipython
from typing import List, Dict, Optional, Sequence, Any, Callable
import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import pubchempy as pcp
from pandas.api.types import is_float_dtype, is_categorical_dtype
from anndata import AnnData

import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro

from helpers import *
from standartization import *

class Pseudobulk:
    '''
    A class to process raw scRNA-seq AnnData and
    aggregate counts to pseudobulk.

    Parameters:
    -----------
    files_input : Dict[str, Optional[str]]
        A path to the input file containing raw data.
    file_output : str
        A path to the output file containing pseudobulk data.
    groupby_fields : List[str]
        A list of columns to construct sample_id based on them.
    filter_cells_params : Optional[Dict[str, int]]
        A dictionary containing the params for cell filtering.
    filter_genes_params : Optional[Dict[str, int]]
        A dictionary containing the params for gene filtering.

    Attributes:
    -----------
    file_input : str
        A path to the input file containing raw data.
    file_output : str
        A path to the output file containing pseudobulk data.
    groupby_fields : List[str]
        A list of columns to construct sample_id based on them.
    filter_cells_params : Optional[Dict[str, int]]
        A dictionary containing the params for cell filtering.
    filter_genes_params : Optional[Dict[str, int]]
        A dictionary containing the params for gene filtering.
    padata : AnnData
        A pseudobulk dataset.
    '''
    def __init__(
                 self,
                 files_input: Dict[str, Optional[str]],
                 file_output: str,
                 filter_malat1: bool,
                 filter_tahoe: bool,
                 filter_nans: bool,
                 filter_singlets: bool,
                 groupby_fields: List[str],
                 dataset_standardization_: Dict[str, Dict[str, Callable]],
                 filter_cells_params: Optional[Dict[str, int]] = None,
                 filter_genes_params: Optional[Dict[str, int]] = None,
                 nworkers: int = 4,
                 ):

        self.files_input = files_input

        dir_output = os.path.dirname(file_output)
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)
        self.file_output = file_output

        self.filter_malat1 = filter_malat1
        self.filter_tahoe = filter_tahoe
        self.filter_nans = filter_nans
        self.filter_singlets = filter_singlets

        if filter_cells_params is None:
            self.filter_cells_params = {}
        else:
            self.filter_cells_params = filter_cells_params

        if filter_genes_params is None:
            self.filter_genes_params = {}
        else:
            self.filter_genes_params = filter_genes_params

        self.groupby_fields = groupby_fields
        self.standartize_obs = dataset_standardization_['standardize_obs']
        self.standartize_var = dataset_standardization_['standardize_var']
        self.pbulk_columns = ['plate', 'well', 'pert_time', 'cell_type', 'perturbagen', 'pert_type',
                              'pert_dose', 'tissue_type', 'tissue', 'disease', 'is_control',
                              'library', 'stimulation', 'guide', 'dataset', 'psbulk_cells',
                              'psbulk_counts', 'pubchem_cid'
                              ]
        self.nworkers = nworkers
        self.padata = None
        
    def determine_groupby_fields(self, adata):
        groupby_fields = np.array(self.groupby_fields)
        ignore_fields = groupby_fields[adata.obs[groupby_fields].isna().all()]
        filterd_fields = [item.item() for item in groupby_fields if item not in ignore_fields]
        self.groupby_fields = filterd_fields

    def run_pseudobulking(self, ) -> None:
        '''
        Run pseudobulking pipeline.
        '''
        adata = self.get_adata()
        logger.info('Standartize anndata')
        adata = self.standartize(adata)
        logger.info('Determine groupby fields')
        self.determine_groupby_fields(adata)
        adata = self.prefilter(adata)
        logger.info("Build pseudobulk")
        self.padata = self.build_pseudoulk(adata)

    def get_adata(self, ):
        logger.info("Read data")
        adata, obs, var = self.read_data()
        logger.info("Merge datasets")
        adata = self.merge_data(adata.copy(), obs, var)
        return adata

    def read_data(self, ):
        if self.files_input.get('path2adata'):
            adata = sc.read_h5ad(self.files_input.get('path2adata'))
        else:
            raise Exception('No path to the adata file')

        if self.files_input.get('path2obs'):
            obs = pd.read_parquet(self.files_input.get('path2obs'))
        else:
            obs = None

        if self.files_input.get('path2var'):
            var = pd.read_parquet(self.files_input.get('path2var'))
        else:
            var = None

        return adata, obs, var

    def merge_data(self, adata, obs, var):
        
        if (obs is not None) and (not adata.obs.equals(obs)):
            if adata.obs.index.name is None:
                adata.obs.index.name = 'index'
            if obs.index.name is None:
                obs.index.name = 'index'
            col_index = adata.obs.index.name
            if (adata.obs.index.name != obs.index.name):
                obs.set_index(col_index, inplace=True)
            if len(adata.obs.index.intersection(obs.index)) == 0:
                    raise Exception('The indices in andata.obs and obs are different!')
            adata.obs = adata.obs.merge(obs, how='left', on=col_index, suffixes=('_old', ''))
    
        if (var is not None) and (not adata.var.equals(var)):
            if adata.var.index.name is None:
                adata.var.index.name = 'index'
            if var.index.name is None:
                var.index.name = 'index'
            col_index = adata.var.index.name
            if adata.var.index.name != var.index.name:
                var.set_index(col_index, inplace=True)
            if len(adata.var.index.intersection(var.index)) == 0:
                    raise Exception('The indices in andata.var and var are different!')
            adata.var = adata.var.merge(var, how='left', on=col_index, suffixes=('_old', ''))
        return adata
    
    def standartize(self, adata):
        adata = self.standartize_obs(adata.copy())
        adata = self.standartize_var(adata.copy())
        adata.layers['counts'] = adata.X.copy()
        return adata


    def prefilter(self, adata: AnnData) -> AnnData:
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
        is_outlier = None
        if self.filter_tahoe:
            logger.info('Filter by counts')
            is_outlier = self.filter_by_counts(adata)
        if self.filter_nans:
            logger.info('Filter by nans')
            is_outlier = is_outlier | self.filter_by_nans(adata)
        if self.filter_malat1:
            logger.info('Filter by humanMALAT1')
            is_outlier = is_outlier | self.filter_by_humanMALAT1(adata)
        if self.filter_singlets:
            logger.info('Filter by Singlets')
            is_outlier = is_outlier | self.filter_by_singlets(adata, nworkers=self.nworkers)
        if is_outlier is not None:
            adata = adata[~is_outlier].copy()
        self.filter_cells_(adata)
        self.filter_genes_(adata)
        if adata.obs.empty:
            logger.warning("The AnnData object is empty after pre-filtering!")
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

    
    def process_pseudobulk(self) -> None:
        '''
        Select columns in padata.obs.
        '''
        if not self.padata is None:
            if len(self.pbulk_columns) > 0:
                self.padata.obs = self.padata.obs[self.pbulk_columns].copy()
        else:
            raise Exception("The pseudobulk dataset is None")

    
    def filter_by_nans(self, adata):
        is_outlier = adata.obs[self.groupby_fields].isna().any(axis=1)
        return is_outlier

    def filter_by_counts(self, adata, min_ngenes=250, min_ncounts=700, max_pcnt_mito=0.2):
        is_outlier = ~((adata.obs['ngenes'] >= min_ngenes)\
                   & (adata.obs['ncounts'] >= min_ncounts)\
                   & (adata.obs['pcnt_mito'] <= max_pcnt_mito))
        return is_outlier
    
    def filter_by_humanMALAT1(self, adata, ens_id='ENSG00000251562', scaling=10000, threshold=3.5):
        is_malat1 = adata.var_names.str.startswith(ens_id)
        fraction_counts_malat1 = adata[:, is_malat1].X.toarray().sum(1)/adata.obs['ncounts'].values
        norm_malat1 = np.log1p(fraction_counts_malat1 * scaling)
        is_outlier = norm_malat1 < threshold
        return is_outlier

    def filter_by_singlets(self, adata, nworkers=4, seed=42, singlet=1):
        X = adata.X.T
        r_script = f'''
        library(scDblFinder)
        library(BiocParallel)
        
        
        set.seed({seed})
        
        sce = scDblFinder(
            SingleCellExperiment(
                list(counts=X),
            ),
            BPPARAM = MulticoreParam(workers = {nworkers})
        )
        
        
        doublet_score = sce$scDblFinder.score
        doublet_class = sce$scDblFinder.class
        '''
        with localconverter(ro.default_converter + anndata2ri.converter):
            X_r = ro.conversion.py2rpy(X)
        ro.globalenv['X'] = X_r
        ro.r(r_script)

        with localconverter(ro.default_converter + anndata2ri.converter):
            doublet_class = ro.conversion.rpy2py(ro.globalenv['doublet_class'])
        is_outlier = (doublet_class != 'singlet')
        return is_outlier

    def filter_cells_(self, adata: AnnData) -> None:
        '''
        A function to filter cell outliers by specified params
        (min_counts, min_genes, max_counts, max_genes).

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        if len(self.filter_cells_params.keys()) != 0:
            logger.info('Filter cells')
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
        if len(self.filter_genes_params.keys()) != 0:
            logger.info('Filter genes')

        for key in self.filter_genes_params.keys():
            kwarg = dict({key: self.filter_genes_params[key]})
            sc.pp.filter_genes(adata, **kwarg, inplace=True)

    def define_sample_group_cols(self, adata: AnnData) -> None:
        '''
        A function to construct 'sample_id' column for pseudobulking 
        procedure based on the list of chosen columns (self.groupby_fields).

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        def f_to_numeric(df_s):
            try:
                pd.to_numeric(df_s, downcast='float')
            except Exception:
                return df_s

        for f in self.groupby_fields:
            if f not in adata.obs:
                raise KeyError(f"adata.obs lacks '{f}'")
            if not (adata.obs[f].isnull().values.any()):
                if is_float_dtype(f_to_numeric(adata.obs[f])):
                    adata.obs[f] = adata.obs[f].astype(int)

        adata.obs['sample_id'] = adata.obs[self.groupby_fields]\
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
                                   drug_col: str = 'perturbagen'
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
    def check_sub_arg_values(d_args, arg_names):
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
    parser.add_argument('--filter_tahoe', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_nans', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_singlets', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_cells_params', nargs='+', action=ParseKW)
    parser.add_argument('--filter_genes_params', nargs='+', action=ParseKW)
    parser.add_argument('--groupby_fields', nargs='+', default=['plate', 'well', 'perturbagen', 'cell_type', 'guide'])
    parser.add_argument('--drug_col', type=str, default='perturbagen')
    parser.add_argument('--nworkers', type=int, default=4)
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
                        filter_malat1=d_args['filter_malat1'],
                        filter_tahoe=d_args['filter_tahoe'],
                        filter_nans=d_args['filter_nans'],
                        filter_singlets=d_args['filter_singlets'],
                        filter_cells_params=d_args['filter_cells_params'],
                        filter_genes_params=d_args['filter_genes_params'],
                        groupby_fields=d_args['groupby_fields'],
                        dataset_standardization_=dataset_standardization[d_args['dataset_name']],
                        nworkers=d_args['nworkers']
                        )
    pseudo.run_pseudobulking()
    logger.info("Mapping drugs to PubChem idx")
    pseudo.add_pubchem_cids_to_padata({}, drug_col=d_args['drug_col'])
    logger.info("Process pseudobulk")
    pseudo.process_pseudobulk()
    logger.info("Save dataset")
    pseudo.save()



if __name__ == "__main__":
    main()
