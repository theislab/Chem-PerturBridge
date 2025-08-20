import os
import json
import argparse
from IPython import get_ipython
from typing import Tuple, List, Dict, Optional, Any, Callable
import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import pubchempy as pcp
from pandas.api.types import is_float_dtype, is_categorical_dtype
from anndata import AnnData
from pandas import DataFrame, Series

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
        Paths to the input files containing raw data.
    file_output : str
        A path to the output file containing pseudobulk data.
    groupby_fields : List[str]
        A list of column names to construct sample_id based on them.
    dataset_standardization_ : Dict[str, Callable]
        A dictionary containing the dataset-specific functions
        for its processing/standardization.
    filter_malat1 : bool
        A flag, True value means to filter observations by expression 
        of MALAT1
    filter_low_counts : bool
        A flag, True value means to filter observations which have
        low values of counts
    filter_nans : bool
        A flag, True value means to filter observations which have
        missed information
    filter_singlets : bool
        A flag, True value means to find doublets in the observations
        and filter them, (currently is not available)
    filter_cells_params : Optional[Dict[str, int]]
        A dictionary containing the params for cell filtering.
    filter_genes_params : Optional[Dict[str, int]]
        A dictionary containing the params for gene filtering.
    ignore_cell_lines : List[str]
        A list of cellosaurus ids to ignore during MALAT1
        filtering. (To exclude cells without a nucleus from
        filtration procedure)
    nworkers : int
        A number of workers to parallel doublets filtering

    Attributes:
    -----------
    files_input : Dict[str, Optional[str]]
        Paths to the input files containing raw data (AnnData as .h5ad, 
        .obs and .var as .parquet).
    file_output : str
        A path to the output file containing pseudobulk data.
    groupby_fields : List[str]
        A list of column names to construct sample_id based on them.
    standartize_obs : Callable
        A dataset-specific function for processing .obs dataframe.
    standartize_var : Callable
        A dataset-specific function for processing .var dataframe.
    filter_malat1 : bool
        A flag, True value means to filter observations by expression 
        of MALAT1
    filter_low_counts : bool
        A flag, True value means to filter observations which have
        low values of counts
    filter_nans : bool
        A flag, True value means to filter observations which have
        missed information
    filter_singlets : bool
        A flag, True value means to find doublets in the observations
        and filter them, (currently is not available)
    filter_cells_params : Optional[Dict[str, int]]
        A dictionary containing the params for cell filtering.
    filter_genes_params : Optional[Dict[str, int]]
        A dictionary containing the params for gene filtering.
    ignore_cell_lines : List[str]
        A list of cellosaurus ids to ignore during MALAT1
        filtering. (To exclude cells without a nucleus from
        filtration procedure)
    pbulk_columns : List[str]
        A list of column names which the final pseudobulk dataset
        should contain.
    nworkers : int
        A number of workers to parallel doublets filtering
    padata : AnnData
        A pseudobulk dataset.
    '''
    def __init__(
                 self,
                 files_input: Dict[str, Optional[str]],
                 file_output: str,
                 groupby_fields: List[str],
                 dataset_standardization_: Dict[str, Dict[str, Callable]],
                 filter_malat1: bool = False,
                 filter_low_counts: bool = False,
                 filter_nans: bool = False, 
                 filter_singlets: bool = False,
                 filter_cells_params: Optional[Dict[str, int]] = None,
                 filter_genes_params: Optional[Dict[str, int]] = None,
                 ignore_cell_lines: List[str] = [],
                 nworkers: int = 4,
                 ):

        self.files_input = files_input

        dir_output = os.path.dirname(file_output)
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)
        self.file_output = file_output

        self.groupby_fields = groupby_fields
        self.standartize_obs = dataset_standardization_['standardize_obs']
        self.standartize_var = dataset_standardization_['standardize_var']

        self.filter_malat1 = filter_malat1
        self.filter_low_counts = filter_low_counts
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

        self.ignore_cell_lines = ignore_cell_lines

        self.pbulk_columns = [
                'plate', 'well', 'cell_type', 'perturbagen', 
                'pert_type', 'is_control', 'pert_dose', 'pert_dose_unit', 
                'pert_time', 'pert_time_unit', 'suspension_type', 'tissue', 
                'tissue_type', 'disease', 'library', 'stimulation', 'guide', 
                'dataset', 'assay', 'development_stage', 'organism', 'sex', 
                'self_reported_ethnicity', 'psbulk_cells', 'psbulk_counts', 'pubchem_cid']
        self.nworkers = nworkers
        self.padata = None
        
    def determine_groupby_fields(self, adata: AnnData) -> None:
        '''
        A function to determine the groupby column names
        to use them for pseudobulking.

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
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
        self.padata = self.build_pseudoulk(adata.copy())

    def get_adata(self, ) -> AnnData:
        '''
        A function to
        - load the datasets
        - integrate the separately stored .obs and .var into 
          the AnnData object.Be aware that the .obs and .var 
          contained within the AnnData object may differ from 
          the updated versions saved externally.
        - return AnnData object.
        '''
        logger.info("Read data")
        adata, obs, var = self.read_data()
        logger.info("Merge datasets")
        adata = self.merge_data(adata, obs, var)
        return adata

    def read_data(self, ) -> Tuple[AnnData, 
                                   DataFrame, 
                                   DataFrame]:
        '''
        A function to read datasets:
        - the AnnData object stored as h5ad
        - observations and variables stored as
          parquet files.
        '''
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

    def merge_data(self, 
                   adata: AnnData, 
                   obs: DataFrame, 
                   var: DataFrame):
        '''
        A function to merge AnnData object
        with observations  and
        variables DataFrames

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        obs : DataFrame
            A curated observation dataset.
        var : DataFrame
            A curated variable dataset.
        '''
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
    
    def standartize(self, adata: AnnData) -> AnnData:
        '''
        Standartize an input AnnData object
        to a specified format:
        - run standardization on obs
        - run standardization on var

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        adata = self.standartize_obs(adata)
        adata = self.standartize_var(adata)
        return adata


    def prefilter(self, adata: AnnData) -> AnnData:
        '''
        A function to preprocess a single cell dataset:
        - if filter_low_counts True, remove observations 
          based on the counts criteria
        - if filter_nans True, remove observations 
          containing nan values in the specified columns
        - if filter_malat1 True, remove observations
          which have low counts of MALAT1.
        - if filter_singlets True, remove observations
          which are classified as doublets
        - filter cell/gene outliers
        - return the filtered dataset

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        logger.info(f'Initial adata size: {adata.shape}')
        is_outlier = pd.Series(np.zeros(adata.obs.shape[0], dtype=bool), index=adata.obs.index)
        if self.filter_low_counts:
            logger.info('Filter by counts')
            n_obs_prev = adata[~is_outlier].shape[0]
            is_outlier = self.filter_by_counts(adata)
            logger.info(f'n_obs: {n_obs_prev} --> {adata[~is_outlier].shape[0]}')
        if self.filter_nans:
            logger.info('Filter by nans')
            n_obs_prev = adata[~is_outlier].shape[0]
            is_outlier = is_outlier | self.filter_by_nans(adata)
            logger.info(f'n_obs: {n_obs_prev} --> {adata[~is_outlier].shape[0]}')
        if self.filter_malat1:
            logger.info('Filter by humanMALAT1')
            n_obs_prev = adata[~is_outlier].shape[0]
            is_outlier = is_outlier | self.filter_by_humanMALAT1(adata)
            logger.info(f'n_obs: {n_obs_prev} --> {adata[~is_outlier].shape[0]}')
        #if self.filter_singlets:
        #    logger.info('Filter by Singlets')
        #    n_obs_prev = adata[~is_outlier].shape[0]
        #    is_outlier = is_outlier | self.filter_by_singlets(adata, nworkers=self.nworkers)
        #    logger.info(f'n_obs: {n_obs_prev} --> {adata[~is_outlier].shape[0]}')
        adata = adata[~is_outlier]
        self.filter_cells_(adata)
        self.filter_genes_(adata)
        logger.info(f'Final adata size: {adata.shape}')
        if adata.obs.empty:
            logger.warning("The AnnData object is empty after pre-filtering!")
        return adata


    def build_pseudoulk(self, adata: AnnData) -> AnnData:
        '''
        A function
        - to define the sample_id which is used for
          counts aggregation while running pseudobulking.
        - to construct a pseudobulk dataset from
          scRNA-seq data by applying the Decoupler package.

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        self.define_sample_group_cols(adata)
        padata = dc.pp.pseudobulk(adata,
                                  sample_col='sample_id',
                                  groups_col='pseudo_group',
                                  )
        self.set_sample_idx(padata)
        del padata.layers["psbulk_props"]
        return padata

    
    def process_pseudobulk(self) -> None:
        '''
        A function to process a pseudobulk dataset:
        - to separate values in the column containing the information
          about doses and dose units (such as uM).
        - to separate values in the column containing the information
          about time and time units (such as hours).
        - to specify the types for columns
        - select the predetermined columns in padata.obs.
        '''
        if (self.padata is not None) and (self.padata.shape[0] != 0):
            self.padata.obs[['pert_dose', 'pert_dose_unit']] = self.padata.obs\
                    .apply(lambda x: split_val_units(x.pert_dose), axis=1, result_type='expand')
            self.padata.obs[['pert_dose', 'pert_dose_unit']] = self.padata.obs\
                    .apply(lambda x: standardize_dose(x.pert_dose, x.pert_dose_unit), axis=1, result_type='expand')
            self.padata.obs[['pert_time', 'pert_time_unit']] = self.padata.obs\
                    .apply(lambda x: split_val_units(x.pert_time), axis=1, result_type='expand')
            self.padata.obs[['pert_dose']] = self.padata.obs[['pert_dose']].astype(float)
            self.padata.obs[['pert_dose_unit']] = self.padata.obs[['pert_dose_unit']].astype('category')
            self.padata.obs[['pert_time']] = self.padata.obs[['pert_time']].astype(float)
            self.padata.obs[['pert_time_unit']] = self.padata.obs[['pert_time_unit']].astype('category')
            self.padata.obs[['psbulk_cells']] = self.padata.obs[['psbulk_cells']].astype(int)
            self.padata.obs[['psbulk_counts']] = self.padata.obs[['psbulk_counts']].astype(int)
            if len(self.pbulk_columns) > 0:
                self.padata.obs = self.padata.obs[self.pbulk_columns]
        else:
            raise Exception('The pseudobulk dataset is None')

    
    def filter_by_nans(self, adata: AnnData) -> Series:
        '''
        Filter observations with missed values in the groupby_fields.
        
        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        is_outlier = adata.obs[self.groupby_fields].isna().any(axis=1)
        return is_outlier

    def filter_by_counts(self, 
                         adata: AnnData, 
                         min_ngenes: int = 250, 
                         min_ncounts: int = 700, 
                         max_pcnt_mito: float = 0.2) -> Series:
        '''
        Filter observations with low counts and number of genes.
        The filtration parameters are chosen according to the Tahoe's
        paper: https://doi.org/10.1101/2025.02.20.639398
        
        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        min_ngenes : int
            A threshold for min number of genes detected in a observation. 
        min_ncounts : int
            A threshold for min number of counts detected in a obsevation.
        max_pcnt_mito : float
            A threshold for the estimated fraction of mitochondrial reads.
        '''
        is_outlier = ~((adata.obs['ngenes'] >= min_ngenes)\
                   & (adata.obs['ncounts'] >= min_ncounts)\
                   & (adata.obs['pcnt_mito'] <= max_pcnt_mito))
        return is_outlier
    
    def filter_by_humanMALAT1(self, 
                              adata: AnnData, 
                              ens_id: str = 'ENSG00000251562', 
                              scaling: int = 10000, 
                              threshold: float = 3.5) -> Series:
        '''
        Filter observations with low values in MALAT1 expression.
        The reference paper: https://doi.org/10.1186/s12864-024-11015-5
        The scaling and threshold values are taken from the paper.

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        ens_id : str
            MALAT1 ENSG_id
        scaling : int
            A scaling factor
        threshold : float
            A threshold to filter observations based on transformed
            count values
        '''
        if ens_id in adata.var_names:
            is_malat1 = adata.var_names.str.startswith(ens_id)
            fraction_counts_malat1 = adata[:, is_malat1].X.toarray().sum(1)/adata.obs['ncounts'].values
            norm_malat1 = np.log1p(fraction_counts_malat1 * scaling)
            is_outlier = (norm_malat1 < threshold) & (~adata.obs['cell_type'].isin(self.ignore_cell_lines))
        else:
            is_outlier = pd.Series(np.zeros(adata.obs.shape[0], dtype=bool), index=adata.obs.index)
            logger.warning('There is no information about MALAT1. There no filtration occured.')
        return is_outlier

    #def filter_by_singlets(self, 
    #                       adata: AnnData, 
    #                       nworkers: int = 4, 
    #                       seed: int = 42) -> Series:
    #'''
    #Find doublets utilizing expression data by 
    #running R-based script and using the
    #scDblFinder package and to filter them.
    #scDblFinder uses default parameters
    #(such as Expected doublet rate) for datasets 
    #obtained from 10x Genomics experiments.
    #
    #The parameters might vary with the datasets.
    #The decision for applying filter_by_singlets
    #should be made after careful consideration
    #of the experimental methods.
    #
    #The example of the usage of scDblFinder is mentioned here:
    #https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
    #
    #Parameters:
    #-----------
    #adata : AnnData
    #   An input scRNA-seq dataset.
    #nworkers : int
    #   The number of worker to parallel
    #seed : int
    #   Predefined seed to reproduce calculations
    #'''
    #TODO
    #    X = adata.X.T
    #    r_script = f'''
    #    library(scDblFinder)
    #    library(BiocParallel)
    #    set.seed({seed})
    #    sce = scDblFinder(
    #        SingleCellExperiment(
    #            list(counts=X),
    #        ),
    #        BPPARAM = MulticoreParam(workers = {nworkers})
    #    )
    #    
    #    doublet_score = sce$scDblFinder.score
    #    doublet_class = sce$scDblFinder.class
    #    '''
    #    with localconverter(ro.default_converter + anndata2ri.converter):
    #        X_r = ro.conversion.py2rpy(X)
    #    ro.globalenv['X'] = X_r
    #    ro.r(r_script)
    #
    #    with localconverter(ro.default_converter + anndata2ri.converter):
    #        doublet_class = ro.conversion.rpy2py(ro.globalenv['doublet_class'])
    #    is_outlier = (doublet_class != 'singlet')
    #    return is_outlier

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
        procedure based on the list of chosen column names (self.groupby_fields).

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        def f_to_numeric(df_s):
            '''
            Try to convert pandas Series column to the float
            type if it is possible.

            Parameters:
            -----------
            df_s : Series
                An input Series column
            '''
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
        if (not self.padata is None) and (self.padata.shape[0] != 0):
            self.padata.write_h5ad(self.file_output, compression="gzip")
        else:
            raise Exception("The pseudobulk dataset is None or empty")

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
    #parser.add_argument('--filter_singlets', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_cells_params', nargs='+', action=ParseKW)
    parser.add_argument('--filter_genes_params', nargs='+', action=ParseKW)
    parser.add_argument('--groupby_fields', nargs='+', default=['plate', 'well', 'perturbagen', 'cell_type', 'guide'])
    parser.add_argument('--drug_col', type=str, default='perturbagen')
    parser.add_argument('--ignore_cell_lines', nargs='+', default=[])
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
                        groupby_fields=d_args['groupby_fields'],
                        dataset_standardization_=dataset_standardization[d_args['dataset_name']],
                        filter_malat1=d_args['filter_malat1'],
                        filter_low_counts=d_args['filter_low_counts'],
                        filter_nans=d_args['filter_nans'],
                        filter_cells_params=d_args['filter_cells_params'],
                        filter_genes_params=d_args['filter_genes_params'],
                        ignore_cell_lines=d_args['ignore_cell_lines'],
                        nworkers=d_args['nworkers'],
                        #filter_singlets=d_args['filter_singlets'],
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
