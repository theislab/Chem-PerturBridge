import time
import os
import json
import argparse
from typing import Tuple, List, Dict, Optional, Any, Callable
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import decoupler as dc
import pubchempy as pcp
from pandas.api.types import is_float_dtype, is_categorical_dtype
from anndata import AnnData
from pandas import DataFrame, Series
from scipy.sparse import csr_matrix

from src.utils.parsing_utils import *
from .standardization import *
from .pubchem_imputation import *

class Pseudobulk:
    '''
    A class to process raw scRNA-seq AnnData and
    aggregate counts to pseudobulk.

    Parameters:
    -----------
    dataset_name : str
        A name of the dataset to process
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
    filter_cells_params : Optional[Dict[str, int]]
        A dictionary containing the params for cell filtering.
    filter_genes_params : Optional[Dict[str, int]]
        A dictionary containing the params for gene filtering.
    ignore_cell_lines : List[str]
        A list of cellosaurus ids to ignore during MALAT1
        filtering. (To exclude cells without a nucleus from
        filtration procedure)

    Attributes:
    -----------
    dataset_name : str
        A name of the dataset to process
    files_input : Dict[str, Optional[str]]
        Paths to the input files containing raw data (AnnData as .h5ad, 
        .obs and .var as .parquet).
    file_output : str
        A path to the output file containing pseudobulk data.
    groupby_fields : List[str]
        A list of column names to construct sample_id based on them.
    standardize_obs : Callable
        A dataset-specific function for processing .obs dataframe.
    standardize_var : Callable
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
    padata : AnnData
        A pseudobulk dataset.
    '''
    def __init__(
                 self,
                 dataset_name: str,
                 files_input: Dict[str, Optional[str]],
                 file_output: str,
                 groupby_fields: List[str],
                 dataset_standardization_: Dict[str, Dict[str, Callable]],
                 sm2pubchem_: Optional[Dict[str, int]] = None,
                 filter_malat1: bool = False,
                 filter_low_counts: bool = False,
                 filter_nans: bool = False, 
                 filter_cells_params: Optional[Dict[str, int]] = None,
                 filter_genes_params: Optional[Dict[str, int]] = None,
                 ignore_cell_lines: List[str] = [],
                 standard_dose_units = 'uM',
                 standard_time_units = 'h',
                 ):

        self.dataset_name = dataset_name
        self.files_input = files_input

        dir_output = os.path.dirname(file_output)
        if not os.path.exists(dir_output):
            try:
                os.makedirs(dir_output)
            except FileExistsError as e:
                logger.warning('%s', str(e))
        self.file_output = file_output

        self.groupby_fields = groupby_fields
        self.standardize_obs = dataset_standardization_['standardize_obs']
        self.standardize_var = dataset_standardization_['standardize_var']
        self.sm2pubchem = sm2pubchem_

        self.filter_malat1 = filter_malat1
        self.filter_low_counts = filter_low_counts
        self.filter_nans = filter_nans

        if filter_cells_params is None:
            self.filter_cells_params = {}
        else:
            self.filter_cells_params = filter_cells_params

        if filter_genes_params is None:
            self.filter_genes_params = {}
        else:
            self.filter_genes_params = filter_genes_params

        self.ignore_cell_lines = ignore_cell_lines
        self.standard_dose_units = standard_dose_units
        self.standard_time_units = standard_time_units

        self.pbulk_columns = [
                'plate', 'well', 'cell_type', 'perturbagen', 
                'pert_type', 'is_control', f'pert_dose_{self.standard_dose_units}', 
                f'pert_time_{self.standard_time_units}', 'suspension_type', 'tissue', 
                'tissue_type', 'disease', 'library', 'stimulation', 'guide', 
                'dataset', 'assay', 'development_stage', 'organism', 'sex', 
                'self_reported_ethnicity', 'psbulk_cells', 'psbulk_counts', 'pubchem_cid']
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
        logger.info('Standardize anndata')
        adata = self.standardize(adata)
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
            logger.info("Read data from %s", self.files_input.get('path2adata'))
            adata = sc.read_h5ad(self.files_input.get('path2adata'))
        else:
            raise Exception('No path to the adata file')

        if self.files_input.get('path2obs'):
            logger.info("Read data from %s", self.files_input.get('path2obs'))
            obs = pd.read_parquet(self.files_input.get('path2obs'))
        else:
            obs = None

        if self.files_input.get('path2var'):
            logger.info("Read data from %s", self.files_input.get('path2var'))
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
    
    def standardize(self, adata: AnnData) -> AnnData:
        '''
        Standardize an input AnnData object
        to a specified format:
        - run standardization on obs
        - run standardization on var

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        adata = self.standardize_obs(adata)
        adata = self.standardize_var(adata)
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
        adata = adata[~is_outlier]
        self.filter_cells_(adata)
        self.filter_genes_(adata)
        logger.info(f'Final adata size: {adata.shape}')
        if adata.obs.empty:
            logger.warning("The AnnData object is empty after pre-filtering!")
        return adata

    def chunking(self, adata: AnnData, processing: Callable):
        '''
        A function to chunk the large AnnData objects to avoid
        OOM errors.

        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        processing : Callable
            A processing function which might require extra
            memory.
        '''
        n_chunks = 1
        idx = adata.obs['sample_id'].unique()
        max_n_chunks = len(idx)
        n_idx = len(idx)
        while n_chunks <= max_n_chunks:
            try:
                if n_chunks == 1:
                    return processing(adata)
                else:
                    size = math.ceil(n_idx / n_chunks)
                    results = []
                    for i in range(0, n_idx, size):
                        chunk = adata[adata.obs['sample_id'].isin(idx[i:i+size])]
                        results.append(processing(chunk))
                    return ad.concat(results, merge='same')
            except MemoryError as e:
                n_chunks *= 2


        raise RuntimeError('Processing failed: the maximum chunk size reached')
    
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
        pseudo = lambda x: dc.pp.pseudobulk(x,
                                  sample_col='sample_id',
                                  groups_col='pseudo_group',
                                  )

        padata = self.chunking(adata, pseudo)
        self.set_sample_idx(padata)
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
        logger.info("Process pseudobulk")
        if (self.padata is not None) and (self.padata.shape[0] != 0):
            self.padata.obs[f'pert_dose_{self.standard_dose_units}'] = self.padata.obs\
                    .apply(lambda x: standardize_units(x.pert_dose, self.standard_dose_units), axis=1)
            self.padata.obs[f'pert_time_{self.standard_time_units}'] = self.padata.obs\
                    .apply(lambda x: standardize_units(x.pert_time, self.standard_time_units), axis=1)
            self.padata.obs[[f'pert_dose_{self.standard_dose_units}']] = self.padata.obs[[f'pert_dose_{self.standard_dose_units}']].astype(float)
            self.padata.obs[[f'pert_time_{self.standard_time_units}']] = self.padata.obs[[f'pert_time_{self.standard_time_units}']].astype(float)
            self.padata.obs[['psbulk_cells']] = self.padata.obs[['psbulk_cells']].astype(int)
            self.padata.obs[['psbulk_counts']] = self.padata.obs[['psbulk_counts']].astype(int)
            if len(self.pbulk_columns) > 0:
                self.padata.obs = self.padata.obs[self.pbulk_columns]
            self.padata.X = csr_matrix(self.padata.X.astype(int))
            self.padata.layers['psbulk_props'] = csr_matrix(self.padata.layers['psbulk_props'])
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
            logger.info("Save data to %s", self.file_output)
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
                            cache: Dict[str, Optional[int]],
                            n_retries: int = 5,
                            ) -> Optional[int]:
            '''
            Fetch PubChem CID for a given drug name, using cache to skip repeat lookups.
            
            Parameters:
            -----------
            drug_name : str
                A drug name for mapping to PubChem CID.
            cache : Dict[str, Optional[int]]
                A dictionary storing the mapping of drugs to PubChem CIDs.
            n_retries : int
                The number of retries when connecting to PubChem
                goes wrong.
            '''
            if pd.isna(drug_name) or not drug_name:
                return None

            if drug_name in cache:
                return cache[drug_name]
            
            cnt = 0
            while cnt < n_retries:
                try:
                    compounds = pcp.get_compounds(drug_name, 'name')
                    cid = compounds[0].cid if compounds else None
                    logger.debug("CID for '%s': %d", drug_name, cid)
                    break
                except Exception as e:
                    if (isinstance(e, pcp.PubChemHTTPError) \
                            or isinstance(e, pcp.TimeoutError) \
                            or isinstance(e, pcp.ServerError)) \
                            and ((e.code == 503) or (e.code == 504)):
                        logger.warning("PubChem lookup failed for '%s': %s. Retry %d.", drug_name, str(e), cnt)
                        cnt += 1
                        time.sleep(5)
                    else:
                        logger.warning("PubChem lookup failed for '%s': %s", drug_name, str(e))
                        cid = None
                        break

            cache[drug_name] = cid
            return cid

        logger.info("Mapping drugs to PubChem idx")
        if not self.padata is None:
            unique_drugs = self.padata.obs[drug_col].dropna().unique().tolist()
            to_fetch = [d for d in unique_drugs if d not in cache]
            if to_fetch:
                logger.info("Looking up %d new perturbations on PubChem...", len(to_fetch))
                for drug in to_fetch:
                    get_pubchem_cid(drug, cache)

            self.padata.obs['pubchem_cid'] = self.padata.obs[drug_col].map(cache)
            if self.sm2pubchem:
                self.padata.obs['pubchem_cid'] = self.padata.obs['pubchem_cid']\
                                                     .fillna(self.padata.obs[drug_col].map(self.sm2pubchem))
            self.padata.obs['pubchem_cid'] = self.padata.obs['pubchem_cid'].astype('Int64').astype('category')
        else:
            raise Exception("The pseudobulk dataset is empty")


