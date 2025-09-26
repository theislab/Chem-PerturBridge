import time
import os
import math
import re
from typing import Tuple, List, Dict, Optional, Callable
from string import ascii_lowercase
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import decoupler as dc
import pubchempy as pcp
from pandas.api.types import is_float_dtype
from pandas import DataFrame, Series
from anndata import AnnData
from scipy.sparse import csr_matrix
from magnitude import mg, new_mag, Magnitude

from src.utils.parsing_utils import *

# Initialize new units: moles per liter)
new_mag('l', Magnitude(0.001, m=3))
new_mag('M', mg(1, 'mol/l'))

# Define column names for .obs standardization
STANDARD_ADATA_OBS_COLS = [
    'plate',
    'well',
    'ngenes',
    'ncounts',
    'pcnt_mito',
    'cell_type',
    'perturbagen',
    'pert_type',
    'is_control',
    'pert_dose',
    'pert_time',
    'suspension_type',
    'tissue',
    'tissue_type',
    'disease',
    'library',
    'stimulation',
    'guide',
    'dataset',
    'assay',
    'development_stage',
    'organism',
    'sex',
    'self_reported_ethnicity',
    'pubchem_cid',
]

# Define categorical columns
CAT_COLS = [
    'plate',
    'well',
    'cell_type',
    'perturbagen',
    'pert_type',
    'suspension_type',
    'tissue',
    'tissue_type',
    'disease',
    'library',
    'stimulation',
    'guide',
    'dataset',
    'development_stage',
    'organism',
    'sex',
    'self_reported_ethnicity'
]

# Define column names for .var standardization
STANDARD_VAR_COLS = ['symbol']



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
    standardize_dataset : Callable
        A dataset-specific function for processing AnnData dataframe.
    sm2pubchem : Dict[str, int]
        A dictionary for manual mapping of drugs to pubchem cids.
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
    drug_col : Optional[str]
        A name of column containing drugs to map them to
        pubchem cids.
    standard_dose_units : str
        Common units of dose concentration.
    standard_time_units : str
        Common units for time of perturbations.

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
    standardize_dataset : Callable
        A dataset-specific function for processing AnnData dataframe.
    sm2pubchem : Dict[str, int]
        A dictionary for manual mapping of drugs to pubchem cids.
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
    drug_col : Optional[str]
        A name of column containing drugs to map them to
        pubchem cids.
    standard_dose_units : str
        Common units of dose concentration.
    standard_time_units : str
        Common units for time of perturbations.
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
        standardize_dataset: Dict[str, Callable],
        sm2pubchem: Optional[Dict[str, int]] = None,
        filter_malat1: bool = False,
        filter_low_counts: bool = False,
        filter_nans: bool = False,
        filter_cells_params: Optional[Dict[str, int]] = None,
        filter_genes_params: Optional[Dict[str, int]] = None,
        ignore_cell_lines: List[str] = [],
        drug_col: Optional[str] = 'perturbagen',
        standard_dose_units: str = 'uM',
        standard_time_units: str = 'h',
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
        self.standardize_dataset = standardize_dataset
        self.sm2pubchem = sm2pubchem
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
        self.drug_col = drug_col
        self.standard_dose_units = standard_dose_units
        self.standard_time_units = standard_time_units

        self.pbulk_columns = [
            'plate', 'well', 'cell_type', 'perturbagen',
            'pert_type', 'is_control', f'pert_dose_{self.standard_dose_units}',
            f'pert_time_{self.standard_time_units}', 'suspension_type', 'tissue',
            'tissue_type', 'disease', 'library', 'stimulation', 'guide',
            'dataset', 'assay', 'development_stage', 'organism', 'sex',
            'self_reported_ethnicity', 'pubchem_cid', 'psbulk_cells', 'psbulk_counts']
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
        filtered_fields = [item.item()
                          for item in groupby_fields if item not in ignore_fields]
        self.groupby_fields = filtered_fields

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
        logger.info('Build pseudobulk')
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
        logger.info('Merge datasets')
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
            logger.info("Read data from %s",
                        self.files_input.get('path2adata'))
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
            if adata.obs.index.name != obs.index.name:
                obs.set_index(col_index, inplace=True)
            if len(adata.obs.index.intersection(obs.index)) == 0:
                raise Exception(
                    'The indices in andata.obs and obs are different!')
            adata.obs = adata.obs.merge(
                obs, how='left', on=col_index, suffixes=('_old', ''))

        if (var is not None) and (not adata.var.equals(var)):
            if adata.var.index.name is None:
                adata.var.index.name = 'index'
            if var.index.name is None:
                var.index.name = 'index'
            col_index = adata.var.index.name
            if adata.var.index.name != var.index.name:
                var.set_index(col_index, inplace=True)
            if len(adata.var.index.intersection(var.index)) == 0:
                raise Exception(
                    'The indices in andata.var and var are different!')
            adata.var = adata.var.merge(
                var, how='left', on=col_index, suffixes=('_old', ''))
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
        adata = self.standardize_dataset(adata)
        if not self.drug_col is None:
            adata = self.add_pubchem_cids(adata, 
                                          cache={}, 
                                          drug_col=self.drug_col)
        adata = self.standardize_common(adata)
        return adata

    def standardize_obs_common(self, adata: AnnData):
        '''
        Convert obs from AnnData object to the standard format.
        Select the pre-defined columns.
    
        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        categories = {c: 'category' for c in CAT_COLS}
        adata.obs = adata.obs[STANDARD_ADATA_OBS_COLS]
        adata.obs = adata.obs.astype(categories)
        return adata

    def standardize_var_common(self, adata: AnnData) -> AnnData:
        '''
        Convert variables from AnnData object to the standard format.
        Select the pre-defined columns.
        
        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        categories = {c: 'category' for c in STANDARD_VAR_COLS}
        adata.var = adata.var[STANDARD_VAR_COLS]
        adata.var = adata.var.astype(categories)
        return adata

    def standardize_common(self, adata: AnnData) -> AnnData:
        '''
        Convert AnnData object to the standard format.
        Select the pre-defined columns.
        
        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        '''
        adata = self.standardize_obs_common(adata)
        adata = self.standardize_var_common(adata)
        return adata

    def standardize_units(self, s: str, standard_units: str = 'uM') -> Optional[float]:
        '''
        Split combined value and unit strings to store them separately in columns.
        Convert dose/time values to the standard units (uM/hours)
        
        Parameters:
        -----------
        s : str
            Value-unit string to process.
        standard_unit : str
            Standard unit
        '''
        
        r = r'[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?'
        try:
            val = float(re.findall(r, s)[0])
        except Exception as e:
            logger.warning('%s', str(e))
            val = None
        try:
            units = re.split(r, s)[-1]
        except Exception as e:
            logger.warning('%s', str(e))
            units = None
        try:
            val_ = mg(val, units).toval(standard_units)
        except Exception as e:
            logger.warning('%s', str(e))
            val_ = None

        return val_

    def log_filtering_stats(self,
                        n_obs_prev: int,
                        n_var_prev: int,
                        n_ct_prev: int,
                        n_pa_prev: int,
                        n_obs: int,
                        n_var: int,
                        n_ct: int,
                        n_pa: int
                       ) -> None:

        '''
        A function to log the statistics on
        the number of cells (n_obs),
        the number of genes (n_var), the number of perturbations
        (n_perturb)

        Parameters:
        -----------
        n_obs_prev : int
            The number of cells before filtering
        n_var_prev : int
            The number of genes before filtering
        n_ct_prev : int
            The number of unique cell types before filtering
        n_pa_prev : int
            The number of unique perturbagens before filtering
        n_obs: int
            The number of cells after filtering
        n_var: int
            The number of genes after filtering
        n_ct: int
            The number of unique cell types after filtering
        n_pa : int
            The number of unique perturbagens after filtering
        '''
        logger.info('    n_obs: %d --> %d', n_obs_prev, n_obs)
        logger.info('    n_var: %d --> %d', n_var_prev, n_var)
        logger.info('    n_ct: %d --> %d', n_ct_prev, n_ct)
        logger.info('    n_pa: %d --> %d', n_pa_prev, n_pa)

    def compute_filtering_stats(self,
                        adata: AnnData,
                        perturbation_cols: List = [
                                            'cell_type',
                                            'perturbagen'
                                            ]
                                    ) -> None:
        '''
        A function to obtain the statistics after each step 
        of filtering: the number of cells (n_obs), 
        the number of genes (n_var), the number of cell types (n_ct), 
        the number of perturbagens (n_pa) 
        by default.
        
        Parameters:
        -----------
        adata : AnnData
            An input scRNA-seq dataset.
        perturbation_cols : List
            Columns to take into account while computing 
            perturbation statistics.
        '''

        n_obs, n_var = adata.shape
        n_ct = adata.obs.value_counts(perturbation_cols[0]).shape[0]
        n_pa = adata.obs.value_counts(perturbation_cols[1]).shape[0]
        return n_obs, n_var, n_ct, n_pa

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

        logger.info('Initial adata size: %s', str(adata.shape))
        is_outlier = pd.Series(
            np.zeros(adata.obs.shape[0], dtype=bool), index=adata.obs.index)
        if self.filter_low_counts:
            logger.info('Filter by counts')
            n_obs_prev, n_var_prev, n_ct_prev, n_pa_prev = self.compute_filtering_stats(adata[~is_outlier])
            is_outlier = self.filter_by_counts(adata)
            n_obs, n_var, n_ct, n_pa = self.compute_filtering_stats(adata[~is_outlier])
            self.log_filtering_stats(
                        n_obs_prev, 
                        n_var_prev,
                        n_ct_prev,
                        n_pa_prev,
                        n_obs,
                        n_var,
                        n_ct,
                        n_pa)
        if self.filter_nans:
            logger.info('Filter by nans')
            n_obs_prev, n_var_prev, n_ct_prev, n_pa_prev = self.compute_filtering_stats(adata[~is_outlier])
            is_outlier = is_outlier | self.filter_by_nans(adata)
            n_obs, n_var, n_ct, n_pa = self.compute_filtering_stats(adata[~is_outlier])
            self.log_filtering_stats(
                        n_obs_prev,
                        n_var_prev,
                        n_ct_prev,
                        n_pa_prev,
                        n_obs,
                        n_var,
                        n_ct,
                        n_pa)
        if self.filter_malat1:
            logger.info('Filter by humanMALAT1')
            n_obs_prev, n_var_prev, n_ct_prev, n_pa_prev = self.compute_filtering_stats(adata[~is_outlier])
            is_outlier = is_outlier | self.filter_by_humanMALAT1(adata)
            n_obs, n_var, n_ct, n_pa = self.compute_filtering_stats(adata[~is_outlier])
            self.log_filtering_stats(
                        n_obs_prev,
                        n_var_prev,
                        n_ct_prev,
                        n_pa_prev,
                        n_obs,
                        n_var,
                        n_ct,
                        n_pa)
        adata = adata[~is_outlier]
        self.filter_cells_(adata)
        self.filter_genes_(adata)
        logger.info('Final adata size: %s', str(adata.shape))
        if adata.obs.empty:
            logger.warning('The AnnData object is empty after pre-filtering!')
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
                        chunk = adata[adata.obs['sample_id'].isin(
                            idx[i:i + size])]
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

        def pseudo(x): return dc.pp.pseudobulk(x,
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
                .apply(lambda x: self.standardize_units(x.pert_dose, self.standard_dose_units), axis=1)
            self.padata.obs[f'pert_time_{self.standard_time_units}'] = self.padata.obs\
                .apply(lambda x: self.standardize_units(x.pert_time, self.standard_time_units), axis=1)
            self.padata.obs[[f'pert_dose_{self.standard_dose_units}']] = self.padata.obs[[
                f'pert_dose_{self.standard_dose_units}']].astype(float)
            self.padata.obs[[f'pert_time_{self.standard_time_units}']] = self.padata.obs[[
                f'pert_time_{self.standard_time_units}']].astype(float)
            self.padata.obs[['psbulk_cells']
                            ] = self.padata.obs[['psbulk_cells']].astype(int)
            self.padata.obs[['psbulk_counts']
                            ] = self.padata.obs[['psbulk_counts']].astype(int)
            if len(self.pbulk_columns) > 0:
                self.padata.obs = self.padata.obs[self.pbulk_columns]
            self.padata.X = csr_matrix(self.padata.X.astype(int))
            self.padata.layers['psbulk_props'] = csr_matrix(
                self.padata.layers['psbulk_props'])
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
        is_outlier = ~((adata.obs['ngenes'] >= min_ngenes)
                       & (adata.obs['ncounts'] >= min_ncounts)
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
            fraction_counts_malat1 = adata[:, is_malat1].X.toarray().sum(
                1) / adata.obs['ncounts'].values
            norm_malat1 = np.log1p(fraction_counts_malat1 * scaling)
            is_outlier = (norm_malat1 < threshold) & (
                ~adata.obs['cell_type'].isin(self.ignore_cell_lines))
        else:
            is_outlier = pd.Series(
                np.zeros(adata.obs.shape[0], dtype=bool), index=adata.obs.index)
            logger.warning(
                'There is no information about MALAT1. There no filtration occured.')
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
            n_obs_prev, n_var_prev, n_ct_prev, n_pa_prev = self.compute_filtering_stats(adata)
            for key in self.filter_cells_params.keys():
                kwarg = dict({key: self.filter_cells_params[key]})
                sc.pp.filter_cells(adata, **kwarg, inplace=True)
            n_obs, n_var, n_ct, n_pa = self.compute_filtering_stats(adata)
            self.log_filtering_stats(
                        n_obs_prev, 
                        n_var_prev,
                        n_ct_prev,
                        n_pa_prev,
                        n_obs,
                        n_var,
                        n_ct,
                        n_pa)

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
            n_obs_prev, n_var_prev, n_ct_prev, n_pa_prev = self.compute_filtering_stats(adata)
            for key in self.filter_genes_params.keys():
                kwarg = dict({key: self.filter_genes_params[key]})
                sc.pp.filter_genes(adata, **kwarg, inplace=True)
            n_obs, n_var, n_ct, n_pa = self.compute_filtering_stats(adata)
            self.log_filtering_stats(
                        n_obs_prev,
                        n_var_prev,
                        n_ct_prev,
                        n_pa_prev,
                        n_obs,
                        n_var,
                        n_ct,
                        n_pa)

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
                return pd.to_numeric(df_s, downcast='float')
            except Exception:
                return df_s

        for f in self.groupby_fields:
            if f not in adata.obs:
                raise KeyError(f"adata.obs lacks '{f}'")
            if not adata.obs[f].isnull().values.any():
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
            logger.info('Save data to %s', self.file_output)
            self.padata.write_h5ad(self.file_output, compression="gzip")
        else:
            raise Exception("The pseudobulk dataset is None or empty")

    def add_pubchem_cids(
        self,
        adata,
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
            cid = None
            while cnt < n_retries:
                try:
                    compounds = pcp.get_compounds(drug_name, 'name')
                    cid = compounds[0].cid if compounds else None
                    logger.debug("CID for '%s': %d", drug_name, cid)
                    break
                except Exception as e:
                    if (isinstance(e, pcp.PubChemHTTPError)
                            or isinstance(e, pcp.TimeoutError)
                            or isinstance(e, pcp.ServerError)
                            or isinstance(e, pcp.ServerBusyError)):
                        logger.warning(
                            "PubChem lookup failed for '%s': %s. Retry %d.", drug_name, str(e), cnt)
                        cnt += 1
                        time.sleep(5)
                    else:
                        logger.warning(
                            "PubChem lookup failed for '%s': %s", drug_name, str(e))
                        break

            cache[drug_name] = cid
            return cid

        logger.info("Mapping drugs to PubChem idx")
        if not adata is None:
            unique_drugs = adata.obs[drug_col].dropna().unique().tolist()
            to_fetch = [d for d in unique_drugs if d not in cache]
            if to_fetch:
                logger.info(
                    "Looking up %d new perturbations on PubChem...", len(to_fetch))
                for drug in to_fetch:
                    get_pubchem_cid(drug, cache)

            adata.obs['pubchem_cid'] = adata.obs[drug_col].map(
                cache)
            if self.sm2pubchem:
                adata.obs['pubchem_cid'] = adata.obs['pubchem_cid']\
                                                     .fillna(adata.obs[drug_col].map(self.sm2pubchem))
            adata.obs['pubchem_cid'] = adata.obs['pubchem_cid'].astype(
                'Int64').astype('category')
            return adata
        else:
            raise Exception("The pseudobulk dataset is empty")
