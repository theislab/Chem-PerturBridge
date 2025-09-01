from string import ascii_lowercase
import re
from typing import List, Any, Dict, Optional
import scanpy as sc
from anndata import AnnData
from magnitude import mg, new_mag, Magnitude

from src.utils.parsing_utils import *

#Initialize new units: moles per liter)
new_mag('l', Magnitude(0.001, m=3))
new_mag('M', mg(1, 'mol/l'))

#Define column names for .obs standardization
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
                 ]

#Define categorical columns
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

#Define column names for .var standardization
STANDARD_VAR_COLS = ['symbol']

def standardize_units(s: str, standard_units: str = 'uM') -> Optional[float]:
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


def standardize_obs_tahoe(adata: AnnData) -> AnnData:
    '''
    Convert observations from AnnData object to the standard format.
    For Tahoe100 dataset.

    Parameters:
    -----------
    adata : AnnData
       An input scRNA-seq dataset. 
    '''
    rename_cols = {'gene_count': 'ngenes', 
                   'tscp_count': 'ncounts',
                   'cell_line_ontology_id': 'cell_type',
                   'pert_compound': 'perturbagen',
                  }
    categories = {c: 'category' for c in CAT_COLS}
    wellnum2id_mapping = wellnum2id()
    adata.obs['well'] = adata.obs.index.to_series().str.split('_').str[0].astype(int).map(wellnum2id_mapping)
    adata.obs.rename(columns=rename_cols, inplace=True)
    adata.obs['is_control'] = adata.obs['perturbagen'] == 'DMSO_TF'
    adata.obs['library'] = None
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'tahoe100'
    adata.obs = adata.obs[STANDARD_ADATA_OBS_COLS]
    adata.obs = adata.obs.astype(categories)
    return adata


def standardize_var_tahoe(adata: AnnData) -> AnnData:
    '''
    Convert variables from AnnData object to the standard format.
    For Tahoe100 dataset.

    Parameters:
    -----------
    adata : AnnData
       An input scRNA-seq dataset.
    '''
    rename_cols = {'ensembl_gene_id': 'ensembl_id'}
    categories = {c: 'category' for c in STANDARD_VAR_COLS}
    
    adata.var.rename(columns=rename_cols, inplace=True)
    adata.var.set_index('ensembl_id', inplace=True)
    adata.var = adata.var[STANDARD_VAR_COLS]
    adata.var = adata.var.astype(categories)
    return adata

def wellnum2id(nwells: int = 96, 
               ncols: int = 12) -> Dict[int, str]:
    '''
    Convert well numeric ids to letter-based ids.

    Parameters:
    -----------
    nwells : int
        Total number of wells on the experimental plate.
    ncols : int
        The number of columns on the experimental plate.
    '''
    mapping = {}
    for i in range(nwells):
        mapping[i + 1] = ascii_lowercase[i // ncols].upper() + str(i % ncols + 1)
    return mapping

def standardize_obs_sciplex(adata: AnnData) -> AnnData:
    '''
    Convert observations from AnnData object to the standard format.
    For Sci-plex3 dataset.

    Parameters:
    -----------
    adata : AnnData
       An input scRNA-seq dataset.
    '''
    rename_cols = {'n_genes': 'ngenes', 
                   'pct_counts_mt': 'pcnt_mito', 
                   'cell_type': 'cell_type_', 
                   'cell_line': 'cell_type',
                   'pert_compound': 'perturbagen'
                  }
    categories = {c: 'category' for c in CAT_COLS}

    map_clnames2ids = {'A549': 'CVCL_0023',
                       'K-562': 'CVCL_0004',
                       'MCF7': 'CVCL_0031',
                       }
    
    sc.pp.filter_cells(adata, min_genes=0)
    adata.var["mt"] = adata.var['symbol'].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
    adata.obs['pct_counts_mt'] /= 100
    adata.obs.rename(columns=rename_cols, inplace=True)
    adata.obs['cell_type'] = adata.obs['cell_type'].cat.rename_categories(map_clnames2ids)
    adata.obs['is_control'] = adata.obs['perturbagen'] == 'control'
    adata.obs['perturbagen'] = adata.obs['perturbagen'].cat.add_categories(['DMSO'])
    adata.obs.loc[adata.obs['perturbagen'] == 'control', 'perturbagen'] = 'DMSO'
    adata.obs['perturbagen'] = adata.obs['perturbagen'].cat.remove_categories(['control'])
    adata.obs['library'] = None
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'srivatsan20_sciplex3'
    adata.obs = adata.obs[STANDARD_ADATA_OBS_COLS]
    adata.obs = adata.obs.astype(categories)
    return adata

def standardize_var_sciplex(adata: AnnData) -> AnnData:
    '''
    Convert variables from AnnData object to the standard format.
    For Sci-plex3 dataset.

    Parameters:
    -----------
    adata : AnnData
       An input scRNA-seq dataset.
    '''
    categories = {c: 'category' for c in STANDARD_VAR_COLS}
    adata.var = adata.var[STANDARD_VAR_COLS]
    adata.var = adata.var.astype(categories)
    return adata

# Define the dictionary with the dataset-specific standardization functions
dataset_standardization = {'tahoe': {'standardize_obs': standardize_obs_tahoe,
                                     'standardize_var': standardize_var_tahoe},
                           'sciplex': {'standardize_obs': standardize_obs_sciplex,
                                     'standardize_var': standardize_var_sciplex}
                          }
