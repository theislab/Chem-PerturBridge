from string import ascii_lowercase
import re
from typing import List, Any, Dict, Optional
import scanpy as sc
from anndata import AnnData
import pandas as pd

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
    map_clnames2ids = {'A549': 'CVCL_0023',
                       'K-562': 'CVCL_0004',
                       'MCF7': 'CVCL_0031',
                       }

    sc.pp.filter_cells(adata, min_genes=0)
    adata.var["mt"] = adata.var['symbol'].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
    adata.obs['pct_counts_mt'] /= 100
    adata.obs.rename(columns=rename_cols, inplace=True)
    adata.obs['cell_type'] = adata.obs['cell_type'].cat.rename_categories(
        map_clnames2ids)
    adata.obs['is_control'] = adata.obs['perturbagen'] == 'control'
    adata.obs['perturbagen'] = adata.obs['perturbagen'].cat.add_categories([
                                                                           'DMSO'])
    adata.obs.loc[adata.obs['perturbagen']
                  == 'control', 'perturbagen'] = 'DMSO'
    adata.obs['perturbagen'] = adata.obs['perturbagen'].cat.remove_categories([
                                                                              'control'])
    adata.obs = adata.obs.drop(columns=['pert_time'])
    adata.obs['pert_time'] = adata.obs['time'].apply(lambda x: str(x) + 'h' if not pd.isna(x) else x)
    adata.obs['library'] = None
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'srivatsan20_sciplex3'
    return adata

def standardize_var_sciplex(adata: AnnData) -> AnnData:
    '''
    Convert variables from AnnData object to the standard format.
    For Sci-plex3 dataset. To maintain the same set of 
    standardization functions.

    Parameters:
    -----------
    adata : AnnData
       An input scRNA-seq dataset.
    '''
    return adata

def standardize_sciplex(adata: AnnData) -> AnnData:
    '''
    Convert AnnData object to the standard format.
    For Sci-plex3 dataset.

    Parameters:
    -----------
    adata : AnnData
       An input scRNA-seq dataset.
    '''
    adata = standardize_obs_sciplex(adata)
    adata = standardize_var_sciplex(adata)
    return adata
