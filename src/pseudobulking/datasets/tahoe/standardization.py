from string import ascii_lowercase
import re
from typing import List, Any, Dict, Optional
import scanpy as sc
from anndata import AnnData

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
        mapping[i + 1] = ascii_lowercase[i // ncols].upper() + str(i %
                                                                   ncols + 1)
    return mapping

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
    wellnum2id_mapping = wellnum2id()
    adata.obs['well'] = adata.obs.index.to_series().str.split(
        '_').str[0].astype(int).map(wellnum2id_mapping)
    adata.obs.rename(columns=rename_cols, inplace=True)
    adata.obs['is_control'] = adata.obs['perturbagen'] == 'DMSO_TF'
    adata.obs['library'] = None
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'tahoe100'
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

    adata.var.rename(columns=rename_cols, inplace=True)
    adata.var.set_index('ensembl_id', inplace=True)
    return adata

def standardize_tahoe(adata: AnnData) -> AnnData:
    '''
    Convert AnnData object to the standard format.
    For Tahoe dataset.

    Parameters:
    -----------
    adata : AnnData
       An input scRNA-seq dataset.
    '''
    adata = standardize_obs_tahoe(adata)
    adata = standardize_var_tahoe(adata)
    return adata
