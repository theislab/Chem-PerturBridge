"""
Dataset-specific processing functions for Tahoe dataset.
"""
from typing import Optional
import os
import pandas as pd
import anndata as ad
from datasets import load_dataset

from src.utils.parsing_utils import logger


# Constants
TAHOE_HF_DATASET = "vevotx/Tahoe-100M"



def process_tahoe_dataset(padata: ad.AnnData,
                         var_lamindb_path: Optional[str] = None) -> ad.AnnData:
    """
    Process Tahoe dataset by remapping Ensembl IDs to HuggingFace version.
    
    This function merges gene information from three sources:
    1. Current AnnData var (gene annotations in the pseudobulked data)
    2. var.parquet file containing LaminDB Ensembl IDs (preserving in current AnnData var) 
       and the original Ensembl ID mappings
    3. HuggingFace Tahoe-100M dataset (canonical gene metadata with Ensembl IDs)
    
    Parameters
    ----------
    padata : ad.AnnData
        Input AnnData object with pseudobulked data.
    var_lamindb_path : Optional[str], default=None
        Path to the var.parquet file from LaminDB containing gene annotations.
        Must be provided (cannot be None).
    
    Returns
    -------
    ad.AnnData
        Processed AnnData object with updated gene annotations in the .var attribute.
        The .var dataframe will have 'ensembl_id' as index and 'symbol' column.
    
    Raises
    ------
    ValueError
        If var_lamindb_path is None.
    FileNotFoundError
        If the file at var_lamindb_path does not exist.
    
    Notes
    -----
    - The function creates a copy of the input AnnData to avoid modifying the original.
    - Gene metadata is downloaded from HuggingFace dataset 'vevotx/Tahoe-100M'.
    - The merge process maps genes through LaminDB to HuggingFace Ensembl IDs.
    """
    
    # Validate var_lamindb_path
    if var_lamindb_path is None:
        raise ValueError("var_lamindb_path must be provided")
    
    if not os.path.exists(var_lamindb_path):
        raise FileNotFoundError(f'var.parquet file not found at: {var_lamindb_path}')
    
    # Copy to avoid modifying original
    padata_processed = padata.copy()
    
    # Get current var dataframe from AnnData
    var_current = padata_processed.var.copy().reset_index()
    
    logger.info(f'  Loading var.parquet from: {var_lamindb_path}')
    var_lamindb = pd.read_parquet(var_lamindb_path)
    
    # Ensure gene_name is the index in var_lamindb
    if 'gene_name' in var_lamindb.columns:
        var_lamindb = var_lamindb.set_index('gene_name')
    elif var_lamindb.index.name != 'gene_name':
        var_lamindb.index.name = 'gene_name'
    
    logger.info(f'  Loaded {len(var_lamindb)} genes from var.parquet')
    
    # Download gene metadata from HuggingFace
    logger.info('  Downloading gene metadata from HuggingFace')
    gene_metadata = load_dataset(TAHOE_HF_DATASET, name="gene_metadata", split="train")
    
    # Create vocabulary DataFrame
    gene_vocab = {
        'gene_symbol': [entry["gene_symbol"] for entry in gene_metadata],
        'ensembl_id': [entry["ensembl_id"] for entry in gene_metadata],
    }
    
    var_hf = pd.DataFrame(gene_vocab)
    # Set gene_symbol as index for merging
    var_hf = var_hf.set_index('gene_symbol')

    # Merge var_lamindb with HuggingFace metadata
    var_lamindb_merged_hf = var_lamindb.merge(var_hf, left_index=True, right_index=True, how='left')
    
    # Merge with current var
    var_current_updated = var_current.merge(
        var_lamindb_merged_hf, 
        left_on='ensembl_id', 
        right_on='ensembl_gene_id', 
        suffixes=('', '_hf')
    )
    var_current_updated = var_current_updated[['ensembl_id_hf', 'symbol']].rename(
        columns={'ensembl_id_hf': 'ensembl_id'}
    ).set_index('ensembl_id')

    padata_processed.var = var_current_updated.copy()
    
    return padata_processed

