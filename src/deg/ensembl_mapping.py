import requests
import time
import pandas as pd
import numpy as np
from collections import Counter
from typing import Dict, List, Optional, Tuple, Any
from requests.exceptions import HTTPError
from anndata import AnnData
from scipy import sparse

from src.utils.parsing_utils import *


class EnsemblArchiveClient:
    """
    Simple client for Ensembl archive ID lookups with rate limiting.
    
    The Ensembl Archive API allows querying archive information for Ensembl IDs
    to determine if they are current or deprecated, and find replacement IDs.

    https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
    """
    
    def __init__(self, server: str = 'https://rest.ensembl.org', 
                       reqs_per_sec: int = 10,
                       n_retries: int = 10,
                       sleep: int = 10):
        """
        Initialize the client.
        
        Parameters:
        -----------
        server : str, default='https://rest.ensembl.org'
            Ensembl REST API server URL
        reqs_per_sec : int, default=10
            Maximum requests per second
        n_retries : int, default=10
            Number of retry attempts for failed requests
        sleep : int, default=10
            Sleep time in seconds between retries
        """
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.n_retries = n_retries
        self.sleep = sleep
    
    def _rate_limit(self):
        """Apply rate limiting to avoid overwhelming the server."""
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
    
    def get_archive_ids(self, ensembl_ids: List[str]) -> List[Dict[str, Any]]:
        """
        Look up archive information for Ensembl IDs.
        
        Parameters:
        -----------
        ensembl_ids : List[str]
            List of Ensembl gene IDs (e.g., ['ENSG00000251678', ...])
        
        Returns:
        --------
        List[Dict[str, Any]]
            List of dictionaries with archive information for each ID
        """
        self._rate_limit()
        
        url = f"{self.server}/archive/id"
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        payload = {"id": ensembl_ids}
        
        for i in range(self.n_retries):
            try:
                response = requests.post(url, headers=headers, json=payload, timeout=30)
                response.raise_for_status()
                self.req_count += 1
                return response.json()
            except HTTPError as e:
                logger.warning(f"Error fetching archive IDs: {e}")
                
                code = e.response.status_code
                if code == 429:
                    retry_after = int(response.headers.get('Retry-After', self.sleep))
                else:
                    retry_after = self.sleep
                logger.info(f"Retrying {i + 1} out of {self.n_retries}")
                time.sleep(retry_after)
            except Exception as e:
                logger.warning(f"Unexpected error: {e}")
                if i < self.n_retries - 1:
                    time.sleep(self.sleep)
                else:
                    raise
        
        logger.error(f"Failed to fetch archive IDs after {self.n_retries} attempts")
        return []
    
    def get_archive_ids_batch(self, ensembl_ids: List[str], query_size: int = 1000, verbose: bool = True) -> List[Dict[str, Any]]:
        """
        Look up archive information in batches (useful for large lists).
        
        Parameters:
        -----------
        ensembl_ids : List[str]
            List of Ensembl gene IDs
        query_size : int, default=1000
            Number of IDs to query per request
        verbose : bool, default=True
            Show progress with logger
        
        Returns:
        --------
        List[Dict[str, Any]]
            Combined list of dictionaries with archive information
        """
        all_results = []
        
        # Create queries
        queries = [ensembl_ids[i:i + query_size] for i in range(0, len(ensembl_ids), query_size)]
        n_queries = len(queries)
        
        if verbose:
            logger.info(f'Processing {len(ensembl_ids)} IDs in {n_queries} queries')
        
        for query_idx, query in enumerate(queries):
            if verbose:
                logger.info(f'  Query {query_idx + 1}/{n_queries}: processing {len(query)} IDs')
            results = self.get_archive_ids(query)
            all_results.extend(results)
        
        return all_results


def summarize_results(results: List[Dict[str, Any]], 
                     ens_ids: Optional[List[str]] = None, 
                     verbose: bool = True) -> Tuple[Dict[str, str], List[List[str]]]:
    """
    Summarize Ensembl archive API results by creating a mapping dictionary
    and counting current, deprecated, and replaced gene IDs.
    
    Parameters:
    -----------
    results : List[Dict[str, Any]]
        List of dictionaries from Ensembl archive API results
    ens_ids : List[str], optional
        Original list of Ensembl IDs (for comparison)
    verbose : bool, default=True
        If True, print warnings and summary statistics
        
    Returns:
    --------
    Tuple[Dict[str, str], List[List[str]]]
        - mapping: Dictionary mapping old IDs to current/replacement IDs
        - repl: List of lists containing multiple replacements for genes
    """
    if not results:
        if verbose:
            logger.warning("Empty results list provided")
        return {}, []
    
    mapping = {}
    cnt_current = 0
    cnt_deprecated = 0
    cnt_replaced_single = 0
    repl = []
    
    for item in results:
        if item['is_current'] == '1':
            if len(item['possible_replacement']) != 0:
                if verbose:
                    logger.warning(f"Warning, there are replacements {item['possible_replacement']} for existing ENS_ID {item['id']}")
            mapping[item['id']] = item['id']
            cnt_current += 1
        else:
            reps = item['possible_replacement']
            if not reps:
                cnt_deprecated += 1
            elif len(reps) > 1:
                if verbose:
                    logger.warning(f"Warning, the gene {item['id']} has more than one replacements: {reps}")
                repl.append([r['stable_id'] for r in reps])
            else:
                mapping[item['id']] = reps[0]['stable_id']
                cnt_replaced_single += 1
    
    cnt_total = cnt_current + cnt_deprecated + cnt_replaced_single + len(repl)
    assert cnt_total == len(results), f"Count mismatch: {cnt_total} != {len(results)}"
    
    # Check for duplicate mappings (multiple old IDs mapping to the same new ID)
    duplicates = {k: v for k, v in Counter(mapping.values()).items() if v > 1}
    
    if verbose:
        logger.info(f"Summary:")
        if ens_ids is not None:
            logger.info(f"  Total results: {cnt_total} out of {len(ens_ids)}")
        else:
            logger.info(f"  Total results: {cnt_total}")
        logger.info(f"  Current IDs: {cnt_current}")
        logger.info(f"  Deprecated (no replacement): {cnt_deprecated}")
        logger.info(f"  Single replacements: {cnt_replaced_single}")
        logger.info(f"  Multiple replacements: {len(repl)}")
        logger.info(f"  Unique target IDs with duplicates: {len(duplicates)}")
        logger.info(f"  Mapping size: {len(mapping)}")
    
    return mapping, repl


def _aggregate_matrix(matrix, groups, new_id_to_index, var_index):
    """
    Aggregate matrix (X or layer) for duplicated genes using vectorized matrix multiplication.
    
    Parameters:
    -----------
    matrix : scipy.sparse or np.ndarray
        Gene expression matrix (n_obs x n_genes_old)
    groups : pd.core.groupby.DataFrameGroupBy
        Grouped dataframe by new_id
    new_id_to_index : dict
        Mapping from new gene ID to column index
    var_index : pd.Index
        Original var index with old gene IDs
        
    Returns:
    --------
    scipy.sparse or np.ndarray
        Aggregated matrix (n_obs x n_genes_new)
    """
    # Build aggregation matrix: maps old gene indices to new gene indices
    # Shape: (n_genes_old, n_genes_new)
    # Each row corresponds to an old gene, each column to a new gene
    # Value is 1 if old gene maps to new gene, 0 otherwise
    old_cols = []
    new_cols = []
    
    for new_id, group in groups:
        if pd.isna(new_id):
            continue
        new_col = new_id_to_index[new_id]
        for idx in group.index:
            old_cols.append(var_index.get_loc(idx))
            new_cols.append(new_col)
    
    # Create sparse aggregation matrix
    # Use matrix dtype for data to preserve precision
    data = np.ones(len(old_cols), dtype=matrix.dtype if hasattr(matrix, 'dtype') else np.float32)
    n_genes_old = matrix.shape[1]
    n_genes_new = len(new_id_to_index)
    
    aggregation_matrix = sparse.csr_matrix(
        (data, (old_cols, new_cols)),
        shape=(n_genes_old, n_genes_new),
        dtype=matrix.dtype if hasattr(matrix, 'dtype') else np.float32
    )
    
    # Matrix multiplication: matrix @ aggregation_matrix
    # This sums columns (old genes) that map to the same new gene in a single operation
    agg_matrix = matrix @ aggregation_matrix
    
    return agg_matrix


def map_and_aggregate_duplicated_genes(adata: AnnData, 
                                      mapping: Dict[str, str],
                                      var_index_col: Optional[str] = None,
                                      keep_unmapped: bool = False, 
                                      verbose: bool = True) -> AnnData:
    """
    Map gene IDs and aggregate counts for duplicated genes in AnnData.
    
    When multiple old gene IDs map to the same new gene ID, this function aggregates
    their expression counts using vectorized matrix multiplication for efficiency.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object to map and aggregate
    mapping : Dict[str, str]
        Dictionary mapping old gene IDs to new gene IDs
    var_index_col : str, optional
        Column name in var if IDs are not in index. If None, uses var.index
    keep_unmapped : bool, default=False
        If True, keep genes not in mapping (retain original IDs). If False, remove them.
    verbose : bool, default=True
        Show progress and summary statistics
        
    Returns:
    --------
    AnnData
        AnnData object with mapped and aggregated gene IDs
    """
    if var_index_col and var_index_col not in adata.var.columns:
        raise ValueError(f"Column '{var_index_col}' not found in adata.var")
    
    gene_ids = adata.var[var_index_col] if var_index_col else adata.var.index
    mapped_ids = np.array([mapping.get(gid, gid if keep_unmapped else None) for gid in gene_ids], dtype=object)
    
    if not keep_unmapped:
        mask = mapped_ids != None
        adata = adata[:, mask].copy()
        mapped_ids = mapped_ids[mask]
        gene_ids = gene_ids[mask]
    
    mapping_df = pd.DataFrame({'old_id': gene_ids.values, 'new_id': mapped_ids}, index=adata.var.index)
    groups = mapping_df.groupby('new_id')
    unique_new_ids = mapping_df['new_id'].dropna().unique()
    
    
    new_id_to_index = {gid: idx for idx, gid in enumerate(unique_new_ids)}
    X_agg = _aggregate_matrix(adata.X, groups, new_id_to_index, adata.var.index)
    
    # Create new var dataframe with first occurrence of each new_id group
    new_var = pd.DataFrame({
        col: [adata.var.loc[groups.get_group(new_id).index[0], col] for new_id in unique_new_ids]
        for col in adata.var.columns if col != var_index_col
    }, index=unique_new_ids)
    
    # Create new AnnData with aggregated data
    adata_mapped = AnnData(
        X=X_agg, 
        obs=adata.obs.copy(), 
        var=new_var,
        uns=getattr(adata, 'uns', {}).copy(),
        obsm=getattr(adata, 'obsm', {}).copy(),
        varm=getattr(adata, 'varm', {}).copy()
    )
    
    
    return adata_mapped


def map_ensembl_ids_for_dataset(padata: AnnData,
                                ensembl_id_col: str = 'ensembl_id',
                                query_size: int = 1000,
                                verbose: bool = True) -> AnnData:
    """
    Map Ensembl gene IDs in AnnData object to current/replacement IDs.
    
    This function:
    1. Extracts Ensembl IDs from var.index or var column
    2. Queries Ensembl Archive API to find current/replacement IDs
    3. Aggregates counts for duplicated genes
    4. Filters out genes that couldn't be mapped (deprecated with no replacement or multiple replacements)
    
    Parameters:
    -----------
    padata : AnnData
        AnnData object with Ensembl IDs in var.index or var[ensembl_id_col]
    ensembl_id_col : str, default='ensembl_id'
        Column name in var if IDs are not in index
    query_size : int, default=1000
        Number of IDs to query per batch
    verbose : bool, default=True
        Show progress and summary statistics
        
    Returns:
    --------
    AnnData
        AnnData object with updated Ensembl IDs (current/replacement IDs)
    """
    logger.info('Mapping Ensembl IDs to current/replacement IDs')
    
    # Extract current Ensembl IDs
    if ensembl_id_col in padata.var.columns:
        current_ids = padata.var[ensembl_id_col].values.tolist()
        id_source = 'column'
    else:
        current_ids = padata.var_names.tolist()
        id_source = 'index'
    
    logger.info(f'  Found {len(current_ids)} genes with Ensembl IDs in var.{id_source}')
    
    # Create client and get archive information
    client = EnsemblArchiveClient()
    results = client.get_archive_ids_batch(
        ensembl_ids=current_ids,
        query_size=query_size,
        verbose=verbose
    )
    
    # Summarize results and create mapping
    mapping, repl = summarize_results(results, ens_ids=current_ids, verbose=verbose)
    
    # Count unmapped genes (deprecated with no replacement + genes with multiple replacements)
    unmapped_count = len(current_ids) - len(mapping)
    
    if unmapped_count > 0:
        logger.warning(f'Filtering out {unmapped_count} genes that could not be mapped (deprecated with no replacement or multiple replacements)')
    
    # Check for duplicates (multiple old IDs mapping to same new ID)
    duplicates = {k: v for k, v in Counter(mapping.values()).items() if v > 1}
    var_index_col = None if id_source == 'index' else ensembl_id_col
    
    if len(duplicates) > 0:
        logger.warning(f'Found {len(duplicates)} target IDs with duplicates - will aggregate counts')
        logger.info(f'Aggregating running...')
        logger.info(f"  The layer containing 'psbulk_props' will be deprecated.")
        padata_mapped = map_and_aggregate_duplicated_genes(
            padata, mapping, var_index_col=var_index_col, keep_unmapped=False, verbose=verbose)
        if id_source != 'index' and padata_mapped.var.index.name != ensembl_id_col:
            padata_mapped.var.set_index(ensembl_id_col, inplace=True)
        # Add is_merged column: True for genes that resulted from merging multiple old IDs
        padata_mapped.var['is_merged'] = padata_mapped.var.index.isin(duplicates.keys())
        
    else:
        # No duplicates - simple mapping without aggregation
        padata_mapped = padata.copy()
        if id_source == 'index':
            mapped_ids = [old_id for old_id in padata_mapped.var_names if old_id in mapping]
            padata_mapped = padata_mapped[:, mapped_ids].copy()
            padata_mapped.var.index = pd.Index([mapping[old_id] for old_id in padata_mapped.var_names], name='ensembl_id')
        else:
            padata_mapped.var[ensembl_id_col] = padata_mapped.var[ensembl_id_col].map(mapping)
            padata_mapped = padata_mapped[:, padata_mapped.var[ensembl_id_col].notna()].copy()
            if padata_mapped.var.index.name != ensembl_id_col:
                padata_mapped.var.set_index(ensembl_id_col, inplace=True)
        # Add is_merged column (all False since no duplicates)
        padata_mapped.var['is_merged'] = False
    if verbose:
        logger.info(f"  Result: {padata.n_vars} → {padata_mapped.n_vars} genes (reduced by {padata.n_vars - padata_mapped.n_vars})")
    
    return padata_mapped
