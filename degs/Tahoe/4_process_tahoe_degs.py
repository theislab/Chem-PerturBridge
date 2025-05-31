#!/usr/bin/env python3
"""
Process differential expression results for all cell lines in the Tahoe dataset.
Converts DEG results to AnnData format with proper metadata annotations.
"""

import pandas as pd
import anndata
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Optional
import pubchempy as pcp
from tqdm import tqdm

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def get_cell_line_names(data_dir: Path) -> List[str]:
    """Extract cell line names from available parquet files."""
    parquet_files = data_dir.glob("*_differential_expression_results.parquet")
    cell_lines = []
    
    for file in parquet_files:
        # Extract cell line name from filename
        cell_line = file.stem.replace("_differential_expression_results", "")
        cell_lines.append(cell_line)
    
    return sorted(cell_lines)


def load_experiment_dict(dict_path: Path) -> Dict[str, str]:
    """Load the experiment dictionary mapping conditions to drug info."""
    df = pd.read_csv(dict_path)
    return dict(zip(df['Unnamed: 0'], df['0']))


def process_deg_results(deg_results: pd.DataFrame, cell_line_name: str) -> pd.DataFrame:
    """Process differential expression results for a single cell line."""
    # Filter out plate-related conditions
    deg_results = deg_results[~deg_results['condition'].str.startswith('plate')]
    
    # Add cell line to condition name
    deg_results["condition"] = deg_results["condition"] + "_" + deg_results["cellline"]
    
    return deg_results


def create_pivoted_matrices(deg_results: pd.DataFrame) -> tuple:
    """Create pivoted matrices for each metric."""
    metrics = {
        'logFC':     'logfc',
        'AveExpr':   'ave_expr',
        'P.Value':   'p_value',
        'adj.P.Val': 'adj_p_value',
        'B':         'b'
    }
    
    pivoted = {}
    for orig_col, layer_name in metrics.items():
        pivoted[layer_name] = deg_results.pivot(
            index='condition',
            columns='gene',
            values=orig_col
        )
    
    # Capture ordering
    conditions = pivoted['logfc'].index
    genes = pivoted['logfc'].columns
    
    return pivoted, conditions, genes


def create_obs_var(deg_results: pd.DataFrame, conditions: pd.Index, genes: pd.Index) -> tuple:
    """Create observation and variable metadata."""
    # obs: one row per condition
    obs = (
        deg_results
        .drop_duplicates(subset=['condition'])
        .set_index('condition')
        [['cellline']]
        .loc[conditions]  # align order
    )
    
    # var: one row per gene
    var = pd.DataFrame(index=genes)
    var['gene_symbol'] = genes
    
    return obs, var


def annotate_drug_info(adata: anndata.AnnData, experiment_dict: Dict[str, str]) -> anndata.AnnData:
    """Add drug name and dosage information to AnnData object."""
    # Extract base condition name (remove cell line suffix)
    adata.obs['condition_base'] = adata.obs.index.str.replace(r'_[^_]+$', '', regex=True)
    
    # Map to drug info
    adata.obs['drug_name'] = adata.obs['condition_base'].map(experiment_dict)
    
    # Extract dosage and clean drug name
    adata.obs['dosage_uM'] = adata.obs['drug_name'].apply(
        lambda x: eval(x)[0][1] if isinstance(x, str) else None
    )
    adata.obs['drug_name'] = adata.obs['drug_name'].apply(
        lambda x: eval(x)[0][0] if isinstance(x, str) else None
    )
    
    return adata


def get_pubchem_cid(drug_name: str, drug_to_cid_cache: Dict[str, Optional[int]]) -> Optional[int]:
    """
    Fetch PubChem CID for a given drug name, using cache to avoid repeated lookups.
    
    Args:
        drug_name: Name of the drug to look up
        drug_to_cid_cache: Dictionary cache of previously fetched CIDs
        
    Returns:
        PubChem CID or None if not found
    """
    if pd.isna(drug_name) or drug_name is None:
        return None
    
    # Check cache first
    if drug_name in drug_to_cid_cache:
        return drug_to_cid_cache[drug_name]
    
    # Not in cache, fetch from PubChem
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        cid = compounds[0].cid if compounds else None
        logger.debug(f"Fetched CID for '{drug_name}': {cid}")
    except Exception as e:
        logger.warning(f"Error fetching CID for {drug_name}: {str(e)}")
        cid = None
    
    # Cache the result (including None values to avoid repeated failed lookups)
    drug_to_cid_cache[drug_name] = cid
    return cid


def add_pubchem_cids_to_adata(adata: anndata.AnnData, drug_to_cid_cache: Dict[str, Optional[int]]) -> None:
    """
    Add PubChem CID column to AnnData object, fetching CIDs as needed.
    
    Updates the cache with any new drug names encountered.
    """
    # Get unique drug names in this dataset
    unique_drugs = adata.obs['drug_name'].dropna().unique()
    
    # Fetch CIDs for any new drugs not in cache
    new_drugs = [drug for drug in unique_drugs if drug not in drug_to_cid_cache]
    if new_drugs:
        logger.info(f"Fetching CIDs for {len(new_drugs)} new drugs...")
        for drug in tqdm(new_drugs, desc="Fetching new drug CIDs"):
            get_pubchem_cid(drug, drug_to_cid_cache)
    
    # Map all drugs to their CIDs
    adata.obs['pubchem_cid'] = adata.obs['drug_name'].map(drug_to_cid_cache)


def process_cell_line(cell_line_name: str, data_dir: Path, experiment_dict: Dict[str, str],
                     output_dir: Path, drug_to_cid_cache: Dict[str, Optional[int]]) -> Optional[anndata.AnnData]:
    """Process a single cell line's differential expression results."""
    logger.info(f"Processing cell line: {cell_line_name}")
    
    try:
        # Load data
        deg_file = data_dir / f"{cell_line_name}_differential_expression_results.parquet"
        deg_results = pd.read_parquet(deg_file)
        
        # Process DEG results
        deg_results = process_deg_results(deg_results, cell_line_name)
        
        # Create pivoted matrices
        pivoted, conditions, genes = create_pivoted_matrices(deg_results)
        
        # Create obs and var
        obs, var = create_obs_var(deg_results, conditions, genes)
        
        # Create AnnData object
        # Note: Using empty X matrix since we'll store everything in layers
        adata = anndata.AnnData(
            X   = np.empty((len(conditions), len(genes))),
            obs = obs,
            var = var
        )
        
        # Add all metrics as layers
        for layer_name, df in pivoted.items():
            adata.layers[layer_name] = df.loc[conditions, genes].values
        
        # Remove empty X matrix
        del adata.X
        
        # Annotate with drug information
        adata = annotate_drug_info(adata, experiment_dict)
        
        # Add PubChem CIDs
        add_pubchem_cids_to_adata(adata, drug_to_cid_cache)
        
        # Save
        output_file = output_dir / f"tahoe_{cell_line_name.lower()}_deg_adata.h5ad"
        adata.write_h5ad(output_file, compression="gzip")
        logger.info(f"Saved: {output_file}")
        
        return adata
        
    except Exception as e:
        logger.error(f"Error processing {cell_line_name}: {str(e)}")
        return None


def main():
    """Main processing function."""
    # Set up paths
    raw_data_dir = Path("../../data/raw/tahoe")
    output_dir = Path("../../data/degs")
    experiment_dict_path = Path("../Tahoe/Experiment_dict_Tahoe100m.csv")
    
    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load experiment dictionary
    logger.info("Loading experiment dictionary...")
    experiment_dict = load_experiment_dict(experiment_dict_path)
    
    # Get all cell line names
    cell_lines = get_cell_line_names(raw_data_dir)
    logger.info(f"Found {len(cell_lines)} cell lines: {', '.join(cell_lines)}")
    
    # Initialize drug to CID cache
    drug_to_cid_cache = {}
    
    # Process each cell line
    successful = 0
    failed = []
    
    for cell_line in cell_lines:
        result = process_cell_line(cell_line, raw_data_dir, experiment_dict, 
                                 output_dir, drug_to_cid_cache)
        if result is not None:
            successful += 1
        else:
            failed.append(cell_line)
    
    # Summary
    logger.info("\n" + "="*60)
    logger.info("Processing complete!")
    logger.info(f"Successfully processed: {successful}/{len(cell_lines)} cell lines")
    if failed:
        logger.warning(f"Failed cell lines: {', '.join(failed)}")
    
    # Report CID cache statistics
    total_drugs = len(drug_to_cid_cache)
    drugs_with_cid = sum(1 for cid in drug_to_cid_cache.values() if cid is not None)
    logger.info(f"\nPubChem CID Statistics:")
    logger.info(f"Total unique drugs encountered: {total_drugs}")
    logger.info(f"Drugs with valid CIDs: {drugs_with_cid} ({drugs_with_cid/total_drugs*100:.1f}%)")
    logger.info(f"Drugs without CIDs: {total_drugs - drugs_with_cid}")
    
    # Optional: Save the drug-to-CID mapping for future use
    if drug_to_cid_cache:
        cache_file = output_dir / "drug_to_pubchem_cid_mapping.csv"
        cache_df = pd.DataFrame(
            [(drug, cid) for drug, cid in drug_to_cid_cache.items()],
            columns=['drug_name', 'pubchem_cid']
        ).sort_values('drug_name')
        cache_df.to_csv(cache_file, index=False)
        logger.info(f"\nSaved drug-to-CID mapping to: {cache_file}")


if __name__ == "__main__":
    main()