#!/usr/bin/env python3
"""
Process differential expression results from Sciplex data for multiple cell lines.
This script processes the DEG results for A549, K562, and MCF7 cell lines and saves
them as AnnData objects.
"""

import pandas as pd
import anndata
import mygene
from pathlib import Path
import logging
from typing import Dict, List, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
CELL_LINES = ['A549', 'K562', 'MCF7']
METRICS = {
    'logFC': 'logfc',
    'AveExpr': 'ave_expr',
    'P.Value': 'p_value',
    'adj.P.Val': 'adj_p_value',
    'B': 'b'
}

def load_experiment_dict() -> Dict[str, str]:
    """Load the experiment dictionary mapping condition names to drug names."""
    experiment_dict = pd.read_csv(
        Path(__file__).parent.parent / "Sciplex/Experiment_dict_Sciplex.csv",
        index_col='Unnamed: 0'
    )['0'].to_dict()
    return experiment_dict

def process_deg_results(deg_results: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], List[str]]:
    """
    Process DEG results by pivoting metrics and creating obs/var dataframes.
    
    Args:
        deg_results: DataFrame containing differential expression results
        
    Returns:
        Tuple of (pivoted metrics dict, obs dataframe, var dataframe, conditions list, genes list)
    """
    # Filter out plate conditions
    deg_results = deg_results[~deg_results['condition'].str.startswith('plate')]
    
    # Add cell-line into condition
    deg_results["condition"] = deg_results["condition"] + "_" + deg_results["cellline"]
    
    # Pivot each metric
    pivoted = {}
    for orig_col, layer_name in METRICS.items():
        pivoted[layer_name] = deg_results.pivot(
            index='condition',
            columns='gene',
            values=orig_col
        )
    
    # Capture the ordering
    conditions = pivoted['logfc'].index
    genes = pivoted['logfc'].columns
    
    # Build obs (one row per condition)
    obs = (
        deg_results
        .drop_duplicates(subset=['condition'])
        .set_index('condition')
        [['cellline']]
        .loc[conditions]
    )
    
    # Build var (one row per gene)
    var = pd.DataFrame(index=genes)
    var['gene_name'] = genes
    
    return pivoted, obs, var, conditions, genes

def add_gene_symbols(var_df: pd.DataFrame) -> pd.DataFrame:
    """Add gene symbols to the var dataframe using MyGene."""
    mg = mygene.MyGeneInfo()
    
    # Split human vs mouse IDs
    human_ids = [i for i in var_df['gene_name'] if i.startswith('ENSG')]
    mouse_ids = [i for i in var_df['gene_name'] if i.startswith('ENSMUSG')]
    
    # Query gene symbols
    human_q = mg.querymany(
        human_ids,
        scopes='ensembl.gene',
        fields='symbol',
        species='human',
        as_dataframe=True
    ) if human_ids else pd.DataFrame()
    
    mouse_q = mg.querymany(
        mouse_ids,
        scopes='ensembl.gene',
        fields='symbol',
        species='mouse',
        as_dataframe=True
    ) if mouse_ids else pd.DataFrame()
    
    # Combine results
    mapping = {}
    if not human_q.empty:
        mapping.update(human_q['symbol'].to_dict())
    if not mouse_q.empty:
        mapping.update(mouse_q['symbol'].to_dict())
    
    # Add gene symbols
    var_df['gene_symbol'] = var_df['gene_name'].map(mapping)
    return var_df

def create_anndata(pivoted: Dict[str, pd.DataFrame], 
                  obs: pd.DataFrame, 
                  var: pd.DataFrame,
                  conditions: List[str],
                  genes: List[str],
                  experiment_dict: Dict[str, str]) -> anndata.AnnData:
    """Create and process AnnData object."""
    # Create AnnData with logFC as X
    adata = anndata.AnnData(
        X=pivoted['logfc'].loc[conditions, genes].values,
        obs=obs,
        var=var
    )
    
    # Add other stats as layers
    for layer_name in ['ave_expr', 'p_value', 'adj_p_value', 'b']:
        adata.layers[layer_name] = pivoted[layer_name].loc[conditions, genes].values
    
    adata.layers["logfc"] = adata.X
    del adata.X
    
    # Process observation metadata
    adata.obs = adata.obs.reset_index()
    adata.obs['drug'] = adata.obs['condition'].str.extract(
        r'(?:factor_valid_drug_names_)?(.*?)_\d+_0_'
    )[0]
    adata.obs['dosage_uM'] = adata.obs.condition.str.extract(
        r'_(\d+)_0_[^_]+$'
    )[0].astype(float) / 1000
    
    # Add drug names
    adata.obs['condition_base'] = adata.obs['condition'].str.replace(r'_[^_]+$', '', regex=True)
    adata.obs['drug_name'] = adata.obs['condition_base'].map(experiment_dict)
    adata.obs['drug_name'] = adata.obs['drug_name'].str.replace(r'_[^_]+$', '', regex=True)
    
    # Set index back to condition
    adata.obs = adata.obs.set_index('condition')
    
    return adata

def process_cell_line(cell_line: str, experiment_dict: Dict[str, str], output_dir: Path) -> None:
    """Process a single cell line's DEG results."""
    logger.info(f"Processing {cell_line}...")
    
    # Load data
    input_path = Path(f"../../data/raw/sciplex/{cell_line}_differential_expression_results.parquet")
    deg_results = pd.read_parquet(input_path)
    
    # Process results
    pivoted, obs, var, conditions, genes = process_deg_results(deg_results)
    
    # Add gene symbols
    var = add_gene_symbols(var)
    
    # Create AnnData
    adata = create_anndata(pivoted, obs, var, conditions, genes, experiment_dict)
    
    # Save results
    output_path = output_dir / f"sciplex3_{cell_line.lower()}_deg_adata.h5ad"
    adata.write_h5ad(output_path, compression="gzip")
    logger.info(f"Saved results to {output_path}")

def main():
    """Main function to process all cell lines."""
    # Create output directory if it doesn't exist
    output_dir = Path("../../data/degs")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load experiment dictionary
    experiment_dict = load_experiment_dict()
    
    # Process each cell line
    for cell_line in CELL_LINES:
        try:
            process_cell_line(cell_line, experiment_dict, output_dir)
        except Exception as e:
            logger.error(f"Error processing {cell_line}: {str(e)}")
            continue

if __name__ == "__main__":
    main() 