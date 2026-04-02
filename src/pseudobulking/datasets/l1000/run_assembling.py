"""
L1000 dataset assembly script.

This script assembles the L1000 Level 3 dataset from GCTX files into a standardized AnnData format.
"""
import os
import argparse
import anndata as ad

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.l1000.assembling import assemble_l1000_dataset


def main():
    """
    Assemble L1000 dataset from GCTX files into AnnData format.
    
    This is the main entry point for the L1000 assembly pipeline. It:
    1. Parses command-line arguments for configuration
    2. Checks if output file already exists (skips if present)
    3. Builds configuration from command-line arguments
    4. Calls the assembly function to process GCTX files
    5. Saves the assembled dataset as an H5AD file
    
    Command-line Arguments
    ----------------------
    --output_file : str, required
        Output path for the assembled H5AD file
    --data_root : str, default='./lincs_data'
        Root directory containing L1000 GCTX and metadata files
    --perturbation_types_to_keep : list of str, default=['trt_cp', 'ctl_vehicle']
        List of perturbation types to include in the dataset
    --control : list of str, default=['DMSO']
        List of control perturbagen names
    --full_gene_matrix : bool, default=False
        If set, uses full gene matrix instead of landmark genes only
    --subsampling : bool, default=False
        If set, uses subsampling of the dataset
    --no-download : bool, default=False
        If set, disables automatic downloading of missing data files
    --annotate-pubchem : bool, default=False
        If set, annotates perturbations with PubChem CIDs using multiple lookup strategies.
        This may involve API calls to PubChem and can be time-consuming.
    --dataset : str, required
        Dataset name (l1000_phase1 or l1000_phase2)
    """
    parser = argparse.ArgumentParser(
        description='Assemble L1000 Level 3 dataset from GCTX files'
    )
    parser.add_argument('--output_file', type=str, required=True,
                        help='Output path for assembled H5AD file')
    parser.add_argument('--data_root', type=str, default='./lincs_data',
                        help='Root directory containing L1000 GCTX files')
    parser.add_argument('--perturbation_types_to_keep', nargs='*', 
                        default=['trt_cp', 'ctl_vehicle'],
                        help='List of perturbation types to keep')
    parser.add_argument('--control', nargs='*', default=['DMSO'],
                        help='List of control perturbagen names')
    parser.add_argument('--full_gene_matrix', action='store_true',
                        help='Use full gene matrix instead of landmark genes only')
    parser.add_argument('--subsampling', action='store_true',
                        help='Use subsampling of the dataset')
    parser.add_argument('--no-download', action='store_true',
                        help='Disable automatic downloading of missing data files')
    parser.add_argument('--annotate-pubchem', action='store_true',
                        help='Annotate perturbations with PubChem CIDs (may involve API calls)')
    parser.add_argument('--dataset', type=str, required=True,
                        choices=['l1000_phase1', 'l1000_phase2'],
                        help='Dataset name (l1000_phase1 or l1000_phase2)')
    
    args = parser.parse_args()
    
    # Check if output already exists
    if os.path.isfile(args.output_file):
        logger.info(f"Assembled dataset already exists: {args.output_file}")
        return
    
    # Build configuration
    config = {
        'perturbation_types_to_keep': set(args.perturbation_types_to_keep),
        'control': set(args.control),
        'full_gene_matrix': args.full_gene_matrix,
        'subsampling': args.subsampling,
        'annotate_pubchem': args.annotate_pubchem,
        'download_if_missing': not args.no_download,
        'dataset': args.dataset,
    }
    
    # Assemble the dataset
    logger.info(f"Assembling {args.dataset} dataset")
    assembled_adata = assemble_l1000_dataset(
        data_root=args.data_root,
        config=config
    )
    
    # Save the assembled dataset
    logger.info(f"Saving assembled dataset to: {args.output_file}")
    assembled_adata.write_h5ad(args.output_file, compression='gzip')
    logger.info(f"{args.dataset} assembly complete")


if __name__ == "__main__":
    main()
