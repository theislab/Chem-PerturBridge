"""
OP3 dataset standardization script.

This script standardizes the OP3 NeurIPS 2023 dataset to match the unified pseudobulk schema.
"""
import os
import argparse

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.op3.standardization import standardize_op3_dataset


def main():
    """
    Standardize OP3 dataset to unified pseudobulk schema.
    
    This is the main entry point for the OP3 standardization pipeline. It:
    1. Parses command-line arguments for configuration
    2. Checks if output file already exists (skips if present)
    3. Downloads data from S3 if missing (unless --no-download is set)
    4. Calls the standardization function to process the dataset
    5. Saves the standardized dataset as an H5AD file
    
    Command-line Arguments
    ----------------------
    --output_file : str, required
        Output path for the standardized H5AD file
    --data_root : str, default='./op3_data'
        Root directory containing OP3 data files
    --file : str, default='pseudobulk_filtered_with_uns.h5ad'
        Input filename for the raw OP3 data
    --no-download : bool, default=False
        If set, disables automatic downloading of missing data files
    --annotate-pubchem : bool, default=False
        If set, annotates compounds with PubChem CIDs (may involve API calls)
    """
    parser = argparse.ArgumentParser(
        description='Standardize OP3 NeurIPS 2023 dataset to unified pseudobulk schema'
    )
    parser.add_argument('--output_file', type=str, required=True,
                        help='Output path for standardized H5AD file')
    parser.add_argument('--data_root', type=str, default='./op3_data',
                        help='Root directory for OP3 data files')
    parser.add_argument('--file', type=str, default='pseudobulk_filtered_with_uns.h5ad',
                        help='Input filename')
    parser.add_argument('--no-download', action='store_true',
                        help='Disable automatic downloading of missing data files')
    parser.add_argument('--annotate-pubchem', action='store_true',
                        help='Annotate compounds with PubChem CIDs (may involve API calls)')
    
    args = parser.parse_args()
    
    # Check if output already exists
    if os.path.isfile(args.output_file):
        logger.info(f"Standardized dataset already exists: {args.output_file}")
        return
    
    # Standardize the dataset
    logger.info("Standardizing OP3 dataset")
    adata_standardized = standardize_op3_dataset(
        data_root=args.data_root,
        file=args.file,
        download_if_missing=not args.no_download,
        annotate_pubchem=args.annotate_pubchem
    )
    
    # Save the standardized dataset
    logger.info(f"Saving standardized dataset to: {args.output_file}")
    adata_standardized.write_h5ad(args.output_file, compression='gzip')
    logger.info("OP3 standardization complete")


if __name__ == "__main__":
    main()
