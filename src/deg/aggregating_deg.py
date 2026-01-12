import os
import re
import argparse
import json
import anndata as ad
import h5py
from typing import List
from os.path import join

from src.utils.parsing_utils import *


def rewrite_h5ad(path2file: str) -> None:
    '''
    Rewrite h5ad file by completely removing X matrix.
    
    Completely removes the X matrix from an h5ad file using direct HDF5 access.
    This is equivalent to the R rewrite_h5ad function and ensures older
    Python anndata versions can read the file.
    
    Parameters:
    -----------
    path2file : str
        Path to the h5ad file
        
    Returns:
    --------
    None
        Modifies file in place
    '''
    if not os.path.exists(path2file):
        raise FileNotFoundError(f"File not found: {path2file}")
    
    logger.info(f'    Rewriting h5ad file: {path2file}')
    
    # Open file in read-write mode
    with h5py.File(path2file, 'r+') as f:
        # Check if X group exists and remove it completely
        if 'X' in f.keys():
            del f['X']
            logger.info('    Removed X matrix')
        else:
            logger.info('    X matrix not found or already removed')
    
    logger.info('    Rewriting h5ad file is done!')


def get_batch_files(cell_type_dir: str, cell_type: str) -> List[str]:
    '''
    Get all batch files for a cell type, sorted by filename.
    
    Parameters:
    -----------
    cell_type_dir : str
        Directory containing batch files for the cell type
    cell_type : str
        Cell type name (not used, kept for API compatibility)
        
    Returns:
    --------
    List[str]
        Sorted list of batch file paths
    '''
    # List all .h5ad files in the directory
    batch_files = []
    for filename in os.listdir(cell_type_dir):
        if filename.endswith('.h5ad'):
            batch_files.append(join(cell_type_dir, filename))
    
    # Sort by batch number (numeric order: 1, 2, 3, ..., 10, 11, 12, ...)
    def get_batch_num(filepath):
        """Extract batch number from filename for numeric sorting using regex."""
        basename = os.path.basename(filepath)
        match = re.search(r'_batch_(\d+)', basename)
        return int(match.group(1)) if match else float('inf')
    
    return sorted(batch_files, key=get_batch_num)


def aggregate_cell_type_batches(input_dir: str, 
                                cell_type: str,
                                output_dir: str) -> None:
    '''
    Aggregate all batch files for a cell type into a single file.
    
    Adds a 'processing_batch' column to .obs indicating which batch each observation
    came from (e.g., 'batch_1', 'batch_2', etc.).
    
    Parameters:
    -----------
    input_dir : str
        Base input directory containing cell type subdirectories
    cell_type : str
        Cell type name (subdirectory name)
    output_dir : str
        Output directory where aggregated file will be saved
        
    Returns:
    --------
    None
        Saves aggregated file to disk with batch information in .obs['processing_batch']
    '''
    cell_type_dir = join(input_dir, cell_type)
    
    if not os.path.isdir(cell_type_dir):
        logger.warning(f'Cell type directory not found: {cell_type_dir}')
        return
    
    batch_files = get_batch_files(cell_type_dir, cell_type)
    
    if len(batch_files) == 0:
        logger.warning(f'No batch files found for {cell_type} in {cell_type_dir}')
        return
    
    logger.info(f'Aggregating {len(batch_files)} batch files for {cell_type}')
    
    # Read all batch files and add batch identifier
    adatas = []
    for batch_idx, batch_file in enumerate(batch_files):
        logger.info(f'  Reading: {os.path.basename(batch_file)}')
        adata = ad.read_h5ad(batch_file)
        
        # Extract batch number from filename
        basename = os.path.basename(batch_file)
        try:
            batch_num = basename.split('_batch_')[1].split('_de.h5ad')[0]
            batch_id = f"batch_{batch_num}"
        except (IndexError, ValueError):
            batch_id = f"batch_{batch_idx + 1}"
        
        # Add batch column to obs
        adata.obs['processing_batch'] = batch_id
        adatas.append(adata)
    
    # Concatenate all batches
    logger.info(f'  Concatenating {len(adatas)} batches...')
    aggregated = ad.concat(adatas, merge='same', join='outer', uns_merge='same')
    
    # Log batch information
    if 'processing_batch' in aggregated.obs.columns:
        logger.info(f"  Batch column added with {aggregated.obs['processing_batch'].nunique()}")
    
    # Save aggregated file
    os.makedirs(output_dir, exist_ok=True)
    output_file = join(output_dir, f"{cell_type}_de.h5ad")
    
    # Remove X matrix from aggregated data before saving (to avoid encoding issues)
    if aggregated.X is not None:
        logger.info('  Removing X matrix from aggregated data before saving')
        aggregated.X = None
    
    aggregated.write_h5ad(output_file, compression='gzip')
    
    # Remove X matrix from file using HDF5 access (for compatibility with older anndata versions)
    rewrite_h5ad(output_file)
    
    logger.info(f'  Saved aggregated file: {output_file} ({aggregated.n_obs} obs, {aggregated.n_vars} vars)')


def aggregate_all_cell_types(input_dir: str, output_dir: str) -> None:
    '''
    Aggregate batch files for all cell types in the input directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing cell type subdirectories with batch files
    output_dir : str
        Output directory where aggregated files will be saved
        
    Returns:
    --------
    None
        Saves aggregated files to disk
    '''
    if not os.path.isdir(input_dir):
        logger.warning(f'Input directory not found: {input_dir}')
        return
    
    # Get all cell type directories
    cell_type_dirs = [d for d in os.listdir(input_dir) 
                     if os.path.isdir(join(input_dir, d))]
    
    if len(cell_type_dirs) == 0:
        logger.warning(f'No cell type directories found in {input_dir}')
        return
    
    logger.info(f'Found {len(cell_type_dirs)} cell type(s) to aggregate')
    
    for cell_type in sorted(cell_type_dirs):
        aggregate_cell_type_batches(input_dir, cell_type, output_dir)
    
    logger.info('Aggregation completed')


def main():
    '''
    Main function to aggregate batch files from input directory.
    
    Command Line Arguments:
    -----------------------
    --config : str, optional
        Path to JSON configuration file. Arguments in config file are merged
        with command line arguments (command line takes precedence).
    --input_dir : str, required
        Path to input directory containing cell type subdirectories
    --output_dir : str, required
        Output directory where aggregated files will be saved
        
    Raises:
    -------
    Exception
        If required arguments are not set
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str)
    parser.add_argument('--input_dir', type=str)
    parser.add_argument('--output_dir', type=str)
    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']
    
    if not args.config is None:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)
    
    for key in ['input_dir', 'output_dir']:
        if not d_args.get(key):
            raise Exception(f'The argument {key} is not set')
    
    aggregate_all_cell_types(d_args['input_dir'], d_args['output_dir'])


if __name__ == '__main__':
    main()


