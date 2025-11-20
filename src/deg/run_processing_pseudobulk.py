import os
import json
import time
import fcntl
import argparse
import math
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import Optional, Tuple

from src.utils.parsing_utils import *

#MAX_ELEMENTS 2^31 - 1 from https://svn.r-project.org/R/trunk/src/library/base/R/qr.R
MAX_ELEMENTS = 2147483647
RATIO_MAX_ELEMENTS = 0.9
RANDOM_SEED = 0

def create_perturbation_label(is_control: bool, 
                              pert: str, 
                              dose: float, 
                              time: float,
                              well: str,
                              plate: str,
                              design_param: str) -> str:
    '''
    Create a perturbation label based on design parameters.
    
    Parameters:
    -----------
    is_control : bool
        Whether the perturbation is a control condition
    pert : str
        Name of the perturbagen
    dose : float
        Dose concentration in uM
    time : float
        Time point in hours
    well : str
        Well identifier
    plate : str
        Plate identifier
    design_param : str
        Design parameter determining label format ('group_all_replicates' or 'separate_replicates')
        
    Returns:
    --------
    str
        Formatted perturbation label
    '''
    if is_control:
        return str(pert) \
                + '_' + str(time) + 'h'
    else:
        if design_param == 'group_all_replicates':
            return str(pert) \
                + '_' + str(dose) \
                + 'uM_' + str(time) + 'h'
        elif design_param == 'separate_replicates':
            return str(pert) \
                + '_' + str(dose) \
                + 'uM_' + str(time) + 'h' \
                + '_' + str(well) \
                + '_' + str(plate)


def create_dir_if_not_exists(file_output: str) -> None:
    '''
    Create directory for output file if it doesn't exist.
    
    Parameters:
    -----------
    file_output : str
        Path to the output file
    '''
    dir_output = os.path.dirname(file_output)
    os.makedirs(dir_output, exist_ok=True)

def save_read(file_input: str, 
              n_retries: int = 10,
              delay: float = 5.0) -> ad.AnnData:
    '''
    Read AnnData file with file locking and retry mechanism.
    
    Parameters:
    -----------
    file_input : str
        Path to the input h5ad file
    n_retries : int, default=10
        Number of retry attempts if file is locked
    delay : float, default=5.0
        Delay in seconds between retry attempts
        
    Returns:
    --------
    ad.AnnData
        Loaded AnnData object
        
    Raises:
    -------
    RuntimeError
        If file cannot be read after all retry attempts
    '''

    for attempt in range(n_retries):
        try:
            with open(file_input, 'r') as f:
                fcntl.flock(f, fcntl.LOCK_SH | fcntl.LOCK_NB)
                padata = ad.read_h5ad(file_input)
                fcntl.flock(f, fcntl.LOCK_UN)
                return padata
        except BlockingIOError:
            print(f'Attempt {attempt + 1}: file locked, retrying...')
            time.sleep(delay)
    raise RuntimeError(f'Could not read {file_input} after {retries} attempts.')

def get_output_path_combined(file_input: str,
                        dir_output: str) -> str:
    '''
    Generate output path for processed pseudobulk file.
    
    Parameters:
    -----------
    file_input : str
        Path to the input file
    dir_output : str
        Output directory path
        
    Returns:
    --------
    str
        Full path to the output file
    '''
    dataset_name = os.path.splitext(os.path.basename(file_input))[0]
    return os.path.join(dir_output, f'{dataset_name}_processed.h5ad')


def chunk_celltype_by_matrix_size(adata_ct: ad.AnnData,
                                   dir_output: str,
                                   ct_clean: str) -> None:
    '''
    Check matrix size and split cell type dataset into chunks if needed.
    
    This function checks if the matrix size (n_obs * n_perturbs) exceeds MAX_ELEMENTS.
    If it does, the dataset is recursively split into two parts based on unique
    perturbagens. Controls are excluded from the split and are included in all
    partitions. The function returns the total number of chunks created.
    
    Parameters:
    -----------
    adata_ct : ad.AnnData
        AnnData object for a single cell type
    dir_output : str
        Output directory path
    ct_clean : str
        Sanitized cell type name for filenames
        
    Returns:
    --------
    None
    '''
    import math
    from collections import deque
    
    # Set random seed for reproducible splitting
    np.random.seed(RANDOM_SEED)
    
    # Initialize queue with the full dataset
    # Each element is a tuple: (adata_subset, is_original)
    queue = deque([(adata_ct, True)])
    chunk_counter = 0

    # Separate controls and non-controls for this subset
    if 'is_control' in adata_ct.obs.columns:
        controls_mask = adata_ct.obs['is_control'] == True
        controls = adata_ct[controls_mask].copy() if controls_mask.any() else None
    else:
        raise ValueError(f'  No is_control column found for {ct_clean}')
    
    # Process queue until empty
    while queue:
        adata_subset, is_original = queue.popleft()
        
        n_obs_subset = adata_subset.n_obs
        n_perturbs_subset = adata_subset.obs['perturbation_label'].nunique()
        matrix_elements_subset = n_obs_subset * n_perturbs_subset
        
        # If matrix size is within limits, save as single file
        if matrix_elements_subset <= round(RATIO_MAX_ELEMENTS * MAX_ELEMENTS):
            if not is_original:
                chunk_counter += 1
            save_chunk(adata_subset, dir_output, ct_clean, is_original, chunk_counter, n_obs_subset, n_perturbs_subset)
            continue
        
        # Matrix is too large, need to split into two parts
        logger.info(f'  Splitting {"original" if is_original else "chunk"}: {n_obs_subset} obs × {n_perturbs_subset} perturbs = {matrix_elements_subset} elements (exceeds {MAX_ELEMENTS})')
        
        # Separate controls and non-controls for this subset
        if 'is_control' in adata_subset.obs.columns:
            controls_mask_subset = adata_subset.obs['is_control'] == True
            non_controls_subset = adata_subset[~controls_mask_subset].copy()
        else:
            raise ValueError(f'  No is_control column found for {ct_clean}')
        
        # Get unique perturbagens from non-controls
        unique_perturbagens = non_controls_subset.obs['perturbagen'].unique()
        n_perturbagens = len(unique_perturbagens)
        
        if n_perturbagens < 2:
            # Can't split further by perturbagens, save as is
            logger.warning(f'  Cannot split further by perturbagens, saving as is')
            if not is_original:
                chunk_counter += 1
            save_chunk(adata_subset, dir_output, ct_clean, is_original, chunk_counter, n_obs_subset, n_perturbs_subset)
            continue
        
        # Split perturbagens into two groups and create two parts
        part1, part2 = split_by_perturbagens(
            non_controls_subset,
            unique_perturbagens,
            controls=controls,
        )
        
        # Add both parts back to queue for further processing (both are now split parts, not original)
        queue.append((part1, False))
        queue.append((part2, False))
    
    total_chunks = chunk_counter if chunk_counter > 0 else 1
    logger.info(f'  Total chunks created for {ct_clean}: {total_chunks}')
    return

def save_chunk(adata_subset: ad.AnnData,
               dir_output: str,
               ct_clean: str,
               is_original: bool,
               chunk_counter: int,
               n_obs_subset: int,
               n_perturbs_subset: int) -> None:
    '''
    Save a chunk of the dataset.
    
    Parameters:
    -----------
    adata_subset : ad.AnnData
        Subset of the dataset to save
    dir_output : str
        Output directory path
    ct_clean : str
        Sanitized cell type name for filenames
    is_original : bool
        Whether the dataset is the original dataset
    chunk_counter : int
        Counter for the chunk
    n_obs_subset : int
        Number of observations in the subset
    n_perturbs_subset : int
        Number of perturbagens in the subset

    Returns:
    --------
    None
    '''
    if is_original:
        outfile = os.path.join(dir_output, f"{ct_clean}_processed.h5ad")
    else:
        outfile = os.path.join(dir_output, f"{ct_clean}_processed_chunk_{chunk_counter}.h5ad")
    adata_subset.write_h5ad(outfile, compression='gzip')
    logger.info(f'  Saved chunk {chunk_counter if not is_original else 1}: {outfile}\
     \ ({n_obs_subset} obs, {n_perturbs_subset} perturbs)')


    
def split_by_perturbagens(non_controls: ad.AnnData,
                          unique_perturbagens: np.ndarray,
                          controls: Optional[ad.AnnData] = None,
                          ) -> Tuple[ad.AnnData, ad.AnnData]:
    '''
    Split dataset into two parts based on perturbagens.
    
    Randomly shuffles perturbagens and divides them into two groups,
    then creates two subsets of non-controls based on these groups.
    Adds controls to each part if available.
    
    Parameters:
    -----------
    non_controls : ad.AnnData
        Non-control observations to split
    unique_perturbagens : np.ndarray
        Array of unique perturbagen names
    controls : ad.AnnData, optional
        Original controls to add to each part
        
    Returns:
    --------
    Tuple[ad.AnnData, ad.AnnData]
        Two parts of the dataset, each with controls added if available
    '''
    
    # Split perturbagens into two groups
    shuffled_perturbagens = unique_perturbagens.copy()
    np.random.seed(RANDOM_SEED)
    np.random.shuffle(shuffled_perturbagens)
    
    mid_point = math.ceil(len(shuffled_perturbagens) / 2)
    perturbagens_part1 = shuffled_perturbagens[:mid_point]
    perturbagens_part2 = shuffled_perturbagens[mid_point:]
    
    # Create two subsets based on perturbagens
    part1_non_controls = non_controls[non_controls.obs['perturbagen'].isin(perturbagens_part1)].copy()
    part2_non_controls = non_controls[non_controls.obs['perturbagen'].isin(perturbagens_part2)].copy()
    
    # Add controls to each part
    if controls is not None and controls.n_obs > 0:
        part1 = ad.concat([part1_non_controls, controls], merge='same')
        part2 = ad.concat([part2_non_controls, controls], merge='same')
    else:
        part1 = part1_non_controls
        part2 = part2_non_controls
    
    return part1, part2    


def save_by_celltype(padata: ad.AnnData,
                     dir_output: str) -> None:
    '''
    Save pseudobulk data split by cell type.
    
    Parameters:
    -----------
    padata : ad.AnnData
        Pseudobulk AnnData object
    dir_output : str
        Output directory path
        
    Raises:
    -------
    ValueError
        If no valid cell types are found for splitting
    '''

    celltype_dir = os.path.join(dir_output, 'by_celltype')
    os.makedirs(celltype_dir, exist_ok=True)
    
    # Get unique cell types and filter out None/NaN
    cell_types = padata.obs['cell_type'].unique()
    cell_types = [ct for ct in cell_types if ct is not None and (isinstance(ct, str) or not np.isnan(ct))]
    cell_types = sorted(cell_types)
    
    if len(cell_types) == 0:
        raise ValueError('No valid cell types found for splitting')
    
    logger.info(f'Splitting into {len(cell_types)} cell types')
    
    for ct in cell_types:
        ct_clean = sanitize_celltype_name(ct)
        adata_ct = padata[padata.obs['cell_type'] == ct].copy()
        chunk_celltype_by_matrix_size(adata_ct, celltype_dir, ct_clean)



def add_perturbation_label_to_padata(file_input: str,
                                     dir_output: str,
                                     design_param: str,
                                     split_by_celltype: bool = False,
                                     ) -> None:
    '''
    Add perturbation labels to pseudobulk data and save processed files.
    
    Parameters:
    -----------
    file_input : str
        Path to the input pseudobulk h5ad file
    dir_output : str
        Output directory path
    design_param : str
        Design parameter for perturbation labeling ('group_all_replicates' or 'separate_replicates')
    split_by_celltype : bool, default=False
        Whether to split output by cell type
    '''
    
    logger.info('Read pseudobulk file')

    padata = save_read(file_input)
    obs = padata.obs.copy()
    obs['pert_dose_uM'] = obs['pert_dose_uM'].apply(lambda x: format(x, ".15g"))
    obs['pert_time_h'] = obs['pert_time_h'].apply(lambda x: format(x, ".15g"))
    
    logger.info('Add perturbation label')
    padata.obs['perturbation_label'] =  obs.apply(lambda x: create_perturbation_label(x.is_control,
                                                x.perturbagen,
                                                x.pert_dose_uM,
                                                x.pert_time_h,
                                                x.well,
                                                x.plate,
                                                design_param), axis=1).astype("category")
    
    logger.info('Save data')
    file_output = get_output_path_combined(file_input, 
                                           dir_output)
    create_dir_if_not_exists(file_output)
    padata.write_h5ad(file_output, compression='gzip')
    
    # Optionally split by cell type
    if split_by_celltype:
        save_by_celltype(padata, dir_output)

    


def sanitize_celltype_name(celltype: str) -> str:
    '''
    Sanitize cell type name for use in filenames.
    
    Parameters:
    -----------
    celltype : str
        Original cell type name
        
    Returns:
    --------
    str
        Sanitized cell type name safe for filenames
    '''
    return celltype.replace(' ', '_').replace('/', '-').replace('(', '').replace(')', '')



def main():
    '''
    Main function to process pseudobulk data and add perturbation labels.
    
    Processes command line arguments and configuration files to run the
    pseudobulk processing pipeline with perturbation labeling.
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str)
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--output_dir', type=str)
    parser.add_argument('--design_param', choices=['group_all_replicates',
                                                   'separate_replicates'])
    parser.add_argument('--split_by_celltype', action='store_true', default=False)
    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']

    if not args.config is None:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)

    for key in ['input_file', 'output_dir', 'design_param']:
        if not d_args.get(key):
            raise Exception(f'The argument {key} is not set')

    
    add_perturbation_label_to_padata(d_args['input_file'],
                                     d_args['output_dir'],
                                     d_args['design_param'],
                                     d_args.get('split_by_celltype', False))



if __name__ == '__main__':
    main()
