import os
import json
import time
import fcntl
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import Optional, Tuple, Union

from src.utils.parsing_utils import *

#MAX_ELEMENTS 2^31 - 1 from https://svn.r-project.org/R/trunk/src/library/base/R/qr.R
MAX_ELEMENTS = 2147483647
RATIO_MAX_ELEMENTS = 0.9
RATIO_MIN_ELEMENTS = 0.85
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
            with open(file_input, 'rb') as f:
                fcntl.flock(f, fcntl.LOCK_SH | fcntl.LOCK_NB)
                padata = ad.read_h5ad(file_input)
                fcntl.flock(f, fcntl.LOCK_UN)
                return padata
        except BlockingIOError:
            logger.warning(f'Attempt {attempt + 1}: file locked, retrying...')
            time.sleep(delay)
    raise RuntimeError(f'Could not read {file_input} after {n_retries} attempts.')

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


def calculate_perturbagen_info(non_controls: ad.AnnData,
                               unique_perturbagens: np.ndarray) -> dict:
    '''
    Calculate properties for each perturbagen.
    
    Parameters:
    -----------
    non_controls : ad.AnnData
        Non-control observations
    unique_perturbagens : np.ndarray
        Array of unique perturbagen names
        
    Returns:
    --------
    dict
        Dictionary mapping perturbagen names to their properties:
        - 'n_obs': number of observations
        - 'n_labels': number of unique perturbation labels
        - 'data': AnnData subset for this perturbagen
    '''
    perturbagen_info = {}
    for pert in unique_perturbagens:
        pert_data = non_controls[non_controls.obs['perturbagen'] == pert]
        perturbagen_info[pert] = {
            'n_obs': pert_data.n_obs,
            'n_labels': pert_data.obs['perturbation_label'].nunique(),
            'data': pert_data
        }
    return perturbagen_info


def get_control_info(controls: Optional[ad.AnnData]) -> Tuple[int, int]:
    '''
    Extract control observation and label counts.
    
    Parameters:
    -----------
    controls : Optional[ad.AnnData]
        Control observations
        
    Returns:
    --------
    Tuple[int, int]
        (ctl_n_obs, ctl_n_labels)
    '''
    if (controls is not None) and (controls.n_obs > 0):
        return controls.n_obs, controls.obs['perturbation_label'].nunique()
    return 0, 0


def get_target_sizes() -> Tuple[int, int]:
    '''
    Calculate target matrix sizes.
    
    Returns:
    --------
    Tuple[int, int]
        (target_min, target_max)
    '''
    return (round(RATIO_MIN_ELEMENTS * MAX_ELEMENTS),
            round(RATIO_MAX_ELEMENTS * MAX_ELEMENTS))


def calculate_matrix_size(n_obs: int, n_labels: int,
                         ctl_n_obs: int = 0, ctl_n_labels: int = 0) -> int:
    '''
    Calculate matrix size including controls.
    
    Parameters:
    -----------
    n_obs : int
        Number of observations (non-controls)
    n_labels : int
        Number of unique perturbation labels (non-controls)
    ctl_n_obs : int
        Number of control observations
    ctl_n_labels : int
        Number of unique control labels
        
    Returns:
    --------
    int
        Matrix size: (n_obs + ctl_n_obs) * (n_labels + ctl_n_labels)
    '''
    return (n_obs + ctl_n_obs) * (n_labels + ctl_n_labels)


def separate_controls(adata: ad.AnnData, ct_clean: str = None) -> Tuple[Optional[ad.AnnData], ad.AnnData]:
    '''
    Separate controls and non-controls from AnnData.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object with 'is_control' column in obs
    ct_clean : str, optional
        Cell type name for error messages
        
    Returns:
    --------
    Tuple[Optional[ad.AnnData], ad.AnnData]
        (controls, non_controls)
        
    Raises:
    -------
    ValueError
        If 'is_control' column is not found
    '''
    if 'is_control' not in adata.obs.columns:
        error_msg = f'  No is_control column found'
        if ct_clean:
            error_msg += f' for {ct_clean}'
        raise ValueError(error_msg)
    
    controls_mask = adata.obs['is_control'] == True
    controls = adata[controls_mask].copy() if controls_mask.any() else None
    non_controls = adata[~controls_mask].copy()
    return controls, non_controls


def get_perturbagen_properties(perturbagen_info: dict, pert: str) -> Tuple[int, int, ad.AnnData]:
    '''
    Extract properties for a perturbagen from perturbagen_info.
    
    Parameters:
    -----------
    perturbagen_info : dict
        Dictionary mapping perturbagen names to their properties
    pert : str
        Perturbagen name
        
    Returns:
    --------
    Tuple[int, int, ad.AnnData]
        (pert_n_obs, pert_n_labels, pert_data)
    '''
    pert_info = perturbagen_info[pert]
    return pert_info['n_obs'], pert_info['n_labels'], pert_info['data']


def calculate_new_batch_size(batch: dict, pert_n_obs: int, pert_n_labels: int,
                            ctl_n_obs: int, ctl_n_labels: int) -> Tuple[int, int, int]:
    '''
    Calculate what the batch size would be if we add a perturbagen (without modifying batch).
    
    Parameters:
    -----------
    batch : dict
        Batch dictionary with 'n_obs' and 'n_labels'
    pert_n_obs : int
        Number of observations in the perturbagen
    pert_n_labels : int
        Number of unique labels in the perturbagen
    ctl_n_obs : int
        Number of control observations
    ctl_n_labels : int
        Number of unique control labels
        
    Returns:
    --------
    Tuple[int, int, int]
        (new_n_obs, new_n_labels, new_matrix_size)
    '''
    new_n_obs = batch['n_obs'] + pert_n_obs
    new_n_labels = batch['n_labels'] + pert_n_labels
    new_matrix_size = calculate_matrix_size(new_n_obs, 
                                            new_n_labels, 
                                            ctl_n_obs=ctl_n_obs, 
                                            ctl_n_labels=ctl_n_labels)
    return new_n_obs, new_n_labels, new_matrix_size


def add_data_to_batch(batch: dict, n_obs: int, n_labels: int, indices: list, pert: Optional[Union[str, set, list]] = None) -> None:
    '''
    Add observations to a batch (modifies batch in place).
    
    Can be used for both perturbagens and controls. If pert is provided,
    it will be added to batch['perturbagens'] set. Can be a single name (str)
    or multiple names (set/list) for controls.
    
    Parameters:
    -----------
    batch : dict
        Batch dictionary to modify
    n_obs : int
        Number of observations to add
    n_labels : int
        Number of unique perturbation labels to add
    indices : list
        List of observation indices to add to batch['obs_indices']
    pert : str, set, or list, optional
        Perturbagen/control name(s). If provided, adds to batch['perturbagens'] set.
        Can be a single string, or a set/list of strings for multiple names.
        If None, no names are added to batch['perturbagens'].
    '''
    if pert is not None:
        if isinstance(pert, str):
            batch['perturbagens'].add(pert)
        else:
            batch['perturbagens'].update(pert)
    batch['obs_indices'].extend(indices)
    batch['n_obs'] += n_obs
    batch['n_labels'] += n_labels


def create_batches(non_controls: ad.AnnData,
                    unique_perturbagens: np.ndarray,
                    controls: Optional[ad.AnnData] = None,
                    ) -> list:
    '''
    Create multiple batches from perturbagens using Best Fit algorithm.
    
    Randomly shuffles perturbagens and assigns them to batches using Best Fit algorithm,
    targeting 85-95% of MAX_ELEMENTS per batch. The algorithm selects the batch that
    minimizes waste (gets closest to target_max without exceeding it) for each perturbagen.
    Controls are considered in matrix size calculations but not added to batches yet.
    
    Parameters:
    -----------
    non_controls : ad.AnnData
        Non-control observations to split
    unique_perturbagens : np.ndarray
        Array of unique perturbagen names
    controls : ad.AnnData, optional
        Controls to consider in matrix size calculations (not added to batches here)
        
    Returns:
    --------
    list
        List of batch dictionaries, each containing:
        - 'perturbagens': set of perturbagen names in the batch
        - 'n_obs': number of observations (non-controls only)
        - 'n_labels': number of unique perturbation labels (non-controls only)
        - 'obs_indices': list of observation indices from non_controls
    '''
    
    perturbagen_info = calculate_perturbagen_info(non_controls, unique_perturbagens)
    
    shuffled_perturbagens = list(unique_perturbagens.copy())
    np.random.seed(RANDOM_SEED)
    np.random.shuffle(shuffled_perturbagens)
    
    batches = []
    target_min, target_max = get_target_sizes()
    ctl_n_obs, ctl_n_labels = get_control_info(controls)
    
    for pert in shuffled_perturbagens:
        pert_n_obs, pert_n_labels, pert_data = get_perturbagen_properties(perturbagen_info, pert)
        
        best_batch_idx = None
        best_waste = float('inf')
        
        for batch_idx, batch in enumerate(batches):
            _, _, new_matrix_size = calculate_new_batch_size(batch, pert_n_obs, pert_n_labels, ctl_n_obs, ctl_n_labels)
            
            if new_matrix_size <= target_max:
                waste = target_max - new_matrix_size
                if waste < best_waste:
                    best_waste = waste
                    best_batch_idx = batch_idx
        
        if best_batch_idx is not None:
            batch = batches[best_batch_idx]
            add_data_to_batch(batch, pert_n_obs, pert_n_labels, pert_data.obs.index.tolist(), pert=pert)
        else:
            new_batch = {
                'perturbagens': {pert},
                'n_obs': pert_n_obs,
                'n_labels': pert_n_labels,
                'obs_indices': pert_data.obs.index.tolist()
            }
            batches.append(new_batch)
    

    return batches
    
def pad_batches(batches: list,
                non_controls: ad.AnnData,
                unique_perturbagens: np.ndarray,
                controls: Optional[ad.AnnData] = None,
                ) -> list:
    '''
    Pad smaller batches to balance sizes using perturbagens not in current batch.
    
    For batches with matrix size below target_min, this function iteratively adds
    perturbagens from other batches (not already in the current batch). Perturbagens
    are shuffled and added one at a time, checking that the matrix size doesn't exceed
    target_max. The process stops when the batch reaches target_matrix_size or no more
    suitable perturbagens are available.
    
    Parameters:
    -----------
    batches : list
        List of batch dictionaries from create_batches()
    non_controls : ad.AnnData
        Non-control observations (used to get perturbagen data)
    unique_perturbagens : np.ndarray
        Array of unique perturbagen names
    controls : ad.AnnData, optional
        Controls to consider in matrix size calculations (not added here)
        
    Returns:
    --------
    list
        List of batch dictionaries with updated 'obs_indices', 'n_obs', and 'n_labels'
    '''
    
    perturbagen_info = calculate_perturbagen_info(non_controls, unique_perturbagens)
    target_min, target_max = get_target_sizes()
    ctl_n_obs, ctl_n_labels = get_control_info(controls)
    
    if len(batches) > 1:
        target_matrix_size = target_min
        
        for batch_idx, batch in enumerate(batches):
            matrix_size = calculate_matrix_size(batch['n_obs'], 
                                                batch['n_labels'], 
                                                ctl_n_obs=ctl_n_obs, 
                                                ctl_n_labels=ctl_n_labels)
            if matrix_size < target_matrix_size:
                current_batch_perts = batch['perturbagens']
                all_perts = set(non_controls.obs['perturbagen'].unique())
                available_perts = list(all_perts - current_batch_perts)
                
                if len(available_perts) == 0:
                    continue
                
                np.random.seed(RANDOM_SEED)
                shuffled_perts = available_perts.copy()
                np.random.shuffle(shuffled_perts)
                
                for pert in shuffled_perts:
                    pert_n_obs, pert_n_labels, pert_data = get_perturbagen_properties(perturbagen_info, pert)
                    
                    _, _, new_matrix_size = calculate_new_batch_size(batch, pert_n_obs, pert_n_labels, ctl_n_obs, ctl_n_labels)
                    
                    if new_matrix_size <= target_max:
                        add_data_to_batch(batch, pert_n_obs, pert_n_labels, pert_data.obs.index.tolist(), pert=pert)
                        
                        current_matrix_size = calculate_matrix_size(batch['n_obs'], 
                                                                   batch['n_labels'], 
                                                                   ctl_n_obs=ctl_n_obs, 
                                                                   ctl_n_labels=ctl_n_labels)
                        if current_matrix_size >= target_matrix_size:
                            break
    
    # Log final batch statistics after padding
    if len(batches) > 0:
        batch_sizes = [calculate_matrix_size(batch['n_obs'], 
                                             batch['n_labels'], 
                                             ctl_n_obs=ctl_n_obs, 
                                             ctl_n_labels=ctl_n_labels) for batch in batches]
        logger.info(f'  After padding: {len(batches)} batches, matrix sizes: {[f"{s/1e6:.1f}M" for s in batch_sizes]}')
    
    return batches

def add_controls_to_batches(batches: list,
                            controls: Optional[ad.AnnData],
                            ) -> list:
    '''
    Add control observations to all batches.
    
    Appends control observation indices to each batch's obs_indices list and updates
    the batch's n_obs and n_labels counts. If controls are None or empty, returns
    batches unchanged.
    
    Parameters:
    -----------
    batches : list
        List of batch dictionaries to add controls to
    controls : Optional[ad.AnnData]
        Control observations to add to all batches. If None or empty, batches are
        returned unchanged.
        
    Returns:
    --------
    list
        List of batch dictionaries with controls added to 'obs_indices', 'n_obs', and 'n_labels'
    '''
    if controls is None or controls.n_obs == 0:
        return batches
    ctl_n_obs = controls.n_obs
    ctl_n_labels = controls.obs['perturbation_label'].nunique()
    ctl_indices = controls.obs.index.tolist()
    ctl_names = None
    if 'perturbagen' in controls.obs.columns:
        ctl_names = set(controls.obs['perturbagen'].unique())
    for batch in batches:
        add_data_to_batch(batch, ctl_n_obs, ctl_n_labels, ctl_indices, pert=ctl_names)
    return batches





def save_single_file(adata: ad.AnnData, file_output: str) -> None:
    '''
    Save AnnData as a single processed h5ad file.
    
    Creates the output directory if it doesn't exist.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object to save
    file_output : str
        Full path to the output file
    '''
    create_dir_if_not_exists(file_output)
    adata.write_h5ad(file_output, compression='gzip')
    logger.info(f'  Saved: {file_output} ({adata.n_obs} obs)')


def check_and_save_if_small(adata_ct: ad.AnnData, dir_output: str, ct_clean: str) -> bool:
    '''
    Check if matrix size is within limits and save as single file if so.
    
    Parameters:
    -----------
    adata_ct : ad.AnnData
        AnnData object to check and potentially save
    dir_output : str
        Output directory path
    ct_clean : str
        Sanitized cell type name for filenames
        
    Returns:
    --------
    bool
        True if saved (matrix was small enough), False otherwise
    '''
    n_obs = adata_ct.n_obs
    n_labels = adata_ct.obs['perturbation_label'].nunique()
    matrix_size = calculate_matrix_size(n_obs, n_labels)
    
    if matrix_size <= round(RATIO_MAX_ELEMENTS * MAX_ELEMENTS):
        file_output = os.path.join(dir_output, f"{ct_clean}_processed.h5ad")
        save_single_file(adata_ct, file_output)
        return True
    return False


def get_unique_perturbagens_or_save(non_controls: ad.AnnData, adata_ct: ad.AnnData,
                                    dir_output: str, ct_clean: str) -> Optional[np.ndarray]:
    '''
    Get unique perturbagens from non-controls, or save and return None if none found.
    
    Parameters:
    -----------
    non_controls : ad.AnnData
        Non-control observations
    adata_ct : ad.AnnData
        Full AnnData object (for saving if no perturbagens)
    dir_output : str
        Output directory path
    ct_clean : str
        Sanitized cell type name for filenames
        
    Returns:
    --------
    Optional[np.ndarray]
        Array of unique perturbagens, or None if none found (and saved)
    '''
    unique_perturbagens = non_controls.obs['perturbagen'].unique()
    
    if len(unique_perturbagens) < 1:
        logger.warning(f'  No perturbagens found, saving as is')
        file_output = os.path.join(dir_output, f"{ct_clean}_processed.h5ad")
        save_single_file(adata_ct, file_output)
        return None
    
    return unique_perturbagens


def save_batches(batches: list, adata_ct: ad.AnnData, dir_output: str, ct_clean: str) -> None:
    '''
    Save batches to disk as separate h5ad files.
    
    Parameters:
    -----------
    batches : list
        List of batch dictionaries with 'obs_indices' keys
    adata_ct : ad.AnnData
        Full AnnData object to select observations from
    dir_output : str
        Output directory path
    ct_clean : str
        Sanitized cell type name for filenames
        
    Returns:
    --------
    None
        Saves batch files to disk
    '''
    batch_counter = 0
    for batch_idx, batch in enumerate(batches):
        batch_counter += 1
        batch_data = adata_ct[batch['obs_indices']].copy()
        outfile = os.path.join(dir_output, f"{ct_clean}_processed_batch_{batch_counter}.h5ad")
        save_single_file(batch_data, outfile)
        logger.info(f'  Saved batch {batch_counter}: {outfile} ({batch_data.n_obs} obs, {batch_data.obs["perturbation_label"].nunique()} perturbs)')
    
    logger.info(f'  Total batches created for {ct_clean}: {batch_counter}')


def batch_celltype_by_matrix_size(adata_ct: ad.AnnData,
                                   dir_output: str,
                                   ct_clean: str) -> None:
    '''
    Check matrix size and split cell type dataset into batches if needed.
    
    This function checks if the matrix size (n_obs * n_perturbs) exceeds MAX_ELEMENTS.
    If it does, the dataset is split into multiple batches using Best Fit algorithm
    with random shuffling. Controls are excluded from the split and are included in all
    partitions. Smaller batches are padded to balance sizes.
    
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
    if check_and_save_if_small(adata_ct, dir_output, ct_clean):
        return
    
    # Log splitting information
    n_obs = adata_ct.n_obs
    n_labels = adata_ct.obs['perturbation_label'].nunique()
    matrix_size = calculate_matrix_size(n_obs, n_labels)
    logger.info(f'  Splitting {ct_clean}: {n_obs} obs, {n_labels} perturbation labels, {matrix_size/1e6:.1f}M matrix size')
    logger.warning(f'  Matrix size for {ct_clean} exceeds {round(RATIO_MAX_ELEMENTS * MAX_ELEMENTS)}, creating batches')
    
    controls, non_controls = separate_controls(adata_ct, ct_clean)
    
    unique_perturbagens = get_unique_perturbagens_or_save(non_controls, adata_ct, dir_output, ct_clean)
    if unique_perturbagens is None:
        return
    
    n_perturbagens = len(unique_perturbagens)
    n_non_control_obs = non_controls.n_obs
    n_non_control_labels = non_controls.obs['perturbation_label'].nunique()
    logger.info(f'  After separation: {n_perturbagens} perturbagens, {n_non_control_obs} non-control obs, {n_non_control_labels} non-control labels')
    
    logger.info('Create batches')
    batches = create_batches(non_controls,
                          unique_perturbagens,
                          controls=controls,
                          )

    logger.info('Pad batches')
    batches = pad_batches(batches,
                          non_controls,
                          unique_perturbagens,
                          controls=controls,
                          )
    
    logger.info('Add controls to batches')
    batches = add_controls_to_batches(batches,
                                      controls,
                                      )
    
    logger.info('Save batches')
    save_batches(batches, adata_ct, dir_output, ct_clean)
    return


def save_by_celltype(padata: ad.AnnData,
                     dir_output: str) -> None:
    '''
    Save pseudobulk data split by cell type.
    
    Splits the pseudobulk data by cell type and processes each cell type separately.
    For each cell type, checks if matrix size exceeds limits and batches if necessary.
    Files are saved in dir_output/by_celltype/ directory.
    
    Parameters:
    -----------
    padata : ad.AnnData
        Pseudobulk AnnData object
    dir_output : str
        Output directory path. Cell type files are saved in dir_output/by_celltype/
        
    Returns:
    --------
    None
        Saves processed files to disk. Each cell type is saved as
        '{celltype_clean}_processed.h5ad' or batched files if matrix size is too large.
        
    Raises:
    -------
    ValueError
        If no valid cell types are found for splitting
    '''

    celltype_dir = os.path.join(dir_output, 'by_celltype')
    os.makedirs(celltype_dir, exist_ok=True)
    
    cell_types = padata.obs['cell_type'].unique()
    cell_types = [ct for ct in cell_types if ct is not None and (isinstance(ct, str) or not np.isnan(ct))]
    cell_types = sorted(cell_types)
    
    if len(cell_types) == 0:
        raise ValueError('No valid cell types found for splitting')
    
    logger.info(f'Splitting into {len(cell_types)} cell types')
    
    for ct in cell_types:
        ct_clean = sanitize_celltype_name(ct)
        adata_ct = padata[padata.obs['cell_type'] == ct].copy()
        batch_celltype_by_matrix_size(adata_ct, celltype_dir, ct_clean)



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
        If True, splits output by cell type and processes each separately.
        Large cell types are automatically batched if matrix size exceeds limits.
        
    Returns:
    --------
    None
        Saves processed files to disk. Main output file is saved as
        '{basename}_processed.h5ad' in dir_output. If split_by_celltype is True,
        additional files are saved in dir_output/by_celltype/.
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
    file_output = get_output_path_combined(file_input, dir_output)
    save_single_file(padata, file_output)
    
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
    pseudobulk processing pipeline with perturbation labeling. Supports both
    command line arguments and JSON configuration files.
    
    Command Line Arguments:
    -----------------------
    --config : str, optional
        Path to JSON configuration file. Arguments in config file are merged
        with command line arguments (command line takes precedence).
    --input_file : str, required
        Path to the input pseudobulk h5ad file
    --output_dir : str, required
        Output directory path for processed files
    --design_param : str, required
        Design parameter for perturbation labeling. Must be one of:
        - 'group_all_replicates': Groups all replicates together
        - 'separate_replicates': Keeps replicates separate with well/plate info
    --split_by_celltype : bool, optional
        If True, splits output by cell type and processes each separately.
        Default: False
        
    Raises:
    -------
    Exception
        If required arguments (input_file, output_dir, design_param) are not set
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
