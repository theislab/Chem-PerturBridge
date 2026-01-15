import re
import time
import json
import os
import pandas as pd
from typing import Any, Optional, Dict, Callable
import pubchempy as pcp
from rdkit import Chem

from src.utils.parsing_utils import *


def is_valid_pubchem_cid(cid) -> bool:
    """Check if pubchem_cid is valid (positive integer)."""
    if cid is None:
        return False
    try:
        return int(cid) > 0
    except (ValueError, TypeError):
        return False


def is_valid_inchikey(inchikey: Optional[str]) -> bool:
    """Check if InChIKey follows standard format (XXXXXXXXXXXXXX-XXXXXXXXXX-X)."""
    if not inchikey or not isinstance(inchikey, str):
        return False
    pattern = r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'
    return bool(re.match(pattern, inchikey.strip()))


def is_valid_smiles(smiles: str) -> bool:
    """Basic SMILES validation (non-empty string with common SMILES characters)."""
    return Chem.MolFromSmiles(smiles, sanitize=True) is not None


def load_cache_from_json(cache_path: str) -> Dict[str, Optional[int]]:
    """
    Load cache dictionary from JSON file if it exists.
    
    Parameters:
    -----------
    cache_path : str
        Path to the JSON file containing the cache
        
    Returns:
    --------
    Dict[str, Optional[int]]
        Cache dictionary loaded from file, or empty dict if file doesn't exist
    """
    if os.path.exists(cache_path):
        try:
            with open(cache_path, 'r') as f:
                cache = json.load(f)
            # Convert loaded JSON values to Optional[int] type
            # JSON stores None as null, which becomes None in Python
            # Integer values are stored as numbers in JSON, which become Python int when loaded
            cache_typed = {}
            for key, value in cache.items():
                cache_typed[key] = int(value) if value is not None else None
            logger.info(f"Loaded cache from {cache_path} with {len(cache_typed)} entries")
            return cache_typed
        except Exception as e:
            logger.warning(f"Failed to load cache from {cache_path}: {e}. Starting with empty cache.")
            return {}
    else:
        logger.info(f"Cache file {cache_path} not found. Starting with empty cache.")
        return {}


def save_cache_to_json(cache: Dict[str, Optional[int]], cache_path: str) -> None:
    """
    Save cache dictionary to JSON file.
    
    Parameters:
    -----------
    cache : Dict[str, Optional[int]]
        Cache dictionary to save
    cache_path : str
        Path to save the JSON file
    """
    try:
        # Create directory if it doesn't exist
        cache_dir = os.path.dirname(cache_path)
        if cache_dir:  # Only create directory if path contains a directory
            os.makedirs(cache_dir, exist_ok=True)
        
        # Save cache directly to JSON (cache is already JSON-serializable: Dict[str, Optional[int]])
        # None values are preserved as null, integers are stored as numbers (JSON natively supports integers)
        with open(cache_path, 'w') as f:
            json.dump(cache, f, indent=2)
        logger.debug(f"Saved cache to {cache_path} with {len(cache)} entries")
    except Exception as e:
        logger.warning(f"Failed to save cache to {cache_path}: {e}")





def _fetch_pubchem_cid_with_retry(identifier: str,
                                  lookup_type: str,
                                  cache: Dict[str, Optional[int]],
                                  cache_key: str,
                                  identifier_label: str,
                                  universal_cache_key: Optional[str] = None,
                                  n_retries: int = 5) -> Optional[int]:
    """
    Helper function to fetch PubChem CID with retry logic.
    
    Parameters:
    -----------
    identifier : str
        The identifier to look up (InChIKey, SMILES, or drug name)
    lookup_type : str
        PubChem lookup type: 'inchikey', 'smiles', or 'name'
    cache : Dict[str, Optional[int]]
        Cache dictionary for storing results
    cache_key : str
        Method-specific key to use in cache dictionary
    identifier_label : str
        Label for logging (e.g., 'InChIKey', 'SMILES', 'drug name')
    universal_cache_key : Optional[str]
        Universal cache key (e.g., perturbagen name) to store result under
    n_retries : int
        Number of retries on failure
        
    Returns:
    --------
    Optional[int]
        PubChem CID or None if not found
    """
    cnt = 0
    cid = None
    while cnt < n_retries:
        try:
            compounds = pcp.get_compounds(identifier, lookup_type)
            cid = compounds[0].cid if compounds else None
            logger.debug("CID for %s '%s': %s", identifier_label, identifier, cid)
            break
        except Exception as e:
            if (isinstance(e, pcp.PubChemHTTPError) or isinstance(e, pcp.TimeoutError) or
                isinstance(e, pcp.ServerError) or isinstance(e, pcp.ServerBusyError)):
                logger.warning("PubChem lookup failed for %s '%s': %s. Retry %d.", 
                             identifier_label, identifier, str(e), cnt)
                cnt += 1
                time.sleep(5)
            else:
                logger.warning("PubChem lookup failed for %s '%s': %s", 
                             identifier_label, identifier, str(e))
                break
    
    # Store in both method-specific and universal cache keys (only if cid is not None)
    if cid is not None:
        cache[cache_key] = cid
        if universal_cache_key:
            cache[universal_cache_key] = cid
    
    return cid


def get_pubchem_cid_by_inchikey(inchikey: str, 
                                 cache: Dict[str, Optional[int]], 
                                 universal_cache_key: Optional[str] = None,
                                 n_retries: int = 5) -> Optional[int]:
    """
    Fetch PubChem CID for a given InChIKey, using cache to skip repeat lookups.
    
    Parameters:
    -----------
    inchikey : str
        InChIKey for mapping to PubChem CID
    cache : Dict[str, Optional[int]]
        A dictionary storing the mapping of identifiers to PubChem CIDs
    universal_cache_key : Optional[str]
        Universal cache key (e.g., perturbagen name) to check/store result
    n_retries : int
        The number of retries when connecting to PubChem goes wrong
        
    Returns:
    --------
    Optional[int]
        PubChem CID or None if not found
    """
    if not is_valid_inchikey(inchikey):
        return None
    
    # Check universal cache first
    if universal_cache_key and universal_cache_key in cache:
        return cache[universal_cache_key]
    
    # Check method-specific cache
    cache_key = f"inchikey:{inchikey}"
    if cache_key in cache:
        return cache[cache_key]
    
    return _fetch_pubchem_cid_with_retry(inchikey, 'inchikey', cache, cache_key, 'InChIKey', 
                                        universal_cache_key, n_retries)


def get_pubchem_cid_by_smiles(smiles: str, 
                              cache: Dict[str, Optional[int]], 
                              universal_cache_key: Optional[str] = None,
                              n_retries: int = 5) -> Optional[int]:
    """
    Fetch PubChem CID for a given SMILES string, using cache to skip repeat lookups.
    
    Parameters:
    -----------
    smiles : str
        SMILES string for mapping to PubChem CID
    cache : Dict[str, Optional[int]]
        A dictionary storing the mapping of identifiers to PubChem CIDs
    universal_cache_key : Optional[str]
        Universal cache key (e.g., perturbagen name) to check/store result
    n_retries : int
        The number of retries when connecting to PubChem goes wrong
        
    Returns:
    --------
    Optional[int]
        PubChem CID or None if not found
    """
    if not is_valid_smiles(smiles):
        return None
    
    # Check universal cache first
    if universal_cache_key and universal_cache_key in cache:
        return cache[universal_cache_key]
    
    # Check method-specific cache
    cache_key = f"smiles:{smiles}"
    if cache_key in cache:
        return cache[cache_key]
    
    return _fetch_pubchem_cid_with_retry(smiles, 'smiles', cache, cache_key, 'SMILES', 
                                        universal_cache_key, n_retries)


def get_pubchem_cid_by_name(drug_name: str, 
                            cache: Dict[str, Optional[int]], 
                            universal_cache_key: Optional[str] = None,
                            n_retries: int = 5) -> Optional[int]:
    """
    Fetch PubChem CID for a given drug name, using cache to skip repeat lookups.
    
    Parameters:
    -----------
    drug_name : str
        Drug name for mapping to PubChem CID
    cache : Dict[str, Optional[int]]
        A dictionary storing the mapping of identifiers to PubChem CIDs
    universal_cache_key : Optional[str]
        Universal cache key (e.g., perturbagen name) to check/store result.
        If None, uses drug_name as universal key.
    n_retries : int
        The number of retries when connecting to PubChem goes wrong
        
    Returns:
    --------
    Optional[int]
        PubChem CID or None if not found
    """
    if pd.isna(drug_name) or not drug_name:
        return None
    
    # Use drug_name as universal key if not provided
    if universal_cache_key is None:
        universal_cache_key = drug_name
    
    # Check universal cache first
    if universal_cache_key in cache:
        return cache[universal_cache_key]
    
    # Check method-specific cache
    cache_key = f"name:{drug_name}"
    if cache_key in cache:
        return cache[cache_key]
    
    return _fetch_pubchem_cid_with_retry(drug_name, 'name', cache, cache_key, 'drug name', 
                                        universal_cache_key, n_retries)


def lookup_pubchem_cids(df: pd.DataFrame,
                           cache: Dict[str, Optional[int]],
                           pert_id_col: Optional[str] = 'pert_id',
                           drug_col: str = 'perturbagen',
                           pubchem_cid_col: str = 'pubchem_cid',
                           inchikey_col: str = 'inchi_key',
                           smiles_col: str = 'canonical_smiles',
                           cache_path: Optional[str] = None,
                           manual_mapping_func: Optional[Callable[[], Dict]] = None,
                           dataset_key: Optional[str] = None) -> pd.DataFrame:
    """
    Add or update 'pubchem_cid' column in a dataframe.
    
    Uses multiple strategies in order of preference:
    1. Use existing valid pubchem_cid if present
    2. Lookup by InChIKey if available and valid
    3. Lookup by SMILES if available and valid
    4. Lookup by drug name (perturbagen)
    5. Use manual mapping from manual_mapping_func (if provided)
    
    Uses universal cache keys (pert_id or perturbagen) to avoid redundant lookups
    when the same compound is identified by different methods.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame to add pubchem_cid column to
    cache : Dict[str, Optional[int]]
        A dictionary storing the mapping of identifiers to PubChem CIDs.
        Will be updated with new entries during processing.
    pert_id_col : Optional[str], default='pert_id'
        Column name for perturbation ID to use as universal cache key.
        If None, uses drug_col (perturbagen) as universal cache key.
    drug_col : str, default='perturbagen'
        Column name containing drug names (used as universal cache key if pert_id_col is None)
    pubchem_cid_col : str, default='pubchem_cid'
        Column name for PubChem CID (will be created/updated)
    inchikey_col : str, default='inchi_key'
        Column name for InChIKey
    smiles_col : str, default='canonical_smiles'
        Column name for SMILES
    cache_path : Optional[str], default=None
        Path to JSON file for persistent cache storage.
        If provided, cache will be loaded from this file at start and saved periodically.
    manual_mapping_func : Optional[Callable[[], Dict]], default=None
        Function that returns manual PubChem CID mappings. If provided, should return
        either a dict directly (e.g., {'drug_name': cid}) or a dict with dataset keys
        (e.g., {'dataset_name': {'drug_name': cid}}). If None, no manual mapping is used.
    dataset_key : Optional[str], default=None
        Key to extract from the dict returned by manual_mapping_func if it returns
        a nested dict structure. If None and manual_mapping_func returns a nested dict,
        defaults to 'l1000' for backward compatibility.
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with updated pubchem_cid column
    """
    df = df.copy()

    if df is None:
        raise Exception("The pseudobulk dataset is empty")
    
    
    # Load cache from file if path is provided
    if cache_path:
        file_cache = load_cache_from_json(cache_path)
        # Merge file cache with provided cache (provided cache takes precedence)
        cache.update(file_cache)
    
    # Initialize pubchem_cid column if it doesn't exist
    if pubchem_cid_col not in df.columns:
        df[pubchem_cid_col] = None
    
    # Get manual mappings
    if manual_mapping_func is not None:
        sm2pubchem = manual_mapping_func()
        # Handle both flat dicts and nested dicts
        if isinstance(sm2pubchem, dict):
            if dataset_key is not None:
                # Extract specific key from nested dict
                manual_mapping = sm2pubchem.get(dataset_key, {})
            elif len(sm2pubchem) == 1:
                # If single key, use it automatically
                manual_mapping = list(sm2pubchem.values())[0]
            else:
                # Check if it's a nested dict (values are dicts) or flat dict (values are ints)
                first_value = list(sm2pubchem.values())[0] if sm2pubchem else None
                if isinstance(first_value, dict):
                    # Nested dict but no dataset_key specified - use first key as fallback
                    manual_mapping = first_value
                else:
                    # Flat dict with drug_name: cid mappings
                    manual_mapping = sm2pubchem
        else:
            manual_mapping = {}
    else:
        manual_mapping = {}
    
    # Process each row to determine best CID source
    cids = []
    iteration_count = 0
    n_compounds = len(df)
    logger.info(f"Processing {n_compounds} compounds")
    for idx, row in df.iterrows():
        iteration_count += 1
        cid = None
        
        # Determine universal cache key (pert_id or perturbagen)
        if pert_id_col and pert_id_col in row and pd.notna(row[pert_id_col]):
            universal_key = str(row[pert_id_col]).strip()
        elif drug_col in row and pd.notna(row[drug_col]):
            universal_key = str(row[drug_col]).strip()
        else:
            universal_key = None
        
        # Check universal cache first
        if universal_key and universal_key in cache:
            cid = cache[universal_key]
        
        # Strategy 1: Use existing valid pubchem_cid
        if cid is None and pubchem_cid_col in row and is_valid_pubchem_cid(row[pubchem_cid_col]):
            cid = int(row[pubchem_cid_col])
            # Store in universal cache
            if universal_key:
                cache[universal_key] = cid
        
        # Strategy 2: Lookup by InChIKey (if CID not found yet)
        if cid is None and inchikey_col in row and pd.notna(row[inchikey_col]):
            inchikey = str(row[inchikey_col]).strip()
            if is_valid_inchikey(inchikey):
                cid = get_pubchem_cid_by_inchikey(inchikey, cache, universal_key)
        
        # Strategy 3: Lookup by SMILES (if CID not found yet)
        if cid is None and smiles_col in row and pd.notna(row[smiles_col]):
            smiles = str(row[smiles_col]).strip()
            if is_valid_smiles(smiles):
                cid = get_pubchem_cid_by_smiles(smiles, cache, universal_key)
        
        # Strategy 4: Lookup by drug name (if CID not found yet)
        if cid is None and drug_col in row and pd.notna(row[drug_col]):
            drug_name = str(row[drug_col]).strip()
            cid = get_pubchem_cid_by_name(drug_name, cache, universal_key)
        
        # Strategy 5: Manual mapping (if CID not found yet)
        if cid is None and drug_col in row and pd.notna(row[drug_col]):
            drug_name = str(row[drug_col]).strip()
            if drug_name in manual_mapping:
                cid = manual_mapping[drug_name]
                # Store in universal cache
                if universal_key:
                    cache[universal_key] = cid
        
        cids.append(cid)
        
        # Log progress every 50 compounds
        if iteration_count % 50 == 0:
            n_mapped_so_far = sum(1 for c in cids if c is not None)
            logger.info(f"Processed {iteration_count}/{n_compounds} compounds ({n_mapped_so_far} mapped so far)")
        
        # Save cache every 500 iterations if cache_path is provided
        if cache_path and iteration_count % 500 == 0:
            save_cache_to_json(cache, cache_path)
    
    # Final save of cache if cache_path is provided
    if cache_path:
        save_cache_to_json(cache, cache_path)
    
    # Update pubchem_cid column
    df_updated = df.copy()
    df_updated[pubchem_cid_col] = cids
    
    n_mapped = df_updated[pubchem_cid_col].notna().sum()
    logger.info(f"Mapped {n_mapped} out of {len(df)} compounds to PubChem CIDs")
    
    return df_updated
