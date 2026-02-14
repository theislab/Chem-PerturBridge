from __future__ import annotations

from typing import Optional, Union
from pathlib import Path
import os
import numpy as np
import pandas as pd
import anndata as ad
import s3fs

from src.utils.parsing_utils import *
from src.pseudobulking.common.pubchem import lookup_pubchem_cids
from src.pseudobulking.datasets.op3.pubchem_imputation import pubchem_mapping_op3
from src.pseudobulking.datasets.op3.gene_annotation import fetch_ensg_ids_from_symbols


def download_op3_pseudobulk(file_path: Union[str, Path]) -> None:
    """
    Download OP3 pseudobulk data from S3.
    
    Parameters
    ----------
    file_path : str or Path
        Full path where file will be downloaded
    """
    # Extract directory from file_path and create it
    if not isinstance(file_path, Path):
        file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    
    url = "s3://openproblems-data/resources/task_perturbation_prediction/datasets/neurips-2023-data/pseudobulk_filtered_with_uns.h5ad"
    
    logger.info(f"Downloading OP3 data from S3 to {file_path}")
    fs = s3fs.S3FileSystem(anon=True)
    fs.get(url, str(file_path))
    logger.info("Download complete")


def define_obs_schema() -> list:
    """
    Define the strict obs schema for pseudobulk data.
    
    Returns
    -------
    list
        List of tuples defining the schema: (column_name, dtype, description)
    """
    return [
        ("sample_id", "category", "ID of the observation: plate + well + cell_type + perturbagen"),
        ("plate", "category", "Assay detection plate identifier"),
        ("well", "category", "Well ID on the plate"),
        ("cell_type", "category", "Cell type with ontology ID"),
        ("perturbagen", "category", "Perturbagen name"),
        ("pert_type", "category", "Perturbation type"),
        ("is_control", "category", "True/False for controls"),
        ("pert_dose_uM", "float64", "Dose in micromolar"),
        ("pert_time_h", "float64", "Exposure time in hours"),
        ("suspension_type", "category", "Type of biological material that was isolated into suspension and used for profiling"),
        ("tissue", "category", "Primary tissue/site"),
        ("tissue_type", "category", "Type of tissue: tissue, cell culture or organoid"),
        ("disease", "category", "Disease/subtype"),
        ("library", "category", "Library ID"),
        ("stimulation", "category", "High-level stimulus"),
        ("guide", "category", "A guide RNA directs the CRISPR system"),
        ("dataset", "category", "Dataset label"),
        ("assay", "category", "Assay label"),
        ("development_stage", "category", "Derived from donor age"),
        ("organism", "category", "Organism"),
        ("sex", "category", "Donor sex"),
        ("self_reported_ethnicity", "category", "Donor ethnicity"),
        ("pubchem_cid", "category", "PubChem CID"),
        ("psbulk_cells", "int64", "Total #cells contributing"),
        ("psbulk_counts", "int64", "Total #counts contributing"),
        ("split", "category", "Extra column for defining a split of data"),
    ]


def build_development_stage(row: pd.Series, 
                            col: str = "age") -> str:
    """
    Build development stage label from donor age.
    
    Creates a standardized development stage label in the format
    "{age}-year-old stage" from numeric age values.
    
    Parameters
    ----------
    row : pd.Series
        Row of metadata
    col : str, default="age"
        Column name containing age information
        
    Returns
    -------
    str
        Development stage label (e.g., "45-year-old stage") or "unknown"
    """
    if col in row.index and pd.notna(row[col]):
        try:
            age = float(row[col])
            if np.isfinite(age) and age > 0:
                return f"{int(age)}-year-old stage"
        except ValueError:
            pass
    return "unknown"


def get_cell_type_ontology_map() -> dict:
    """
    Get mapping from cell type names to Cell Ontology IDs.
    
    Returns
    -------
    dict
        Dictionary mapping cell type names to CL (Cell Ontology) IDs
    """
    return {
        'NK cells': 'CL:0000623',
        'T cells': 'CL:0000084',
        'Myeloid cells': 'CL:0000763',
        'B cells': 'CL:0000236'
    }


def get_donor_metadata() -> pd.DataFrame:
    """
    Get donor metadata for OP3 dataset.
    
    Returns
    -------
    pd.DataFrame
        Donor metadata with columns: donor_id, age, sex, self_reported_ethnicity
    """
    data = {
        "donor_id": ["Donor 1", "Donor 2", "Donor 3"],
        "age": [45, 52, 45],
        "sex": ["female", "male", "male"],
        "self_reported_ethnicity": ["White", "White", "White"]
    }
    
    df = pd.DataFrame(data)
    df["development_stage"] = df.apply(lambda x: build_development_stage(x, col='age'), axis=1)
    
    return df


def get_column_rename_map() -> dict:
    """
    Get mapping from OP3 column names to standardized schema names.
    
    Returns
    -------
    dict
        Dictionary mapping original column names to schema column names
    """
    return {
        "plate_name": "plate",
        "dose_uM": "pert_dose_uM",
        "sm_name": "perturbagen",
        "control": "is_control",
        "timepoint_hr": "pert_time_h",
        "library_id": "library",
        "cell_count_by_well_celltype": "psbulk_cells",
    }


def materialize_string_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert string columns to object type for AnnData compatibility.
    
    AnnData requires string data to be stored as object dtype rather than
    pandas StringDtype. This function converts all string columns and index
    to object type.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with string-typed columns
        
    Returns
    -------
    pd.DataFrame
        Dataframe with string columns converted to object dtype
    """
    for col in df.select_dtypes(include="string").columns:
        df[col] = df[col].astype(object).where(pd.notna(df[col]), np.nan)
    
    if pd.api.types.is_string_dtype(df.index.dtype):
        df.index = df.index.astype(object)
    return df


def add_fixed_metadata_columns(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Add fixed metadata columns and sample_id specific to OP3 dataset.
    
    Parameters
    ----------
    obs : pd.DataFrame
        Observations dataframe
        
    Returns
    -------
    pd.DataFrame
        Observations with added fixed columns and sample_id
    """
    obs = obs.copy()
    
    obs["disease"] = "normal"
    obs["pert_type"] = "compound"
    obs["organism"] = "human"
    obs["suspension_type"] = "cell"
    obs["tissue"] = "blood"
    obs["tissue_type"] = "cell culture"
    obs["stimulation"] = None
    obs["guide"] = None
    obs["dataset"] = "NeurIPS2023 scPerturb DGE"
    obs["assay"] = "10x 3' v3.1"
    
    # Initialize pubchem_cid if not present
    if 'pubchem_cid' not in obs.columns:
        obs["pubchem_cid"] = None

    if 'pubchem_cid' in obs.columns:
        obs['pubchem_cid'] = pd.to_numeric(obs['pubchem_cid'], errors='coerce').fillna(-666).astype('int64')
    
    # Create composite sample_id
    obs["sample_id"] = (
        obs["plate"].astype(str).str.replace(" ", "", regex=False) + "_" +
        obs["well"].astype(str).str.replace(" ", "", regex=False) + "_" +
        obs["perturbagen"].astype(str).str.replace(" ", "_", regex=False) + "_" +
        obs["cell_type"].astype(str).str.replace(" ", "_", regex=False)
    )
    
    return obs


def process_cell_types(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Map cell type names to Cell Ontology IDs.
    
    Parameters
    ----------
    obs : pd.DataFrame
        Observations dataframe with cell_type column
        
    Returns
    -------
    pd.DataFrame
        Observations with cell_type mapped to ontology IDs
    """
    obs = obs.copy()
    cell_type_map = get_cell_type_ontology_map()
    
    # Save original cell_type values
    original_cell_type = obs["cell_type"].astype(str)
    
    # Drop the cell_type column
    obs = obs.drop(columns=["cell_type"])
    
    # Map cell type names to ontology IDs and add back
    obs["cell_type"] = original_cell_type.map(cell_type_map).astype("category")
    
    return obs


def extract_compound_metadata(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Extract unique compound metadata from observations.
    
    Parameters
    ----------
    obs : pd.DataFrame
        Observations dataframe with sm_lincs_id, sm_name, smiles columns
        
    Returns
    -------
    pd.DataFrame
        Unique compound metadata with sm_name as key
    """
    # Extract compound-related columns
    compound_cols = ['sm_lincs_id', 'sm_name', 'SMILES']
    available_cols = [col for col in compound_cols if col in obs.columns]
    
    compound_metadata = obs[available_cols].drop_duplicates('sm_name').copy()
    logger.info(f"  Extracted {len(compound_metadata):,} unique compounds")
    
    return compound_metadata


def annotate_pubchem_cids(compound_metadata: pd.DataFrame, paths: dict) -> pd.DataFrame:
    """
    Annotate unique compounds with PubChem CIDs.
    
    This function looks up PubChem CIDs for compounds using:
    1. Manual mapping from pubchem_imputation.py
    2. Automated lookup by drug name using PubChem API
    
    The results are cached to avoid redundant API calls.
    
    Parameters
    ----------
    compound_metadata : pd.DataFrame
        Unique compound metadata
    paths : dict
        Dictionary of file paths, must include 'pubchem_cache'
        
    Returns
    -------
    pd.DataFrame
        Compound metadata with pubchem_cid column populated
    """
    compound_metadata = compound_metadata.copy()
    
    logger.info("Annotating compounds with PubChem CIDs")
    
    # Initialize pubchem_cid column if not present
    if 'pubchem_cid' not in compound_metadata.columns:
        compound_metadata['pubchem_cid'] = None
    
    # Create manual mapping function
    def create_pubchem_mapping_func():
        manual_mapping = pubchem_mapping_op3().get('op3', {})
        return manual_mapping
    
    # Set up cache
    pubchem_cache = {}
    cache_path = str(paths["pubchem_cache"])
    
    # Look up PubChem CIDs
    compound_metadata = lookup_pubchem_cids(
        compound_metadata,
        cache=pubchem_cache,
        pert_id_col=None,  # OP3 doesn't have a pert_id column
        drug_col='sm_name',
        smiles_col='SMILES',
        cache_path=cache_path,
        manual_mapping_func=create_pubchem_mapping_func,
        manual_mapping_by_drug_name=True,
        dataset_key='op3'
    )
    
    # Convert to nullable integer
    compound_metadata['pubchem_cid'] = pd.to_numeric(compound_metadata['pubchem_cid'], errors='coerce').astype('Int64')
    
    return compound_metadata


def calculate_psbulk_counts(adata: ad.AnnData) -> np.ndarray:
    """
    Calculate total counts per observation from expression matrix.
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object with expression matrix
        
    Returns
    -------
    np.ndarray
        Total counts per observation
    """
    psbulk_counts = adata.X.sum(axis=1)
    if hasattr(psbulk_counts, 'A1'):  # If sparse matrix
        psbulk_counts = psbulk_counts.A1
    return np.array(psbulk_counts, dtype=int)


def enforce_obs_schema(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Enforce strict obs schema on the observations dataframe.
    
    This function ensures that the obs dataframe follows the strict schema
    defined in define_obs_schema(). It:
    1. Adds missing columns as NaN if needed
    2. Casts columns to the correct dtypes (category, float64, int64)
    3. Selects only the columns in the schema
    4. Materializes string columns for AnnData compatibility
    
    Parameters
    ----------
    obs : pd.DataFrame
        Observations dataframe to enforce schema on
        
    Returns
    -------
    pd.DataFrame
        Observations dataframe conforming to the strict schema
    """
    # Get obs schema
    obs_schema = define_obs_schema()
    
    # Create dtype map from schema
    dtype_map = {col: dtype for col, dtype, _ in obs_schema}
    
    # Create schema dataframe for reference
    obs_schema_df = pd.DataFrame(obs_schema, columns=["column", "dtype", "description"]).set_index("column")
    
    obs_for_schema = obs.copy()
    
    # Ensure all schema columns exist and cast to appropriate dtypes
    for col, dtype in dtype_map.items():
        if col not in obs_for_schema.columns:
            obs_for_schema[col] = np.nan
        
        # Cast to appropriate dtype
        if dtype == "category":
            obs_for_schema[col] = (
                obs_for_schema[col]
                .astype(object)
                .astype("category")
            )
        elif dtype == "float64":
            obs_for_schema[col] = pd.to_numeric(obs_for_schema[col], errors="coerce")
        elif dtype == "int64":
            obs_for_schema[col] = pd.to_numeric(obs_for_schema[col], errors="coerce").astype("int64")
    
    # Select only columns in schema
    obs_for_schema = obs_for_schema[obs_schema_df.index.tolist()]
    
    # Materialize string columns for AnnData compatibility
    obs_for_schema = materialize_string_columns(obs_for_schema)
    
    return obs_for_schema


def process_obs_dataframe(
    adata: ad.AnnData,
    paths: dict,
    annotate_pubchem: bool = False
) -> pd.DataFrame:
    """
    Process OP3 observations dataframe to match standardized schema.
    
    This function:
    1. Extracts compound metadata from observations
    2. Optionally annotates compounds with PubChem CIDs (on unique compounds only)
    3. Renames columns to match schema
    4. Adds fixed metadata columns (disease, organism, sample_id, etc.)
    5. Maps cell types to ontology IDs
    6. Merges compound metadata with PubChem CIDs
    7. Merges donor information
    8. Calculates pseudobulk counts
    9. Enforces strict schema
    
    Parameters
    ----------
    adata : ad.AnnData
        OP3 AnnData object with raw observations
    paths : dict
        Dictionary of file paths including pubchem_cache
    annotate_pubchem : bool, default=False
        If True, annotate compounds with PubChem CIDs
        
    Returns
    -------
    pd.DataFrame
        Processed observations dataframe conforming to schema
    """
    logger.info("Processing OP3 observations dataframe")
    
    obs = adata.obs.copy()
    
    # Extract compound metadata early
    logger.info("  Extracting compound metadata")
    compound_metadata = extract_compound_metadata(obs)
    
    # Annotate compounds with PubChem CIDs if requested
    if annotate_pubchem:
        compound_metadata = annotate_pubchem_cids(compound_metadata, paths)

    # Merge compound metadata (with PubChem CIDs if annotated)
    if annotate_pubchem and 'pubchem_cid' in compound_metadata.columns:
        logger.info("  Merging PubChem annotations")
        obs = obs.drop(columns=['pubchem_cid'], errors='ignore')
        obs = obs.merge(
            compound_metadata[['sm_name', 'pubchem_cid']],
            on='sm_name',
            how='left'
        )
    
    # Rename columns to match schema
    logger.info("  Renaming columns")
    col_rename_map = get_column_rename_map()
    obs = obs.rename(columns=col_rename_map)
    
    # Add fixed metadata columns and sample_id
    logger.info("  Adding fixed metadata columns and sample_id")
    obs = add_fixed_metadata_columns(obs)
    
    # Map cell types to ontology IDs
    logger.info("  Mapping cell types to ontology IDs")
    obs = process_cell_types(obs)
    
    
    
    # Merge donor metadata
    logger.info("  Merging donor metadata")
    donor_df = get_donor_metadata()
    obs = obs.merge(donor_df, on='donor_id', how='left')
    
    # Calculate pseudobulk counts
    logger.info("  Calculating pseudobulk counts")
    obs['psbulk_counts'] = calculate_psbulk_counts(adata)
    
    # Enforce schema
    logger.info("  Enforcing strict obs schema")
    obs_final = enforce_obs_schema(obs)
    
    logger.info(f"  Processed {len(obs_final):,} observations")
    
    return obs_final


def process_gene_annotations(
    var: pd.DataFrame,
    paths: dict,
    annotate_genes: bool = True
) -> pd.DataFrame:
    """
    Process gene annotations and optionally map to Ensembl IDs.
    
    Parameters
    ----------
    var : pd.DataFrame
        Gene/variable annotations with gene symbols as index
    paths : dict
        Dictionary with paths including 'ensembl_cache'
    annotate_genes : bool, default=True
        If True, annotate genes with Ensembl IDs using API lookups
        
    Returns
    -------
    pd.DataFrame
        Processed var dataframe with ensembl_id as index and symbol column
    """
    logger.info("Processing gene annotations")
    
    # Create empty dataframe with gene symbols as index
    var_df = pd.DataFrame(index=var.index)
    
    # Store original gene symbols
    var_df['symbol'] = var_df.index.astype(object).astype("category")
    
    # Get gene symbols
    gene_symbols = var.index.tolist()
    
    if annotate_genes:
        # Fetch Ensembl IDs from symbols
        cache_path = str(paths["ensembl_cache"])
        ensg_ids = fetch_ensg_ids_from_symbols(gene_symbols, cache_path=cache_path)
        var_df['ensembl_id'] = ensg_ids
    else:
        # Initialize with None
        var_df['ensembl_id'] = None
    
    # Fill missing Ensembl IDs with gene symbol
    var_df["ensembl_id"] = var_df["ensembl_id"].fillna(var_df.index.to_series()).astype(object).astype("category")
    
    # Handle duplicate ensembl IDs
    duplicate_mask = var_df["ensembl_id"].duplicated(keep=False)
    if duplicate_mask.any():
        logger.warning(f"  Found {duplicate_mask.sum()} duplicate ensembl IDs")
    
    # Set ensembl_id as index
    var_df = var_df.set_index("ensembl_id", drop=True)
    
    var_df = materialize_string_columns(var_df)
    
    return var_df


def standardize_op3_dataset(
    data_root: Optional[str] = None,
    file: str = "pseudobulk_filtered_with_uns.h5ad",
    download_if_missing: bool = True,
    annotate_pubchem: bool = False
) -> ad.AnnData:
    """
    Standardize OP3 dataset to match pseudobulk schema.
    
    This function:
    1. Downloads data from S3 if missing
    2. Loads AnnData object
    3. Processes observations to match schema
    4. Annotates genes with Ensembl IDs
    5. Optionally annotates compounds with PubChem CIDs
    6. Returns standardized AnnData
    
    Parameters
    ----------
    data_root : str, optional
        Root directory for data files. Defaults to './op3_data'
    file : str, default="pseudobulk_filtered_with_uns.h5ad"
        Filename for the h5ad file
    download_if_missing : bool, default=True
        If True, download data if not found locally
    annotate_pubchem : bool, default=False
        If True, annotate compounds with PubChem CIDs using API lookups.
        This may involve API calls to PubChem and can be time-consuming.
        
    Returns
    -------
    ad.AnnData
        Standardized OP3 AnnData object with schema-compliant obs and var
    """
    logger.info("Applying OP3-specific processing")
    
    if data_root is None:
        data_root = './op3_data'
    
    data_root = Path(data_root)
    data_root.mkdir(parents=True, exist_ok=True)
    
    # Define paths
    processed_dir = data_root / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)
    
    paths = {
        "pubchem_cache": processed_dir / "pubchem_cache.json",
        "ensembl_cache": processed_dir / "ensembl_cache.json",
    }
    
    file_path = data_root / file
    
    # Download if missing
    if not file_path.exists() and download_if_missing:
        logger.info(f"Data file not found at {file_path}")
        download_op3_pseudobulk(file_path)
    elif not file_path.exists():
        raise FileNotFoundError(f"Data file not found at {file_path}. Set download_if_missing=True to download.")
    
    # Load AnnData
    logger.info(f"Loading OP3 data from {file_path}")
    adata = ad.read_h5ad(file_path)
    logger.info(f"  Loaded AnnData: {adata.n_obs:,} × {adata.n_vars:,}")
    
    # Process observations
    obs_standardized = process_obs_dataframe(adata, paths, annotate_pubchem=annotate_pubchem)
    
    # Process gene annotations
    var_standardized = process_gene_annotations(adata.var, paths, annotate_genes=True)
    
    # Create new AnnData with standardized obs and var
    logger.info("Creating standardized AnnData object")
    adata_standardized = ad.AnnData(
        X=adata.X,
        obs=obs_standardized.set_index("sample_id"),
        var=var_standardized
    )
    
    logger.info(f"Standardized AnnData: {adata_standardized.n_obs:,} × {adata_standardized.n_vars:,}")
    logger.info("OP3-specific processing completed")
    
    return adata_standardized
