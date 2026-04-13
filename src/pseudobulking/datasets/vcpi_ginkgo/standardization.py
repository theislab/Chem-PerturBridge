"""
VCPI Ginkgo bulk RNA dataset standardization.

Converts the raw VCPI AnnData to the unified pseudobulk schema.
"""
from __future__ import annotations

import json
import time
import urllib.error
from pathlib import Path
from typing import Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import pubchempy as pcp

from src.pseudobulking.common.pubchem import is_valid_pubchem_cid, lookup_pubchem_cids
from src.pseudobulking.datasets.vcpi_ginkgo.pubchem_imputation import pubchem_mapping_vcpi_ginkgo
from src.pseudobulking.datasets.vcpi_ginkgo.ensembl_symbol_lookup import (
    fetch_symbols_for_ensembl_ids,
)
from src.utils.parsing_utils import logger


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .var standardization
# ---------------------------------------------------------------------------

def add_symbols_to_adata_var(
    adata: ad.AnnData,
    symbol_map: Dict[str, Optional[str]],
    symbol_col: str = "symbol",
    fallback_to_ensembl_id: bool = False,
) -> ad.AnnData:
    """
    Add a ``symbol`` column to ``adata.var`` from a pre-fetched symbol map.

    Parameters
    ----------
    adata:
        AnnData whose ``var.index`` contains Ensembl gene IDs.
    symbol_map:
        dict mapping Ensembl ID → gene symbol (as returned by
        ``fetch_symbols_for_ensembl_ids``).
    symbol_col:
        Name of the output column in ``adata.var``.
    fallback_to_ensembl_id:
        If True, fill missing symbols with the Ensembl ID itself.

    Returns
    -------
    AnnData with ``adata.var[symbol_col]`` populated (dtype: category).
    """
    adata = adata.copy()
    symbols = pd.Series(adata.var_names, index=adata.var_names).map(symbol_map)

    n_mapped   = symbols.notna().sum()
    n_unmapped = symbols.isna().sum()
    logger.info(
        "Symbol mapping: %d mapped, %d unmapped out of %d genes",
        n_mapped, n_unmapped, len(symbols),
    )

    if fallback_to_ensembl_id and n_unmapped > 0:
        logger.info("Filling %d unmapped genes with their Ensembl ID", n_unmapped)
        symbols = symbols.fillna(pd.Series(adata.var_names, index=adata.var_names))

    adata.var[symbol_col] = symbols.astype("category")
    return adata


def standardize_var(
    adata: ad.AnnData,
    ensembl_symbol_cache: Optional[str] = None,
) -> ad.AnnData:
    """
    Drop ERCC spike-ins and set the unified var schema.

    Output conforms to the pseudobulk format spec:
      - var.index : ensembl_id <object>
      - var['symbol'] : gene symbol <category>

    Parameters
    ----------
    adata:
        AnnData with var.index = Ensembl gene IDs (possibly with version suffix).
    ensembl_symbol_cache:
        Optional path to a JSON cache file for Ensembl → symbol lookups.
        Results are cached for reuse across runs. When None, lookups are still
        performed but results are not persisted to disk.

    Returns
    -------
    AnnData with ERCC genes removed and var.index renamed to ensembl_id.
    """
    var = adata.var.copy()
    var.index = var.index.astype(str)
    var.index.name = "ensembl_id"

    ercc = var.index.str.startswith("ERCC-", na=False)
    n_ercc = int(ercc.sum())
    if n_ercc:
        logger.info("Dropping %d ERCC spike-in genes", n_ercc)

    adata = adata[:, ~ercc].copy()
    adata.var.index.name = "ensembl_id"

    logger.info("Fetching gene symbols from Ensembl REST API ...")
    symbol_map = fetch_symbols_for_ensembl_ids(
        list(adata.var_names),
        cache_path=ensembl_symbol_cache,
    )
    adata = add_symbols_to_adata_var(adata, symbol_map)

    logger.info("Genes after ERCC removal: %d", adata.n_vars)
    return adata


# ---------------------------------------------------------------------------
# Sample filtering
# ---------------------------------------------------------------------------

def filter_samples(adata: ad.AnnData) -> ad.AnnData:
    """
    Filter samples from the VCPI Ginkgo dataset.

    Removes:
    1. Samples with total counts ``adata.X.sum(axis=1) <= 0`` (empty library).
    2. Positive control samples: ``is_control == True`` but
       ``user_compound_id != "DMSO"`` (non-DMSO controls).

    Parameters
    ----------
    adata:
        Raw VCPI Ginkgo AnnData object (obs must still contain ``user_compound_id``).

    Returns
    -------
    Filtered AnnData.
    """
    logger.info("Original: %d samples", adata.n_obs)
    is_outlier = pd.Series(False, index=adata.obs.index)

    lib_sum = np.asarray(adata.X.sum(axis=1)).ravel()
    mask_zero_library = pd.Series(lib_sum <= 0, index=adata.obs.index)
    logger.info("Zero or negative library size (sum X <= 0): %d samples", mask_zero_library.sum())
    is_outlier = is_outlier | mask_zero_library

    mask_pos_control = (
        (adata.obs["is_control"] == True)
        & (adata.obs["user_compound_id"] != "DMSO")
    )
    logger.info("Positive controls (non-DMSO is_control): %d samples", mask_pos_control.sum())
    is_outlier = is_outlier | mask_pos_control

    logger.info("Total outliers: %d / %d", is_outlier.sum(), adata.n_obs)
    adata = adata[~is_outlier].copy()
    logger.info("Remaining: %d samples", adata.n_obs)
    return adata


# ---------------------------------------------------------------------------
# Compound annotation
# ---------------------------------------------------------------------------

def annotate_vcpi_ginkgo_pubchem_cids(
    compound_df: pd.DataFrame,
    cache_path: str = './vcpi_ginkgo_pubchem_cache.json',
) -> pd.DataFrame:
    """
    Annotate VCPI Ginkgo compound metadata with PubChem CIDs.

    Lookup order for each compound:
    1. Universal cache keyed by compound
    2. Lookup by inchi_key
    3. Lookup by canonical_smiles
    4. Lookup by compound name

    Parameters
    ----------
    compound_df:
        DataFrame with at minimum the columns
        ``user_compound_id``, ``compound``, ``inchi_key``, ``canonical_smiles``.
    cache_path:
        Path to a JSON file used for persistent caching of CID lookups.

    Returns
    -------
    compound_df with an added / updated ``pubchem_cid`` column (category, None for not found).
    """
    compound_df = compound_df.copy()
    compound_df['compound_name'] = compound_df['user_compound_id'].apply(
        lambda v: None if _is_numeric_string(v) else v
    )

    cache = {}
    annotated = lookup_pubchem_cids(
        compound_df,
        cache=cache,
        pert_id_col='compound',
        drug_col='compound_name',
        inchikey_col='inchi_key',
        smiles_col='canonical_smiles',
        cache_path=cache_path,
        manual_mapping_func=pubchem_mapping_vcpi_ginkgo,
        manual_mapping_by_drug_name=True,
        dataset_key='vcpi_ginkgo',
    )
    annotated['pubchem_cid'] = pd.to_numeric(annotated['pubchem_cid'], errors='coerce').fillna(-666).astype('int64')
    annotated['pubchem_cid'] = annotated['pubchem_cid'].replace({-666: None, '-666': None, 'None': None, 'nan': None, '<NA>': None})
    annotated['pubchem_cid'] = annotated['pubchem_cid'].astype(object).astype('category')
    return annotated


def _fetch_drug_name_by_cid(
    cid: int,
    cache: Dict[int, Optional[str]],
    n_retries: int = 5,
) -> Optional[str]:
    """Fetch the first synonym for a PubChem CID with retry logic.

    Uses compound.synonyms (common/trade names) rather than compound.iupac_name.
    The first synonym is typically the most widely-used name.
    """
    cnt = 0
    while cnt < n_retries:
        try:
            compounds = pcp.get_compounds(cid, 'cid')
            synonyms = compounds[0].synonyms if compounds else []
            name = synonyms[0] if synonyms else None
            cache[cid] = name
            return name
        except Exception as e:
            if isinstance(
                e,
                (
                    pcp.PubChemHTTPError,
                    pcp.TimeoutError,
                    pcp.ServerError,
                    pcp.ServerBusyError,
                    ConnectionError,
                    urllib.error.URLError,
                ),
            ):
                logger.warning('PubChem lookup failed for CID %d: %s. Retry %d.', cid, e, cnt)
                cnt += 1
                time.sleep(5)
            else:
                logger.warning('PubChem lookup failed for CID %d: %s', cid, e)
                break
    cache[cid] = None
    return None


def get_pubchem_name_overrides() -> Dict[int, str]:
    """
    Manual CID → name mappings for compounds where automatic lookup returns
    a poor synonym or nothing at all.

    Returns
    -------
    dict mapping int PubChem CID to the preferred name string.
    """
    return {
        679: 'DMSO',
    }


def lookup_drug_names_by_cid(
    df: pd.DataFrame,
    pubchem_cid_col: str = 'pubchem_cid',
    drug_name_col: str = 'pubchem_name',
    cache_path: Optional[str] = None,
    n_retries: int = 5,
    request_delay_s: float = 0.2,
) -> pd.DataFrame:
    """
    Add a drug name column to a DataFrame by reverse-looking up PubChem CIDs.

    For each row with a valid pubchem_cid the first synonym (common/trade name)
    is fetched from PubChem. Results are cached to avoid redundant API calls.

    Parameters
    ----------
    df:
        DataFrame containing a ``pubchem_cid`` column.
    pubchem_cid_col:
        Name of the column holding PubChem CIDs.
    drug_name_col:
        Name of the output column that will receive the fetched names.
    cache_path:
        Path to a JSON file used for persistent caching of name lookups.
    n_retries:
        Number of retries on transient PubChem errors.
    request_delay_s:
        Delay in seconds between API calls (fair-use throttling).

    Returns
    -------
    df with an added / updated ``drug_name_col`` column containing the first
    PubChem synonym (common/trade name, not IUPAC systematic name).
    """
    df = df.copy()

    cache: Dict[int, Optional[str]] = {}
    if cache_path:
        try:
            with open(cache_path) as f:
                raw = json.load(f)
            for k, v in raw.items():
                try:
                    cache[int(k)] = v
                except (ValueError, TypeError):
                    pass  # skip keys that are not integer CIDs
        except FileNotFoundError:
            pass

    # manual overrides take priority over both cache and API results
    cache.update(get_pubchem_name_overrides())

    unique_cids = [
        int(cid) for cid in df[pubchem_cid_col].dropna().unique()
        if is_valid_pubchem_cid(cid)
    ]
    to_fetch = [cid for cid in unique_cids if cid not in cache]
    logger.info('Unique valid CIDs: %d, to fetch: %d', len(unique_cids), len(to_fetch))

    for i, cid in enumerate(to_fetch, 1):
        _fetch_drug_name_by_cid(cid, cache, n_retries)
        if i % 50 == 0:
            logger.info('Fetched %d/%d', i, len(to_fetch))
        if cache_path and i % 500 == 0:
            with open(cache_path, 'w') as f:
                json.dump({str(k): v for k, v in cache.items()}, f)
        if request_delay_s > 0:
            time.sleep(request_delay_s)

    if cache_path:
        with open(cache_path, 'w') as f:
            json.dump({str(k): v for k, v in cache.items()}, f)

    df[drug_name_col] = df[pubchem_cid_col].apply(
        lambda cid: cache.get(int(cid)) if is_valid_pubchem_cid(cid) else None
    )
    n_mapped = df[drug_name_col].notna().sum()
    logger.info('Mapped %d/%d rows to drug names', n_mapped, len(df))
    return df


def process_compound_df(compound_df: pd.DataFrame) -> pd.DataFrame:
    """
    Enrich compound_df with a perturbagen column.

    Assumes lookup_drug_names_by_cid was already called before this function.

    Steps:
    1. If pubchem_name column is absent, fill it with None.
    2. Build perturbagen from pubchem_name; fall back to user_compound_id for NaN values
       or for pubchem_name values duplicated across different compounds.

    Parameters
    ----------
    compound_df:
        Compound annotation table with at least user_compound_id and pubchem_cid.

    Returns
    -------
    compound_df with added/updated pubchem_name and perturbagen columns.
    """
    df = compound_df.copy()

    if 'pubchem_name' not in df.columns:
        df['pubchem_name'] = None

    perturbagen = df['pubchem_name'].copy()
    nan_mask = perturbagen.isna()
    dup_mask = perturbagen.duplicated(keep=False) & ~nan_mask
    perturbagen[nan_mask | dup_mask] = df.loc[nan_mask | dup_mask, 'user_compound_id']

    df['perturbagen'] = perturbagen

    n_nan = nan_mask.sum()
    n_dup = dup_mask.sum()
    logger.info('perturbagen: %d from pubchem_name, %d NaN -> user_compound_id, %d duplicated -> user_compound_id',
                (~nan_mask & ~dup_mask).sum(), n_nan, n_dup)
    return df


# ---------------------------------------------------------------------------
# .obs standardization
# ---------------------------------------------------------------------------

def define_obs_schema() -> list:
    """Define the strict obs schema for pseudobulk data."""
    return [
        ("sample_id",              "category", "ID of the observation: plate + well + cell_type + perturbagen"),
        ("plate",                  "category", "Assay detection plate identifier"),
        ("well",                   "category", "Well ID on the plate"),
        ("cell_type",              "category", "Cell type with ontology ID"),
        ("perturbagen",            "category", "Perturbagen name"),
        ("pert_type",              "category", "Perturbation type"),
        ("is_control",             "category", "bool categories: [False, True]"),
        ("pert_dose_uM",           "float64",  "Dose in micromolar"),
        ("pert_time_h",            "float64",  "Exposure time in hours"),
        ("suspension_type",        "category", "Type of biological material that was isolated into suspension and used for profiling"),
        ("tissue",                 "category", "Primary tissue/site"),
        ("tissue_type",            "category", "Type of tissue: tissue, cell culture or organoid"),
        ("disease",                "category", "Disease/subtype"),
        ("library",                "category", "Library ID"),
        ("stimulation",            "category", "High-level stimulus"),
        ("guide",                  "category", "A guide RNA directs the CRISPR system"),
        ("dataset",                "category", "Dataset label"),
        ("assay",                  "category", "Assay label"),
        ("development_stage",      "category", "Derived from donor age"),
        ("organism",               "category", "Organism"),
        ("sex",                    "category", "Donor sex"),
        ("self_reported_ethnicity","category", "Donor ethnicity"),
        ("pubchem_cid",            "category", "PubChem CID"),
        ("psbulk_cells",           "int64",    "Total #cells contributing"),
        ("psbulk_counts",          "int64",    "Total #counts contributing"),
    ]


def get_cell_type_ontology_map() -> dict:
    """Map cell type names to Cell Ontology / Cellosaurus IDs."""
    return {
        'THP-1': 'CVCL_0006',
    }


def get_column_rename_map() -> dict:
    """Map raw VCPI column names to the unified pseudobulk schema names."""
    return {
        'container_id': 'plate',
    }


def add_fixed_metadata_columns(obs: pd.DataFrame, dataset_title: str) -> pd.DataFrame:
    """
    Add VCPI-specific fixed metadata and build sample_id.

    - is_control: VCPI ``is_control`` flag AND ``user_compound_id`` == DMSO
    - pert_dose_uM set to 0 for controls
    - well derived from row_id / column_id
    - Fixed dataset-level fields: organism, tissue, assay, etc.
    - sample_id: plate_well_perturbagen_cell_type
    """
    obs = obs.copy()
    cell_type_map = get_cell_type_ontology_map()

    obs["well"] = "row" + obs["row_id"].astype(str) + "_col" + obs["column_id"].astype(str)

    obs["pert_dose_uM"] = _apply_concentration_unit_to_um(
        obs["compound_concentration"],
        obs["compound_concentration_unit"],
    )
    obs["pert_time_h"] = _parse_timepoint_hours(obs["timepoint"])

    obs["cell_type"]  = obs["cell_line"].astype(str).str.strip().map(cell_type_map)
    obs["is_control"] = obs["is_control"].astype(bool).astype("category")

    non_dmso = (obs["is_control"] == True) & (obs["user_compound_id"] != "DMSO")
    if non_dmso.any():
        logger.warning("%d control samples have non-DMSO user_compound_id — check filter_samples", non_dmso.sum())

    obs.loc[obs["is_control"] == True, "pert_dose_uM"] = 0.0

    obs["pert_type"]               = "compound"
    obs["organism"]                = "human"
    obs["suspension_type"]         = "cell"
    obs["tissue"]                  = "peripheral blood"
    obs["tissue_type"]             = "cell culture"
    obs["disease"]                 = "childhood acute monocytic leukemia"
    obs["library"]                 = None
    obs["stimulation"]             = None
    obs["guide"]                   = None
    obs["dataset"]                 = dataset_title
    obs["assay"]                   = "DRUG-seq"
    obs["sex"]                     = "male"
    obs["self_reported_ethnicity"] = "Japanese"
    obs["development_stage"]       = "1-year-old stage"

    obs["plate"] = obs["plate"].astype(str).str.replace(" ", "", regex=False).str.replace("-", "_", regex=False)

    if "pubchem_cid" not in obs.columns:
        obs["pubchem_cid"] = None
    obs["pubchem_cid"] = pd.to_numeric(obs["pubchem_cid"], errors="coerce").fillna(-666).astype("int64")
    obs['pubchem_cid'] = obs['pubchem_cid'].replace({-666: None, '-666': None, 'None': None, 'nan': None, '<NA>': None})
    obs['pubchem_cid'] = obs['pubchem_cid'].astype(object).astype('category')

    obs["sample_id"] = (
        obs["plate"].astype(str).str.replace(" ", "", regex=False) + "_" +
        obs["well"].astype(str).str.replace(" ", "", regex=False) + "_" +
        obs["perturbagen"].astype(str).str.replace(" ", "_", regex=False) + "_" +
        obs["cell_type"].astype(str).str.replace(" ", "_", regex=False)
    )
    return obs


def calculate_psbulk_counts(adata: ad.AnnData) -> np.ndarray:
    """Sum counts per observation from the expression matrix."""
    psbulk_counts = adata.X.sum(axis=1)
    if hasattr(psbulk_counts, "A1"):
        psbulk_counts = psbulk_counts.A1
    return np.array(psbulk_counts, dtype=int)


def enforce_obs_schema(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Enforce strict obs schema: add missing columns, cast dtypes, select schema columns.
    """
    obs_schema = define_obs_schema()
    dtype_map  = {col: dtype for col, dtype, _ in obs_schema}
    schema_cols = [col for col, _, _ in obs_schema]

    obs_out = obs.copy()
    for col, dtype in dtype_map.items():
        if col not in obs_out.columns:
            obs_out[col] = np.nan
        if col == "is_control":
            # preserve bool category: Categories (2, bool): [False, True]
            obs_out[col] = obs_out[col].astype(bool).astype("category")
        elif dtype == "category":
            obs_out[col] = obs_out[col].astype(object).astype("category")
        elif dtype == "float64":
            obs_out[col] = pd.to_numeric(obs_out[col], errors="coerce").astype("float64")
        elif dtype == "int64":
            obs_out[col] = pd.to_numeric(obs_out[col], errors="coerce").astype("int64")

    obs_out = obs_out[schema_cols]
    return obs_out


def process_obs_dataframe(adata: ad.AnnData, compound_df: pd.DataFrame, dataset_title: str) -> pd.DataFrame:
    """
    Process VCPI observations to match the unified pseudobulk schema.

    Steps:
    1. Merge perturbagen and pubchem_cid from compound_df (on compound)
    2. Rename columns to schema names
    3. Add fixed metadata columns and sample_id
    4. Set psbulk_cells = -666 (not available for bulk); compute psbulk_counts from X
    5. Enforce strict obs schema

    Parameters
    ----------
    adata:
        VCPI AnnData (obs indexed by sequenced_id).
    compound_df:
        Compound annotation table with at least compound, perturbagen, pubchem_cid.
        user_compound_id is used internally for perturbagen fallback but not as the merge key.
    dataset_title:
        Value for the ``dataset`` obs column (varies per VCPI experiment).

    Returns
    -------
    pd.DataFrame conforming to the unified obs schema (sample_id column, not index).
    """
    obs = adata.obs.reset_index()
    if obs.columns[0] != "sequenced_id":
        obs = obs.rename(columns={obs.columns[0]: "sequenced_id"})

    # bring processed perturbagen and pubchem_cid from compound_df
    # merge on 'compound'; user_compound_id stays in compound_df for intermediate use only
    merge_cols = [c for c in ["compound", "user_compound_id", "perturbagen", "pubchem_cid"]
                  if c in compound_df.columns]
    drop_cols  = [c for c in ["perturbagen", "pubchem_cid", "user_compound_id"] if c in obs.columns]
    obs = obs.drop(columns=drop_cols, errors="ignore")
    obs = obs.merge(
        compound_df[merge_cols].drop_duplicates("compound"),
        on="compound", how="left",
    )

    obs = obs.rename(columns=get_column_rename_map())
    obs = add_fixed_metadata_columns(obs, dataset_title)

    obs["psbulk_cells"]  = -666   # not available for bulk RNA
    obs["psbulk_counts"] = calculate_psbulk_counts(adata)

    obs_final = enforce_obs_schema(obs)
    logger.info(f"  Processed {len(obs_final):,} observations")
    return obs_final


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def standardize_vcpi_ginkgo_dataset(
    paths: dict,
    annotate_pubchem: bool = False,
    annotate_pubchem_names: bool = False,
) -> ad.AnnData:
    """
    Full standardization pipeline for the VCPI Ginkgo bulk RNA dataset.

    Steps
    -----
    1. Load raw AnnData and compound metadata.
    2. Standardize .var (drop ERCC, set ensembl_id index, fetch gene symbols).
    3. Filter outlier samples (empty libraries, positive controls).
    4. Optionally annotate compounds with PubChem CIDs (annotate_pubchem).
    5. Optionally look up drug synonym names by CID (annotate_pubchem_names).
    6. Build unified .obs conforming to the pseudobulk schema.
    7. Convert .X to int64.

    Parameters
    ----------
    paths:
        Dictionary of file paths as returned by get_vcpi_ginkgo_paths(). Expected keys:
        raw_h5ad, compound_csv, dataset_title, pubchem_cid_cache, pubchem_names_cache,
        ensembl_symbol_cache.
    annotate_pubchem:
        If True, run PubChem CID lookups (requires network access).
    annotate_pubchem_names:
        If True, look up drug synonym names by PubChem CID (requires network access).
        Requires annotate_pubchem=True or a pre-existing pubchem_cid column.

    Returns
    -------
    Standardized AnnData with obs.index = sample_id and var.index = ensembl_id.
    """
    logger.info("Loading raw data ...")
    adata       = ad.read_h5ad(paths["raw_h5ad"])
    compound_df = pd.read_csv(paths["compound_csv"])

    adata.obs["is_control"] = adata.obs["is_control"].astype(str).str.strip().str.lower() == "true"

    logger.info("Standardizing .var ...")
    adata = standardize_var(adata, ensembl_symbol_cache=str(paths["ensembl_symbol_cache"]))

    logger.info("Filtering samples ...")
    adata = filter_samples(adata)

    if annotate_pubchem:
        logger.info("Annotating compounds with PubChem CIDs ...")
        compound_df = annotate_vcpi_ginkgo_pubchem_cids(
            compound_df,
            cache_path=str(paths["pubchem_cid_cache"]),
        )

    if annotate_pubchem_names:
        _has_cids = (
            'pubchem_cid' in compound_df.columns
            and compound_df['pubchem_cid'].astype(object).notna().any()
        )
        if not _has_cids:
            logger.warning(
                "annotate_pubchem_names=True but compound_df has no valid pubchem_cid values. "
                "Run with annotate_pubchem=True first or provide a pre-annotated compound table."
            )
        else:
            logger.info("Looking up drug synonym names by PubChem CID ...")
            compound_df = lookup_drug_names_by_cid(
                compound_df,
                cache_path=str(paths["pubchem_names_cache"]),
            )

    compound_df = process_compound_df(compound_df)

    logger.info("Standardizing .obs ...")
    obs_standardized = process_obs_dataframe(adata, compound_df, paths["dataset_title"]).set_index('sample_id')

    var_standardized = adata.var.copy()
    var_standardized.index = var_standardized.index.astype('category')

    logger.info("Creating standardized AnnData object ...")
    adata_standardized = ad.AnnData(
        X=adata.X.astype("int64"),
        obs=obs_standardized,
        var=var_standardized,
    )

    logger.info(
        "Standardization complete: %d observations x %d genes",
        adata_standardized.n_obs, adata_standardized.n_vars,
    )
    return adata_standardized


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _is_numeric_string(value: object) -> bool:
    """Return True if *value* can be interpreted as an integer (e.g. '12345')."""
    try:
        int(str(value))
        return True
    except (ValueError, TypeError):
        return False



def _parse_timepoint_hours(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    return pd.to_numeric(s.str.extract(r"(?i)^([\d.]+)\s*h$", expand=False), errors="coerce")


def _apply_concentration_unit_to_um(dose: pd.Series, unit: pd.Series) -> pd.Series:
    out = pd.to_numeric(dose, errors="coerce")
    u   = unit.astype(str).str.strip().str.lower()
    nm  = u == "nm"
    if not nm.any():
        return out
    out = out.astype("float64").copy()
    out.loc[nm] = out.loc[nm] / 1000.0
    return out
