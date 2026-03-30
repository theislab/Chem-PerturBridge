"""
Novartis MoABox DRUG-seq dataset standardization.

Converts the raw Novartis AnnData to the unified pseudobulk schema.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import pubchempy as pcp

from src.deg.ensembl_mapping import _aggregate_matrix
from src.pseudobulking.common.pubchem import is_valid_pubchem_cid, lookup_pubchem_cids
from src.utils.parsing_utils import logger


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

POS_CTL = ['BD-11-DV28', 'EA-18-FP00', 'SE-15-AV21']


# ---------------------------------------------------------------------------
# .var standardization
# ---------------------------------------------------------------------------

def filter_var(adata: ad.AnnData, genes: pd.DataFrame) -> pd.DataFrame:
    """
    Build a clean geneID → ensembl_gene_id mapping from adata.var and the gene annotation table.

    Removes ensembl_gene_ids that belong to any geneID mapping to more than one ensembl_gene_id
    (split genes), along with all geneIDs linked to those contaminated ensembl_gene_ids.

    Parameters
    ----------
    adata:
        AnnData with var.index = geneID.
    genes:
        Gene annotation table with columns 'geneID' and 'ensembl_gene_id'.

    Returns
    -------
    Filtered flat DataFrame with columns 'geneID', 'ensembl_gene_id' (and any other gene columns),
    containing only unambiguous 1-to-1 or many-to-1 (ensembl_gene_id side) mappings.
    """
    idx_col = adata.var.index.name or 'index'
    var = adata.var.reset_index().rename(columns={idx_col: 'geneID'})
    var = var.merge(genes, how='left')
    var_by_geneid = var.groupby('geneID')['ensembl_gene_id'].agg(list).reset_index()
    ens_to_filter = (
        var_by_geneid[var_by_geneid['ensembl_gene_id'].str.len() > 1][['ensembl_gene_id']]
        .explode('ensembl_gene_id')
    )
    var_filtered = var[~var['ensembl_gene_id'].isin(ens_to_filter['ensembl_gene_id'].unique())]

    n_no_annotation = var_filtered['ensembl_gene_id'].isna().sum()
    if n_no_annotation > 0:
        logger.warning('Genes with no annotation (ensembl_gene_id=NaN), dropping: %d', n_no_annotation)
        var_filtered = var_filtered.dropna(subset=['ensembl_gene_id'])

    n_genes = adata.n_vars
    logger.info('Genes before filtration: %d', n_genes)
    logger.info('Genes after filtration:  %d (%d removed)',
                var_filtered["geneID"].nunique(), n_genes - var_filtered["geneID"].nunique())
    return var_filtered


def standardize_var(adata: ad.AnnData, var_filtered: pd.DataFrame) -> ad.AnnData:
    """
    Filter adata to clean geneIDs and aggregate X counts for geneIDs
    that share the same ensembl_gene_id (many-to-one, split-free mappings).

    Output conforms to the pseudobulk format spec (docs/Format_Pseudobulk.ipynb):
      - var.index : ensembl_id <object>
      - var['symbol'] : gene symbol <category>

    Parameters
    ----------
    adata:
        AnnData with var.index = geneID.
    var_filtered:
        Gene mapping table (columns: geneID, ensembl_gene_id, ...)
        with contaminated genes removed.

    Returns
    -------
    AnnData with var.index = ensembl_id and summed X counts.
    """
    assert var_filtered['geneID'].nunique() == var_filtered.shape[0], \
        "var_filtered contains duplicate geneIDs — filtration may be incomplete"
    geneID_to_ens = var_filtered.set_index('geneID')['ensembl_gene_id']

    adata = adata[:, adata.var.index.isin(geneID_to_ens.index)].copy()
    logger.info("AnnData subset to filtered geneIDs: %d vars", adata.n_vars)

    var_by_ens = var_filtered.groupby('ensembl_gene_id')['geneID'].agg(list).reset_index()
    to_merge = var_by_ens[var_by_ens['geneID'].str.len() > 1]
    logger.info('ensembl_gene_ids with >1 geneID (will be merged): %d', len(to_merge))

    def _build_var(ens_ids: pd.Index) -> pd.DataFrame:
        """Build var DataFrame with symbol column, indexed by ensembl_id.

        geneID has the format 'SYMBOL,chr' — split on ',' to extract the base symbol.
        Use plain symbol when unique across all ensembl_ids in the result;
        fall back to the full geneID (symbol,chr) for duplicated base symbols.

        drop_duplicates('ensembl_gene_id') is safe here: after the aggregation step
        each ensembl_id maps to exactly one geneID, so no information is lost.
        """
        geneID = (
            var_filtered.drop_duplicates('ensembl_gene_id')
                        .set_index('ensembl_gene_id')['geneID']
                        .reindex(ens_ids)
        )
        base_symbol = geneID.str.split(',').str[0]
        duplicated_mask = base_symbol.duplicated(keep=False)
        symbol = base_symbol.copy()
        symbol[duplicated_mask] = geneID[duplicated_mask]
        return pd.DataFrame({'symbol': symbol.astype('category')}, index=ens_ids.rename('ensembl_id'))

    if len(to_merge) == 0:
        logger.info('No aggregation needed — reindexing var to ensembl_id')
        new_index = pd.Index(adata.var.index.map(geneID_to_ens), name='ensembl_id')
        adata.var = _build_var(new_index)
        return adata

    mapping_df = pd.DataFrame(
        {'old_id': adata.var.index, 'new_id': adata.var.index.map(geneID_to_ens)},
        index=adata.var.index,
    )

    groups          = mapping_df.groupby('new_id')
    unique_ens_ids  = pd.Index(mapping_df['new_id'].dropna().unique())
    new_id_to_index = {eid: idx for idx, eid in enumerate(unique_ens_ids)}

    X_agg = _aggregate_matrix(adata.X, groups, new_id_to_index, adata.var.index)

    adata = ad.AnnData(
        X   = X_agg,
        obs = adata.obs.copy(),
        var = _build_var(unique_ens_ids),
    )
    logger.info('Genes after aggregation (ensembl_id): %d', adata.n_vars)
    return adata


# ---------------------------------------------------------------------------
# Sample filtering
# ---------------------------------------------------------------------------

def filter_samples(adata: ad.AnnData, robust_dmso: pd.DataFrame) -> ad.AnnData:
    """
    Filter cells from the Novartis dataset.

    Removes:
    1. Cells with an empty well_type.
    2. RC wells whose external_biosample_id is not listed in robust_dmso.
    3. SA wells whose cmpd_sample_id is a positive control (POS_CTL).

    Parameters
    ----------
    adata:
        Raw Novartis AnnData object.
    robust_dmso:
        Robust Reference Control (RC) DMSO wells.

    Returns
    -------
    Filtered AnnData.
    """
    logger.info('Original: %d cells', adata.shape[0])

    robust_ext_biosample_ids = robust_dmso['external_biosample_id'].astype(str)
    is_outlier = pd.Series(np.zeros(adata.obs.shape[0], dtype=bool), index=adata.obs.index)

    mask_empty_well_type = adata.obs['well_type'] == 'EMPTY'
    logger.info('Empty well_type: %d cells', mask_empty_well_type.sum())
    is_outlier = is_outlier | mask_empty_well_type

    mask_rc_not_in_dmso = (
        (adata.obs['well_type'] == 'RC')
        & ~adata.obs['external_biosample_id'].astype(str).isin(robust_ext_biosample_ids)
    )
    logger.info('RC not in robust_dmso: %d cells', mask_rc_not_in_dmso.sum())
    is_outlier = is_outlier | mask_rc_not_in_dmso

    mask_sa_poscon = (
        (adata.obs['well_type'] == 'SA')
        & adata.obs['cmpd_sample_id'].isin(POS_CTL)
    )
    logger.info('SA in POS_CTL: %d cells', mask_sa_poscon.sum())
    is_outlier = is_outlier | mask_sa_poscon

    logger.info('Total outliers: %d / %d', is_outlier.sum(), adata.obs.shape[0])
    adata = adata[~is_outlier]
    logger.info('Remaining: %d cells', adata.shape[0])
    return adata


# ---------------------------------------------------------------------------
# Compound annotation
# ---------------------------------------------------------------------------

def annotate_novartis_pubchem_cids(
    compound_df: pd.DataFrame,
    cache_path: str = './novartis_pubchem_cache.json',
) -> pd.DataFrame:
    """
    Annotate Novartis compound metadata with PubChem CIDs.

    Lookup order for each compound:
    1. Universal cache keyed by cmpd_sample_id
    2. Lookup by inchi_key
    3. Lookup by smiles
    4. Lookup by cas_number (used as drug name)

    Parameters
    ----------
    compound_df:
        DataFrame with at minimum the columns
        ``cmpd_sample_id``, ``inchi_key``, ``smiles``, ``cas_number``.
    cache_path:
        Path to a JSON file used for persistent caching of CID lookups.

    Returns
    -------
    compound_df with an added / updated ``pubchem_cid`` column (category, None for not found).
    """
    cache = {}
    annotated = lookup_pubchem_cids(
        compound_df,
        cache=cache,
        pert_id_col='cmpd_sample_id',
        drug_col='cas_number',
        inchikey_col='inchi_key',
        smiles_col='smiles',
        cache_path=cache_path,
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
            if isinstance(e, (pcp.PubChemHTTPError, pcp.TimeoutError,
                               pcp.ServerError, pcp.ServerBusyError)):
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
    2. Build perturbagen from pubchem_name; fall back to cmpd_sample_id for NaN values
       or for pubchem_name values duplicated across different compounds.

    Parameters
    ----------
    compound_df:
        Compound annotation table with at least cmpd_sample_id and pubchem_cid.

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
    perturbagen[nan_mask | dup_mask] = df.loc[nan_mask | dup_mask, 'cmpd_sample_id']

    df['perturbagen'] = perturbagen

    n_nan = nan_mask.sum()
    n_dup = dup_mask.sum()
    logger.info('perturbagen: %d from pubchem_name, %d NaN -> cmpd_sample_id, %d duplicated -> cmpd_sample_id',
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
        ("is_control",             "category", "True/False for controls"),
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
        'U-2 OS': 'CVCL_0042'
    }


def get_column_rename_map() -> dict:
    """Map raw Novartis column names to the unified pseudobulk schema names."""
    return {
        'plate_barcode':        'plate',
        'well_id':              'well',
        'concentration':        'pert_dose_uM',
        'hours_post_treatment': 'pert_time_h',
    }



def add_fixed_metadata_columns(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Add Novartis-specific fixed metadata and build sample_id.

    - is_control derived from well_type (RC = reference control)
    - pert_dose_uM set to 0 for controls
    - Donor info fixed: U-2 OS was established from a 15-year-old female (Caucasian)
    - Fixed dataset-level fields: organism, tissue, assay, etc.
    - sample_id: plate_well_perturbagen_cell_type
    """
    obs = obs.copy()
    cell_type_map = get_cell_type_ontology_map()

    obs["is_control"] = obs["well_type"] == "RC"
    obs.loc[obs["is_control"], "pert_dose_uM"] = 0.0

    obs["cell_type"]               = obs["cell_line_name"].map(cell_type_map)
    obs["disease"]                 = "osteosarcoma"
    obs["pert_type"]               = "compound"
    obs["organism"]                = "human"
    obs["suspension_type"]         = "cell"
    obs["tissue"]                  = "bone"
    obs["tissue_type"]             = "cell culture"
    obs["library"]                 = None
    obs["stimulation"]             = None
    obs["guide"]                   = None
    obs["dataset"]                 = "Novartis MoABox DRUG-seq"
    obs["assay"]                   = "DRUG-seq"

    # Donor metadata: U-2 OS established from a 15-year-old Caucasian female
    obs["sex"]                     = "female"
    obs["self_reported_ethnicity"] = "Caucasian"
    obs["development_stage"]       = "15-year-old stage"

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
        if dtype == "category":
            obs_out[col] = obs_out[col].astype(object).astype("category")
        elif dtype == "float64":
            obs_out[col] = pd.to_numeric(obs_out[col], errors="coerce")
        elif dtype == "int64":
            obs_out[col] = pd.to_numeric(obs_out[col], errors="coerce").astype("int64")

    obs_out = obs_out[schema_cols]
    return obs_out


def process_obs_dataframe(adata: ad.AnnData, compound_df: pd.DataFrame) -> pd.DataFrame:
    """
    Process Novartis observations to match the unified pseudobulk schema.

    Steps:
    1. Merge perturbagen and pubchem_cid from compound_df
    2. Rename columns to schema names
    3. Add fixed metadata columns and sample_id (including donor info and is_control)
    4. Set psbulk_cells = -666 (not available for DRUG-seq); compute psbulk_counts from X
    5. Enforce strict obs schema

    Parameters
    ----------
    adata:
        Filtered Novartis AnnData.
    compound_df:
        Compound annotation table with at least cmpd_sample_id and pubchem_cid.

    Returns
    -------
    pd.DataFrame conforming to the unified obs schema (sample_id column, not index).
    """
    obs = adata.obs.copy()

    # bring processed perturbagen and pubchem_cid from compound_df
    merge_cols = [c for c in ["cmpd_sample_id", "perturbagen", "pubchem_cid"]
                  if c in compound_df.columns]
    drop_cols  = [c for c in ["perturbagen", "pubchem_cid"] if c in obs.columns]
    obs = obs.drop(columns=drop_cols, errors="ignore")
    obs = obs.merge(
        compound_df[merge_cols].drop_duplicates("cmpd_sample_id"),
        on="cmpd_sample_id", how="left",
    )

    obs = obs.rename(columns=get_column_rename_map())
    obs = add_fixed_metadata_columns(obs)

    obs["psbulk_cells"]  = -666   # not available for DRUG-seq
    obs["psbulk_counts"] = calculate_psbulk_counts(adata)

    obs_final = enforce_obs_schema(obs)
    logger.info(f"  Processed {len(obs_final):,} observations")
    return obs_final


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def standardize_novartis_dataset(
    paths: dict,
    annotate_pubchem: bool = False,
    annotate_pubchem_names: bool = False,
) -> ad.AnnData:
    """
    Full standardization pipeline for the Novartis MoABox DRUG-seq dataset.

    Steps
    -----
    1. Load raw AnnData, gene annotation, compound metadata, robust DMSO table.
    2. Filter and map .var geneIDs → ensembl_ids; aggregate counts where needed.
    3. Filter outlier cells (empty wells, non-robust RC wells, positive controls).
    4. Optionally annotate compounds with PubChem CIDs (annotate_pubchem).
    5. Optionally look up drug synonym names by CID (annotate_pubchem_names).
    6. Build unified .obs conforming to the pseudobulk schema.
    7. Convert .X to int64.

    Parameters
    ----------
    paths:
        Dictionary of file paths as returned by get_novartis_paths(). Expected keys:
        raw_h5ad, genes_csv, compound_tsv, robust_dmso,
        pubchem_cid_cache, pubchem_names_cache.
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
    genes       = pd.read_csv(paths["genes_csv"])
    compound_df = pd.read_csv(paths["compound_tsv"], sep='\t')
    robust_dmso = pd.read_csv(paths["robust_dmso"], sep='\t')

    logger.info("Standardizing .var ...")
    var_filtered = filter_var(adata, genes)
    adata = standardize_var(adata, var_filtered)

    logger.info("Filtering cells ...")
    adata = filter_samples(adata, robust_dmso)

    if annotate_pubchem:
        logger.info("Annotating compounds with PubChem CIDs ...")
        compound_df = annotate_novartis_pubchem_cids(
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
    obs_standardized = process_obs_dataframe(adata, compound_df).set_index('sample_id')

    var_standardized = adata.var.copy()
    var_standardized.index = var_standardized.index.astype('category')

    logger.info("Creating standardized AnnData object ...")
    adata_standardized = ad.AnnData(
        X=adata.X.astype(np.int64),
        obs=obs_standardized,
        var=var_standardized,
    )

    logger.info(
        "Standardization complete: %d observations x %d genes",
        adata_standardized.n_obs, adata_standardized.n_vars,
    )
    return adata_standardized
