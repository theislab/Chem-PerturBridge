"""
DILImap training dataset standardization.

Converts the DILImap AnnData to the common pseudobulk schema used by
the rest of the pipeline, mirroring standardize_op3_dataset in structure.

Schema reference: docs/Format_Pseudobulk.ipynb

Output layout
-------------
adata.obs  : DataFrame, index = sample_id (category), columns per define_obs_schema()
adata.var  : DataFrame, index = ensembl_id (category), columns = ["symbol"]
adata.X    : int64 counts
adata.uns  : {}
"""
from pathlib import Path
from typing import Callable, Dict, Optional

from anndata import AnnData
import anndata as ad
import dilimap as dmap
import numpy as np
import pandas as pd
import scipy.sparse

from src.utils.parsing_utils import logger
from src.pseudobulking.common.pubchem import lookup_pubchem_cids
from src.pseudobulking.datasets.dilimap.pubchem_imputation import pubchem_mapping_dilimap


# Compounds that are vehicle controls
_CONTROL_COMPOUNDS = {"DMSO"}


# DILImap uses primary human hepatocytes (PHH)
# Source: https://www.nature.com/articles/s41467-025-65690-3
# Cell culture: cryopreserved primary human hepatocytes from multiple adult donors
_CELL_TYPE = "CL_0000182"  # hepatocyte (Cell Ontology)
_DEVELOPMENT_STAGE = "37-year-old stage"
_SEX = "female"
_ETHNICITY = "Caucasian"


# ── Schema ────────────────────────────────────────────────────────────────────

def define_obs_schema() -> list:
    """
    Define the strict obs schema for pseudobulk data.

    Returns
    -------
    list
        List of tuples: (column_name, dtype, description)
    """
    return [
        ("sample_id",               "category", "ID of the observation: plate + well + perturbagen + dose"),
        ("plate",                   "category", "Assay detection plate identifier"),
        ("well",                    "category", "Well ID on the plate"),
        ("cell_type",               "category", "Cell type with ontology ID"),
        ("perturbagen",             "category", "Perturbagen name"),
        ("pert_type",               "category", "Perturbation type"),
        ("is_control",              "category", "True/False for controls"),
        ("pert_dose_uM",            "float64",  "Dose in micromolar"),
        ("pert_time_h",             "float64",  "Exposure time in hours"),
        ("suspension_type",         "category", "Type of biological material isolated into suspension"),
        ("tissue",                  "category", "Primary tissue/site"),
        ("tissue_type",             "category", "Type of tissue: tissue, cell culture or organoid"),
        ("disease",                 "category", "Disease/subtype"),
        ("library",                 "category", "Library ID"),
        ("stimulation",             "category", "High-level stimulus"),
        ("guide",                   "category", "A guide RNA directs the CRISPR system"),
        ("dataset",                 "category", "Dataset label"),
        ("assay",                   "category", "Assay label"),
        ("development_stage",       "category", "Derived from donor age"),
        ("organism",                "category", "Organism"),
        ("sex",                     "category", "Donor sex"),
        ("self_reported_ethnicity", "category", "Donor ethnicity"),
        ("pubchem_cid",             "category", "PubChem CID"),
        ("psbulk_cells",            "int64",    "Total #cells contributing"),
        ("psbulk_counts",           "int64",    "Total #counts contributing"),
        ("split",                   "category", "Extra column for defining a split of data"),
    ]


def materialize_string_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert pandas StringDtype columns to object dtype for AnnData compatibility.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    pd.DataFrame
    """
    for col in df.select_dtypes(include="string").columns:
        df[col] = df[col].astype(object).where(pd.notna(df[col]), np.nan)
    if pd.api.types.is_string_dtype(df.index.dtype):
        df.index = df.index.astype(object)
    return df


def enforce_obs_schema(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Enforce the strict obs schema defined in :func:`define_obs_schema`.

    Missing columns are added as NaN; all columns are cast to the correct
    dtype; extra columns are dropped; string columns are materialised for
    AnnData compatibility.

    Parameters
    ----------
    obs : pd.DataFrame

    Returns
    -------
    pd.DataFrame
        Schema-conformant observations dataframe.
    """
    schema = define_obs_schema()
    dtype_map = {col: dtype for col, dtype, _ in schema}
    col_order = [col for col, _, _ in schema]

    obs = obs.copy()

    for col, dtype in dtype_map.items():
        if col not in obs.columns:
            obs[col] = np.nan
        if dtype == "category":
            obs[col] = obs[col].astype(object).astype("category")
        elif dtype == "float64":
            obs[col] = pd.to_numeric(obs[col], errors="coerce").astype("float64")
        elif dtype == "int64":
            obs[col] = pd.to_numeric(obs[col], errors="coerce").astype("int64")

    obs = obs[col_order]
    obs = materialize_string_columns(obs)
    return obs


# ── QC filter (run on raw AnnData, before renaming) ──────────────────────────

def filter_samples(adata: AnnData) -> AnnData:
    """
    Apply QC filters to the raw DILImap AnnData before standardization.

    Retains only samples that pass all three QC checks:

    1. **LDH viability QC** — ``LDH_QC`` flag is empty and the compound is
       not ``DMSO_replaced``.
    2. **RNA QC** — log-total-RNA is above the per-plate lower threshold
       (median − 2.5 SD) and mitochondrial RNA fraction is below 9 %.
    3. **Cross-replicate correlation QC** — samples cluster well with their
       within-compound replicates.

    Uses raw (pre-rename) column names.
    Reference: https://www.dilimap.org/tutorials/1_Compute_Pathway_Signatures.html#2.-Quality-control

    Parameters
    ----------
    adata : AnnData
        Raw DILImap dataset.

    Returns
    -------
    AnnData
        Filtered dataset.
    """
    adata = adata.copy()
    n_initial = adata.n_obs
    logger.info(f"Filtering samples: {n_initial:,} initial observations")

    # QC flags from viability screen
    adata.obs["ldh_qc_pass"] = (adata.obs["LDH_QC"] == "") & (
        adata.obs["COMPOUND"] != "DMSO_replaced"
    )
    n_ldh_fail = (~adata.obs["ldh_qc_pass"]).sum()
    logger.info(f"  LDH QC: {n_ldh_fail:,} samples removed (failed LDH or DMSO_replaced)")

    # QC of total counts and mtRNA
    dmap.pp.qc_metrics(adata)

    thresh_counts = (
        np.median(adata.obs["log_totalRNA"])
        - 2.5 * adata.obs["log_totalRNA"].std()
    )

    adata.obs["rna_qc_pass"] = (
        adata.obs["log_totalRNA"] > thresh_counts
    ) & (adata.obs["pct_mtRNA"] < 9)
    n_rna_fail = (~adata.obs["rna_qc_pass"]).sum()
    logger.info(f"  RNA QC: {n_rna_fail:,} samples removed (low counts or high mtRNA; threshold log_totalRNA > {thresh_counts:.3f})")

    # QC (cross-replicate correlation)
    dmap.pp.qc_cross_rep_correlation(
        adata, group_key="CMPD_DOSE", plate_key="PLATE_NAME"
    )
    n_rep_fail = (~adata.obs["rep_corr_qc_pass"]).sum()
    logger.info(f"  Cross-replicate correlation QC: {n_rep_fail:,} samples removed")

    adata = adata[
        adata.obs["ldh_qc_pass"]
        & adata.obs["rna_qc_pass"]
        & adata.obs["rep_corr_qc_pass"]
    ].copy()

    n_final = adata.n_obs
    logger.info(f"  Retained: {n_final:,} / {n_initial:,} samples ({n_initial - n_final:,} removed in total)")

    return adata


# ── obs / var processing ──────────────────────────────────────────────────────

def standardize_obs_dilimap(adata: AnnData) -> pd.DataFrame:
    """
    Process DILImap observations into a schema-conformant DataFrame.

    Parameters
    ----------
    adata : AnnData
        Raw (or QC-filtered) DILImap dataset.

    Returns
    -------
    pd.DataFrame
        Observations conforming to :func:`define_obs_schema`, with
        ``sample_id`` as a regular column (not the index).
    """
    obs = adata.obs.copy()

    # ── rename to standard names ──────────────────────────────────────
    rename_cols = {
        "COMPOUND":          "perturbagen",
        "CONCENTRATION_UM":  "pert_dose_uM",
        "TIMEPOINT_HOURS":   "pert_time_h",
        "PLATE_NAME":        "plate",
        "WELL_ID":           "well",
        "LIBRARY_ID":        "library",
        "SPECIES":           "organism",
        "BATCH_ID":          "batch",
    }
    obs.rename(columns=rename_cols, inplace=True)

    # ── controls ──────────────────────────────────────────────────────
    obs["is_control"] = obs["perturbagen"].isin(_CONTROL_COMPOUNDS) | (
        obs["DOSE_LEVEL"] == "Control"
    )
    obs.loc[obs["is_control"], "pert_dose_uM"] = 0.0

    # ── organism label ────────────────────────────────────────────────
    obs["organism"] = obs["organism"].str.strip().str.lower()

    # ── fixed metadata columns ────────────────────────────────────────
    obs["cell_type"]               = _CELL_TYPE
    obs["pert_type"]               = "compound"
    obs["suspension_type"]         = "cell"
    obs["tissue"]                  = "liver"
    obs["tissue_type"]             = "cell culture"
    obs["disease"]                 = "normal"
    obs["stimulation"]             = None
    obs["guide"]                   = None
    obs["assay"]                   = "SMARTSeq bulk RNA-seq"
    obs["development_stage"]       = _DEVELOPMENT_STAGE
    obs["sex"]                     = _SEX
    obs["self_reported_ethnicity"] = _ETHNICITY
    if "pubchem_cid" not in obs.columns:
        obs["pubchem_cid"]         = None
    obs["split"]                   = obs["SPLIT"]
    obs["psbulk_cells"]            = -666  # sentinel: bulk data, cell count not available
    obs["psbulk_counts"]           = (
        np.round(np.asarray(adata.X.sum(axis=1)).flatten()).astype(int)
    )


    # ── pubchem_cid cleaning ──────────────────────────────────────────
    obs["pubchem_cid"] = pd.to_numeric(obs["pubchem_cid"], errors="coerce").fillna(-666).astype("int64")
    obs["pubchem_cid"] = obs["pubchem_cid"].replace({-666: None, "-666": None, "None": None, "nan": None, "<NA>": None})
    obs["pubchem_cid"] = obs["pubchem_cid"].astype(object).astype("category")

    # ── sample_id (composite key) ─────────────────────────────────────
    obs["sample_id"] = (
        obs["plate"].astype(str) + "_"
        + obs["well"].astype(str) + "_"
        + obs["perturbagen"].astype(str).str.replace(" ", "_", regex=False) + "_"
        + obs["pert_dose_uM"].astype(str)
    )

    return enforce_obs_schema(obs)


def standardize_var_dilimap(adata: AnnData) -> pd.DataFrame:
    """
    Process DILImap variables into a schema-conformant DataFrame.

    Parameters
    ----------
    adata : AnnData
        DILImap dataset (obs standardized or not).

    Returns
    -------
    pd.DataFrame
        var DataFrame with ``ensembl_id`` as a categorical index and a
        single ``symbol`` column.
    """
    var = adata.var.copy()

    # var.index = gene symbols, var['gene_id'] = Ensembl IDs
    var["symbol"] = var.index.astype(str).astype("category")
    var.index = var["gene_id"].astype(str)
    var.index.name = "ensembl_id"
    var.index = var.index.astype("category")

    return var[["symbol"]]



def annotate_pubchem_cids(
    obs: pd.DataFrame,
    cache_path: str,
    manual_mapping_func: Callable[[], Dict] = pubchem_mapping_dilimap,
    drug_col: str = "COMPOUND",
) -> pd.DataFrame:
    """
    Look up PubChem CIDs for unique compounds and merge back into obs.

    Steps
    -----
    1. Extract the deduplicated set of compound names from ``obs[drug_col]``.
    2. Left-join with :func:`dmap.datasets.compound_DILI_labels` to enrich
       each compound with its SMILES string (more reliable for PubChem matching
       than compound name alone).
    3. Run :func:`~src.pseudobulking.common.pubchem.lookup_pubchem_cids` using
       both SMILES and compound name as lookup strategies.
    4. Merge the resulting ``pubchem_cid`` back into obs by ``drug_col``.

    Should be called on the **raw** AnnData obs (before column renaming) so
    that ``drug_col`` still holds the original compound name (``"COMPOUND"``).

    Parameters
    ----------
    obs : pd.DataFrame
        Raw observations dataframe.
    cache_path : str
        Path to a JSON file for persistent PubChem lookup caching. Results are
        loaded from and saved to this file to avoid redundant API calls.
    manual_mapping_func : callable, default :func:`pubchem_mapping_dilimap`
        Zero-argument function returning a ``{"dilimap": {drug_name: cid}}``
        mapping dict.
    drug_col : str, default ``"COMPOUND"``
        Column in ``obs`` containing compound names.

    Returns
    -------
    pd.DataFrame
        obs with ``pubchem_cid`` column populated (int64 dtype).
    """
    obs = obs.copy()
    logger.info("Annotating compounds with PubChem CIDs")

    # 1. Unique compounds from obs
    unique_compounds = (
        obs[[drug_col]]
        .drop_duplicates(drug_col)
        .copy()
        .reset_index(drop=True)
    )
    logger.info(f"  Unique compounds: {len(unique_compounds):,}")

    # 2. Enrich with SMILES from DILImap compound database
    dilimap_compounds = dmap.datasets.compound_DILI_labels()[["compound_name", "smiles"]]
    unique_compounds = unique_compounds.merge(
        dilimap_compounds,
        left_on=drug_col,
        right_on="compound_name",
        how="left",
    ).drop(columns=["compound_name"], errors="ignore")

    n_with_smiles = unique_compounds["smiles"].notna().sum()
    logger.info(f"  Compounds matched with SMILES: {n_with_smiles:,} / {len(unique_compounds):,}")

    # 3. Look up PubChem CIDs using SMILES (preferred) and compound name (fallback)
    logger.info(f"  Using PubChem cache: {cache_path}")
    unique_compounds = lookup_pubchem_cids(
        unique_compounds,
        cache={},
        pert_id_col=drug_col,
        drug_col=drug_col,
        smiles_col="smiles",
        cache_path=cache_path,
        manual_mapping_func=manual_mapping_func,
        manual_mapping_by_drug_name=True,
        dataset_key="dilimap",
    )

    n_mapped = unique_compounds["pubchem_cid"].notna().sum()
    n_missing = unique_compounds["pubchem_cid"].isna().sum()
    logger.info(f"  PubChem CIDs found: {n_mapped:,} / {len(unique_compounds):,}")
    if n_missing:
        missing = unique_compounds.loc[unique_compounds["pubchem_cid"].isna(), drug_col].tolist()
        logger.warning(f"  {n_missing} compound(s) without PubChem CID — consider adding to pubchem_imputation.py: {missing}")

    # 4. Merge pubchem_cid back into obs
    obs = obs.drop(columns=["pubchem_cid"], errors="ignore")
    obs = obs.merge(
        unique_compounds[[drug_col, "pubchem_cid"]],
        on=drug_col,
        how="left",
    )
    return obs


def standardize_dilimap(
    adata: AnnData,
    dataset: str,
    paths: dict,
    manual_mapping_func: Callable[[], Dict] = pubchem_mapping_dilimap,
    annotate_pubchem: bool = False,
) -> AnnData:
    """
    Full standardization of a DILImap AnnData object.

    Returns a fresh AnnData with:
    - ``obs`` indexed by ``sample_id`` (category)
    - ``var`` indexed by ``ensembl_id`` (category), ``symbol`` column only
    - ``X`` as int64 counts
    - ``uns`` cleared to ``{}``

    Parameters
    ----------
    adata : AnnData
        Raw DILImap dataset.
    dataset : str
        Dataset label written to ``obs["dataset"]`` (e.g. ``"dilimap_train"``).
    paths : dict
        Dictionary of file paths as returned by :func:`get_dilimap_paths`.
        Expected keys: ``pubchem_cache``.
    manual_mapping_func : callable, default :func:`pubchem_mapping_dilimap`
        Zero-argument function returning a ``{"dilimap": {drug_name: cid}}``
        manual mapping dict passed to :func:`annotate_pubchem_cids`.
    annotate_pubchem : bool, default False
        If True, enrich compounds with SMILES from the DILImap compound
        database and look up PubChem CIDs via the PubChem API. This involves
        network calls and can be slow; disable for quick reruns.

    Returns
    -------
    AnnData
        Standardized dataset.
    """
    adata = filter_samples(adata)

    if annotate_pubchem:
        Path(paths["pubchem_cache"]).parent.mkdir(parents=True, exist_ok=True)
        adata.obs = annotate_pubchem_cids(
            adata.obs,
            cache_path=str(paths["pubchem_cache"]),
            manual_mapping_func=manual_mapping_func,
            drug_col="COMPOUND",
        )

    obs = standardize_obs_dilimap(adata)
    var = standardize_var_dilimap(adata)

    obs["dataset"] = pd.Categorical([dataset] * len(obs))

    # Cast X to int64
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.astype(np.int64)
    else:
        X = np.array(X, dtype=np.int64)

    adata_out = ad.AnnData(
        X=X,
        obs=obs.set_index("sample_id"),
        var=var,
    )
    adata_out.obs.index = adata_out.obs.index.astype("category")
    adata_out.uns = {}
    return adata_out


