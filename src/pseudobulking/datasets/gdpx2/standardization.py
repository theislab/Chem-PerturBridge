"""
Ginkgo GDPx2 bulk RNA dataset standardization.

GDPx2 is a transcriptomic dataset from The Virtual Cell Pharmacology Initiative (VCPI)
available on LaminLabs (laminlabs/pertdata, key prefix: ginkgo-datapoints/vcpi/).

Unlike vcpi-0001/0002 (THP-1 only), GDPx2 spans four primary cell types:
  - aortic smooth muscle cell  (CL_0002539)
  - epithelial melanocyte      (CL_0002484)
  - fibroblast of dermis       (CL_0002551)
  - skeletal muscle myoblast   (CL_0000515)

Lamindb curation notebook already annotates obs with:
  pert_name, pert_type, pert_compound, pert_dose (string, e.g. "1000.0nM"),
  pert_time (string, e.g. "24.0h"), organism, tissue_type, assay, suspension_type.

This module only adds the columns absent from the curation:
  plate, well, perturbagen (= pert_name), is_control, pert_dose_uM, pert_time_h,
  percent_volume_dmso, tissue, disease, sex, development_stage, self_reported_ethnicity,
  pubchem_cid, library, stimulation, guide, dataset, psbulk_cells, psbulk_counts.

Note: perturbagen is set directly from lamindb-curated pert_name.
No pubchem_name reverse-lookup is needed because pert_name is already
a human-readable compound name.

Data download keys (laminlabs/pertdata):
  obs : ginkgo-datapoints/vcpi/obs.parquet
  var : ginkgo-datapoints/vcpi/var.parquet
  X   : ginkgo-datapoints/vcpi/X.h5ad
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from src.pseudobulking.common.pubchem import is_valid_pubchem_cid, lookup_pubchem_cids
from src.pseudobulking.datasets.gdpx2.pubchem_imputation import pubchem_mapping_gdpx2
from src.utils.parsing_utils import logger


# ---------------------------------------------------------------------------
# Cell-type ontology mapping
# ---------------------------------------------------------------------------

def get_cell_type_ontology_map() -> dict:
    """
    Map lamindb-curated cell_type text labels to Cell Ontology IDs.

    CL IDs use underscore format (CL_XXXXXXX) consistent with the rest of the pipeline.
    """
    return {
        "aortic smooth muscle cell": "CL_0002539",
        "epithelial melanocyte":     "CL_0002484",
        "fibroblast of dermis":      "CL_0002551",
        "skeletal muscle myoblast":  "CL_0000515",
    }


# ---------------------------------------------------------------------------
# .var standardization
# ---------------------------------------------------------------------------

def standardize_var(adata: ad.AnnData) -> ad.AnnData:
    """
    Drop ERCC spike-ins and set the unified var schema.

    Output:
      - var.index : ensembl_id <category>  — stripped of version suffix;
                    NA entries are filled with the symbol value.
      - var['symbol'] : gene symbol <category>

    Input var is expected to have 'ensembl_gene_id' and 'symbol' columns (from lamindb).
    All extra var columns are dropped; only 'symbol' is kept.
    """
    var = adata.var[["symbol"]].copy()
    ensembl_ids = adata.var["ensembl_gene_id"]
    na_mask = ensembl_ids.isna()
    var.index = ensembl_ids.where(~na_mask).astype(str).str.replace(r"\.\d+$", "", regex=True)
    var.index.name = "ensembl_id"

    ercc = var.index.str.startswith("ERCC-", na=False) | var["symbol"].astype(str).str.startswith("ERCC-", na=False)
    if ercc.any():
        logger.info("Dropping %d ERCC spike-in genes", ercc.sum())

    keep = ~ercc.values
    adata = adata[:, keep].copy()
    var = var.iloc[keep]
    na_mask = na_mask.values[keep]

    if na_mask.any():
        logger.info("Filling %d NA index values with symbol", na_mask.sum())
        var.index = pd.Index(np.where(na_mask, var["symbol"].astype(str), var.index), name="ensembl_id")

    dup_mask = var.index.duplicated(keep=False)
    if dup_mask.any():
        logger.info("Filling %d duplicate index values with symbol", dup_mask.sum())
        var.index = pd.Index(np.where(dup_mask, var["symbol"].astype(str), var.index), name="ensembl_id")

    var["symbol"] = var["symbol"].astype("category")
    var.index = pd.Index(var.index.astype(str), name="ensembl_id")
    adata.var = var
    logger.info("Genes after ERCC removal: %d", adata.n_vars)
    return adata


# ---------------------------------------------------------------------------
# Sample filtering
# ---------------------------------------------------------------------------

def filter_samples(adata: ad.AnnData) -> ad.AnnData:
    """
    Filter samples from the GDPx2 dataset.

    Removes:
    1. Samples with total counts <= 0 (empty library).
    2. Positive controls: ``sample_type`` starts with ``'Ginkgo Pos Control'``.
    """
    logger.info("Original: %d samples", adata.n_obs)
    is_outlier = pd.Series(False, index=adata.obs.index)

    lib_sum = np.asarray(adata.X.sum(axis=1)).ravel()
    mask_zero = pd.Series(lib_sum <= 0, index=adata.obs.index)
    logger.info("Zero or negative library size: %d samples", mask_zero.sum())
    is_outlier = is_outlier | mask_zero

    if "sample_type" in adata.obs.columns:
        mask_pos = adata.obs["sample_type"].str.startswith("Ginkgo Pos Control")
        logger.info("Positive controls (sample_type starts with 'Ginkgo Pos Control'): %d samples", mask_pos.sum())
        is_outlier = is_outlier | mask_pos

    logger.info("Total outliers: %d / %d", is_outlier.sum(), adata.n_obs)
    adata = adata[~is_outlier].copy()
    logger.info("Remaining: %d samples", adata.n_obs)
    return adata


# ---------------------------------------------------------------------------
# Compound annotation (PubChem)
# ---------------------------------------------------------------------------

def annotate_gdpx2_pubchem_cids(
    compound_df: pd.DataFrame,
    cache_path: str = "./gdpx2_pubchem_cache.json",
) -> pd.DataFrame:
    """
    Annotate GDPx2 compound metadata with PubChem CIDs.

    compound_df is keyed on 'compound' (internal UUID from obs/CSV).
    pert_name is used for name-based lookups; perturbagen is taken from obs directly.
    Lookup order: inchi_key → smiles → drug name (pert_name) → manual mapping.
    """
    compound_df = compound_df.copy()

    inchikey_col = next(
        (c for c in ("inchi_key", "inchikey", "InChIKey") if c in compound_df.columns), None
    )
    smiles_col = next(
        (c for c in ("smiles", "canonical_smiles", "SMILES") if c in compound_df.columns), None
    )

    cache: dict = {}
    annotated = lookup_pubchem_cids(
        compound_df,
        cache=cache,
        pert_id_col="compound",
        drug_col="pert_name",
        inchikey_col=inchikey_col,
        smiles_col=smiles_col,
        cache_path=cache_path,
        manual_mapping_func=pubchem_mapping_gdpx2,
        manual_mapping_by_drug_name=True,
        dataset_key="gdpx2",
    )
    annotated["pubchem_cid"] = (
        pd.to_numeric(annotated["pubchem_cid"], errors="coerce")
        .fillna(-666)
        .astype("int64")
    )
    annotated["pubchem_cid"] = annotated["pubchem_cid"].replace(
        {-666: None, "-666": None, "None": None, "nan": None, "<NA>": None}
    )
    annotated["pubchem_cid"] = annotated["pubchem_cid"].astype(object).astype("category")
    return annotated


# ---------------------------------------------------------------------------
# .obs standardization helpers
# ---------------------------------------------------------------------------

def define_obs_schema() -> list:
    """Define the strict obs schema for pseudobulk data."""
    return [
        ("sample_id",               "category", "ID: plate_well_perturbagen_cell_type"),
        ("plate",                   "category", "Assay container / plate identifier"),
        ("well",                    "category", "Well position (row{i}_col{j})"),
        ("cell_type",               "category", "Cell type — Cell Ontology ID"),
        ("perturbagen",             "category", "Perturbagen name"),
        ("pert_type",               "category", "Perturbation type"),
        ("is_control",              "category", "bool categories: [False, True]"),
        ("pert_dose_uM",            "float64",  "Dose in µM"),
        ("pert_time_h",             "float64",  "Exposure time in hours"),
        ("percent_volume_dmso",     "float64",  "DMSO vehicle concentration (% v/v)"),
        ("suspension_type",         "category", "Biological material isolation type"),
        ("tissue",                  "category", "Primary tissue/site"),
        ("tissue_type",             "category", "Tissue, cell culture, or organoid"),
        ("disease",                 "category", "Disease / subtype"),
        ("library",                 "category", "Library ID"),
        ("stimulation",             "category", "High-level stimulus"),
        ("guide",                   "category", "CRISPR guide RNA"),
        ("dataset",                 "category", "Dataset label"),
        ("assay",                   "category", "Assay label"),
        ("development_stage",       "category", "Derived from donor age"),
        ("organism",                "category", "Organism"),
        ("sex",                     "category", "Donor sex"),
        ("self_reported_ethnicity", "category", "Donor ethnicity"),
        ("pubchem_cid",             "category", "PubChem CID"),
        ("psbulk_cells",            "int64",    "Total #cells contributing"),
        ("psbulk_counts",           "int64",    "Total #counts contributing"),
    ]


def _parse_pert_dose_uM(pert_dose: pd.Series) -> pd.Series:
    """
    Parse the lamindb pert_dose string column (e.g. '1000.0nM', '3.0uM', '0.15%')
    into a float µM value.

    Rules:
      - nM → divide by 1000
      - uM / µM / um → as-is
      - % → 0.0 (DMSO vehicle; overridden to 0 later anyway)
      - plain number → assumed µM
    """
    def _convert(val):
        if pd.isna(val):
            return np.nan
        s = str(val).strip().lower()
        # percentage (DMSO vehicle)
        if s.endswith("%"):
            return 0.0
        for suffix, factor in [("nm", 1e-3), ("µm", 1.0), ("um", 1.0), ("mm", 1e3)]:
            if s.endswith(suffix):
                try:
                    return float(s[: -len(suffix)].strip()) * factor
                except ValueError:
                    break
        try:
            return float(s)
        except ValueError:
            logger.warning("Could not parse pert_dose: %r", val)
            return np.nan

    return pert_dose.apply(_convert)


def _parse_pert_time_h(pert_time: pd.Series) -> pd.Series:
    """
    Parse the lamindb pert_time string column (e.g. '24.0h') into float hours.
    """
    return pd.to_numeric(
        pert_time.astype(str).str.strip().str.extract(r"(?i)^([\d.]+)\s*h$", expand=False),
        errors="coerce",
    )


def build_obs_columns(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Add the columns that are NOT provided by the lamindb curation notebook.

    Already present from lamindb (left untouched):
      pert_name, pert_type, pert_compound, pert_dose, pert_time,
      organism, tissue_type, assay, suspension_type, cell_type (text).

    Added here:
      plate          — from container_id
      well           — from row_id + column_id
      perturbagen    — directly from pert_name (already curated by lamindb)
      is_control     — pert_name == 'DMSO'
      pert_dose_uM          — parsed from pert_dose string
      pert_time_h           — parsed from pert_time string
      percent_volume_dmso   — passed through from obs (lamindb column)
      cell_type      — remapped from text to CL ontology IDs
      tissue         — 'unknown' (not curated)
      disease        — 'unknown' (not curated)
      sex            — 'unknown' (not curated)
      development_stage       — 'unknown' (not curated)
      self_reported_ethnicity — 'unknown' (not curated)
      library, stimulation, guide — None
      dataset        — 'Ginkgo GDPx2'
      pubchem_cid    — filled later by PubChem annotation step
      sample_id      — plate_well_perturbagen_cell_type
    """
    obs = obs.copy()

    # --- plate ---
    if "container_id" in obs.columns:
        obs["plate"] = obs["container_id"].astype(str)
    elif "plate" not in obs.columns:
        logger.warning("Neither 'container_id' nor 'plate' found — using placeholder")
        obs["plate"] = "plate1"
    obs["plate"] = (
        obs["plate"].astype(str)
        .str.replace(" ", "", regex=False)
        .str.replace("-", "_", regex=False)
    )

    # --- well ---
    if "row_id" in obs.columns and "column_id" in obs.columns:
        obs["well"] = "row" + obs["row_id"].astype(str) + "_col" + obs["column_id"].astype(str)
    elif "well" not in obs.columns:
        raise ValueError(
            "Cannot build 'well': neither 'row_id'/'column_id' nor 'well' found in obs columns. "
            f"Available columns: {list(obs.columns)}"
        )

    # --- perturbagen: lamindb pert_name is already the curated compound name ---
    obs["perturbagen"] = obs["pert_name"].astype(str)

    # --- is_control ---
    obs["is_control"] = (obs["sample_type"] == "Ginkgo Neg Control - DMSO").astype("category")

    non_dmso = (obs["is_control"] == True) & (obs["perturbagen"] != "DMSO")
    if non_dmso.any():
        logger.warning("%d control samples have non-DMSO perturbagen — check filter_samples", non_dmso.sum())

    # --- pert_dose_uM (parse lamindb pert_dose string) ---
    if "pert_dose" in obs.columns:
        obs["pert_dose_uM"] = _parse_pert_dose_uM(obs["pert_dose"])
    else:
        logger.warning("'pert_dose' column not found — filling pert_dose_uM with NaN")
        obs["pert_dose_uM"] = np.nan
    obs.loc[obs["is_control"] == True, "pert_dose_uM"] = 0.0

    # --- pert_time_h (parse lamindb pert_time string) ---
    if "pert_time" in obs.columns:
        obs["pert_time_h"] = _parse_pert_time_h(obs["pert_time"])
    else:
        logger.warning("'pert_time' column not found — filling pert_time_h with NaN")
        obs["pert_time_h"] = np.nan

    # --- percent_volume_dmso: pass through from lamindb obs ---
    obs["percent_volume_dmso"] = pd.to_numeric(obs["percent_volume_dmso"], errors="coerce")
    
    # --- cell_type: text → CL ontology IDs ---
    cell_type_map = get_cell_type_ontology_map()
    obs["cell_type"] = obs["cell_type"].astype(str).str.strip().map(cell_type_map)
    n_unmapped = obs["cell_type"].isna().sum()
    if n_unmapped > 0:
        logger.warning(
            "%d obs have unmapped cell_type — add them to get_cell_type_ontology_map()",
            n_unmapped,
        )

    # --- per-cell-type donor metadata (from paper Methods + ThermoFisher product pages) ---
    #
    # Sources:
    #   CL_0002539  HASMC          Cat. C0075C  lot 2164581 — 27-year-old male donor
    #   CL_0002484  HEMn-LP        Cat. C0025C  — neonatal foreskin → male, newborn
    #   CL_0002551  HDFa           Cat. C0135C  — 45-year-old female donor
    #   CL_0000515  HSkMM          Cat. A11440  — Caucasian male donor
    #
    # All cells are from healthy donors → disease = "normal" for all cell types.
    # Tissues use UBERON IDs (underscore format); ethnicity unknown for all.

    _tissue_map = {
        "CL_0002539": "aorta",                  # UBERON_0000947
        "CL_0002484": "skin",                   # UBERON_0002097
        "CL_0002551": "skin",                   # UBERON_0002067
        "CL_0000515": "skeletal muscle",        # UBERON_0001134
    }
    _sex_map = {
        "CL_0002539": "male",    # CoA C0075C lot 2164581 — 27-year-old male donor
        "CL_0002484": "male",    # neonatal foreskin donor
        "CL_0002551": "female",  # HDFa Cat. C0135C — 45-year-old female donor
        "CL_0000515": "male",   # Caucasian male donor
    }
    _dev_stage_map = {
        "CL_0002539": "27-year-old stage",  # CoA C0075C lot 2164581
        "CL_0002484": "0-year-old stage",   # neonatal foreskin donor
        "CL_0002551": "45-year-old stage",  # HDFa Cat. C0135C
        "CL_0000515": "unknown",            # HSkMM Cat. A11440 — age not specified
    }

    obs["tissue"]                  = obs["cell_type"].map(_tissue_map).fillna("unknown")
    obs["sex"]                     = obs["cell_type"].map(_sex_map).fillna("unknown")
    _ethnicity_map = {
        "CL_0002539": "unknown",
        "CL_0002484": "unknown",
        "CL_0002551": "unknown",
        "CL_0000515": "Caucasian",  # HSkMM Cat. A11440
    }

    obs["development_stage"]       = obs["cell_type"].map(_dev_stage_map).fillna("unknown")
    obs["disease"]                 = "normal"
    obs["self_reported_ethnicity"] = obs["cell_type"].map(_ethnicity_map).fillna("unknown")
    obs["library"]                 = None
    obs["stimulation"]             = None
    obs["guide"]                   = None
    obs["dataset"]                 = "Ginkgo GDPx2"

    if "pubchem_cid" not in obs.columns:
        obs["pubchem_cid"] = None
    obs["pubchem_cid"] = (
        pd.to_numeric(obs["pubchem_cid"], errors="coerce")
        .fillna(-666)
        .astype("int64")
    )
    obs["pubchem_cid"] = obs["pubchem_cid"].replace(
        {-666: None, "-666": None, "None": None, "nan": None, "<NA>": None}
    )
    obs["pubchem_cid"] = obs["pubchem_cid"].astype(object).astype("category")

    # --- sample_id ---
    obs["sample_id"] = (
        obs["plate"].astype(str).str.replace(" ", "", regex=False) + "_"
        + obs["well"].astype(str).str.replace(" ", "", regex=False) + "_"
        + obs["perturbagen"].str.replace(" ", "_", regex=False) + "_"
        + obs["cell_type"].astype(str).str.replace(" ", "_", regex=False)
    )
    return obs


def calculate_psbulk_counts(adata: ad.AnnData) -> np.ndarray:
    """Sum counts per observation from the expression matrix."""
    psbulk_counts = adata.X.sum(axis=1)
    if hasattr(psbulk_counts, "A1"):
        psbulk_counts = psbulk_counts.A1
    return np.array(psbulk_counts, dtype=int)


def enforce_obs_schema(obs: pd.DataFrame) -> pd.DataFrame:
    """Enforce strict obs schema: add missing columns, cast dtypes, select schema columns."""
    obs_schema = define_obs_schema()
    dtype_map   = {col: dtype for col, dtype, _ in obs_schema}
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
            obs_out[col] = pd.to_numeric(obs_out[col], errors="coerce")
        elif dtype == "int64":
            obs_out[col] = pd.to_numeric(obs_out[col], errors="coerce").astype("int64")

    obs_out = obs_out[schema_cols]
    return obs_out


def process_obs_dataframe(adata: ad.AnnData, compound_df: pd.DataFrame) -> pd.DataFrame:
    """
    Process GDPx2 observations to the unified pseudobulk schema.

    Steps:
    1. Merge pubchem_cid from compound_df (keyed on compound UUID).
    2. Build all missing columns (plate, well, perturbagen=pert_name, is_control,
       pert_dose_uM, pert_time_h, cell_type CL IDs, tissue/disease/sex/etc.).
    3. Set psbulk_cells = -666 (bulk); compute psbulk_counts from X.
    4. Enforce strict obs schema.
    """
    obs = adata.obs.reset_index()
    if obs.columns[0] != "sequenced_id":
        obs = obs.rename(columns={obs.columns[0]: "sequenced_id"})

    # Bring pubchem_cid from compound_df (keyed on compound UUID)
    if "pubchem_cid" in compound_df.columns:
        obs = obs.drop(columns=["pubchem_cid"], errors="ignore")
        obs = obs.merge(
            compound_df[["compound", "pubchem_cid"]].drop_duplicates("compound"),
            on="compound",
            how="left",
        )

    obs = build_obs_columns(obs)

    obs["psbulk_cells"]  = -666
    obs["psbulk_counts"] = calculate_psbulk_counts(adata)

    obs_final = enforce_obs_schema(obs)
    logger.info("Processed %d observations", len(obs_final))
    return obs_final


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def standardize_gdpx2_dataset(
    paths: dict,
    annotate_pubchem: bool = False,
) -> ad.AnnData:
    """
    Full standardization pipeline for the Ginkgo GDPx2 bulk RNA dataset.

    Steps
    -----
    1. Load X.h5ad and obs/var parquets from the LaminLabs download.
    2. Attach the curated obs/var from parquets (overriding the embedded ones).
    3. Standardize .var (drop ERCC, set ensembl_id index).
    4. Filter outlier samples (empty libraries, positive controls).
    5. Load compound_df from compound_csv (downloaded from LaminLabs); falls
       back to deriving it from obs if the file is absent.
    6. Optionally annotate PubChem CIDs (annotate_pubchem).
    7. Build unified .obs conforming to the pseudobulk schema.
       perturbagen is set directly from lamindb-curated pert_name.
    8. Convert .X to int64.

    Parameters
    ----------
    paths:
        Dictionary from get_gdpx2_paths(). Expected keys:
        raw_h5ad, raw_obs, raw_var, pubchem_cid_cache.
    annotate_pubchem:
        If True, run PubChem CID lookups (requires network access).

    Returns
    -------
    Standardized AnnData with obs.index = sample_id and var.index = ensembl_id.
    """
    logger.info("Loading count matrix from %s", paths["raw_h5ad"])
    adata = ad.read_h5ad(paths["raw_h5ad"])

    logger.info("Loading obs from %s", paths["raw_obs"])
    obs_raw = pd.read_parquet(paths["raw_obs"])

    logger.info("Loading var from %s", paths["raw_var"])
    var_raw = pd.read_parquet(paths["raw_var"])

    # Use lamindb-curated obs/var (they contain pert_name, pert_dose, etc.)
    obs_raw.index = obs_raw.index.astype(str)
    var_raw.index = var_raw.index.astype(str)
    adata.obs = obs_raw.reindex(adata.obs_names.astype(str))
    adata.var = var_raw.reindex(adata.var_names.astype(str))

    if "is_control" in adata.obs.columns:
        adata.obs["is_control"] = adata.obs["is_control"].astype(bool)

    logger.info("Standardizing .var ...")
    adata = standardize_var(adata)

    logger.info("Filtering samples ...")
    adata = filter_samples(adata)

    # Load compound_df from the dedicated compound CSV, then restrict it to the
    # compounds present in adata.obs (left join: obs compounds → CSV rows).
    # Falls back to deriving compound_df from obs if the file was not downloaded.
    if paths.get("compound_csv") and Path(paths["compound_csv"]).is_file():
        logger.info("Loading compound metadata from %s", paths["compound_csv"])
        compound_df_full = pd.read_csv(paths["compound_csv"])
        logger.info("Compound CSV: %d rows total", len(compound_df_full))
        obs_compound_map = (
            adata.obs[["compound", "pert_name"]]
            .drop_duplicates("compound")
            .reset_index(drop=True)
        )
        compound_df = compound_df_full.merge(obs_compound_map, on="compound", how="left")
        logger.info(
            "Compounds after filtering to obs: %d / %d",
            len(compound_df), len(compound_df_full),
        )
    else:
        raise FileNotFoundError(
            f"Compound CSV not found: {paths.get('compound_csv')}. "
            "Run download_compound_csv() before standardization."
        )

    if annotate_pubchem:
        logger.info("Annotating compounds with PubChem CIDs ...")
        compound_df = annotate_gdpx2_pubchem_cids(
            compound_df,
            cache_path=str(paths["pubchem_cid_cache"]),
        )

    logger.info("Standardizing .obs ...")
    obs_standardized = process_obs_dataframe(adata, compound_df).set_index("sample_id")

    var_standardized = adata.var.copy()
    var_standardized.index = var_standardized.index.astype("category")

    logger.info("Creating standardized AnnData ...")
    adata_standardized = ad.AnnData(
        X=adata.X.astype("int64"),
        obs=obs_standardized,
        var=var_standardized,
    )

    logger.info(
        "Standardization complete: %d observations × %d genes",
        adata_standardized.n_obs, adata_standardized.n_vars,
    )
    return adata_standardized
