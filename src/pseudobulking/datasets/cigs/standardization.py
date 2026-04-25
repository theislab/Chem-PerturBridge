"""CIGS (Chemical-Induced Gene Signature) dataset standardization."""
import re
import shutil
import string
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Literal, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd

from src.pseudobulking.common.pubchem import lookup_pubchem_cids
from src.pseudobulking.datasets.cigs.pubchem_imputation import pubchem_mapping_cigs
from src.pseudobulking.datasets.cigs.ensembl_symbol_lookup import (
    fetch_ensembl_ids_for_symbols,
)
from src.utils.parsing_utils import *


# ---------------------------------------------------------------------------
# Cell line → Cellosaurus / tissue / donor metadata
# ---------------------------------------------------------------------------

# Maps every raw cell_id variant observed in CIGS metadata to a canonical cell
# line name. MCE uses "cell_293T" / "cell_231" / "cell_MDA_MB_231"; TCM uses
# "293T" / "231" / "MDA_MB_231" (sometimes as integers). Keys are the string
# representation of the raw value.
_CELL_ID_MAP: Dict[str, str] = {
    "cell_293T":       "HEK293T",
    "293T":            "HEK293T",
    "293":             "HEK293T",
    "cell_MDA_MB_231": "MDA_MB_231",
    "cell_231":        "MDA_MB_231",
    "MDA_MB_231":      "MDA_MB_231",
    "231":             "MDA_MB_231",
}

_CELL_LINE_META: Dict[str, dict] = {
    "HEK293T": {
        "cell_type":               "CVCL_0063",
        "tissue":                  "kidney",
        "disease":                 "normal",
        "sex":                     "female",
        "development_stage":       "fetal stage",
        "self_reported_ethnicity": "unknown",
    },
    "MDA_MB_231": {
        "cell_type":               "CVCL_0062",
        "tissue":                  "breast",
        "disease":                 "breast adenocarcinoma",
        "sex":                     "female",
        "development_stage":       "51-year-old stage",
        "self_reported_ethnicity": "Caucasian",
    },
}


# ---------------------------------------------------------------------------
# Subset-level constants (dose / time are fixed per subset_key by design)
# ---------------------------------------------------------------------------

# (dose_uM, time_h) for every CIGS subset_key. These are design parameters of
# the library, not parsed from individual metadata rows — using a lookup keeps
# the numeric values robust to placeholder strings like '-666' or 'NA' and
# keeps DMSO wells numerically coherent with their plate.
_SUBSET_DOSE_TIME: Dict[str, Dict[str, float]] = {
    "mce_hek293t_10uM":    {"dose_uM": 10.0, "time_h": 24.0},
    "mce_mda_mb_231_10uM": {"dose_uM": 10.0, "time_h": 24.0},
    "tcm_hek293t_10":      {"dose_uM": 10.0, "time_h": 24.0},
    "tcm_hek293t_20":      {"dose_uM": 20.0, "time_h": 24.0},
    "tcm_mda_mb_231_10":   {"dose_uM": 10.0, "time_h": 24.0},
    "tcm_mda_mb_231_20":   {"dose_uM": 20.0, "time_h": 24.0},
}


# Grabs the numeric part of strings like "10.0 uM", "20uM", "24 h", "24h".
_NUMBER_RE = re.compile(r"([\d.]+)")

# Target field in _SUBSET_DOSE_TIME → raw per-sample columns to parse from.
# MCE uses ``pert_idose`` / ``pert_itime``; TCM uses ``Dose`` / ``Time``.
_DOSE_TIME_SOURCES: Dict[str, Tuple[str, ...]] = {
    "dose_uM": ("pert_idose", "Dose"),
    "time_h":  ("pert_itime", "Time"),
}


def _assert_dose_time_consistent(obs: pd.DataFrame) -> None:
    """Raise if raw per-sample dose/time disagrees with ``_SUBSET_DOSE_TIME``."""
    if "_subset_key" not in obs.columns:
        return

    def _parse(*cols: str) -> pd.Series:
        out = pd.Series(np.nan, index=obs.index, dtype=float)
        for c in cols:
            if c in obs.columns:
                raw = obs[c].astype(str).str.strip()
                # CIGS uses "-666" as a "not applicable" sentinel for control
                # wells (Blank / DMSO / RNA) in both pert_idose / pert_itime
                # (MCE) and Dose / Time (TCM). Mask it out so the regex below
                # doesn't parse it as a literal 666.
                raw = raw.mask(raw.eq("-666"))
                nums = raw.str.extract(_NUMBER_RE, expand=False)
                out = out.combine_first(nums.astype(float))
        return out

    for field, sources in _DOSE_TIME_SOURCES.items():
        raw = _parse(*sources)
        tbl = obs["_subset_key"].map(
            lambda k, f=field: _SUBSET_DOSE_TIME.get(k, {}).get(f, np.nan)
        )
        mask = raw.notna() & tbl.notna() & (raw != tbl)
        if mask.any():
            raise ValueError(
                f"{field}: raw metadata disagrees with _SUBSET_DOSE_TIME for "
                f"{int(mask.sum()):,} samples "
                f"(subset_key={obs.loc[mask, '_subset_key'].iloc[0]!r}, "
                f"raw={raw[mask].iloc[0]}, table={tbl[mask].iloc[0]}). "
                f"Update _SUBSET_DOSE_TIME."
            )


# ---------------------------------------------------------------------------
# String column materialization
# ---------------------------------------------------------------------------

def materialize_string_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Convert pandas StringDtype columns to object dtype for AnnData compatibility."""
    for col in df.select_dtypes(include="string").columns:
        df[col] = df[col].astype(object).where(pd.notna(df[col]), np.nan)

    if pd.api.types.is_string_dtype(df.index.dtype):
        df.index = df.index.astype(object)
    return df


# ---------------------------------------------------------------------------
# Column normalization: MCE- and TCM-specific obs normalizers
# ---------------------------------------------------------------------------

def normalize_mce_obs(obs: pd.DataFrame) -> pd.DataFrame:
    """Map MCE obs columns to unified names: treatment→catalog_id,
    sample_plate→plate (kept with ``_X\\d+`` replicate tag),
    sample_row+sample_column→well (e.g. row 4, col 10 → ``D10``)."""
    obs = obs.copy()
    obs["catalog_id"] = obs["treatment"]
    # Keep the full plate ID including the ``_X1`` / ``_X2`` replicate suffix
    # so each physical replicate plate is identifiable downstream.
    obs["plate"] = obs["sample_plate"]

    def _well(row):
        try:
            letter = string.ascii_uppercase[int(row["sample_row"]) - 1]
            return f"{letter}{int(row['sample_column'])}"
        except (ValueError, IndexError, KeyError):
            return np.nan

    obs["well"] = obs.apply(_well, axis=1)
    return obs


def normalize_tcm_obs(obs: pd.DataFrame) -> pd.DataFrame:
    """Map TCM obs columns to unified names: Treat→catalog_id (``Compd``→``Cpd``
    prefix rewrite to match ``compounds_tcm.csv``), Plate→plate and
    Sample_unique_id→sample_unique_id (both prefixed with ``{subset_key}.``
    since the 4 TCM subsets reuse plate/sample IDs). ``well`` is left NaN
    (HiMAP-seq is well-less)."""
    obs = obs.copy()
    obs["catalog_id"] = obs["Treat"].str.replace(r"^Compd", "Cpd", regex=True)
    obs["cell_id"]    = obs["Cell"]

    subset_prefix = obs["_subset_key"].astype(str) + "."

    obs["plate"] = (
        subset_prefix + obs["Plate"].astype(str)
        if "Plate" in obs.columns else np.nan
    )
    obs["well"] = np.nan
    # obs.index is already "{subset_key}.{Sample_unique_id}" upstream; only
    # the legacy column path still needs the manual prefix.
    if obs.index.name == "Sample_unique_id":
        obs["sample_unique_id"] = obs.index.astype(str)
    elif "Sample_unique_id" in obs.columns:
        obs["sample_unique_id"] = subset_prefix + obs["Sample_unique_id"].astype(str)
    return obs


# ---------------------------------------------------------------------------
# .var standardization (gene symbol → Ensembl ID)
# ---------------------------------------------------------------------------

# Fix Excel/R-mangled symbols to canonical HGNC form so MCE and TCM share an
# identical gene axis.
_GENE_SYMBOL_FIXES: Dict[str, str] = {
    "2024-09-07 00:00:00": "SEPT7",
    "X7.Sep":               "SEPT7",
    "HLA.A":   "HLA-A",
    "HLA.B":   "HLA-B",
    "HLA.C":   "HLA-C",
    "HLA.DMA": "HLA-DMA",
    "HLA.DPA1":"HLA-DPA1",
    "HLA.DRA": "HLA-DRA",
    "HLA.E":   "HLA-E",
    "NKX2.1":  "NKX2-1",
    "NKX3.1":  "NKX3-1",
}


def _canonicalize_gene_symbols(symbols: list) -> Tuple[list, int]:
    """Apply ``_GENE_SYMBOL_FIXES`` to ``symbols``; return (fixed, n_changed)."""
    fixed = [_GENE_SYMBOL_FIXES.get(s, s) for s in symbols]
    n_fixed = sum(1 for a, b in zip(symbols, fixed) if a != b)
    return fixed, n_fixed


def process_gene_annotations(
    var: pd.DataFrame,
    paths: dict,
    annotate_genes: bool = True,
) -> pd.DataFrame:
    """Convert var from HGNC symbols to ``var.index = ensembl_id``,
    ``var["symbol"] = symbol`` via the Ensembl REST API. Unresolved symbols
    keep the symbol as their index placeholder."""
    logger.info("Processing gene annotations")

    raw_symbols = var.index.astype(str).tolist()
    symbols, n_fixed = _canonicalize_gene_symbols(raw_symbols)
    if n_fixed:
        logger.info(f"  Canonicalized {n_fixed} Excel/R-mangled gene symbols")
    logger.info(f"  Resolving Ensembl IDs for {len(symbols):,} gene symbols ...")

    if annotate_genes:
        ensembl_map = fetch_ensembl_ids_for_symbols(
            symbols,
            cache_path=str(paths["ensembl_cache"]),
        )
    else:
        ensembl_map = {s: None for s in symbols}

    ensembl_series = pd.Series(symbols, index=symbols).map(ensembl_map)

    n_mapped   = ensembl_series.notna().sum()
    n_unmapped = ensembl_series.isna().sum()
    logger.info(f"  Symbol→Ensembl mapping: {n_mapped:,} mapped, {n_unmapped:,} unmapped")

    if n_unmapped:
        logger.warning(
            f"  {n_unmapped:,} symbols did not resolve; their symbol string will "
            "be used as the Ensembl ID placeholder"
        )
        ensembl_series = ensembl_series.fillna(pd.Series(symbols, index=symbols))

    # Keep the ``var`` row order aligned with the incoming var.
    ensembl_ids = ensembl_series.values

    var_df = pd.DataFrame(
        {"symbol": pd.Categorical(symbols)},
        index=pd.Index(ensembl_ids, name="ensembl_id"),
    )

    dupes = var_df.index.duplicated().sum()
    if dupes:
        logger.warning(f"  {dupes:,} duplicated Ensembl IDs after mapping — keeping all rows")

    var_df = materialize_string_columns(var_df)

    logger.info(f"  Genes after annotation: {len(var_df):,}")
    return var_df


# ---------------------------------------------------------------------------
# Sample filtering (delegated to the paper's R code)
# ---------------------------------------------------------------------------

_QC_RSCRIPT_PATH = Path(__file__).parent / "qc_filter.R"


def _check_rscript(rscript: str) -> None:
    """Verify that ``Rscript`` is available on PATH."""
    if shutil.which(rscript) is None:
        raise RuntimeError(
            f"'{rscript}' not found on PATH. "
            "Make sure R is installed and the conda environment is active."
        )


def _raw_treatment(obs: pd.DataFrame) -> pd.Series:
    """Combine MCE ``treatment`` and TCM ``Treat`` columns into one Series."""
    treat = pd.Series(pd.NA, index=obs.index, dtype=object)
    if "treatment" in obs.columns:
        treat = treat.combine_first(obs["treatment"])
    if "Treat" in obs.columns:
        treat = treat.combine_first(obs["Treat"])
    return treat


def filter_samples(
    adata: ad.AnnData,
    rscript: str = "Rscript",
    work_dir: Optional[Path] = None,
) -> ad.AnnData:
    """Quality-filter CIGS samples by delegating to :mod:`qc_filter.R` (BLANK/RNA
    drop + per-subset IQR outlier rule, mirroring the paper's R code).
    Python side drops zero-library-size samples first, then exchanges
    ``(sample_id, subset_key, total_reads, treatment, sample_plate)`` via
    TSV and keeps only ``status == "pass"``."""
    logger.info(f"Original: {adata.n_obs:,} samples")

    lib_sum = np.asarray(adata.X.sum(axis=1)).ravel()
    mask_zero = lib_sum <= 0
    n_zero = int(mask_zero.sum())
    if n_zero:
        logger.warning(f"Removing {n_zero:,} samples with zero or negative library size")
        adata = adata[~mask_zero].copy()
        lib_sum = lib_sum[~mask_zero]

    _check_rscript(rscript)

    # sample_plate is only present on MCE (e.g. "MCE4_293T_24H_X1"); the R
    # side uses it to reproduce the paper's hard-coded plate exclusions.
    # Defaults to "" for subsets that don't carry the column (TCM).
    sample_plate = (
        adata.obs["sample_plate"].astype(str).fillna("")
        if "sample_plate" in adata.obs.columns
        else pd.Series("", index=adata.obs_names)
    )

    qc_df = pd.DataFrame({
        "sample_id":    adata.obs_names.astype(str),
        "subset_key":   adata.obs["_subset_key"].values,
        "total_reads":  lib_sum.astype(np.int64),
        "treatment":    _raw_treatment(adata.obs).values,
        "sample_plate": sample_plate.values,
    })

    ctx = tempfile.TemporaryDirectory() if work_dir is None else None
    try:
        tmp_root = Path(ctx.name) if ctx is not None else Path(work_dir)
        tmp_root.mkdir(parents=True, exist_ok=True)
        input_tsv  = tmp_root / "cigs_qc_input.tsv"
        output_tsv = tmp_root / "cigs_qc_status.tsv"

        qc_df.to_csv(input_tsv, sep="\t", index=False)
        logger.info(f"Running qc_filter.R on {len(qc_df):,} samples ...")
        subprocess.run(
            [
                rscript, str(_QC_RSCRIPT_PATH),
                "--input_tsv",  str(input_tsv),
                "--output_tsv", str(output_tsv),
            ],
            check=True,
            text=True,
        )

        status = pd.read_csv(output_tsv, sep="\t", dtype={"sample_id": str, "status": str})
    finally:
        if ctx is not None:
            ctx.cleanup()

    status = status.set_index("sample_id").reindex(adata.obs_names.astype(str))
    if status["status"].isna().any():
        missing = int(status["status"].isna().sum())
        raise RuntimeError(
            f"qc_filter.R returned no status for {missing:,} samples "
            "(sample_id mismatch between Python and R)."
        )

    counts = status["status"].value_counts().to_dict()
    logger.info("QC status: " + ", ".join(f"{k}={v:,}" for k, v in sorted(counts.items())))

    keep = (status["status"].values == "pass")
    adata = adata[keep].copy()

    logger.info(f"Remaining: {adata.n_obs:,} samples")
    return adata


# ---------------------------------------------------------------------------
# Compound metadata: merging the paper's supplementary table into obs
# ---------------------------------------------------------------------------

def load_compound_table(csv_path: Path) -> pd.DataFrame:
    """Load the paper's compound annotation CSV with columns
    ``compound_name``, ``catalog_number``, ``cas_number`` (+ MCE extras
    ``moa``, ``clinical_info``, ``approved_type``)."""
    if not Path(csv_path).is_file():
        raise FileNotFoundError(
            f"Compound annotation CSV not found: {csv_path}. "
            "Run extract_compound_tables() from the download step first."
        )
    df = pd.read_csv(csv_path, dtype=str)
    df = df.dropna(subset=["catalog_number"]).drop_duplicates("catalog_number")
    logger.info(f"  Loaded compound table {csv_path.name}: {len(df):,} compounds")
    return df


def merge_compound_metadata(
    obs: pd.DataFrame,
    compound_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge the paper's compound annotations (``compound_name``, ``cas_number``)
    onto obs by ``catalog_id``. For rows with no table match (controls,
    unlisted entries), the catalog ID itself is used as ``compound_name``
    so downstream perturbagens stay dense."""
    obs = obs.copy()
    keep_cols = [c for c in ("catalog_number", "compound_name", "cas_number")
                 if c in compound_df.columns]
    merged = obs.merge(
        compound_df[keep_cols],
        how="left",
        left_on="catalog_id",
        right_on="catalog_number",
    )
    if "catalog_number" in merged.columns:
        merged = merged.drop(columns=["catalog_number"])

    # Fall back to the catalog ID when the compound table does not cover the
    # value (controls, rare library entries); this keeps compound_name dense.
    missing = merged["compound_name"].isna()
    if missing.any():
        merged.loc[missing, "compound_name"] = merged.loc[missing, "catalog_id"]
        logger.info(
            f"  {int(missing.sum()):,} rows without a compound-table hit "
            "(controls / unlisted) — using catalog_id as compound_name"
        )

    return merged


def extract_compound_df(obs: pd.DataFrame) -> pd.DataFrame:
    """Build a unique-compound table for PubChem lookup, deduped by
    ``catalog_id``. ``compound_name`` and ``cas_number`` are kept as the
    actual lookup identifiers."""
    cols = [c for c in ("catalog_id", "compound_name", "cas_number") if c in obs.columns]
    df = obs[cols].drop_duplicates("catalog_id").copy()
    logger.info(f"  Extracted {len(df):,} unique compounds from obs")
    return df


# ---------------------------------------------------------------------------
# PubChem annotation
# ---------------------------------------------------------------------------

def annotate_pubchem_cids(
    compound_df: pd.DataFrame,
    cache_path: str = "./cigs_pubchem_cache.json",
) -> pd.DataFrame:
    """Annotate CIGS compounds with PubChem CIDs via :func:`lookup_pubchem_cids`.

    Lookup order: manual override → CAS → compound_name. Universal cache
    key is ``catalog_id``; ``cas:<cas>`` and ``name:<name>`` namespaces
    are also populated so CAS-resolved CIDs are reused across catalog IDs.
    """
    cache: dict = {}
    annotated = lookup_pubchem_cids(
        compound_df,
        cache=cache,
        pert_id_col="catalog_id",
        drug_col="compound_name",
        inchikey_col=None,
        smiles_col=None,
        cas_col="cas_number" if "cas_number" in compound_df.columns else None,
        cache_path=cache_path,
        manual_mapping_func=pubchem_mapping_cigs,
        manual_mapping_by_drug_name=False,
        dataset_key="cigs",
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
# .obs standardization
# ---------------------------------------------------------------------------

def define_obs_schema() -> list:
    """Define the strict obs schema for pseudobulk data."""
    return [
        ("sample_id",               "category", "Unique sample identifier"),
        ("plate",                   "category", "Assay detection plate identifier"),
        ("well",                    "category", "Well ID on the plate"),
        ("cell_type",               "category", "Cell type ontology ID (Cellosaurus CVCL)"),
        ("perturbagen",             "category", "Perturbagen name"),
        ("pert_type",               "category", "Perturbation type"),
        ("is_control",              "category", "True for DMSO / vehicle controls"),
        ("pert_dose_uM",            "float64",  "Dose in µM (MCE and TCM are both dosed in µM)"),
        ("pert_time_h",             "float64",  "Exposure time in hours"),
        ("suspension_type",         "category", "Type of biological material"),
        ("tissue",                  "category", "Primary tissue name"),
        ("tissue_type",             "category", "Type of tissue"),
        ("disease",                 "category", "Disease or 'normal'"),
        ("library",                 "category", "Library ID"),
        ("stimulation",             "category", "High-level stimulus"),
        ("guide",                   "category", "CRISPR guide RNA"),
        ("dataset",                 "category", "Dataset label"),
        ("assay",                   "category", "Assay platform (HTS2 for MCE, HiMAP-seq for TCM)"),
        ("development_stage",       "category", "Development stage"),
        ("organism",                "category", "Organism"),
        ("sex",                     "category", "Donor sex"),
        ("self_reported_ethnicity", "category", "Donor ethnicity"),
        ("pubchem_cid",             "category", "PubChem CID"),
        ("psbulk_cells",            "int64",    "Total cells contributing (−666 for bulk)"),
        ("psbulk_counts",           "int64",    "Total counts per sample"),
    ]


def define_tcm_obs_schema() -> list:
    """TCM-only extra columns appended to the strict obs schema."""
    return [
        ("sample_unique_id", "category", "Upstream-assigned running sample counter"),
    ]


def add_fixed_metadata_columns(obs: pd.DataFrame, source: Literal["mce", "tcm"]) -> pd.DataFrame:
    """Add fixed CIGS metadata (cell-line ontology, assay, dose/time from
    ``_SUBSET_DOSE_TIME``, controls) and build ``sample_id`` from plate/
    well/perturbagen/cell_type."""
    obs = obs.copy()

    cell_line = obs["cell_id"].astype(str).str.strip().map(_CELL_ID_MAP)
    unknown_cells = obs.loc[cell_line.isna(), "cell_id"].astype(str).unique()
    if len(unknown_cells):
        logger.warning(f"Unknown cell_id values (will produce NaN metadata): {unknown_cells}")

    for col in ["cell_type", "tissue", "disease", "sex", "development_stage",
                "self_reported_ethnicity"]:
        obs[col] = cell_line.map({k: v[col] for k, v in _CELL_LINE_META.items()})

    # CIGS doesn't fetch PubChem canonical names, so perturbagen is just
    # compound_name. Rename rather than duplicate the column.
    obs = obs.rename(columns={"compound_name": "perturbagen"})
    
    obs["organism"]        = "human"
    obs["assay"]           = "HTS2" if source == "mce" else "HiMAP-seq"
    obs["tissue_type"]     = "cell culture"
    obs["suspension_type"] = "cell"
    obs["pert_type"]       = "compound"
    obs["dataset"]         = "CIGS MCE" if source == "mce" else "CIGS TCM"
    obs["library"]         = None
    obs["stimulation"]     = None
    obs["guide"]           = None

    # Dose (µM) and time (h) are fixed per (library × cell line × dose) subset.
    subset_keys = obs["_subset_key"]
    dose = subset_keys.map(lambda k: _SUBSET_DOSE_TIME.get(k, {}).get("dose_uM", np.nan))
    time = subset_keys.map(lambda k: _SUBSET_DOSE_TIME.get(k, {}).get("time_h", np.nan))

    unknown = subset_keys[~subset_keys.isin(_SUBSET_DOSE_TIME)].unique()
    if len(unknown):
        logger.warning(f"Unknown _subset_key values (dose/time = NaN): {list(unknown)}")

    

    # DMSO is the only vehicle control in CIGS; BLANK/RNA wells are dropped
    # upstream by qc_filter.R, so they should never reach this function.
    is_ctrl = obs["perturbagen"].str.upper() == "DMSO"
    dose = dose.mask(is_ctrl, 0.0)

    obs["pert_dose_uM"] = dose.astype(float)
    obs["pert_time_h"]  = time.astype(float)
    obs["is_control"]   = is_ctrl.astype(bool).astype("category")

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


    # MCE: {plate}_{well}_{perturbagen}_{cell_type}
    # TCM: {sample_unique_id}.{raw_plate}_{perturbagen}_{cell_type}  (no wells;
    # `sample_unique_id` already carries the `{subset_key}.` prefix, so we
    # strip it from `plate` to avoid duplicating it).
    if "sample_unique_id" in obs.columns:
        raw_plate = obs["plate"].astype(str).str.split(".", n=1).str[1]
        prefix = obs["sample_unique_id"].astype("object").fillna("") + "." + raw_plate.fillna("")
    else:
        # .astype("object") avoids Categorical fillna TypeError.
        plate_str = obs["plate"].astype("object").fillna("")
        well_str  = obs["well"].astype("object").fillna("")
        prefix = plate_str + "_" + well_str

    obs["sample_id"] = (
        prefix + "_" +
        obs["perturbagen"].fillna("").str.replace(" ", "_", regex=False) + "_" +
        obs["cell_type"].fillna("").str.replace(" ", "_", regex=False)
    )

    dups = obs["sample_id"].duplicated(keep=False)
    if dups.any():
        logger.warning(
            f"{int(dups.sum())} duplicate sample_ids detected. "
            "Consider adding more distinguishing info (e.g. dose or batch) to sample_id."
        )

    return obs


def calculate_psbulk_counts(adata: ad.AnnData) -> np.ndarray:
    """Sum counts per observation from the expression matrix."""
    counts = adata.X.sum(axis=1)
    if hasattr(counts, "A1"):
        counts = counts.A1
    return np.asarray(counts, dtype=int)


def enforce_obs_schema(
    obs: pd.DataFrame,
    extra_schema: Optional[list] = None,
) -> pd.DataFrame:
    """Add missing columns, cast dtypes, and select schema columns.
    ``extra_schema`` is an optional per-source schema appended after the
    strict one; extras are only kept if present in ``obs``."""
    obs_schema = define_obs_schema() + [
        t for t in (extra_schema or []) if t[0] in obs.columns
    ]
    dtype_map   = {col: dt for col, dt, _ in obs_schema}
    schema_cols = [col for col, _, _ in obs_schema]

    out = obs.copy()
    for col, dtype in dtype_map.items():
        if col not in out.columns:
            out[col] = np.nan
        if col == "is_control":
            out[col] = out[col].astype(bool).astype("category")
        elif dtype == "category":
            out[col] = out[col].astype(object).astype("category")
        elif dtype == "float64":
            out[col] = pd.to_numeric(out[col], errors="coerce").astype("float64")
        elif dtype == "int64":
            out[col] = pd.to_numeric(out[col], errors="coerce").astype("int64")

    out = out[schema_cols]
    out = materialize_string_columns(out)
    return out


def process_obs_dataframe(
    adata: ad.AnnData,
    paths: dict,
    source: Literal["mce", "tcm"],
    annotate_pubchem: bool = False,
) -> pd.DataFrame:
    """Process one source's obs into the unified schema: normalize columns,
    merge compound table, annotate PubChem (optional), build perturbagen/
    sample_id, set psbulk stats, and enforce the strict schema."""
    if source not in ("mce", "tcm"):
        raise ValueError(f"source must be 'mce' or 'tcm'; got {source!r}")

    logger.info(f"Processing CIGS observations dataframe (source={source})")

    raw_obs = adata.obs.copy()
    if source == "mce":
        obs = normalize_mce_obs(raw_obs)
    else:
        obs = normalize_tcm_obs(raw_obs)

    # preserve _subset_key for add_fixed_metadata_columns (dose/time lookup)
    if "_subset_key" not in obs.columns:
        obs["_subset_key"] = raw_obs["_subset_key"].values

    logger.info("  Merging supplementary compound table onto obs")
    compound_csv = paths[f"compounds_{source}_csv"]
    compound_df = load_compound_table(Path(compound_csv))
    obs = merge_compound_metadata(obs, compound_df)

    if annotate_pubchem:
        logger.info("  Annotating compounds with PubChem CIDs")
        unique_df = annotate_pubchem_cids(
            extract_compound_df(obs),
            cache_path=str(paths["pubchem_cid_cache"]),
        )
        obs = obs.drop(columns=["pubchem_cid"], errors="ignore").merge(
            unique_df[["catalog_id", "pubchem_cid"]].drop_duplicates("catalog_id"),
            on="catalog_id",
            how="left",
        )

    logger.info("  Adding fixed metadata columns and sample_id")
    obs = add_fixed_metadata_columns(obs, source=source)

    obs["psbulk_cells"]  = -666
    obs["psbulk_counts"] = calculate_psbulk_counts(adata)

    logger.info("  Enforcing strict obs schema")
    extra_schema = define_tcm_obs_schema() if source == "tcm" else None
    obs_final = enforce_obs_schema(obs, extra_schema=extra_schema)

    logger.info(f"  Processed {len(obs_final):,} observations")
    return obs_final


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def _subset_adata_by_source(
    adata: ad.AnnData,
    source: Literal["mce", "tcm"],
) -> ad.AnnData:
    """Subset the merged CIGS AnnData to rows belonging to ``source``."""
    is_mce = adata.obs["_subset_key"].str.startswith("mce_").values
    mask = is_mce if source == "mce" else ~is_mce
    n = int(mask.sum())
    if n == 0:
        raise ValueError(f"No samples found for source={source!r} in merged AnnData")
    logger.info(f"  Source '{source}': keeping {n:,} / {adata.n_obs:,} samples")
    return adata[mask].copy()


def _filter_compound_by_name(
    adata: ad.AnnData,
    source: Literal["mce", "tcm"],
    compound_name: str,
) -> ad.AnnData:
    """Drop samples whose raw catalog literal equals ``compound_name``.
    Intended for non-library entries written directly into
    ``treatment``/``Treat`` (``JQ1``, ``DMSO``, ``Blank``, ``RNA``)."""
    catalog_col = "treatment" if source == "mce" else "Treat"
    drop_mask = adata.obs[catalog_col].eq(compound_name).values
    logger.info(f"  Dropping {int(drop_mask.sum()):,} {compound_name} samples")
    return adata[~drop_mask].copy()


def standardize_cigs_dataset(
    paths: dict,
    source: Literal["mce", "tcm"],
    adata_raw: Optional[ad.AnnData] = None,
    annotate_pubchem: bool = False,
) -> ad.AnnData:
    """Full standardization pipeline for one CIGS source: load raw → R-based
    QC filter → drop JQ1 → process obs (incl. optional PubChem) → process
    gene annotations → emit the standardized AnnData (int64 X)."""
    logger.info(f"Applying CIGS-specific processing (source={source})")

    if adata_raw is None:
        raw_h5ad = paths[f"raw_{source}_h5ad"]
        logger.info(f"Loading merged raw {source.upper()} AnnData from {raw_h5ad}")
        adata = ad.read_h5ad(raw_h5ad)
    else:
        adata = adata_raw

    logger.info(f"  Loaded AnnData: {adata.n_obs:,} × {adata.n_vars:,}")

    logger.info("Filtering empty and low-QC samples ...")
    adata = filter_samples(adata)

    _assert_dose_time_consistent(adata.obs)

    # JQ1 is an outlier BET bromodomain inhibitor in CIGS whose downstream
    # effects dominate signatures; drop it to stay consistent with the
    # paper's downstream analyses.
    logger.info("Filtering JQ1 samples ...")
    adata = _filter_compound_by_name(adata, source=source, compound_name="JQ1")

    obs_standardized = process_obs_dataframe(
        adata,
        paths,
        source=source,
        annotate_pubchem=annotate_pubchem,
    )

    var_standardized = process_gene_annotations(adata.var, paths, annotate_genes=True)

    logger.info("Creating standardized AnnData object")
    adata_out = ad.AnnData(
        X=adata.X.astype("int64"),
        obs=obs_standardized.set_index("sample_id"),
        var=var_standardized,
    )

    logger.info(f"Standardized AnnData: {adata_out.n_obs:,} × {adata_out.n_vars:,}")
    logger.info("CIGS-specific processing completed")
    return adata_out
