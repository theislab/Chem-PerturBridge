"""CIGS dataset download and Excel → AnnData conversion."""
from __future__ import annotations

import re
import subprocess
import urllib.parse
from pathlib import Path
from typing import Literal, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.utils.parsing_utils import logger


_BASE_URL = "https://cigs.iomicscloud.com/docs/cigs"
_COMPOUNDS_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41592-025-02781-5/MediaObjects/"
    "41592_2025_2781_MOESM3_ESM.xlsx"
)
_COMPOUNDS_FNAME = "41592_2025_2781_MOESM3_ESM.xlsx"

# Maps subset_key → (counts_filename, metadata_filename)
CIGS_SUBSETS = {
    "mce_hek293t_10uM": (
        "MCE_Bioactive_Compounds_HEK293T_10μM_Counts.xlsx",
        "MCE_Bioactive_Compounds_HEK293T_10μM_MetaData.xlsx",
    ),
    "mce_mda_mb_231_10uM": (
        "MCE_Bioactive_Compounds_MDA_MB_231_10μM_Counts.xlsx",
        "MCE_Bioactive_Compounds_MDA_MB_231_10μM_MetaData.xlsx",
    ),
    "tcm_hek293t_10": (
        "TCM_Compounds_HEK293T_10_Counts.xlsx",
        "TCM_Compounds_HEK293T_10_MetaData.xlsx",
    ),
    "tcm_hek293t_20": (
        "TCM_Compounds_HEK293T_20_Counts.xlsx",
        "TCM_Compounds_HEK293T_20_MetaData.xlsx",
    ),
    "tcm_mda_mb_231_10": (
        "TCM_Compounds_MDA_MB_231_10_Counts.xlsx",
        "TCM_Compounds_MDA_MB_231_10_MetaData.xlsx",
    ),
    "tcm_mda_mb_231_20": (
        "TCM_Compounds_MDA_MB_231_20_Counts.xlsx",
        "TCM_Compounds_MDA_MB_231_20_MetaData.xlsx",
    ),
}


# ---------------------------------------------------------------------------
# Download manifest
# ---------------------------------------------------------------------------

def get_cigs_download_manifest(
    data_root: Union[str, Path],
    source: Optional[Literal["mce", "tcm"]] = None,
) -> pd.DataFrame:
    """
    Return a manifest DataFrame of CIGS files to download.

    Parameters
    ----------
    data_root :
        Directory where downloaded files will be stored.
    source :
        If set, include only subsets whose key starts with ``{source}_``
        (e.g. ``"mce"`` → MCE subsets only). The shared compounds Excel is
        always included.

    Returns
    -------
    pd.DataFrame with columns: file, kind, url, path, notes, curl_example.
    """
    data_root = Path(data_root)

    rows = []
    for subset_key, (counts_fname, meta_fname) in CIGS_SUBSETS.items():
        if source is not None and not subset_key.startswith(f"{source}_"):
            continue
        for fname, kind in [(counts_fname, "counts"), (meta_fname, "metadata")]:
            if fname is None:
                continue
            url = f"{_BASE_URL}/{urllib.parse.quote(fname)}"
            rows.append({
                "file":  fname,
                "kind":  kind,
                "url":   url,
                "path":  data_root / fname,
                "notes": subset_key,
            })

    # Compound annotation table (Supplementary Table 3 from the CIGS paper)
    rows.append({
        "file":  _COMPOUNDS_FNAME,
        "kind":  "compounds",
        "url":   _COMPOUNDS_URL,
        "path":  data_root / _COMPOUNDS_FNAME,
        "notes": "compound annotation (Supplementary Table 3)",
    })

    manifest_df = pd.DataFrame(rows)
    manifest_df["curl_example"] = manifest_df.apply(
        lambda row: f"curl -L '{row['url']}' -o '{row['path']}'",
        axis=1,
    )
    return manifest_df


def extract_compound_tables(
    data_root: Union[str, Path],
    skip_existing: bool = True,
) -> None:
    """
    Extract MCE and TCM compound annotation tables from the supplementary Excel
    file and save them as CSV files.

    Reads Table 3 (MCE compounds) and Table 4 (TCM compounds) from
    ``41592_2025_2781_MOESM3_ESM.xlsx`` and writes:
      - ``compounds_mce.csv``
      - ``compounds_tcm.csv``

    Parameters
    ----------
    data_root :
        Directory containing the downloaded supplementary Excel file.
    skip_existing :
        Skip extraction if both CSV files already exist.
    """
    data_root = Path(data_root)
    xlsx_path  = data_root / _COMPOUNDS_FNAME
    mce_csv    = data_root / "compounds_mce.csv"
    tcm_csv    = data_root / "compounds_tcm.csv"

    if skip_existing and mce_csv.exists() and tcm_csv.exists():
        logger.info("Compound CSVs already exist, skipping extraction")
        return

    if not xlsx_path.exists():
        raise FileNotFoundError(
            f"Compound annotation file not found: {xlsx_path}. "
            "Run download_cigs_files() first."
        )

    logger.info(f"Extracting compound tables from {_COMPOUNDS_FNAME} ...")

    # Row 0 of each sheet is the "Supplementary Table N. ..." title; the real
    # header lives on row 1. We rename columns explicitly because the raw TCM
    # header reuses "Catalog Number" for what is actually the CAS number.
    df_mce = pd.read_excel(xlsx_path, sheet_name="Supplementary Table 3", header=1, engine="calamine")
    df_tcm = pd.read_excel(xlsx_path, sheet_name="Supplementary Table 4", header=1, engine="calamine")

    mce_cols = ["compound_name", "catalog_number", "cas_number",
                "moa", "clinical_info", "approved_type"]
    tcm_cols = ["compound_name", "catalog_number", "cas_number"]

    if df_mce.shape[1] != len(mce_cols):
        logger.warning(f"  MCE table has {df_mce.shape[1]} columns, expected {len(mce_cols)}")
    if df_tcm.shape[1] != len(tcm_cols):
        logger.warning(f"  TCM table has {df_tcm.shape[1]} columns, expected {len(tcm_cols)}")

    df_mce.columns = mce_cols[: df_mce.shape[1]]
    df_tcm.columns = tcm_cols[: df_tcm.shape[1]]

    df_mce = _clean_compound_table(df_mce)
    df_tcm = _clean_compound_table(df_tcm)

    df_mce.to_csv(mce_csv, index=False)
    logger.info(f"Saved MCE compound table → {mce_csv.name} ({len(df_mce):,} rows)")

    df_tcm.to_csv(tcm_csv, index=False)
    logger.info(f"Saved TCM compound table → {tcm_csv.name} ({len(df_tcm):,} rows)")


_WHITESPACE_RE = re.compile(r"\s+")


def _clean_compound_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Strip the text-wrap artefacts Excel leaves in the paper's compound tables.

    The TCM sheet of ``41592_2025_2781_MOESM3_ESM.xlsx`` in particular ships
    with narrow columns, so long compound names get word-wrapped and the wrap
    is stored as a literal ``\\n`` inside the cell (e.g.
    ``"Clematichinenos\\nide AR"``, ``"Ganoderic    acid\\nB"``). We collapse
    all runs of whitespace (incl. ``\\n`` and ``\\r``) to a single space and
    strip leading/trailing whitespace, applied to every string column. Named
    text integrity is not guaranteed when the wrap falls *inside* a word
    (``"Syringaresin ol"`` stays as two tokens), but downstream PubChem
    lookups key off ``cas_number`` so this only affects display strings.
    """
    df = df.copy()
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = (
                df[col]
                .astype("string")
                .str.replace(_WHITESPACE_RE, " ", regex=True)
                .str.strip()
            )
    return df


def download_cigs_files(
    data_root: Union[str, Path],
    skip_existing: bool = True,
    source: Optional[Literal["mce", "tcm"]] = None,
) -> None:
    """
    Download CIGS Excel files listed in the download manifest.

    Parameters
    ----------
    data_root :
        Directory where files will be saved.
    skip_existing :
        Skip files that already exist locally.
    source :
        If set, download only counts/metadata files for the given source
        (``"mce"`` or ``"tcm"``). The shared compounds Excel is always
        downloaded so both ``compounds_mce.csv`` and ``compounds_tcm.csv``
        can be extracted.

    Raises
    ------
    subprocess.CalledProcessError
        If a ``curl`` download command fails.
    """
    data_root = Path(data_root)
    data_root.mkdir(parents=True, exist_ok=True)

    manifest_df = get_cigs_download_manifest(data_root, source=source)
    logger.info("Downloading CIGS files to: %s", data_root)

    for _, row in manifest_df.iterrows():
        file_path = Path(row["path"])
        if skip_existing and file_path.exists():
            logger.info("%s already exists, skipping", row["file"])
            continue
        logger.info("Downloading %s ...", row["file"])
        try:
            subprocess.run(["curl", "-L", row["url"], "-o", str(file_path)], check=True)
        except subprocess.CalledProcessError as e:
            logger.error("  Failed to download %s: %s", row["file"], e)
            raise

    logger.info("All CIGS downloads complete")

    extract_compound_tables(data_root, skip_existing=skip_existing)


# ---------------------------------------------------------------------------
# Excel → AnnData conversion
# ---------------------------------------------------------------------------

_CIGS_EXPECTED_N_GENES = 3407


def _read_counts_excel(path: Path) -> pd.DataFrame:
    """
    Read a CIGS counts Excel file into a (samples × genes) DataFrame.

    Layout verified against all 6 CIGS ``*_Counts.xlsx`` files distributed
    with the paper (Wang-lab302/CIGS):

    * **Row 0** — single license-notice cell
      (``"Note: All these datasets are for non-commercial research use only"``).
      Discarded via ``header=1``.
    * **Row 1** — true header: first cell is ``Sample_id`` (MCE files) or
      ``Sample_unique_id`` (TCM files), followed by 3,407 gene symbols.
    * **Rows 2..N** — sample ID in column A, integer UMI counts in
      columns B..end.

    All 6 files ship as samples × genes, so no orientation auto-detection
    is required; a sanity check is emitted if the gene count departs from
    the expected 3,407.

    Parameters
    ----------
    path :
        Path to the counts ``.xlsx`` file.

    Returns
    -------
    DataFrame with samples as rows (indexed by ``Sample_id`` /
    ``Sample_unique_id``) and 3,407 gene symbols as columns.
    """
    logger.info("Reading counts from %s ...", path.name)

    df = pd.read_excel(
        path,
        engine="calamine",
        header=1,          # skip the licence-notice row, use row 1 as header
        index_col=0,
    )

    logger.info(
        "  Shape: %d samples × %d genes | index='%s' | first genes: %s",
        df.shape[0], df.shape[1], df.index.name, list(df.columns[:3]),
    )

    if df.shape[1] != _CIGS_EXPECTED_N_GENES:
        logger.warning(
            "  Unexpected gene count %d (expected %d) in %s",
            df.shape[1], _CIGS_EXPECTED_N_GENES, path.name,
        )
    if not df.index.is_unique:
        logger.warning(
            "  %d duplicated sample IDs in index of %s",
            int(df.index.duplicated().sum()), path.name,
        )
    if df.index.isna().any():
        logger.warning(
            "  %d NaN sample IDs in index of %s",
            int(df.index.isna().sum()), path.name,
        )

    non_numeric = [c for c, t in df.dtypes.items() if not pd.api.types.is_numeric_dtype(t)]
    if non_numeric:
        logger.warning(
            "  Coercing %d non-numeric columns in %s (examples: %s)",
            len(non_numeric), path.name, non_numeric[:3],
        )
        df[non_numeric] = df[non_numeric].apply(pd.to_numeric, errors="coerce")

    return df


def _read_metadata_excel(path: Path) -> pd.DataFrame:
    """
    Read a CIGS metadata Excel file.

    Some files (MCE subsets) contain a licence/note row before the real
    column header.  The data is loaded without a header assumption; if the
    first cell starts with "Note", that row is dropped and the next row is
    promoted to column names.

    Uses the first column as index if it contains unique string identifiers.

    Parameters
    ----------
    path :
        Path to the metadata .xlsx file.

    Returns
    -------
    DataFrame indexed by sample identifier.
    """
    logger.info("Reading metadata from %s ...", path.name)
    df = pd.read_excel(path, engine="calamine", header=None)

    first_cell = str(df.iloc[0, 0])
    if first_cell.startswith("Note:"):
        logger.info("  Note row detected — dropping and promoting next row to header")
        df.columns = df.iloc[1].values
        df = df.iloc[2:].reset_index(drop=True)
    else:
        df.columns = df.iloc[0].values
        df = df.iloc[1:].reset_index(drop=True)

    logger.info("  Metadata shape: %d rows × %d cols", *df.shape)
    logger.info("  Metadata columns: %s", df.columns.tolist())

    first_col = df.columns[0]
    if df[first_col].nunique() == len(df) and df[first_col].dtype == object:
        df = df.set_index(first_col)
        logger.info("  Set index to column '%s'", first_col)

    return df


def _harmonize_obs_dtypes_for_h5ad(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Coerce each obs column to an AnnData-writable dtype.

    ``calamine`` preserves native Excel dtypes, so a column like
    ``sample_row`` can land as ``object`` with a mix of ``int``, ``NaN`` and
    the occasional string — which breaks anndata's vlen-string writer.
    For each object column we try a numeric cast first (covering well-row /
    well-column / dose-like columns); if that loses information we fall
    back to a pure-string object column with ``NaN`` preserved, which the
    anndata writer handles fine.
    """
    def _coerce(s: pd.Series) -> pd.Series:
        if s.dtype != object:
            return s
        numeric = pd.to_numeric(s, errors="coerce")
        if not (numeric.isna() & s.notna()).any():
            return numeric
        return s.map(lambda x: str(x) if pd.notna(x) else None)

    return obs.apply(_coerce)


def _build_anndata_from_subset(
    counts_df: pd.DataFrame,
    meta_df: Optional[pd.DataFrame],
    subset_key: str,
) -> ad.AnnData:
    """
    Combine a counts DataFrame and optional metadata DataFrame into an AnnData.

    The two DataFrames are aligned on their shared index. The ``_subset_key``
    column is stamped onto obs so that standardization.py can look up
    per-subset metadata (cell line, dose) from its own registry.

    Parameters
    ----------
    counts_df :
        (samples × genes) count matrix.
    meta_df :
        Optional sample metadata; index must align with ``counts_df.index``.
    subset_key :
        Experiment subset identifier (e.g. ``"mce_hek293t_10uM"``).

    Returns
    -------
    AnnData with obs indexed by sample ID and var indexed by gene symbol.
    """
    if meta_df is not None:
        shared = counts_df.index.intersection(meta_df.index)
        n_counts_only = len(counts_df.index) - len(shared)
        n_meta_only   = len(meta_df.index)   - len(shared)
        if n_counts_only or n_meta_only:
            logger.warning(
                "Index alignment mismatch for %s: "
                "%d samples only in counts, %d only in metadata — using intersection (%d)",
                subset_key, n_counts_only, n_meta_only, len(shared),
            )
        counts_df = counts_df.loc[shared]
        obs = meta_df.loc[shared].copy()
    else:
        obs = pd.DataFrame(index=counts_df.index)

    # Prefix obs_names with subset_key so merged TCM has unique sample IDs
    # (the 4 TCM subsets reuse the same plate×well grid).
    if not counts_df.index.equals(obs.index):
        if set(counts_df.index) != set(obs.index):
            raise RuntimeError(
                f"{subset_key}: counts_df and obs have different sample sets"
            )
        logger.warning("%s: reordering counts_df to match obs.index", subset_key)
        counts_df = counts_df.loc[obs.index]
    obs.index = subset_key + "." + obs.index.astype(str)

    obs["_subset_key"] = subset_key
    obs = _harmonize_obs_dtypes_for_h5ad(obs)

    var = pd.DataFrame(index=counts_df.columns)
    var.index.name = "gene_symbol"

    X = sp.csr_matrix(counts_df.values.astype('int64'))

    adata = ad.AnnData(X=X, obs=obs, var=var)
    logger.info(
        "  Built AnnData for %s: %d samples × %d genes",
        subset_key, adata.n_obs, adata.n_vars,
    )
    return adata


def convert_cigs_subset_to_adata(
    subset_key: str,
    raw_dir: Union[str, Path],
    output_h5ad: Union[str, Path],
    force: bool = False,
) -> ad.AnnData:
    """
    Convert a CIGS (counts + metadata) Excel pair to an .h5ad file.

    Parameters
    ----------
    subset_key :
        One of the keys in ``CIGS_SUBSETS``.
    raw_dir :
        Directory where the downloaded Excel files live.
    output_h5ad :
        Destination path for the converted .h5ad file.
    force :
        Re-run conversion even if *output_h5ad* already exists.

    Returns
    -------
    AnnData loaded from *output_h5ad*.
    """
    output_h5ad = Path(output_h5ad)
    if output_h5ad.exists() and not force:
        logger.info("h5ad already exists, skipping: %s", output_h5ad)
        return ad.read_h5ad(output_h5ad)

    if subset_key not in CIGS_SUBSETS:
        raise ValueError(
            f"Unknown subset_key {subset_key!r}. Valid keys: {list(CIGS_SUBSETS)}"
        )

    counts_fname, meta_fname = CIGS_SUBSETS[subset_key]
    raw_dir = Path(raw_dir)

    counts_df = _read_counts_excel(raw_dir / counts_fname)

    meta_df: Optional[pd.DataFrame] = None
    if meta_fname is not None:
        meta_path = raw_dir / meta_fname
        if meta_path.exists():
            meta_df = _read_metadata_excel(meta_path)
        else:
            logger.warning("Metadata file not found: %s — proceeding without it", meta_path)

    adata = _build_anndata_from_subset(counts_df, meta_df, subset_key=subset_key)

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad, compression="gzip")
    logger.info("Saved subset h5ad: %s", output_h5ad)
    return adata


def load_and_merge_cigs_subsets(
    paths: dict,
    source: Optional[Literal["mce", "tcm"]] = None,
    subset_keys: Optional[list] = None,
    force_convert: bool = False,
) -> ad.AnnData:
    """
    Convert CIGS subsets to h5ad (if needed) and concatenate them.

    Each subset is first converted to an individual .h5ad under
    ``paths["raw_dir"] / "subsets" / "{subset_key}.h5ad"``, then the
    selected subsets are concatenated via an inner join on genes. MCE and
    TCM use different gene panels (hyphens vs dots in HLA names, Excel
    MARCH/SEPT date-autocorrupt artefacts on the MCE side), so they must
    never be merged together — pass ``source="mce"`` or ``"tcm"``.

    Parameters
    ----------
    paths :
        Path dict as returned by ``get_cigs_paths()``.
    source :
        Restrict to subsets whose key starts with ``{source}_`` (e.g.
        ``source="mce"`` keeps ``mce_hek293t_10uM`` and ``mce_mda_mb_231_10uM``).
        Mutually exclusive with ``subset_keys``.
    subset_keys :
        Explicit subset keys to include. Defaults to all keys in
        ``CIGS_SUBSETS`` filtered by ``source``.
    force_convert :
        Re-run Excel → h5ad conversion even if per-subset .h5ad already exists.

    Returns
    -------
    Concatenated AnnData covering all requested subsets.
    """
    if subset_keys is None:
        subset_keys = [k for k in CIGS_SUBSETS if source is None or k.startswith(f"{source}_")]
    if not subset_keys:
        raise ValueError(
            f"No CIGS subsets selected (source={source!r}, "
            f"available keys: {list(CIGS_SUBSETS)})"
        )

    adatas = []
    for key in subset_keys:
        h5ad_path = paths["raw_dir"] / "subsets" / f"{key}.h5ad"
        adata = convert_cigs_subset_to_adata(
            subset_key=key,
            raw_dir=paths["raw_dir"],
            output_h5ad=h5ad_path,
            force=force_convert,
        )
        adatas.append(adata)

    if len(adatas) == 1:
        merged = adatas[0].copy()
    else:
        logger.info(f"Concatenating {len(adatas)} subsets (inner join on genes) ...")
        gene_sets = [set(a.var_names) for a in adatas]
        shared = set.intersection(*gene_sets)
        for key, a in zip(subset_keys, adatas):
            dropped = set(a.var_names) - shared
            if dropped:
                logger.info(f"  {key}: {len(dropped):,} genes dropped (not shared across all subsets)")
        logger.info(f"  Shared genes across all subsets: {len(shared):,}")

        merged = ad.concat(adatas, join="inner", merge="same")

    # Re-harmonize obs dtypes after concat: ad.concat can promote columns to
    # object dtype (e.g. when subsets disagree on int vs float or introduce
    # NaN on control wells), and h5py cannot serialize mixed-object columns
    # (observed failure: `Cell` column on TCM, `sample_row` on MCE). This
    # also defends against stale per-subset h5ads on disk that predate the
    # harmonizer. Idempotent when already clean.
    merged.obs = _harmonize_obs_dtypes_for_h5ad(merged.obs)

    logger.info(f"Merged CIGS AnnData: {merged.n_obs:,} samples × {merged.n_vars:,} genes")
    return merged
