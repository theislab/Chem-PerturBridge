"""
Dataset-specific processing functions for L1000 Level 3 assembly.

This module contains functions to assemble L1000 Level 3 data from GCTX files
into a standardized AnnData format matching the pseudobulk schema.
"""
from __future__ import annotations

from typing import Optional, Iterable
from pathlib import Path
import gzip
import shutil
import subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from src.utils.parsing_utils import *
from src.pseudobulking.common.pubchem import lookup_pubchem_cids
from src.pseudobulking.datasets.l1000.pubchem_imputation import pubchem_mapping_l1000

# Try to import cmapPy for GCTX parsing
try:
    from cmapPy.pandasGEXpress import parse
    import h5py
    HAS_CMAPPY = True
except ImportError:
    HAS_CMAPPY = False
    logger.warning("cmapPy not available. GCTX parsing will not work.")

# Global cache for GCTX metadata keyed by file path
_GCTX_COL_CACHE = {}

def _read_table(path: Path, **kwargs) -> pd.DataFrame:
    """
    Read a table file, handling gzipped files automatically.
    
    Attempts to read the file at the specified path. If not found, looks for
    a gzipped version (.gz) and reads it instead.
    
    Parameters
    ----------
    path : Path
        Path to the table file (CSV format)
    **kwargs
        Additional arguments passed to pd.read_csv()
        
    Returns
    -------
    pd.DataFrame
        Loaded table data
        
    Raises
    ------
    FileNotFoundError
        If neither the file nor its gzipped version exists
    """
    if not path.exists():
        gz = path.with_suffix(path.suffix + ".gz")
        if gz.exists():
            with gzip.open(gz, "rt") as fh:
                df = pd.read_csv(fh, **kwargs)
        else:
            raise FileNotFoundError(path)
    else:
        df = pd.read_csv(path, **kwargs)
    df.columns = (
        df.columns
        .str.strip()
        .str.replace(" ", "_", regex=False)
        .str.replace("-", "_", regex=False)
        .str.lower()
    )
    return df


def get_download_manifest(data_root: Path, dataset: str = "l1000_phase1") -> pd.DataFrame:
    """
    Get manifest of L1000 Level 3 files to download.
    
    Returns a DataFrame with download information for all required L1000 Level 3
    files from GEO and CLUE resources.
    
    Parameters
    ----------
    data_root : Path
        Root directory where files will be downloaded
        
    Returns
    -------
    pd.DataFrame
        Download manifest with columns: file, kind, size, url, path, notes, curl_example
    """
    if dataset == "l1000_phase1":
        download_manifest = [
            {
                "file": "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz",
                "kind": "expression",
                "size": "48.8 GB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz",
                "notes": "Raw epsilon (landmark genes)"
            },
            {
                "file": "GSE92742_Broad_LINCS_inst_info.txt.gz",
                "kind": "metadata",
                "size": "~150 MB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_inst_info.txt.gz",
                "notes": "Instance-level annotations"
            },
            {
                "file": "GSE92742_Broad_LINCS_cell_info.txt.gz",
                "kind": "metadata",
                "size": "<10 KB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz",
                "notes": "Cell line annotations"
            },
            {
                "file": "GSE92742_Broad_LINCS_pert_info.txt.gz",
                "kind": "metadata",
                "size": "~5 MB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz",
                "notes": "Perturbagen annotations"
            },
            {
                "file": "GSE92742_Broad_LINCS_gene_info.txt.gz",
                "kind": "metadata",
                "size": "<1 MB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz",
                "notes": "Gene annotations"
            },
            {
                "file": "geneinfo_beta.txt",
                "kind": "metadata",
                "size": "<1 MB",
                "url": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt",
                "notes": "Beta gene info with Ensembl IDs"
            },
            {
                "file": "cellinfo_beta.txt",
                "kind": "metadata",
                "size": "<100 KB",
                "url": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt",
                "notes": "Beta cell info with Cellosaurus IDs"
            }
        ]
    elif dataset == "l1000_phase2":
        download_manifest = [
                {
                "file": "GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz",
                "kind": "expression",
                "size": "12.6 GB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz",
                "notes": "Raw epsilon (landmark genes)"
            },
            {
                "file": "GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz",
                "kind": "metadata",
                "size": "~150 MB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz",
                "notes": "Instance-level annotations"
            },
            {
                "file": "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz",
                "kind": "metadata",
                "size": "<10 KB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz",
                "notes": "Cell line annotations"
            },
            {
                "file": "GSE70138_Broad_LINCS_pert_info.txt.gz",
                "kind": "metadata",
                "size": "~5 MB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_pert_info.txt.gz",
                "notes": "Perturbagen metadata"
            },
            {
                "file": "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
                "kind": "metadata",
                "size": "~210 KB",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
                "notes": "Landmark gene annotations"
            },
            {
                "file": "geneinfo_beta.txt",
                "kind": "metadata",
                "size": "1.09 MB",
                "url": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt",
                "notes": "2020 CLUE gene dictionary (for Ensembl IDs)"
            },
            {
                "file": "cellinfo_beta.txt",
                "kind": "metadata",
                "size": "<100 KB",
                "url": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt",
                "notes": "2020 CLUE Cell line annotations (for Cellosaurus IDs)"
            }
        ]
    else:
        raise ValueError(f"Invalid dataset: {dataset}")
    
    manifest_df = pd.DataFrame(download_manifest)
    manifest_df["path"] = manifest_df["file"].apply(lambda f: data_root / f)
    manifest_df["curl_example"] = manifest_df.apply(
        lambda row: f"curl -L '{row['url']}' -o {row['path']}",
        axis=1
    )
    
    return manifest_df


def download_l1000_files(data_root: Optional[str] = None, dataset: str = "l1000_phase1", skip_existing: bool = True) -> None:
    """
    Download L1000 Level 3 data files from GEO and CLUE.
    
    Downloads all required L1000 Level 3 files including:
    - GCTX expression file (48.8 GB)
    - Instance metadata
    - Cell line metadata
    - Perturbagen metadata
    - Gene annotations
    
    Parameters
    ----------
    data_root : str, optional
        Root directory for data files. If None, uses default from define_paths()
    skip_existing : bool, default=True
        If True, skips downloading files that already exist
        
    Raises
    ------
    subprocess.CalledProcessError
        If download command fails
    """
    if data_root is None:
        paths = define_paths(dataset=dataset)
        data_root = Path(paths["level3_gctx"]).parent
    else:
        data_root = Path(data_root)
    
    data_root.mkdir(parents=True, exist_ok=True)
    
    manifest_df = get_download_manifest(data_root, dataset=dataset)
    
    logger.info(f"Downloading L1000 Level 3 files to: {data_root}")
    
    for _, row in manifest_df.iterrows():
        file_path = row["path"]
        
        if skip_existing and file_path.exists():
            logger.info(f"  {row['file']} already exists, skipping")
            continue
        
        logger.info(f"  Downloading {row['file']} ({row['size']})...")
        cmd = f"curl -L '{row['url']}' -o {file_path}"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            logger.info(f"    Downloaded {row['file']}")
        except subprocess.CalledProcessError as e:
            logger.error(f"    Failed to download {row['file']}: {e}")
            raise
    
    logger.info("All downloads complete")


def decompress_l1000_files(data_root: Optional[str] = None, dataset: str = "l1000_phase1") -> None:
    """
    Decompress gzipped L1000 data files.
    
    Decompresses all .gz files in the data directory, including the large GCTX
    expression file. The GCTX decompression may take several minutes.
    
    Parameters
    ----------
    data_root : str, optional
        Root directory containing compressed files. If None, uses default from define_paths()
    dataset : str, default="l1000_phase1"
        Dataset name: "l1000_phase1" or "l1000_phase2"
        
    Raises
    ------
    FileNotFoundError
        If compressed files are not found
    """
    if data_root is None:
        paths = define_paths(dataset=dataset)
        data_root = Path(paths["level3_gctx"]).parent
    else:
        data_root = Path(data_root)
    
    paths = define_paths(str(data_root), dataset=dataset)
    to_decompress = []
    already_done = []
    
    for key, path in paths.items():
        if key.endswith("_gz"):
            continue
        
        path = Path(path)
        if not path.suffix == ".gz":
            gz_candidate = path.with_suffix(path.suffix + ".gz")
            if gz_candidate.exists() and not path.exists():
                to_decompress.append((key, gz_candidate, path))
            elif path.exists():
                already_done.append((key, path))
    
    if already_done:
        logger.info(f"{len(already_done)} file(s) already decompressed")
        for key, path in already_done:
            logger.info(f"  - {key}: {path.name}")
    
    if to_decompress:
        logger.info(f"Decompressing {len(to_decompress)} file(s)...")
        for key, gz_path, target_path in to_decompress:
            logger.info(f"  - {key}: {gz_path.name} -> {target_path.name}")
            if "gctx" in key.lower():
                logger.info("    (GCTX is large, this may take a few minutes...)")
            
            with gzip.open(gz_path, "rb") as src, open(target_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
            logger.info(f"    Done")
        
        logger.info(f"All {len(to_decompress)} file(s) successfully decompressed")
    else:
        if not already_done:
            logger.warning("No compressed files found. Download them first using download_l1000_files()")
        else:
            logger.info("All files are already decompressed")


def check_l1000_files(data_root: Optional[str] = None, dataset: str = "l1000_phase1") -> dict:
    """
    Check status of L1000 data files.
    
    Checks which files are missing, compressed, or ready to use.
    
    Parameters
    ----------
    data_root : str, optional
        Root directory to check. If None, uses default from define_paths()
    dataset : str, default="l1000_phase1"
        Dataset name: "l1000_phase1" or "l1000_phase2"
        
    Returns
    -------
    dict
        Dictionary with keys 'missing', 'compressed', 'ready' containing lists of
        (key, path) tuples for each category
    """
    if data_root is None:
        paths = define_paths(dataset=dataset)
        data_root = Path(paths["level3_gctx"]).parent
    else:
        data_root = Path(data_root)
    
    paths = define_paths(str(data_root), dataset=dataset)
    
    missing = []
    compressed = []
    ready = []
    
    for key, path in paths.items():
        if key.endswith("_gz"):
            continue
        
        # Skip pubchem cache JSON files
        path_str = str(path)
        if "pubchem_cache" in key.lower() and path_str.endswith(".json"):
            continue
        
        path = Path(path)
        if path.suffix == ".gz":
            if not path.exists():
                missing.append((key, path))
        else:
            gz_candidate = path.with_suffix(path.suffix + ".gz")
            if path.exists():
                ready.append((key, path))
            elif gz_candidate.exists():
                compressed.append((key, gz_candidate))
            else:
                missing.append((key, path))
    
    # Log status
    if missing:
        logger.warning(f"{len(missing)} file(s) missing:")
        for key, path in missing:
            logger.warning(f"  - {key}: {path.name}")
        logger.info("  Run download_l1000_files() to download them")
    else:
        logger.info("All required files are present")
    
    if compressed:
        logger.warning(f"{len(compressed)} file(s) still compressed:")
        for key, path in compressed:
            logger.warning(f"  - {key}: {path.name}")
        logger.info("  Run decompress_l1000_files() to extract them")
    else:
        if not missing:
            logger.info("All files are uncompressed and ready to use")
    
    return {
        "missing": missing,
        "compressed": compressed,
        "ready": ready
    }


def standardize_inst(inst: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize instance metadata.
    
    This function standardizes the instance ID column to 'lincs_inst_id' and adds plate
    information if missing by parsing the instance ID.
    
    Parameters
    ----------
    inst : pd.DataFrame
        Instance metadata with inst_id, sample_id, or distil_id column
        
    Returns
    -------
    pd.DataFrame
        Standardized dataframe with 'lincs_inst_id' as index and column, and 'det_plate' column added if missing
    """
    df = inst.copy()
    
    if "inst_id" not in df.columns:
        raise KeyError("inst info missing inst_id")
    
    # Add plate info if missing
    if 'det_plate' not in df.columns:
        df['det_plate'] = df['inst_id'].str.split(':').str[0]

    df['lincs_inst_id'] = df['inst_id'].astype("string")
    
    return df.set_index("lincs_inst_id", drop=False)





def process_cellinfo(cellinfo: pd.DataFrame,
                     cellinfo_extra: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Process cell info with optional extra metadata.
    
    Merges base cell info with extra metadata with Cellosaurus IDs.
    
    Parameters
    ----------
    cellinfo : pd.DataFrame
        Base cell line information with cell_id column
    cellinfo_extra : pd.DataFrame, optional
        Extra cell line metadata with Cellosaurus IDs
        
    Returns
    -------
    pd.DataFrame
        Processed cell info indexed by cell_id with prefixed column names
    """
    NAME_TO_CVCL = {
        "HA1E": "CVCL_VU89",
        "HEK293T": "CVCL_0063",
        "HS27A": "CVCL_3719",
        "FIBRNPC": "CVCL_UK07",
        "U266": "CVCL_0566",
        "HUES3":   "CVCL_B161",
        "HUVEC":   "CVCL_2959"
    }

    df = cellinfo.copy()
    cellinfo_id = 'cell_id'
    if cellinfo_extra is not None and not cellinfo_extra.empty:
        df_extra = cellinfo_extra.copy()
        df_extra[cellinfo_id] = df_extra['cell_iname'].str.replace('_', '.')
        df = pd.concat([df, df_extra])[df.columns].drop_duplicates(cellinfo_id).copy()
        df = df.merge(df_extra[[cellinfo_id, 'cellosaurus_id']], on=cellinfo_id, how='left')
        df[cellinfo_id + '_mixed'] = df['cellosaurus_id'].fillna(df[cellinfo_id])
    rename_map = {col: f"cellinfo_{col}" for col in df.columns if col != "cell_id"}
    df = df.rename(columns=rename_map).set_index("cell_id").copy()
    if 'cellinfo_cell_id_mixed' in df.columns:
        df['cellinfo_cell_id_mixed'] = df['cellinfo_cell_id_mixed'].replace(NAME_TO_CVCL)
    return df

def process_pert_metadata(pert: pd.DataFrame, inst: pd.DataFrame) -> pd.DataFrame:
    """
    Process perturbation metadata.
    
    Parameters
    ----------
    pert : pd.DataFrame
        Perturbation metadata
    inst : pd.DataFrame
        Instance metadata
        
    Returns
    -------
    pd.DataFrame
        Processed perturbation metadata
    """
    pert = pert.copy()
    inst = inst.copy()

    if 'pubchem_cid' not in pert.columns:
        pert['pubchem_cid'] = None

    return pd.concat([pert, inst]).drop_duplicates('pert_id', keep='first')[pert.columns]


def standardize_dose(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize dose units to micromolar.
    
    Converts dose values from various units (um, nm, mm) to micromolar (uM).
    Handles unit variations including µm and full names.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with 'pert_dose' and 'pert_dose_unit' columns
        
    Returns
    -------
    pd.DataFrame
        Input dataframe with additional 'pert_dose_um' column in micromolar
    """
    out = df.copy()
    dose = pd.to_numeric(out.get("pert_dose"), errors="coerce")
    unit = out.get("pert_dose_unit").astype("string").str.strip().str.lower()
    dose_um = pd.Series(np.nan, index=out.index, dtype=float)
    mask_um = unit.isin({"um", "µm", "micromolar"})
    mask_nm = unit.isin({"nm", "nanomolar"})
    mask_mm = unit.isin({"mm", "millimolar"})
    dose_um[mask_um] = dose[mask_um]
    dose_um[mask_nm] = dose[mask_nm] / 1000.0
    dose_um[mask_mm] = dose[mask_mm] * 1000.0
    out["pert_dose_um"] = dose_um
    return out


def standardize_time(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize time units to hours.
    
    Converts time values from various units (hours, days, minutes) to hours.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with 'pert_time' and 'pert_time_unit' columns
        
    Returns
    -------
    pd.DataFrame
        Input dataframe with additional 'pert_time_h' column in hours
    """
    out = df.copy()
    time_val = pd.to_numeric(out.get("pert_time"), errors="coerce")
    time_unit = out.get("pert_time_unit").astype("string").str.strip().str.lower()
    hours = pd.Series(np.nan, index=out.index, dtype=float)
    hours[time_unit.isin({"h", "hr", "hrs", "hour", "hours"})] = time_val[time_unit.isin({"h", "hr", "hrs", "hour", "hours"})]
    hours[time_unit.isin({"d", "day", "days"})] = time_val[time_unit.isin({"d", "day", "days"})] * 24.0
    hours[time_unit.isin({"m", "min", "mins", "minute", "minutes"})] = time_val[time_unit.isin({"m", "min", "mins", "minute", "minutes"})] / 60.0
    out["pert_time_h"] = hours
    return out


def build_development_stage(row: pd.Series, 
                            col: str = "cellinfo_donor_age") -> str:
    """
    Build development stage label from donor age.
    
    Creates a standardized development stage label in the format
    "{age}-year-old stage" from numeric age values.
    
    Parameters
    ----------
    row : pd.Series
        Row of instance metadata
    col : str, default="cellinfo_donor_age"
        Column name containing age information
        
    Returns
    -------
    str
        Development stage label (e.g., "69-year-old stage") or "unknown"
    """
    if col in row.index and pd.notna(row[col]):
        try:
            age = float(row[col])
            if np.isfinite(age) and age > 0:
                return f"{int(age)}-year-old stage"
        except ValueError:
            pass
    return "unknown"


def process_gene_annotations(geneinfo: pd.DataFrame,
                            geneinfo_beta: Optional[pd.DataFrame] = None,
                            full_gene_matrix: bool = False) -> pd.DataFrame:
    """
    Process gene annotations from L1000 gene info tables.
    
    Parameters
    ----------
    geneinfo : pd.DataFrame
        Level 3 gene info table with L1000 gene annotations
    geneinfo_beta : pd.DataFrame, optional
        Beta gene info table with Ensembl ID mappings
    full_gene_matrix : bool, default=False
        If False, restrict to landmark genes only
        
    Returns
    -------
    pd.DataFrame
        Processed var dataframe with gene annotations
    """
    var = geneinfo.copy()
    
    # Rename columns to standard schema
    var = var.rename(columns={
        "pr_gene_id": "gene_id",
        "pr_gene_symbol": "symbol",
        "pr_gene_title": "gene_title",
        "pr_is_lm": "is_landmark"
    })
    
    var = var.drop_duplicates(subset=["gene_id"]).set_index("gene_id")
    var["symbol"] = var["symbol"].astype("string")
    
    # Map to Ensembl IDs if available
    if geneinfo_beta is not None and {"gene_id", "ensembl_id"}.issubset(geneinfo_beta.columns):
        sym_to_ens = (
            geneinfo_beta[["gene_id", "ensembl_id"]]
            .dropna()
            .drop_duplicates(subset=["gene_id"])
            .set_index("gene_id")
            ["ensembl_id"]
        )
        var["ensembl_id"] = var.index.map(sym_to_ens)
    else:
        var["ensembl_id"] = var.index.astype("string")
    
    # Fill missing Ensembl IDs with gene_id
    var["ensembl_id"] = var["ensembl_id"].fillna(var.index.to_series().astype("string"))
    var = var[["symbol", "ensembl_id", "is_landmark"]]
    
    # Filter to landmark genes if requested
    if not full_gene_matrix:
        landmark_mask = var["is_landmark"].fillna(0).astype(int)
        var = var.loc[landmark_mask == 1]
        logger.info(f'  Restricting to {len(var):,} landmark genes')
    else:
        logger.info(f'  Using all {len(var):,} gene features')
    
    return var


def build_obs_dataframe(inst: pd.DataFrame, dataset: str = "l1000_phase1") -> pd.DataFrame:
    """
    Build standardized .obs dataframe from instance metadata.
    Sticked to https://lamin.ai/laminlabs/pertdata/transform/REAvqqdo3sbH0000
    
    Parameters
    ----------
    inst : pd.DataFrame
        Instance metadata dataframe with LINCS information
        
    Returns
    -------
    pd.DataFrame
        Standardized obs dataframe with pseudobulk schema
    """
    obs = pd.DataFrame(index=inst.index)
    
    # Copy required identifier columns
    obs["inst_id"] = inst["inst_id"].astype("string")
    
    # Map perturbation types to standard schema
    PERT_TYPE_MAP = {
        "trt_cp": "compound",
        "trt_lig": "biologic",
        "trt_sh": "genetic",
        "trt_oe": "genetic",
        "trt_oe.mut": "genetic",
        "trt_xpr": "genetic",
        "ctl_vehicle": "compound",
        "trt_poscon": "compound",
        "ctl_vector": "genetic",
        "ctl_untrt": "biologic"
    }

    
    # Build standard obs columns
    obs["plate"] = inst.get("det_plate", None)
    obs["well"] = inst.get("rna_well", inst.get("det_well", None))
    obs["cell_type"] = inst.get("cellinfo_cell_id_mixed", inst.get("cell_id", None)).fillna(inst.get("cell_id", None))
    obs["perturbagen"] = inst.get("pert_iname", None)
    if 'pert_type' in inst.columns:
        obs["pert_type"] = inst["pert_type"].map(PERT_TYPE_MAP)
        obs["is_control"] = inst["pert_type"].str.startswith("ctl")
    else:
        raise ValueError("pert_type column not found in instance metadata")
    
    obs["pert_dose_uM"] = inst["pert_dose_um"].astype(float)
    obs.loc[obs['is_control'], 'pert_dose_uM'] = 0
    obs["pert_time_h"] = inst["pert_time_h"].astype(float)
    obs["suspension_type"] = "cell"
    obs["tissue"] = inst.get("cellinfo_primary_site", "unknown")
    obs["tissue_type"] = "cell culture"
    obs["disease"] = inst.get("cellinfo_subtype", "unknown")
    obs["library"] = None
    obs["stimulation"] = None
    obs["guide"] = None

    if dataset == "l1000_phase1":
        obs["dataset"] = "LINCS_phase1_level3_epsilon"
    elif dataset == "l1000_phase2":
        obs["dataset"] = "LINCS_phase2_level3"
    else:
        raise ValueError(f"Invalid dataset: {dataset}")

    obs["assay"] = "L1000 mRNA profiling assay"
    obs["development_stage"] = inst.apply(build_development_stage, axis=1)
    obs["organism"] = "human"

    if 'cellinfo_donor_sex' in inst.columns:
        obs["sex"] = inst["cellinfo_donor_sex"].map({"M": "male", "F": "female"})
    else:
        obs["sex"] = "unknown"
        
    obs["self_reported_ethnicity"] = inst.get("cellinfo_donor_ethnicity", "unknown")
    if 'pubchem_cid' in inst.columns:
        inst['pubchem_cid'] = pd.to_numeric(inst['pubchem_cid'], errors='coerce').fillna(-666).astype('int64')
    obs["pubchem_cid"] = inst.get("pubchem_cid", None)
    obs["psbulk_cells"] = None
    obs["psbulk_counts"] = None
    
    # Add metadata columns
    obs["lincs_inst_id"] = inst.index.astype("string")
    obs["source_gctx"] = inst["source_gctx"].astype("string")
    
    # Create composite sample_id
    obs["sample_id"] = (
        obs["plate"].astype(str).str.replace(" ", "", regex=False) + "_" +
        obs["well"].astype(str).str.replace(" ", "", regex=False) + "_" +
        obs["perturbagen"].astype(str).str.replace(" ", "_", regex=False) + "_" +
        obs["cell_type"].astype(str).str.replace(" ", "_", regex=False)
    )
    # Clean up missing values and duplicates
    obs = obs.replace({-666: None, '-666': None, 'None': None, 'nan': None, '<NA>': None})
    obs = obs.set_index("inst_id", drop=True)
    obs = materialize_string_columns(obs)
    return obs


def enforce_obs_schema(obs: pd.DataFrame) -> pd.DataFrame:
    """
    Enforce strict obs schema on the observations dataframe.
    
    This function ensures that the obs dataframe follows the strict 24-column
    schema defined in OBS_SCHEMA. It:
    1. Adds missing columns as NaN
    2. Casts columns to the correct dtypes (category, float64, int64)
    3. Fills missing values appropriately (unknown for certain category columns, -666 for int64)
    4. Selects only the columns in the schema
    5. Materializes string columns for AnnData compatibility
    
    Parameters
    ----------
    obs : pd.DataFrame
        Observations dataframe to enforce schema on
        
    Returns
    -------
    pd.DataFrame
        Observations dataframe conforming to the strict schema
    """
    # Columns that should be filled with "unknown" instead of NaN
    cols_fillna_unknown = [
        'tissue',
        'tissue_type',
        'disease',
        'development_stage',
        'sex',
        'self_reported_ethnicity'
    ]
    
    # Get obs schema
    obs_schema = define_obs_schema()
    
    # Create dtype map from schema
    dtype_map = {col: dtype for col, dtype, _ in obs_schema}
    
    # Create schema dataframe for reference
    obs_schema_df = pd.DataFrame(obs_schema, columns=["column", "dtype", "description"]).set_index("column")
    
    obs_for_schema = obs.copy()
    
    # Ensure all schema columns exist, add as NaN if missing
    for col, dtype in dtype_map.items():
        if col not in obs_for_schema.columns:
            obs_for_schema[col] = np.nan
        
        # Cast to appropriate dtype
        if dtype == "category":
            if col in cols_fillna_unknown:
                obs_for_schema[col] = (
                    obs_for_schema[col]
                    .fillna("unknown")
                    .astype("string")
                    .astype(object)
                    .astype("category")
                )
            else:
                obs_for_schema[col] = (
                    obs_for_schema[col]
                    .astype(object)
                    .astype("category")
                )
        elif dtype == "float64":
            obs_for_schema[col] = pd.to_numeric(obs_for_schema[col], errors="coerce")
        elif dtype == "int64":
            obs_for_schema[col] = (
                pd.to_numeric(obs_for_schema[col], errors="coerce")
                .fillna(-666)
                .astype("int64")
            )
    
    # Select only columns in schema
    obs_for_schema = obs_for_schema[obs_schema_df.index.tolist()]
    
    # Materialize string columns for AnnData compatibility
    obs_for_schema = materialize_string_columns(obs_for_schema)
    
    return obs_for_schema


def build_config(config: Optional[dict] = None) -> dict:
    """
    Build L1000 configuration by merging defaults with provided config.
    
    Parameters
    ----------
    config : dict, optional
        User-provided configuration parameters
        
    Returns
    -------
    dict
        Complete configuration with defaults filled in. Keys include:
        - perturbation_types_to_keep: set of perturbation types to include
        - control: set of control perturbagen names
        - full_gene_matrix: bool, whether to use full gene matrix
        - subsampling: bool, whether to use subsampling of the dataset
        - annotate_pubchem: bool, whether to annotate perturbations with PubChem CIDs
        - download_if_missing: bool, whether to auto-download missing files (default: True)
    """
    # Set defaults
    if config is None:
        config = {}
    
    # Configuration defaults
    default_config = {
        "perturbation_types_to_keep": {"trt_cp", "ctl_vehicle"},
        "control": {"DMSO"},
        "full_gene_matrix": False,
        "subsampling": False,
        "annotate_pubchem": False,
        "download_if_missing": True,
    }
    
    # Merge with provided config
    return {**default_config, **config}


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
        ("plate", "category", "Assay detection plate identifier (det_plate)"),
        ("well", "category", "Well ID on the RNA/detection plate (rna_well or det_well)"),
        ("cell_type", "category", "Cell line / cell_id"),
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
        ("psbulk_cells", "int64", "Total #cells contributing (if no info - then -666)"),
        ("psbulk_counts", "int64", "Total #counts contributing (if no info - then -666)"),
    ]


def define_paths(data_root: Optional[str] = None, dataset: str = "l1000_phase1") -> dict:
    """
    Define file paths for L1000 Level 3 data.
    
    Parameters
    ----------
    data_root : str, optional
        Root directory containing L1000 data files.
        Defaults to './lincs_data'
        
    Returns
    -------
    dict
        Dictionary mapping file identifiers to Path objects
    """
    if data_root is None:
        data_root = './lincs_data'
    
    data_root = Path(data_root)
    processed_dir = data_root / "processed"
    processed_dir.mkdir(exist_ok=True)
    
    if dataset == "l1000_phase1":
        return {
            "level3_gctx": data_root / "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",
            "instinfo": data_root / "GSE92742_Broad_LINCS_inst_info.txt",
            "cellinfo": data_root / "GSE92742_Broad_LINCS_cell_info.txt",
            "pert_info": data_root / "GSE92742_Broad_LINCS_pert_info.txt",
            "geneinfo_level3": data_root / "GSE92742_Broad_LINCS_gene_info.txt",
            "geneinfo_beta": data_root / "geneinfo_beta.txt",
            "cellinfo_beta": data_root / "cellinfo_beta.txt",
            "pubchem_cache": processed_dir / "pubchem_cache.json",
        }
    elif dataset == "l1000_phase2":
        return {
                "level3_gctx": data_root / "GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx",
                "level3_gctx_gz": data_root / "GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz",
                "instinfo": data_root / "GSE70138_Broad_LINCS_inst_info_2017-03-06.txt",
                "cellinfo": data_root / "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                "pert_info": data_root / "GSE70138_Broad_LINCS_pert_info.txt",
                "geneinfo_level3": data_root / "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                "geneinfo_beta": data_root / "geneinfo_beta.txt",
                "cellinfo_beta": data_root / "cellinfo_beta.txt",
                "pubchem_cache": processed_dir / "pubchem_cache.json",
        }
    else:
        raise ValueError(f"Invalid dataset: {dataset}")

def add_alternative_identifiers(inst_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize instance metadata and ensure inst_id/lincs_inst_id are present.
    """
    logger.info('  Processing instance metadata')
    inst = standardize_inst(inst_raw)

    if 'inst_id' not in inst_raw.columns:
        raise ValueError("instance metadata must have inst_id column")

    if 'lincs_inst_id' not in inst.columns:
        inst['lincs_inst_id'] = inst['inst_id'].astype("string")


    return inst


def annotate_pubchem_cids(df: pd.DataFrame, paths: dict, config: dict = None) -> pd.DataFrame:
    """
    Annotate metadata containing perturbation information with PubChem CIDs.
    
    This function adds or updates PubChem CID information for perturbations
    using multiple lookup strategies (InChIKey, SMILES, drug name) with
    persistent caching to avoid redundant API calls.
    
    Only compound perturbations (trt_cp and ctl_vehicle) are annotated,
    as other perturbation types (shRNA, CRISPR, etc.) don't have PubChem CIDs.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Metadata dataframe containing perturbation information (e.g., instance metadata
        enriched with perturbation data via merge). Must contain columns: pert_id, pert_type,
        pert_iname, and optionally inchi_key, canonical_smiles.
    paths : dict
        Dictionary of file paths from define_paths(), must include 'pubchem_cache'
        
    Returns:
    --------
    pd.DataFrame
        Metadata dataframe with updated pubchem_cid column
    """
    def standardize_pubchem_cid(cid):
        try:
            if pd.isna(cid):
                return None
            cid_int = int(cid)
            return cid_int if cid_int > 0 else None
        except Exception as e:
            logger.warning(f"Error standardizing pubchem_cid: {e}")
            return None

    # Filter to only compound perturbations (trt_cp and ctl_vehicle)
    # Other types (shRNA, CRISPR, etc.) don't have PubChem CIDs
    df = df.copy()

    if 'pubchem_cid' not in df.columns:
        df['pubchem_cid'] = None

    df['pubchem_cid'] = df['pubchem_cid'].apply(standardize_pubchem_cid)
    df['pubchem_cid'] = pd.to_numeric(df['pubchem_cid'], errors='coerce').fillna(-666).astype('int64')

    compound_types = {"trt_cp", "ctl_vehicle"}
    df_compounds = df[df["pert_type"].isin(compound_types)].copy()

    if len(df_compounds) == 0:
        logger.warning("No compound perturbations found to annotate")
        return df

    if config is not None and config.get("subsampling"):
        df_compounds = df_compounds.sample(n=1000, random_state=0)
    
    
    
    
    pubchem_cache = {}
    cache_path = str(paths["pubchem_cache"])
    df_compounds = lookup_pubchem_cids(
        df_compounds, 
        cache=pubchem_cache, 
        pert_id_col='pert_id',
        drug_col='pert_iname',
        cache_path=cache_path,
        manual_mapping_func=pubchem_mapping_l1000,
        dataset_key='l1000'
    )
    
    df_compounds['pubchem_cid'] = pd.to_numeric(df_compounds['pubchem_cid'], errors='coerce').fillna(-666).astype("int64")
    
    # Update the original df with annotated CIDs
    df.loc[df_compounds.index, "pubchem_cid"] = df_compounds["pubchem_cid"]
    return df


def enrich_instance_metadata(inst: pd.DataFrame, 
                            cellinfo: pd.DataFrame,
                            pert_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Enrich instance metadata by merging with cell and perturbation info.
    
    Parameters
    ----------
    inst : pd.DataFrame
        Processed instance metadata
    cellinfo : pd.DataFrame
        Cell line information
    pert_raw : pd.DataFrame
        Perturbation information
        
    Returns
    -------
    pd.DataFrame
        Enriched instance metadata with merged information
    """
    pert_cols = ["pert_id", "pert_iname", "pert_type", "pubchem_cid"]
    inst = inst.merge(cellinfo, how="left", left_on="cell_id", right_index=True)    
    inst = inst.merge(pert_raw[pert_cols].drop_duplicates("pert_id"), 
                     how="left", on="pert_id", suffixes=("", "_pert"))
    
    return inst


def process_instance_metadata(inst_raw: pd.DataFrame,
                            cellinfo: pd.DataFrame,
                            pert_raw: pd.DataFrame,
                            config: dict,
                            paths: dict,
                            max_samples: int = 1000) -> pd.DataFrame:
    """
    Process instance metadata through complete pipeline.
    
    This function handles the full instance metadata processing pipeline:
    1. Add alternative identifiers (distil_id, inst_id, sample_id, sig_id)
    2. Enrich with cell and perturbation information
    3. Standardize dose units to micromolar
    4. Standardize time units to hours
    5. Apply filters (perturbation types, controls, subsampling)
    6. Add source file information
    
    Parameters
    ----------
    inst_raw : pd.DataFrame
        Raw instance/sample metadata
    cellinfo : pd.DataFrame
        Cell line information
    pert_raw : pd.DataFrame
        Perturbation information
    config : dict
        Configuration with filter settings
    paths : dict
        Dictionary of file paths
        
    Returns
    -------
    pd.DataFrame
        Fully processed and filtered instance metadata
    """
    inst = add_alternative_identifiers(inst_raw)
    inst = enrich_instance_metadata(inst, cellinfo, pert_raw)
    inst = standardize_dose(inst)
    inst = standardize_time(inst)
    
    # Apply filters
    if config["perturbation_types_to_keep"] is not None:
        inst = inst[inst["pert_type"].isin(config["perturbation_types_to_keep"])]

    if config["control"] is not None:
        inst = inst[(inst['pert_type'].str.startswith('ctl') & inst['pert_iname'].isin(config['control'])) |
                    (~inst['pert_type'].str.startswith('ctl'))]

    cell_ids_with_controls = set(inst[inst['pert_type'].str.startswith('ctl')]['cell_id'].unique())
    cell_ids_with_compounds = set(inst[~inst['pert_type'].str.startswith('ctl')]['cell_id'].unique())
    valid_cell_ids = cell_ids_with_controls & cell_ids_with_compounds
    if len(valid_cell_ids) == 0:
        raise ValueError("No cell lines found with both controls and compounds")

    inst = inst[inst['cell_id'].isin(valid_cell_ids)].copy()

    if config["subsampling"]:
        inst = inst.sample(max_samples, random_state=0)
    
    # Add source file information
    inst["source_gctx"] = str(paths["level3_gctx"])
    
    return inst


def load_metadata_tables(paths: dict) -> tuple:
    """
    Load L1000 metadata tables from files.
    
    Parameters
    ----------
    paths : dict
        Dictionary of file paths from define_paths()
        
    Returns
    -------
    tuple
        (inst_raw, cellinfo_raw, pert_raw, geneinfo, geneinfo_beta, cellinfo_beta)
        - inst_raw: Instance/sample information
        - cellinfo_raw: Cell line information
        - pert_raw: Perturbation information
        - geneinfo: Gene annotations (Level 3)
        - geneinfo_beta: Gene annotations (beta) - optional, None if not found
        - cellinfo_beta: Cell line annotations (beta) - optional, None if not found
    """
    logger.info('  Loading metadata tables')
    
    inst_raw = _read_table(paths["instinfo"], sep="\t", low_memory=False)
    cellinfo_raw = _read_table(paths["cellinfo"], sep="\t")
    pert_raw = _read_table(paths["pert_info"], sep="\t")
    geneinfo = _read_table(paths["geneinfo_level3"], sep="\t")
    geneinfo_beta = _read_table(paths["geneinfo_beta"], sep="\t") if paths["geneinfo_beta"].exists() else None
    cellinfo_beta = _read_table(paths["cellinfo_beta"], sep="\t") if paths["cellinfo_beta"].exists() else None
    
    # Log loaded table sizes
    logger.info(f"    Loaded {len(inst_raw):,} instances")
    logger.info(f"    Loaded {len(cellinfo_raw):,} cell lines")
    logger.info(f"    Loaded {len(pert_raw):,} perturbations")
    logger.info(f"    Loaded {len(geneinfo):,} genes")
    
    return inst_raw, cellinfo_raw, pert_raw, geneinfo, geneinfo_beta, cellinfo_beta


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


def _get_gctx_column_ids(gctx_path: Path) -> set:
    """
    Load and cache GCTX column IDs.
    
    Parameters
    ----------
    gctx_path : Path
        Path to the GCTX file.
        
    Returns
    -------
    set
        Set of GCTX column IDs.
    """
    cache_key = str(gctx_path)
    if cache_key not in _GCTX_COL_CACHE:
        if not gctx_path.exists():
            raise FileNotFoundError(gctx_path)
        with h5py.File(gctx_path, "r") as handle:
            root = handle["0"] if "0" in handle else handle
            ids = root["META"]["COL"]["id"][:]
        col_ids = [
            id_.decode("utf-8") if isinstance(id_, bytes) else str(id_)
            for id_ in ids
        ]
        _GCTX_COL_CACHE[cache_key] = {
            "ids": col_ids,
            "id_map": {id_.lower(): id_ for id_ in col_ids},
        }
        logger.info(f"  Loaded {len(col_ids):,} column IDs from {gctx_path.name}")
    return set(_GCTX_COL_CACHE[cache_key]["ids"])


def _get_gctx_id_map(gctx_path: Path) -> dict:
    """
    Return case-insensitive GCTX ID map for a given file path.
    """
    cache_key = str(gctx_path)
    if cache_key not in _GCTX_COL_CACHE:
        _get_gctx_column_ids(gctx_path)
    return _GCTX_COL_CACHE[cache_key]["id_map"]


def _filter_obs_for_gctx(obs_subset: pd.DataFrame, gctx_path: Path) -> pd.DataFrame:
    """
    Filter obs to only include samples present in GCTX.
    
    Parameters
    ----------
    obs_subset : pd.DataFrame
        Observations dataframe with "gctx_id" column.
    gctx_path : Path
        Path to the GCTX file.
        
    Returns
    -------
    pd.DataFrame
        Filtered obs containing only samples present in GCTX.
        
    Raises
    ------
    ValueError
        If no samples remain after filtering.
    """
    available = _get_gctx_column_ids(gctx_path)
    mask = obs_subset["gctx_id"].isin(available)
    if not mask.all():
        dropped = int((~mask).sum())
        logger.warning(f"  Dropping {dropped} samples not present in the Level 3 matrix")
    obs_filtered = obs_subset.loc[mask].copy()
    if obs_filtered.empty:
        raise ValueError("No samples remain after intersecting with the Level 3 matrix.")
    return obs_filtered


def _normalize_gene_ids(gene_ids: Optional[Iterable[str]]) -> Optional[list]:
    """
    Normalize gene IDs to a deduplicated list.
    
    Parameters
    ----------
    gene_ids : Optional[Iterable[str]]
        Gene IDs to normalize.
        
    Returns
    -------
    Optional[list]
        Deduplicated list of gene IDs, or None if input is None.
    """
    if gene_ids is None:
        return None
    return list(dict.fromkeys(str(gid) for gid in gene_ids))


def build_obs_helpers(obs: pd.DataFrame, gctx_path: Path) -> pd.DataFrame:
    """
    Build helper dataframe with GCTX IDs for obs matching.
    
    Tries exact matching first, then falls back to case-insensitive matching.
    
    Parameters
    ----------
    obs : pd.DataFrame
        Observations dataframe with index to match against GCTX IDs.
    gctx_path : Path
        Path to the GCTX file.
        
    Returns
    -------
    pd.DataFrame
        Subset of obs with "source_gctx" and "gctx_id" columns.
        
    Raises
    ------
    ValueError
        If none of the obs IDs match the GCTX metadata.
    """
    gctx_ids = _get_gctx_column_ids(gctx_path)
    id_map = _get_gctx_id_map(gctx_path)
    index_vals = obs.index.astype("string")
    mask = index_vals.isin(gctx_ids)
    
    if mask.any():
        matched = int(mask.sum())
        logger.info(f"  Using obs index to subset GCTX ({matched}/{len(obs)} samples match)")
        helpers = obs.loc[mask, ["source_gctx"]].copy()
        helpers["gctx_id"] = index_vals[mask]
        return helpers
    
    # Fallback: case-insensitive matching
    lower_vals = index_vals.str.lower()
    mask_lower = lower_vals.isin(id_map)
    if mask_lower.any():
        matched = int(mask_lower.sum())
        logger.info(f"  Using obs index (case-insensitive) to subset GCTX ({matched}/{len(obs)} samples match)")
        helpers = obs.loc[mask_lower, ["source_gctx"]].copy()
        helpers["gctx_id"] = lower_vals[mask_lower].map(id_map)
        return helpers
    
    raise ValueError("None of the obs ID columns match the Level 3 GCTX metadata.")


def assemble_anndata(obs: pd.DataFrame, var: pd.DataFrame, expr: pd.DataFrame) -> ad.AnnData:
    """
    Assemble final AnnData object from processed components.
    
    This function:
    1. Prepares var dataframe with Ensembl IDs as index
    2. Handles duplicate Ensembl IDs by adding suffixes
    3. Reindexes expression to match var
    4. Creates final AnnData object with sparse matrix
    
    Parameters
    ----------
    obs : pd.DataFrame
        Observations dataframe with sample metadata
    var : pd.DataFrame
        Gene annotations with gene_id as index
    expr : pd.DataFrame
        Expression matrix (genes x samples)
        
    Returns
    -------
    ad.AnnData
        Final AnnData object with standardized schema
    """
    # Prepare var dataframe
    var_df = var.reindex(expr.columns.astype(int)).copy()
    var_idx = var_df["ensembl_id"].fillna(var_df.index.to_series().astype("string"))
    var_idx = var_idx.astype(object)
    
    # Handle duplicate ensembl IDs
    duplicate_mask = var_idx.duplicated(keep=False)
    if duplicate_mask.any():
        suffix = (
            var_idx[duplicate_mask]
            .groupby(var_idx[duplicate_mask])
            .cumcount()
            .astype("string")
        )
        var_idx = var_idx.astype("string")
        var_idx[duplicate_mask] = var_idx[duplicate_mask] + "_" + suffix
    var_idx = var_idx.astype(object)
    var_df.index = var_idx
    var_df = var_df[["symbol"]]
    var_df["symbol"] = (
        var_df["symbol"]
        .replace({"": None})
        .astype(object)
        .astype("category")
    )
    var_df = materialize_string_columns(var_df)
    
    
    # Prepare final obs
    obs_final = obs.loc[expr.index].set_index('sample_id')
    
    # Create AnnData
    adata = ad.AnnData(
        X=sp.csr_matrix(expr.to_numpy(np.float32)),
        obs=obs_final,
        var=var_df
    )
    
    return adata


def load_expression(obs_helpers: pd.DataFrame, gctx_path: Path, 
                    gene_ids: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """
    Load expression data from GCTX for given observations and genes.
    
    Parameters
    ----------
    obs_helpers : pd.DataFrame
        Helper dataframe with "gctx_id" column mapping to GCTX column IDs.
    gctx_path : Path
        Path to the GCTX file.
    gene_ids : Optional[Iterable[str]]
        Gene IDs to extract. If None, extracts all genes.
        
    Returns
    -------
    pd.DataFrame
        Expression dataframe with obs_helpers.index as columns and genes as rows,
        reindexed to match the requested gene_ids order.
    """
    obs_filtered = _filter_obs_for_gctx(obs_helpers, gctx_path)
    gene_ids_normalized = _normalize_gene_ids(gene_ids)
    lincs_ids = obs_filtered["gctx_id"].tolist()
    
    gctx = parse.parse(
        str(gctx_path),
        cid=lincs_ids,
        rid=gene_ids_normalized,
    )
    expr = gctx.data_df.loc[:, lincs_ids]
    expr.columns = obs_filtered.index

    expr_df = pd.DataFrame(expr.T.to_numpy(np.float32), 
                           index=expr.columns.to_numpy(str), 
                           columns=expr.index.to_numpy(str)
                           )

    expr_df = expr_df.reindex(obs_filtered.index)
    if gene_ids_normalized is not None:
        expr_df = expr_df.loc[:, gene_ids_normalized]

    expr_df = expr_df.reindex(columns=gene_ids)
    
    return expr_df


def assemble_l1000_dataset(padata: ad.AnnData,
                         data_root: Optional[str] = None,
                         config: Optional[dict] = None,
                         **kwargs) -> ad.AnnData:
    """
    Assemble L1000 Level 3 data from GCTX files into standardized AnnData format.
    
    This function assembles L1000 Level 3 data by:
    1. Checking and downloading required data files if missing
    2. Loading and standardizing metadata
    3. Optionally annotating perturbations with PubChem CIDs
    4. Building standardized .obs dataframe
    5. Processing gene annotations
    6. Assembling final AnnData object
    
    Parameters:
    -----------
    padata : ad.AnnData
        Pseudobulk AnnData object (may be empty or placeholder)
    data_root : str, optional
        Path to L1000 data root directory containing GCTX and metadata files.
        If files are missing, they will be automatically downloaded.
    config : dict, optional
        Configuration dictionary with processing parameters.
        Can include:
        - 'dataset' (str, required) to specify the dataset name (l1000_phase1 or l1000_phase2)
        - 'download_if_missing' (bool, default=True) to control automatic downloads
        - 'annotate_pubchem' (bool, default=False) to annotate perturbations with PubChem CIDs
          (may involve API calls to PubChem and can be time-consuming)
    **kwargs
        Additional processing parameters
        
    Returns:
    --------
    ad.AnnData
        Processed AnnData object with standardized schema
    """
    logger.info('Applying L1000-specific processing')

    if not HAS_CMAPPY:
        raise ImportError("cmapPy is required for L1000 processing. Install via: pip install cmapPy")
    
    # Build configuration
    CONFIG = build_config(config)
    download_if_missing = CONFIG.get("download_if_missing", True)
    
    # Define file paths
    PATHS = define_paths(data_root, dataset=CONFIG.get("dataset"))
    
    # Check data availability and download if needed
    if download_if_missing:
        logger.info("Checking data availability...")
        status = check_l1000_files(data_root, dataset=CONFIG.get("dataset"))
        
        if status['missing']:
            logger.info(f"Downloading {len(status['missing'])} missing file(s)...")
            download_l1000_files(data_root, dataset=CONFIG.get("dataset"), skip_existing=True)
            status = check_l1000_files(data_root, dataset=CONFIG.get("dataset"))
        
        
        if status['compressed']:
            logger.info(f"Decompressing {len(status['compressed'])} compressed file(s)...")
            decompress_l1000_files(data_root, dataset=CONFIG.get("dataset"))
        
        # Verify all files are ready
        final_status = check_l1000_files(data_root, dataset=CONFIG.get("dataset"))
        if final_status['missing'] or final_status['compressed']:
            raise FileNotFoundError(
                f"Required L1000 data files are still missing or compressed after download attempt. "
                f"Missing: {len(final_status['missing'])}, Compressed: {len(final_status['compressed'])}"
            )
        logger.info("All required data files are available")
    
    FULL_GENE_MATRIX = CONFIG["full_gene_matrix"]
    
    # Load metadata tables
    inst_raw, cellinfo_raw, pert_raw, geneinfo, geneinfo_beta, cellinfo_beta = load_metadata_tables(PATHS)
    
    # Process cell info
    cellinfo = process_cellinfo(cellinfo_raw, cellinfo_extra=cellinfo_beta)

    # Process perturbation metadata
    pert_raw = process_pert_metadata(pert_raw, inst_raw)

    # Annotate compounds with PubChem CIDs (optional)
    if CONFIG.get("annotate_pubchem", False):
        logger.info("Mapping compounds to PubChem CIDs")
        pert_raw = annotate_pubchem_cids(pert_raw, PATHS, config=CONFIG)
    
    # Process instance metadata
    inst = process_instance_metadata(inst_raw, cellinfo, pert_raw, CONFIG, PATHS)
    
    logger.info('  Building .obs dataframe')
    obs = build_obs_dataframe(inst, dataset=CONFIG.get("dataset"))
    
    logger.info('  Enforcing strict obs schema')
    obs_for_schema = enforce_obs_schema(obs)
    
    logger.info('  Processing gene annotations')
    var = process_gene_annotations(geneinfo, geneinfo_beta, FULL_GENE_MATRIX)
    
    logger.info('  Matching obs to GCTX column IDs')
    obs_helpers = build_obs_helpers(obs, PATHS["level3_gctx"])
    
    
    logger.info('  Extracting expression from GCTX')
    expr = load_expression(obs_helpers, PATHS["level3_gctx"], gene_ids=var.index.astype(str))
    
    logger.info('  Assembling AnnData')
    padata_processed = assemble_anndata(obs_for_schema, var, expr)
    
    logger.info(f'  Assembled AnnData: {padata_processed.n_obs:,} × {padata_processed.n_vars:,}')
    logger.info('L1000-specific processing completed')
    
    return padata_processed
