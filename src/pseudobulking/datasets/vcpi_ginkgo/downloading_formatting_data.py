from __future__ import annotations

import os
import pickle
from pathlib import Path
from typing import Union

import anndata as ad
import pandas as pd

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.vcpi_ginkgo.pkl_to_adata import (
    vcpi_experiments_to_adata,
    extract_compound_df,
)

# ---------------------------------------------------------------------------
# Step 1 — load or download experiments payload
# ---------------------------------------------------------------------------

def load_vcpi_experiment(
    paths: dict,
    experiment_id: str,
    force: bool = False,
) -> dict:
    """
    Load the VCPI experiment payload from a local pickle or download via the
    VCPI client package.

    The payload is saved to ``raw/experiments.pkl`` on first download so that
    subsequent runs skip the network call.

    Parameters
    ----------
    paths :
        Path dict as returned by ``get_vcpi_ginkgo_paths()``.
    experiment_id :
        VCPI experiment identifier (e.g. ``vcpi-0001``).
    force :
        Re-download even if ``experiments.pkl`` already exists.

    Returns
    -------
    The raw experiments payload dict.
    """
    experiments_pkl = paths["experiments_pkl"]
    experiments_pkl.parent.mkdir(parents=True, exist_ok=True)

    if experiments_pkl.is_file() and not force:
        logger.info("Loading VCPI payload from pickle: %s", experiments_pkl)
        with open(experiments_pkl, "rb") as f:
            return pickle.load(f)

    try:
        import vcpi
    except ImportError as e:
        raise ImportError(
            "Install the VCPI client package and set TVC_TOKEN, or place "
            f"experiments.pkl at {experiments_pkl}"
        ) from e

    if not os.environ.get("TVC_TOKEN"):
        logger.warning("TVC_TOKEN is not set; vcpi.load_experiment may fail authentication")

    logger.info("Downloading experiment via vcpi.load_experiment(%r)", experiment_id)
    experiments = vcpi.load_experiment(experiment_id)
    with open(experiments_pkl, "wb") as f:
        pickle.dump(experiments, f)
    logger.info("Saved experiments pickle: %s", experiments_pkl)
    return experiments


# ---------------------------------------------------------------------------
# Step 2 — convert experiments payload → raw h5ad
# ---------------------------------------------------------------------------

def convert_experiments_to_adata(
    experiments: dict,
    paths: dict,
    force: bool = False,
) -> ad.AnnData:
    """
    Convert the VCPI experiments payload to a raw AnnData and save it as
    ``raw/vcpi_ginkgo_raw.h5ad``.

    Parameters
    ----------
    experiments :
        Payload returned by ``load_vcpi_experiment()``.
    paths :
        Path dict as returned by ``get_vcpi_ginkgo_paths()``.
    force :
        Re-run conversion even if the ``.h5ad`` already exists.

    Returns
    -------
    The raw AnnData object.
    """
    raw_h5ad = paths["raw_h5ad"]
    raw_h5ad.parent.mkdir(parents=True, exist_ok=True)

    if raw_h5ad.is_file() and not force:
        logger.info("Raw AnnData already exists, skipping: %s", raw_h5ad)
        return ad.read_h5ad(raw_h5ad)

    logger.info("Converting VCPI payload to AnnData ...")
    adata = vcpi_experiments_to_adata(experiments, sparse=True)
    adata.write_h5ad(raw_h5ad, compression="gzip")
    logger.info("Saved raw AnnData: %s", raw_h5ad)
    return adata


# ---------------------------------------------------------------------------
# Step 3 — extract and format compound metadata
# ---------------------------------------------------------------------------

def format_compound_df(
    experiments: dict,
    paths: dict,
    force: bool = False,
) -> pd.DataFrame:
    """
    Extract the compound metadata table from ``experiments['chemistry']``,
    normalise column names, fill optional chemistry columns, and deduplicate
    by ``user_compound_id``.

    The result is saved to ``raw/df_compounds.csv``.

    Parameters
    ----------
    experiments :
        Payload returned by ``load_vcpi_experiment()``.
    paths :
        Path dict as returned by ``get_vcpi_ginkgo_paths()``.
    force :
        Re-extract even if ``df_compounds.csv`` already exists.

    Returns
    -------
    De-duplicated, normalised compound DataFrame.
    """
    compound_csv = paths["compound_csv"]
    compound_csv.parent.mkdir(parents=True, exist_ok=True)

    if compound_csv.is_file() and not force:
        logger.info("Compound CSV already exists, skipping: %s", compound_csv)
        return pd.read_csv(compound_csv)

    logger.info("Extracting compound metadata from experiments['chemistry'] ...")
    df = extract_compound_df(experiments)

    df.to_csv(compound_csv, index=False)
    logger.info("Saved compound metadata: %s (%d rows)", compound_csv, len(df))
    return df


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def load_vcpi_ginkgo(
    paths: dict,
    experiment_id: str,
    force_download: bool = False,
    force_convert: bool = False,
    force_compounds: bool = False,
) -> None:
    """
    Full pipeline: load / download experiments payload → convert to h5ad →
    extract and format compound metadata.

    Each step is skipped automatically when its output already exists, unless
    the corresponding ``force_*`` flag is set.

    Parameters
    ----------
    paths :
        Path dict as returned by ``get_vcpi_ginkgo_paths()``.
    experiment_id :
        VCPI experiment identifier (e.g. ``vcpi-0001``).
    force_download :
        Re-download ``experiments.pkl`` even if it already exists.
    force_convert :
        Re-run pkl → h5ad conversion even if the ``.h5ad`` already exists.
    force_compounds :
        Re-extract ``df_compounds.csv`` even if it already exists.
    """
    experiments = load_vcpi_experiment(paths, experiment_id, force=force_download)
    convert_experiments_to_adata(experiments, paths, force=force_convert)
    format_compound_df(experiments, paths, force=force_compounds)
