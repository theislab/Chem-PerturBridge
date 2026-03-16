"""
DILImap training dataset assembly module.

Loads the training-only or training+validation combined DILImap datasets.
Reuses the same AnnData H5AD format as the validation (dilimap) dataset.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import anndata as ad

from src.utils.parsing_utils import logger


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def define_dilimap_train_paths(data_root: Optional[str] = None) -> dict:
    """
    Define file paths for the training-only DILImap dataset.

    Parameters
    ----------
    data_root : str, optional
        Root directory containing the training H5AD file.
        Defaults to ``'./op3_v2/data'``.

    Returns
    -------
    dict
        Dictionary mapping file identifiers to :class:`~pathlib.Path` objects.
    """
    if data_root is None:
        data_root = "./op3_v2/data"
    data_root = Path(data_root)

    return {
        "source_h5ad": data_root / "adata_training_counts.h5ad",
    }


def define_dilimap_train_val_paths(data_root: Optional[str] = None) -> dict:
    """
    Define file paths for the combined training+validation DILImap dataset.

    Parameters
    ----------
    data_root : str, optional
        Root directory containing both H5AD files.
        Defaults to ``'./op3_v2/data'``.

    Returns
    -------
    dict
        Dictionary mapping file identifiers to :class:`~pathlib.Path` objects.
    """
    if data_root is None:
        data_root = "./op3_v2/data"
    data_root = Path(data_root)

    return {
        "train_h5ad": data_root / "adata_training_counts.h5ad",
        "val_h5ad": data_root / "adata_validation_counts.h5ad",
    }


# ---------------------------------------------------------------------------
# File availability
# ---------------------------------------------------------------------------

def check_files(paths: dict) -> dict:
    """
    Check status of data files.

    Parameters
    ----------
    paths : dict
        Paths dict from a ``define_*`` function.

    Returns
    -------
    dict
        Dictionary with keys ``'missing'`` and ``'ready'``.
    """
    missing = []
    ready = []

    for key, path in paths.items():
        path = Path(path)
        if path.exists():
            ready.append((key, path))
        else:
            missing.append((key, path))

    if missing:
        logger.warning(f"{len(missing)} file(s) missing:")
        for key, path in missing:
            logger.warning(f"  - {key}: {path}")
    else:
        logger.info("All required files are present")

    return {"missing": missing, "ready": ready}


def ensure_data_available(paths: dict) -> None:
    """
    Verify that all source H5AD files exist; raise if not.

    Parameters
    ----------
    paths : dict
        Paths dict from a ``define_*`` function.

    Raises
    ------
    FileNotFoundError
        If any source file is missing.
    """
    status = check_files(paths)
    if status["missing"]:
        names = ", ".join(str(p) for _, p in status["missing"])
        raise FileNotFoundError(
            f"Required data files are missing: {names}. "
            "Please ensure the H5AD files are available in the data_root directory."
        )
    logger.info("All required data files are available")


# ---------------------------------------------------------------------------
# Main assembly entry-points
# ---------------------------------------------------------------------------

def assemble_dilimap_train_dataset(
    data_root: Optional[str] = None,
) -> ad.AnnData:
    """
    Load the training-only DILImap dataset.

    Parameters
    ----------
    data_root : str, optional
        Root directory containing the training H5AD file.

    Returns
    -------
    ad.AnnData
        The loaded training AnnData object.
    """
    logger.info("Assembling DILImap training dataset")

    PATHS = define_dilimap_train_paths(data_root)
    ensure_data_available(PATHS)

    source = Path(PATHS["source_h5ad"])
    logger.info(f"  Reading {source.name}")
    adata = ad.read_h5ad(source)
    logger.info(
        f"  Loaded AnnData: {adata.n_obs:,} observations × {adata.n_vars:,} variables"
    )

    logger.info("DILImap training dataset assembly completed")
    return adata


def assemble_dilimap_train_val_dataset(
    data_root: Optional[str] = None,
) -> ad.AnnData:
    """
    Load and concatenate the training + validation DILImap datasets.

    Parameters
    ----------
    data_root : str, optional
        Root directory containing both H5AD files.

    Returns
    -------
    ad.AnnData
        The concatenated AnnData object.
    """
    logger.info("Assembling DILImap training+validation dataset")

    PATHS = define_dilimap_train_val_paths(data_root)
    ensure_data_available(PATHS)

    train_path = Path(PATHS["train_h5ad"])
    val_path = Path(PATHS["val_h5ad"])

    logger.info(f"  Reading {train_path.name}")
    adata_train = ad.read_h5ad(train_path)
    logger.info(
        f"  Training: {adata_train.n_obs:,} observations × {adata_train.n_vars:,} variables"
    )

    logger.info(f"  Reading {val_path.name}")
    adata_val = ad.read_h5ad(val_path)
    logger.info(
        f"  Validation: {adata_val.n_obs:,} observations × {adata_val.n_vars:,} variables"
    )

    logger.info("  Concatenating training + validation")
    adata = ad.concat([adata_train, adata_val], join="inner")
    adata.obs_names_make_unique()
    logger.info(
        f"  Combined: {adata.n_obs:,} observations × {adata.n_vars:,} variables"
    )

    logger.info("DILImap training+validation dataset assembly completed")
    return adata
