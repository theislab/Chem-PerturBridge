"""
DILImap training dataset standardization.

Converts the DILImap AnnData to the common pseudobulk schema used by
the rest of the pipeline, mirroring standardize_sciplex in structure.
"""
from anndata import AnnData
import pandas as pd
import numpy as np


# Compounds that are vehicle controls
_CONTROL_COMPOUNDS = {"DMSO", "DMSO_replaced"}


def standardize_obs_dilimap(adata: AnnData) -> AnnData:
    """
    Standardize DILImap observations to the common schema.

    Renames columns, flags controls, normalises control compound names,
    and adds the bookkeeping columns expected downstream.

    Parameters
    ----------
    adata : AnnData
        Raw DILImap dataset.

    Returns
    -------
    AnnData
        Dataset with standardized ``.obs``.
    """
    adata = adata.copy()

    # ── rename to standard names ──────────────────────────────────────
    rename_cols = {
        "COMPOUND": "perturbagen",
        "CONCENTRATION_UM": "pert_dose_uM",
        "TIMEPOINT_HOURS": "pert_time_h",
        "PLATE_ID": "plate",
        "WELL_ID": "well",
        "LIBRARY_ID": "library_id",
        "SPECIES": "organism",
        "BATCH_ID": "batch",
    }
    adata.obs.rename(columns=rename_cols, inplace=True)

    # ── controls ──────────────────────────────────────────────────────
    adata.obs["is_control"] = (
        adata.obs["perturbagen"].isin(_CONTROL_COMPOUNDS)
        | (adata.obs["DOSE_LEVEL"] == "Control")
    )

    # Normalise control compound name to 'DMSO'
    adata.obs.loc[
        adata.obs["perturbagen"] == "DMSO_replaced", "perturbagen"
    ] = "DMSO"

    # Zero out dose for controls
    adata.obs.loc[adata.obs["is_control"], "pert_dose_uM"] = 0.0

    # ── organism label ────────────────────────────────────────────────
    adata.obs["organism"] = (
        adata.obs["organism"]
        .str.strip()
        .str.lower()
        .replace({"human": "human"})      # already lowercase after strip
    )

    # ── columns expected by the common schema ─────────────────────────
    adata.obs["cell_type"] = "hepatocyte"  # bulk primary hepatocytes
    adata.obs["pert_type"] = np.where(
        adata.obs["is_control"], "ctl_vehicle", "trt_cp"
    )
    adata.obs["suspension_type"] = "cell"
    adata.obs["tissue"] = "liver"
    adata.obs["tissue_type"] = "cell culture"
    adata.obs["disease"] = "normal"
    adata.obs["library"] = pd.Categorical([None] * len(adata))
    adata.obs["stimulation"] = pd.Categorical([None] * len(adata))
    adata.obs["guide"] = pd.Categorical([None] * len(adata))
    adata.obs["dataset"] = "dilimap"
    adata.obs["assay"] = "SMARTSeq bulk RNA-seq"
    adata.obs["development_stage"] = "unknown"
    adata.obs["sex"] = "unknown"
    adata.obs["self_reported_ethnicity"] = "unknown"
    adata.obs["pubchem_cid"] = pd.array([pd.NA] * len(adata), dtype="Int64")
    adata.obs["psbulk_cells"] = 1
    adata.obs["psbulk_counts"] = 1

    # ── sample_id (composite key) ─────────────────────────────────────
    adata.obs["sample_id"] = (
        adata.obs["plate"].astype(str) + "_"
        + adata.obs["well"].astype(str) + "_"
        + adata.obs["perturbagen"].astype(str).str.replace(" ", "_", regex=False) + "_"
        + adata.obs["pert_dose_uM"].astype(str)
    )

    return adata


def standardize_var_dilimap(adata: AnnData) -> AnnData:
    """
    Standardize DILImap variables to the common schema.

    Moves Ensembl IDs from the ``gene_id`` column to the var index and
    keeps gene symbols accessible via a ``symbol`` column.

    Parameters
    ----------
    adata : AnnData
        DILImap dataset (obs already standardized or not).

    Returns
    -------
    AnnData
        Dataset with standardized ``.var``.
    """
    adata = adata.copy()

    # Current state: var.index = gene symbols, var['gene_id'] = Ensembl IDs
    adata.var["symbol"] = adata.var.index.astype(str)
    adata.var["ensembl_id"] = adata.var["gene_id"].astype(str)

    # Use Ensembl IDs as the canonical index
    adata.var.index = adata.var["ensembl_id"].values
    adata.var.index.name = None

    # Keep only the columns downstream expects
    adata.var = adata.var[["symbol", "ensembl_id"]]

    return adata


def standardize_dilimap(adata: AnnData) -> AnnData:
    """
    Full standardization of a DILImap AnnData object.

    Parameters
    ----------
    adata : AnnData
        Raw DILImap dataset.

    Returns
    -------
    AnnData
        Standardized dataset.
    """
    adata = standardize_obs_dilimap(adata)
    adata = standardize_var_dilimap(adata)
    return adata


def standardize_dilimap_train(adata: AnnData) -> AnnData:
    """
    Standardize DILImap training dataset to the common schema.

    Delegates to :func:`standardize_dilimap` and overrides the
    ``dataset`` column to ``"dilimap_train"``.
    """
    adata = standardize_dilimap(adata)
    adata.obs["dataset"] = "dilimap_train"
    return adata


def standardize_dilimap_train_val(adata: AnnData) -> AnnData:
    """
    Standardize DILImap training+validation dataset to the common schema.

    Delegates to :func:`standardize_dilimap` and overrides the
    ``dataset`` column to ``"dilimap_train_val"``.
    """
    adata = standardize_dilimap(adata)
    adata.obs["dataset"] = "dilimap_train_val"
    return adata
