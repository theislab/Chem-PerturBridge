"""
DILImap training dataset standardization.

Converts the DILImap AnnData to the common pseudobulk schema used by
the rest of the pipeline, mirroring standardize_sciplex in structure.

Schema reference: docs/Format_Pseudobulk.ipynb
"""
from anndata import AnnData
import pandas as pd
import numpy as np


# Compounds that are vehicle controls
_CONTROL_COMPOUNDS = {"DMSO"}

# DILImap uses primary human hepatocytes (PHH)
# Source: https://www.nature.com/articles/s41467-025-65690-3
# Cell culture: cryopreserved primary human hepatocytes from multiple adult donors
_CELL_TYPE = "CL:0000182"  # hepatocyte (Cell Ontology)
_DEVELOPMENT_STAGE = "adult"
_SEX = "unknown"  # mixed donors, not specified per sample


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
        "PLATE_NAME": "plate",
        "WELL_ID": "well",
        "LIBRARY_ID": "library",
        "SPECIES": "organism",
        "BATCH_ID": "batch",
    }
    adata.obs.rename(columns=rename_cols, inplace=True)

    # ── filter DMSO_replaced rows (QC flag per DILImap pipeline) ─────
    # DMSO_replaced marks wells that failed viability QC
    # Ref: https://www.dilimap.org/tutorials/1_Compute_Pathway_Signatures.html
    dmso_replaced_mask = adata.obs["perturbagen"] == "DMSO_replaced"
    if dmso_replaced_mask.any():
        n_dropped = dmso_replaced_mask.sum()
        adata = adata[~dmso_replaced_mask].copy()

    # ── controls ──────────────────────────────────────────────────────
    adata.obs["is_control"] = (
        adata.obs["perturbagen"].isin(_CONTROL_COMPOUNDS)
        | (adata.obs["DOSE_LEVEL"] == "Control")
    )

    # Zero out dose for controls
    adata.obs.loc[adata.obs["is_control"], "pert_dose_uM"] = 0.0

    # ── organism label ────────────────────────────────────────────────
    adata.obs["organism"] = (
        adata.obs["organism"]
        .str.strip()
        .str.lower()
    )

    # ── columns expected by the common schema ─────────────────────────
    adata.obs["cell_type"] = _CELL_TYPE
    adata.obs["pert_type"] = "compound"
    adata.obs["suspension_type"] = "cell"
    adata.obs["tissue"] = "liver"
    adata.obs["tissue_type"] = "cell culture"
    adata.obs["disease"] = "normal"
    adata.obs["stimulation"] = None
    adata.obs["guide"] = None
    adata.obs["dataset"] = "dilimap"
    adata.obs["assay"] = "SMARTSeq bulk RNA-seq"
    adata.obs["development_stage"] = _DEVELOPMENT_STAGE
    adata.obs["sex"] = _SEX
    adata.obs["self_reported_ethnicity"] = "unknown"
    adata.obs["pubchem_cid"] = None
    adata.obs["psbulk_cells"] = -666  # sentinel: bulk data, cell count not available
    adata.obs["psbulk_counts"] = np.round(np.asarray(adata.X.sum(axis=1)).flatten()).astype(int)

    # ── sample_id (composite key) ─────────────────────────────────────
    adata.obs["sample_id"] = (
        adata.obs["plate"].astype(str) + "_"
        + adata.obs["well"].astype(str) + "_"
        + adata.obs["perturbagen"].astype(str).str.replace(" ", "_", regex=False) + "_"
        + adata.obs["pert_dose_uM"].astype(str)
    )

    # ── enforce dtypes per Format_Pseudobulk.ipynb schema ─────────────
    _category_cols = [
        "plate", "well", "cell_type", "perturbagen", "pert_type",
        "is_control", "suspension_type", "tissue", "tissue_type",
        "disease", "library", "stimulation", "guide", "dataset",
        "assay", "development_stage", "organism", "sex",
        "self_reported_ethnicity", "pubchem_cid",
    ]
    for col in _category_cols:
        adata.obs[col] = adata.obs[col].astype("category")

    adata.obs["pert_dose_uM"] = adata.obs["pert_dose_uM"].astype("float64")
    adata.obs["pert_time_h"] = adata.obs["pert_time_h"].astype("float64")
    adata.obs["psbulk_cells"] = adata.obs["psbulk_cells"].astype("int64")
    adata.obs["psbulk_counts"] = adata.obs["psbulk_counts"].astype("int64")

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
    adata.obs["dataset"] = adata.obs["dataset"].astype("category")
    return adata


def standardize_dilimap_train_val(adata: AnnData) -> AnnData:
    """
    Standardize DILImap training+validation dataset to the common schema.

    Delegates to :func:`standardize_dilimap` and overrides the
    ``dataset`` column to ``"dilimap_train_val"``.
    """
    adata = standardize_dilimap(adata)
    adata.obs["dataset"] = "dilimap_train_val"
    adata.obs["dataset"] = adata.obs["dataset"].astype("category")
    return adata
