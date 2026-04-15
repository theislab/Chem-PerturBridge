"""
pkl_to_adata.py  —  VCPI Ginkgo pickle → AnnData conversion helpers.

Python equivalent of ``novartis/rdata_to_adata.R``.
Converts the serialised VCPI experiment payload to a raw AnnData and extracts
the compound metadata table.  Import and call from ``downloading_formatting_data.py``.
"""
from __future__ import annotations

from typing import Any, Mapping

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.utils.parsing_utils import logger


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def vcpi_experiments_to_adata(experiments: Mapping[str, Any], *, sparse: bool = True) -> ad.AnnData:
    """
    Convert a VCPI experiment payload (dict with ``metadata`` and ``data`` keys) to AnnData.

    Parameters
    ----------
    experiments:
        Payload returned by ``vcpi.load_experiment`` or loaded from a pickle file.
    sparse:
        If True, store the expression matrix as a CSR sparse matrix.

    Returns
    -------
    AnnData with obs.index = sequenced_id and var.index = gene_id.
    """
    meta = experiments["metadata"]
    expr = experiments["data"]

    obs = meta.to_pandas().astype(str).drop_duplicates().set_index("sequenced_id")

    sample_cols = [c for c in expr.columns if c != "gene_id"]
    common = [s for s in sample_cols if s in obs.index]
    if not common:
        raise ValueError("No overlap between expression sample columns and metadata sequenced_id index")

    mat = expr.select(common).to_numpy()
    x = mat.T.astype(np.int64)
    logger.info("Expression matrix shape (samples × genes): %s", x.shape)

    if sparse:
        x = sp.csr_matrix(x)

    var = expr.select("gene_id").to_pandas().set_index("gene_id")
    obs = obs.loc[common].reindex(common)

    return ad.AnnData(X=x, obs=obs, var=var)


def extract_compound_df(experiments: Mapping[str, Any]) -> pd.DataFrame:
    """
    Extract the compound metadata table from the VCPI payload.

    Parameters
    ----------
    experiments:
        Payload returned by ``vcpi.load_experiment`` or loaded from a pickle file.

    Returns
    -------
    De-duplicated compound DataFrame (one row per compound entry).
    """
    chem = experiments.get("chemistry")
    if chem is None:
        raise ValueError(
            "experiments['chemistry'] is missing; cannot extract compound metadata."
        )
    df = chem.to_pandas() if hasattr(chem, "to_pandas") else pd.DataFrame(chem)
    return df.drop_duplicates().reset_index(drop=True)


