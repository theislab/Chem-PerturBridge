"""Dataset-specific post-processing for CIGS, applied in the DEG pipeline."""
import anndata as ad

from src.utils.parsing_utils import logger


def process_cigs_dataset(padata: ad.AnnData) -> ad.AnnData:
    """Backfill `well` from `sample_unique_id` so TCM `separate_replicates`
    labels stay unique (TCM is HiMAP-seq, no wells)."""
    padata = padata.copy()

    if "sample_unique_id" in padata.obs.columns and "well" in padata.obs.columns:
        missing_well = padata.obs["well"].isna()
        if missing_well.any():
            logger.info(
                "CIGS post_processing: backfilling well from sample_unique_id "
                f"for {int(missing_well.sum())} samples"
            )
            # Cast to object first — Categorical setitem across disjoint categories raises.
            well = padata.obs["well"].astype("object")
            suid = padata.obs["sample_unique_id"].astype("object")
            padata.obs["well"] = well.combine_first(suid).astype("category")

    return padata
