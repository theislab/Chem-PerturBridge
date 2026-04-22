"""Dataset-specific post-processing for CIGS, applied in the DEG pipeline."""
import anndata as ad

from src.utils.parsing_utils import logger


def process_cigs_dataset(padata: ad.AnnData) -> ad.AnnData:
    """Backfill `well` from `sample_unique_id` for TCM samples.

    CIGS TCM uses HiMAP-seq (no wells), so `well` is NaN. DEG builds
    `separate_replicates` labels as `{pert}_{dose}uM_{time}h_{well}_{plate}`,
    which would yield literal "nan" tokens. Using `sample_unique_id` (already
    prefixed with `{subset_key}.`) keeps replicate labels unique and readable.
    MCE rows (which have real wells) are untouched.
    """
    padata = padata.copy()
    obs = padata.obs

    if "sample_unique_id" in obs.columns and "well" in obs.columns:
        missing_well = obs["well"].isna()
        if missing_well.any():
            logger.info(
                "CIGS post_processing: backfilling well from sample_unique_id "
                f"for {int(missing_well.sum())} samples"
            )
            obs.loc[missing_well, "well"] = obs.loc[missing_well, "sample_unique_id"]
            padata.obs = obs

    return padata
