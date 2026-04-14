"""
Ginkgo GDPx2 bulk RNA dataset standardization script.

Standardizes the LaminLabs-downloaded GDPx2 AnnData to the unified pseudobulk schema.

Usage:
  python3 -m src.pseudobulking.datasets.gdpx2.run_standardization \
      --data_root ./data/gdpx2 \
      --output_file ./data/gdpx2/pseudobulk/full/gdpx2_standardized.h5ad \
      --annotate-pubchem
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.gdpx2.standardization import standardize_gdpx2_dataset


DATASET_TITLE = "Ginkgo GDPx2"

# LaminLabs artifact keys (laminlabs/pertdata)
LAMIN_KEY_ADATA      = "ginkgo-datapoints/vcpi/X.h5ad"
LAMIN_KEY_OBS        = "ginkgo-datapoints/vcpi/obs.parquet"
LAMIN_KEY_VAR        = "ginkgo-datapoints/vcpi/var.parquet"
LAMIN_KEY_COMPOUNDS  = "ginkgo-datapoints/vcpi/compounds-GDPx2-2026-02-09.csv"
LAMIN_INSTANCE       = "laminlabs/pertdata"


def get_gdpx2_paths(data_root: str | Path) -> dict:
    """
    Build a dictionary of all file paths required by the GDPx2 pipeline.

    Parameters
    ----------
    data_root:
        Root directory for the GDPx2 dataset (e.g. ./data/gdpx2).

    Returns
    -------
    dict mapping file identifiers to Path objects.
    """
    raw_dir = Path(data_root) / "raw"
    return {
        "raw_h5ad":          raw_dir / "gdpx2_raw.h5ad",
        "raw_obs":           raw_dir / "gdpx2_obs.parquet",
        "raw_var":           raw_dir / "gdpx2_var.parquet",
        "compound_csv":      raw_dir / "gdpx2_compounds.csv",
        "pubchem_cid_cache": raw_dir / "gdpx2_pubchem_cache.json",
        "dataset_title":     DATASET_TITLE,
        "lamin_key_adata":   LAMIN_KEY_ADATA,
        "lamin_key_obs":     LAMIN_KEY_OBS,
        "lamin_key_var":     LAMIN_KEY_VAR,
    }


def download_compound_csv(paths: dict, force: bool = False) -> None:
    """
    Download the GDPx2 compound metadata CSV from LaminLabs and save it locally.

    Skipped automatically if the file already exists, unless force=True.

    Parameters
    ----------
    paths:
        Path dict as returned by get_gdpx2_paths().
    force:
        Re-download even if the local file already exists.
    """
    dest = paths["compound_csv"]
    if dest.is_file() and not force:
        logger.info("Compound CSV already exists, skipping: %s", dest)
        return

    try:
        import lamindb as ln
    except ImportError as e:
        raise ImportError(
            "Install lamindb to download the compound CSV, or place "
            f"gdpx2_compounds.csv at {dest}"
        ) from e

    dest.parent.mkdir(parents=True, exist_ok=True)

    try:
        ln.setup.disconnect()
    except Exception:
        pass

    max_retries, delay = 10, 10
    for attempt in range(max_retries):
        try:
            ln.connect(LAMIN_INSTANCE)
            break
        except Exception as e:
            logger.error("LaminDB connection failed: %s (attempt %d)", e, attempt + 1)
            if attempt + 1 == max_retries:
                raise
            time.sleep(delay)

    logger.info("Downloading compound CSV from LaminLabs: %s", LAMIN_KEY_COMPOUNDS)
    df = ln.Artifact.get(key=LAMIN_KEY_COMPOUNDS).load(mute=True)
    ln.setup.disconnect()

    df.to_csv(dest, index=False)
    logger.info("Saved compound CSV: %s (%d rows)", dest, len(df))


def main():
    parser = argparse.ArgumentParser(
        description="Standardize Ginkgo GDPx2 bulk RNA dataset to the unified pseudobulk schema"
    )
    parser.add_argument(
        "--data_root", type=str, required=True,
        help="Root directory for raw GDPx2 files (expects raw/gdpx2_raw.h5ad etc.)",
    )
    parser.add_argument(
        "--output_file", type=str, required=True,
        help="Output path for the standardized H5AD file",
    )
    parser.add_argument(
        "--annotate-pubchem", action="store_true",
        help="Annotate compounds with PubChem CIDs (requires network access)",
    )
    args = parser.parse_args()

    if os.path.isfile(args.output_file):
        logger.info("Standardized dataset already exists: %s", args.output_file)
        return

    paths = get_gdpx2_paths(args.data_root)

    logger.info("Downloading GDPx2 compound metadata")
    download_compound_csv(paths)

    logger.info("Standardizing GDPx2 dataset")
    adata_standardized = standardize_gdpx2_dataset(
        paths=paths,
        annotate_pubchem=args.annotate_pubchem,
    )

    logger.info("Saving standardized dataset to: %s", args.output_file)
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    adata_standardized.write_h5ad(args.output_file, compression="gzip")
    logger.info("GDPx2 standardization complete")


if __name__ == "__main__":
    main()
