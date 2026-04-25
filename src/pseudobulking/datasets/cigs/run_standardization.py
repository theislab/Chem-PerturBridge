"""
Standardize one CIGS source (MCE or TCM) to the unified pseudobulk schema.

Downloads the Excel files for the requested source from
https://cigs.iomicscloud.com/, converts them to AnnData, runs the
source-specific standardization pass, and writes a gzip-compressed
.h5ad to ``--output_file``.

MCE and TCM use different gene panels (e.g. HLA names differ: hyphens
vs dots, plus Excel MARCH/SEPT date-autocorrupt artefacts on the MCE
side) so they are never merged at the raw level. Each source is a
separate pseudobulk dataset (``cigs_mce`` / ``cigs_tcm``) with its own
data tree, mirroring the ``l1000_phase1``/``l1000_phase2`` layout.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.cigs.downloading_formatting_data import (
    download_cigs_files,
    load_and_merge_cigs_subsets,
)
from src.pseudobulking.datasets.cigs.standardization import standardize_cigs_dataset


def get_cigs_paths(data_root: str | Path, source: str) -> dict:
    """Build path dict for a single CIGS source (``"mce"`` or ``"tcm"``)."""
    data_root = Path(data_root)
    raw_dir   = data_root / "raw"
    processed = raw_dir / "processed"
    return {
        "raw_dir":           raw_dir,
        "raw_mce_h5ad":      raw_dir / "cigs_mce_raw.h5ad",
        "raw_tcm_h5ad":      raw_dir / "cigs_tcm_raw.h5ad",
        "compounds_mce_csv": raw_dir / "compounds_mce.csv",
        "compounds_tcm_csv": raw_dir / "compounds_tcm.csv",
        "ensembl_cache":     processed / "ensembl_cache.json",
        "pubchem_cid_cache": processed / "pubchem_cid_cache.json",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Standardize one CIGS source (MCE or TCM) to the unified pseudobulk schema"
    )
    parser.add_argument(
        "--source", choices=["mce", "tcm"], required=True,
        help="Which CIGS source to process",
    )
    parser.add_argument(
        "--output_file", type=str, required=True,
        help="Path to the output standardized .h5ad file",
    )
    parser.add_argument(
        "--data_root", type=str, required=True,
        help="Per-source data root (e.g. ./data/cigs_mce)",
    )
    parser.add_argument(
        "--annotate-pubchem", action="store_true",
        help="Annotate compounds with PubChem CIDs (may involve API calls)",
    )
    parser.add_argument(
        "--force-download", action="store_true",
        help="Re-download Excel files even if they already exist",
    )
    parser.add_argument(
        "--force-convert", action="store_true",
        help="Re-run Excel → h5ad conversion even if per-subset .h5ad files exist",
    )

    args = parser.parse_args()

    source = args.source
    output_file = Path(args.output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if output_file.is_file() and not args.force_convert:
        logger.info("Standardized %s dataset already exists: %s", source.upper(), output_file)
        return

    paths = get_cigs_paths(args.data_root, source)
    paths["raw_dir"].mkdir(parents=True, exist_ok=True)
    paths["ensembl_cache"].parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading CIGS %s raw Excel files ...", source.upper())
    download_cigs_files(
        data_root=paths["raw_dir"],
        skip_existing=not args.force_download,
        source=source,
    )

    raw_h5ad = paths[f"raw_{source}_h5ad"]
    if raw_h5ad.exists() and not args.force_convert:
        logger.info("Loading merged raw %s AnnData from %s", source.upper(), raw_h5ad)
        adata_raw = ad.read_h5ad(raw_h5ad)
    else:
        logger.info("Converting Excel files to AnnData and merging %s subsets ...", source.upper())
        adata_raw = load_and_merge_cigs_subsets(
            paths=paths,
            source=source,
            force_convert=args.force_convert,
        )
        logger.info("Saving merged raw %s AnnData to %s ...", source.upper(), raw_h5ad)
        adata_raw.write_h5ad(raw_h5ad, compression="gzip")

    logger.info("Standardizing CIGS %s subset ...", source.upper())
    adata_std = standardize_cigs_dataset(
        paths=paths,
        source=source,
        adata_raw=adata_raw,
        annotate_pubchem=args.annotate_pubchem,
    )

    logger.info("Saving standardized %s dataset to %s ...", source.upper(), output_file)
    adata_std.write_h5ad(output_file, compression="gzip")

    logger.info("CIGS %s standardization complete.", source.upper())


if __name__ == "__main__":
    main()
