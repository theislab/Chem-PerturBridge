"""
VCPI Ginkgo bulk RNA dataset standardization script.

This script standardizes the raw VCPI Ginkgo AnnData to the unified pseudobulk schema.
"""
from __future__ import annotations

import os
import argparse
from pathlib import Path

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.vcpi_ginkgo.downloading_formatting_data import load_vcpi_ginkgo
from src.pseudobulking.datasets.vcpi_ginkgo.standardization import standardize_vcpi_ginkgo_dataset


EXPERIMENT_TO_SLUG = {
    "vcpi-0001": "vcpi_0001",
    "vcpi-0002": "vcpi_0002",
    "vcpi-0003": "vcpi_0003",
}

EXPERIMENT_TITLES = {
    "vcpi-0001": "Ginkgo VCPI vcpi-0001 (tvc-bhr-009)",
    "vcpi-0002": "Ginkgo VCPI vcpi-0002 (tvc-kdl-010)",
    "vcpi-0003": "Ginkgo VCPI vcpi-0003 (tvc-qnu-012)",
}

# The VCPI client resolves downloads by opaque job id (tvc-...), not by the
# human-readable experiment name (vcpi-000X), so map the latter to the former.
EXPERIMENT_TO_JOB_ID = {
    "vcpi-0001": "tvc-bhr-009",
    "vcpi-0002": "tvc-kdl-010",
    "vcpi-0003": "tvc-qnu-012",
}


def get_vcpi_ginkgo_paths(data_root: str | Path, experiment_id: str) -> dict:
    """
    Build a dictionary of all file paths required by the VCPI Ginkgo pipeline.

    Parameters
    ----------
    data_root :
        Root directory for the VCPI Ginkgo dataset (e.g. ./data/vcpi_0001).
    experiment_id :
        VCPI experiment identifier (e.g. ``vcpi-0001``).

    Returns
    -------
    dict mapping file identifiers to Path objects (plus ``dataset_title`` string).
    """
    if experiment_id not in EXPERIMENT_TO_SLUG:
        raise ValueError(
            f"Unknown experiment {experiment_id!r}; supported: {sorted(EXPERIMENT_TO_SLUG)}"
        )
    raw_dir = Path(data_root) / "raw"
    return {
        # loading intermediates
        "experiments_pkl":    raw_dir / "experiments.pkl",
        # standardization inputs
        "raw_h5ad":           raw_dir / "vcpi_ginkgo_raw.h5ad",
        "compound_csv":       raw_dir / "df_compounds.csv",
        # PubChem caches
        "pubchem_cid_cache":  raw_dir / "vcpi_ginkgo_pubchem_cache.json",
        "pubchem_names_cache": raw_dir / "vcpi_ginkgo_pubchem_names_cache.json",
        # Ensembl symbol cache
        "ensembl_symbol_cache": raw_dir / "vcpi_ginkgo_ensembl_symbol_cache.json",
        # experiment metadata
        "dataset_title":      EXPERIMENT_TITLES[experiment_id],
        # opaque VCPI job id used for the actual download (vcpi.load_experiment)
        "job_id":             EXPERIMENT_TO_JOB_ID[experiment_id],
    }


def main():
    """
    Standardize a VCPI Ginkgo bulk RNA experiment to the unified pseudobulk schema.

    This is the main entry point for the VCPI Ginkgo standardization pipeline. It:
    1. Parses command-line arguments for configuration
    2. Checks if output file already exists (exits early if present)
    3. Loads / downloads raw files and converts payload → h5ad + compound CSV
       (each step skipped if already done)
    4. Standardizes the dataset
    5. Saves the standardized dataset as an H5AD file

    Command-line Arguments
    ----------------------
    --experiment : str, required
        VCPI experiment id (vcpi-0001, vcpi-0002 or vcpi-0003).
    --output_file : str, required
        Output path for the standardized H5AD file.
    --data_root : str, required
        Root directory for raw VCPI files
        (expects raw/experiments.pkl or network access via TVC_TOKEN).
    --annotate-pubchem : bool, default=False
        If set, annotates compounds with PubChem CIDs (may involve API calls).
    --annotate-pubchem-names : bool, default=False
        If set, additionally looks up drug synonym names by PubChem CID
        (requires --annotate-pubchem or pre-existing pubchem_cid column).
    """
    parser = argparse.ArgumentParser(
        description='Standardize VCPI Ginkgo bulk RNA dataset to unified pseudobulk schema'
    )
    parser.add_argument('--experiment', type=str, required=True,
                        choices=sorted(EXPERIMENT_TO_SLUG),
                        help='VCPI experiment id (e.g. vcpi-0001)')
    parser.add_argument('--output_file', type=str, required=True,
                        help='Output path for standardized H5AD file')
    parser.add_argument('--data_root', type=str, required=True,
                        help='Root directory for raw VCPI Ginkgo files')
    parser.add_argument('--annotate-pubchem', action='store_true',
                        help='Annotate compounds with PubChem CIDs (may involve API calls)')
    parser.add_argument('--annotate-pubchem-names', action='store_true',
                        help='Look up drug synonym names by PubChem CID (requires --annotate-pubchem or pre-existing pubchem_cid column)')
    args = parser.parse_args()

    if os.path.isfile(args.output_file):
        logger.info("Standardized dataset already exists: %s", args.output_file)
        return

    paths = get_vcpi_ginkgo_paths(args.data_root, args.experiment)

    logger.info("Loading and converting VCPI Ginkgo raw data")
    load_vcpi_ginkgo(paths, args.experiment)

    logger.info("Standardizing VCPI Ginkgo dataset")
    adata_standardized = standardize_vcpi_ginkgo_dataset(
        paths=paths,
        annotate_pubchem=args.annotate_pubchem,
        annotate_pubchem_names=args.annotate_pubchem_names,
    )

    logger.info("Saving standardized dataset to: %s", args.output_file)
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    adata_standardized.write_h5ad(args.output_file, compression='gzip')
    logger.info("VCPI Ginkgo standardization complete")


if __name__ == "__main__":
    main()
