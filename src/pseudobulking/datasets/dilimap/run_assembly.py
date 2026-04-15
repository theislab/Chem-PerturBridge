"""
DILImap training dataset assembly script.

Loads the DILImap training (or training+validation) H5AD file(s), applies
standardization, and saves to the pipeline-expected output path.
"""
import os
import argparse
from pathlib import Path

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.dilimap.assembling import (
    assemble_dilimap_train_dataset,
    assemble_dilimap_train_val_dataset,
)
from src.pseudobulking.datasets.dilimap.standardization import standardize_dilimap


def get_dilimap_paths(data_root: str) -> dict:
    """
    Build a dictionary of all file paths required by the DILImap pipeline.

    Parameters
    ----------
    data_root :
        Root directory for the DILImap dataset (e.g. ``./data/dilimap_train/raw``).

    Returns
    -------
    dict mapping file identifiers to Path objects.
    """
    raw_dir = Path(data_root)
    return {
        "pubchem_cache": raw_dir / "dilimap_pubchem_cache.json",
    }


DATASET_LABELS = {
    "train": "dilimap_train",
    "train_val": "dilimap_train_val",
}


def main():
    """
    Assemble DILImap training dataset, standardize, and save as H5AD.

    Command-line Arguments
    ----------------------
    --mode : str, required
        Either 'train' (training only) or 'train_val' (training + validation).
    --output_file : str, required
        Output path for the assembled H5AD file.
    --data_root : str, required
        Root directory containing the DILImap H5AD source file(s).
    """
    parser = argparse.ArgumentParser(
        description="Assemble DILImap training dataset into pipeline-standard H5AD"
    )
    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=["train", "train_val"],
        help="Assembly mode: 'train' for training only, 'train_val' for training+validation",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output path for assembled H5AD file",
    )
    parser.add_argument(
        "--data_root",
        type=str,
        required=True,
        help="Root directory containing the DILImap H5AD file(s) (e.g. ./data/dilimap_train/raw)",
    )
    parser.add_argument(
        "--annotate_pubchem",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Enrich compounds with SMILES and look up PubChem CIDs via API (slow, requires network)",
    )
    args = parser.parse_args()

    # Check if output already exists
    if os.path.isfile(args.output_file):
        logger.info(f"Assembled dataset already exists: {args.output_file}")
        return

    # Assemble (load) the dataset
    dataset_label = DATASET_LABELS[args.mode]
    logger.info(f"Assembling {dataset_label} dataset (mode={args.mode})")

    if args.mode == "train":
        assembled_adata = assemble_dilimap_train_dataset(
            data_root=args.data_root,
        )
    else:
        assembled_adata = assemble_dilimap_train_val_dataset(
            data_root=args.data_root,
        )

    # Build file paths for this dataset run
    paths = get_dilimap_paths(args.data_root)
    logger.info(f"PubChem cache: {paths['pubchem_cache']}")

    # Standardize columns to the common pipeline schema and annotate PubChem CIDs
    logger.info(f"Standardizing {dataset_label} dataset")
    standardized_adata = standardize_dilimap(
        assembled_adata,
        dataset=dataset_label,
        paths=paths,
        annotate_pubchem=args.annotate_pubchem,
    )

    # Save to the pipeline-expected output location
    os.makedirs(os.path.dirname(args.output_file) or ".", exist_ok=True)
    logger.info(f"Saving standardized dataset to: {args.output_file}")
    standardized_adata.write_h5ad(args.output_file, compression="gzip")
    logger.info(f"{dataset_label} assembly complete")


if __name__ == "__main__":
    main()
