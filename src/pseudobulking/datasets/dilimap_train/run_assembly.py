"""
DILImap training dataset assembly script.

Loads the DILImap training (or training+validation) H5AD file(s), applies
standardization, and saves to the pipeline-expected output path.
"""
import os
import argparse

import numpy as np

from src.utils.parsing_utils import logger
from src.pseudobulking.common.pubchem import lookup_pubchem_cids
from src.pseudobulking.datasets.dilimap_train.assembling import (
    assemble_dilimap_train_dataset,
    assemble_dilimap_train_val_dataset,
)
from src.pseudobulking.datasets.dilimap_train.standardization import (
    standardize_dilimap_train,
    standardize_dilimap_train_val,
)
from src.pseudobulking.datasets.dilimap_train.pubchem_imputation import (
    pubchem_mapping_dilimap_train,
    pubchem_mapping_dilimap_train_val,
)


SUBSAMPLE_SIZE = 500

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
    --data_root : str, default='./op3_v2/data'
        Root directory containing the DILImap H5AD source file(s).
    --subsampling : flag
        If set, randomly subsample to a smaller dataset for debugging.
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
        default="./op3_v2/data",
        help="Root directory containing the DILImap H5AD file(s)",
    )
    parser.add_argument(
        "--subsampling",
        type=bool,
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Subsample the dataset for debugging",
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

    # Subsample if requested
    if args.subsampling:
        n_obs = assembled_adata.n_obs
        size = min(SUBSAMPLE_SIZE, n_obs)
        logger.info(f"Subsampling: {n_obs} -> {size} observations")
        idx = np.random.choice(n_obs, replace=False, size=size)
        assembled_adata = assembled_adata[idx].copy()

    # Standardize columns to the common pipeline schema
    logger.info(f"Standardizing {dataset_label} dataset")
    if args.mode == "train":
        assembled_adata = standardize_dilimap_train(assembled_adata)
        manual_mapping_func = pubchem_mapping_dilimap_train
    else:
        assembled_adata = standardize_dilimap_train_val(assembled_adata)
        manual_mapping_func = pubchem_mapping_dilimap_train_val

    # Enrich pubchem_cid via automatic PubChem lookup + manual fallback mappings
    logger.info(f"Looking up PubChem CIDs for {dataset_label}")
    assembled_adata.obs = lookup_pubchem_cids(
        assembled_adata.obs,
        cache={},
        pert_id_col=None,
        drug_col="perturbagen",
        manual_mapping_func=manual_mapping_func,
        dataset_key=dataset_label,
    )

    # Save to the pipeline-expected output location
    os.makedirs(os.path.dirname(args.output_file) or ".", exist_ok=True)
    logger.info(f"Saving assembled dataset to: {args.output_file}")
    assembled_adata.write_h5ad(args.output_file, compression="gzip")
    logger.info(f"{dataset_label} assembly complete")


if __name__ == "__main__":
    main()
