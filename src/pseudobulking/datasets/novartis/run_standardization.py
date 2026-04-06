"""
Novartis MoABox DRUG-seq dataset standardization script.

This script standardizes the raw Novartis AnnData to the unified pseudobulk schema.
"""
from __future__ import annotations

import os
import argparse
from pathlib import Path

from src.utils.parsing_utils import logger
from src.pseudobulking.datasets.novartis.downloading_formatting_data import load_novartis
from src.pseudobulking.datasets.novartis.standardization import standardize_novartis_dataset


_ZENODO = "https://zenodo.org/records/14291446/files"
_GITHUB = (
    "https://raw.githubusercontent.com/Novartis/DRUG-seq/main"
    "/data/Novartis_drugseq_U2OS_MoABox"
)

# (filename, url) pairs for all raw files that need to be downloaded
NOVARTIS_FILES = [
    ("MoABox_compounds_metadata.txt",
     f"{_GITHUB}/MoABox_compounds_metadata.txt"),
    ("robust_RC_ReferenceControl_DMSO_wells.txt",
     f"{_ZENODO}/robust_RC_ReferenceControl_DMSO_wells.txt?download=1"),
    ("drugseq_ensembl_v98_annotation_and_entrez_mapping.RData",
     f"{_GITHUB}/drugseq_ensembl_v98_annotation_and_entrez_mapping.RData"),
    ("Exp_gzip.RData",
     f"{_ZENODO}/Exp_gzip.RData?download=1"),
]


def get_novartis_paths(data_root: str | Path) -> dict:
    """
    Build a dictionary of all file paths required by the Novartis pipeline.

    Parameters
    ----------
    data_root :
        Root directory for the Novartis dataset (e.g. ./data/novartis).

    Returns
    -------
    dict mapping file identifiers to Path objects.
    """
    raw_dir = Path(data_root) / "raw"
    return {
        # download / R-conversion intermediates
        "gzip_rdata":          raw_dir / "Exp_gzip.RData",
        "rdata":               raw_dir / "Exp.RData",
        "annotation_rdata":    raw_dir / "drugseq_ensembl_v98_annotation_and_entrez_mapping.RData",
        # standardization inputs
        "raw_h5ad":            raw_dir / "novartis_raw.h5ad",
        "genes_csv":           raw_dir / "drugseq_ensg_v98.csv",
        "compound_tsv":        raw_dir / "MoABox_compounds_metadata.txt",
        "robust_dmso":         raw_dir / "robust_RC_ReferenceControl_DMSO_wells.txt",
        # PubChem caches
        "pubchem_cid_cache":   raw_dir / "novartis_pubchem_cache.json",
        "pubchem_names_cache":  raw_dir / "novartis_pubchem_names_cache.json",
    }


def main():
    """
    Standardize the Novartis MoABox DRUG-seq dataset to the unified pseudobulk schema.

    This is the main entry point for the Novartis standardization pipeline. It:
    1. Parses command-line arguments for configuration
    2. Checks if output file already exists (exits early if present)
    3. Downloads raw files and converts RData -> h5ad (each step skipped if already done)
    4. Standardizes the dataset
    5. Saves the standardized dataset as an H5AD file

    Command-line Arguments
    ----------------------
    --output_file : str, required
        Output path for the standardized H5AD file.
    --data_root : str, default='./novartis_data'
        Root directory containing raw Novartis files
        (expects raw/novartis_raw.h5ad, raw/drugseq_ensg_v98.csv,
        raw/MoABox_compounds_metadata.txt,
        raw/robust_RC_ReferenceControl_DMSO_wells.txt).
    --annotate-pubchem : bool, default=False
        If set, annotates compounds with PubChem CIDs (may involve API calls).
    --annotate-pubchem-names : bool, default=False
        If set, additionally looks up drug synonym names by PubChem CID
        (requires --annotate-pubchem or pre-existing pubchem_cid column).
    """
    parser = argparse.ArgumentParser(
        description='Standardize Novartis MoABox DRUG-seq dataset to unified pseudobulk schema'
    )
    parser.add_argument('--output_file', type=str, required=True,
                        help='Output path for standardized H5AD file')
    parser.add_argument('--data_root', type=str, default='./novartis_data',
                        help='Root directory for raw Novartis files')
    parser.add_argument('--annotate-pubchem', action='store_true',
                        help='Annotate compounds with PubChem CIDs (may involve API calls)')
    parser.add_argument('--annotate-pubchem-names', action='store_true',
                        help='Look up drug synonym names by PubChem CID (requires --annotate-pubchem or pre-existing pubchem_cid column)')

    args = parser.parse_args()

    if os.path.isfile(args.output_file):
        logger.info("Standardized dataset already exists: %s", args.output_file)
        return

    paths = get_novartis_paths(args.data_root)

    logger.info("Downloading and converting Novartis raw data")
    load_novartis(paths, files=NOVARTIS_FILES)

    logger.info("Standardizing Novartis dataset")
    adata_standardized = standardize_novartis_dataset(
        paths=paths,
        annotate_pubchem=args.annotate_pubchem,
        annotate_pubchem_names=args.annotate_pubchem_names,
    )

    logger.info("Saving standardized dataset to: %s", args.output_file)
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    adata_standardized.write_h5ad(args.output_file, compression='gzip')
    logger.info("Novartis standardization complete")


if __name__ == "__main__":
    main()
