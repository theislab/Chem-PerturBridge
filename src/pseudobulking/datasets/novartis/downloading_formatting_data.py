from __future__ import annotations

import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Union

import anndata as ad

from src.utils.parsing_utils import logger

RSCRIPT_PATH = Path(__file__).parent / "rdata_to_adata.R"


# ---------------------------------------------------------------------------
# Download all files
# ---------------------------------------------------------------------------

def download_novartis_data(
    dest_dir: Union[str, Path],
    files: list,
    force: bool = False,
) -> None:
    """
    Download raw Novartis DRUG-seq files.

    Each file is skipped when it already exists, unless force=True.

    Parameters
    ----------
    dest_dir :
        Directory where all files will be saved.
    files :
        List of (filename, url) pairs to download, as defined in NOVARTIS_FILES
        in run_standardization.py.
    force :
        Re-download files even if they already exist.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    for fname, url in files:
        dest = dest_dir / fname
        if dest.exists() and not force:
            logger.info("Already exists, skipping: %s", dest)
        else:
            logger.info("Downloading %s ...", fname)
            urllib.request.urlretrieve(url, dest)
            logger.info("  -> saved to %s", dest)


# ---------------------------------------------------------------------------
# Step 1 — decompress
# ---------------------------------------------------------------------------

def decompress_rdata(
    compressed_path: Union[str, Path],
    uncompressed_path: Union[str, Path],
    rscript: str = "Rscript",
    force: bool = False,
) -> None:
    """
    Load a gzip-compressed .RData file in R and re-save it uncompressed.

    The downloaded Exp_gzip.RData is gzip-compressed, which makes repeated
    loads very slow.  This step converts it to an uncompressed .RData so
    that rdata_to_adata.R can load it quickly.

    Parameters
    ----------
    compressed_path :
        Path to the downloaded Exp_gzip.RData.
    uncompressed_path :
        Destination path for the uncompressed Exp.RData.
    rscript :
        Name or full path of the Rscript executable.
    force :
        Re-run even if the uncompressed file already exists.
    """
    compressed_path   = Path(compressed_path)
    uncompressed_path = Path(uncompressed_path)
    uncompressed_path.parent.mkdir(parents=True, exist_ok=True)

    if uncompressed_path.exists() and not force:
        logger.info(
            "Uncompressed RData already exists at %s — skipping decompression.",
            uncompressed_path,
        )
        return

    _check_rscript(rscript)

    r_code = (
        f"load('{compressed_path}'); "
        f"save(list = ls(), file = '{uncompressed_path}', compress = FALSE)"
    )

    logger.info(
        "Decompressing %s -> %s (this may take a few minutes) ...",
        compressed_path, uncompressed_path,
    )
    subprocess.run(
        [rscript, "-e", r_code],
        check=True,
        text=True,
    )
    logger.info("Decompression complete.")


# ---------------------------------------------------------------------------
# Step 2 — convert RData -> h5ad
# ---------------------------------------------------------------------------

def rdata_to_adata(
    input_rdata: Union[str, Path],
    output_h5ad: Union[str, Path],
    annotation_rdata: Union[str, Path],
    genes_csv: Union[str, Path],
    rscript: str = "Rscript",
    force: bool = False,
) -> ad.AnnData:
    """
    Run rdata_to_adata.R to convert the (uncompressed) Exp RData to .h5ad
    and export the drugseq_ensg_v98 gene annotation table to CSV.

    Parameters
    ----------
    input_rdata :
        Path to the uncompressed .RData file produced by decompress_rdata().
    output_h5ad :
        Path where the .h5ad will be written.
    annotation_rdata :
        Path to drugseq_ensembl_v98_annotation_and_entrez_mapping.RData.
    genes_csv :
        Path where the drugseq_ensg_v98 gene annotation CSV will be saved.
    rscript :
        Name or full path of the Rscript executable.
    force :
        Re-run the R conversion even if output_h5ad already exists.
    """
    input_rdata      = Path(input_rdata)
    output_h5ad      = Path(output_h5ad)
    annotation_rdata = Path(annotation_rdata)
    genes_csv        = Path(genes_csv)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    if not output_h5ad.exists() or force:
        _check_rscript(rscript)
        logger.info("Running rdata_to_adata.R ...")
        subprocess.run(
            [
                rscript, str(RSCRIPT_PATH),
                "--input_file",      str(input_rdata),
                "--output_file",     str(output_h5ad),
                "--annotation_file", str(annotation_rdata),
                "--genes_csv",       str(genes_csv),
            ],
            check=True,
            text=True,
        )
        logger.info("R conversion finished.")
    else:
        logger.info(
            "h5ad already exists at %s — skipping R conversion.", output_h5ad
        )

    logger.info("Loading %s", output_h5ad)
    return ad.read_h5ad(output_h5ad)


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def load_novartis(
    paths: dict,
    files: list,
    rscript: str = "Rscript",
    force_download: bool = False,
    force_decompress: bool = False,
    force_convert: bool = False,
) -> ad.AnnData:
    """
    Full pipeline: download all files → decompress → convert → load.

    Each step is skipped automatically when its output already exists,
    unless the corresponding force_* flag is set.

    Parameters
    ----------
    paths :
        Dictionary of file paths as returned by get_novartis_paths(). Expected keys:
        gzip_rdata, rdata, annotation_rdata, raw_h5ad, genes_csv.
    files :
        List of (filename, url) pairs to download, as defined in NOVARTIS_FILES
        in run_standardization.py.
    rscript :
        Name or full path of the Rscript executable.
    force_download :
        Re-download all files even if they already exist.
    force_decompress :
        Re-decompress even if Exp.RData already exists.
    force_convert :
        Re-run R conversion even if the .h5ad already exists.
    """
    raw_dir = Path(paths["raw_h5ad"]).parent

    download_novartis_data(raw_dir, files=files, force=force_download)

    decompress_rdata(
        paths["gzip_rdata"], paths["rdata"],
        rscript=rscript, force=force_decompress,
    )

    return rdata_to_adata(
        paths["rdata"], paths["raw_h5ad"],
        annotation_rdata=paths["annotation_rdata"],
        genes_csv=paths["genes_csv"],
        rscript=rscript, force=force_convert,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _check_rscript(rscript: str) -> None:
    if shutil.which(rscript) is None:
        raise RuntimeError(
            f"'{rscript}' not found on PATH. "
            "Make sure R is installed and the conda environment is active."
        )
