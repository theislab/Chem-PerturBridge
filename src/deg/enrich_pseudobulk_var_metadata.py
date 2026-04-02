"""
Left-join gene metadata (`symbol`, `is_merged`) onto DEG result AnnData objects.

Walks ``**/*_de.h5ad`` under a given root directory. Gene annotations are left-joined from a
reference processed pseudobulk ``.h5ad`` (paths via ``--deg_dir`` / ``--reference_h5ad``).

Use when DEG outputs lost `.var` columns (e.g. after `anndata.concat` with `merge='same'`).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, Tuple

import anndata as ad
import pandas as pd

from src.deg.aggregating_deg import rewrite_h5ad
from src.utils.parsing_utils import logger


REQUIRED_REF_COLS = ("symbol", "is_merged")


def needs_enrichment(var: pd.DataFrame) -> bool:
    """True if `symbol` or `is_merged` is missing or entirely null (needs join)."""
    if "symbol" not in var.columns or var["symbol"].isna().all():
        return True
    if "is_merged" not in var.columns:
        return True
    return False


def load_reference_var(reference_path: str) -> pd.DataFrame:
    if not os.path.isfile(reference_path):
        raise FileNotFoundError(
            f"Reference h5ad not found: {reference_path}. "
            "Run pseudobulk preprocessing first or pass --reference_h5ad."
        )
    ref = ad.read_h5ad(reference_path)
    missing = [c for c in REQUIRED_REF_COLS if c not in ref.var.columns]
    if missing:
        raise ValueError(
            f"Reference {reference_path} is missing .var columns: {missing}. "
            "Cannot enrich targets."
        )
    out = ref.var[list(REQUIRED_REF_COLS)].copy()
    out.index = out.index.astype(str)
    return out


def enrich_var_from_reference(adata: ad.AnnData, ref_var: pd.DataFrame) -> ad.AnnData:
    """Left-join reference ``symbol`` / ``is_merged`` onto ``adata.var`` by gene index."""
    adata = adata.copy()
    v = adata.var.copy()
    v.index = v.index.astype(str)

    ref = ref_var[list(REQUIRED_REF_COLS)].copy()
    ref.index = ref.index.astype(str)

    to_drop = [c for c in REQUIRED_REF_COLS if c in v.columns]
    if to_drop:
        v = v.drop(columns=to_drop)

    out = v.join(ref, how="left")
    out["symbol"] = out["symbol"].astype("category")

    adata.var = out
    adata.var.index = adata.var.index.astype(object)
    return adata


def iter_deg_h5ad_files(root: str) -> Iterable[str]:
    """All ``*_de.h5ad`` files under ``root`` (recursive), sorted."""
    root_p = Path(root)
    if not root_p.is_dir():
        return
    for path in sorted(root_p.rglob("*_de.h5ad")):
        yield str(path)


def enrich_tree(
    deg_dir: str,
    reference_path: str,
    *,
    dry_run: bool = False,
) -> Tuple[int, int, int]:
    """
    Walk ``**/*_de.h5ad`` under ``deg_dir``; enrich those that need it.

    Returns
    -------
    (n_scanned, n_need, n_written)
    """
    ref_var = load_reference_var(reference_path)

    n_scanned = 0
    n_need = 0
    n_written = 0

    for path in iter_deg_h5ad_files(deg_dir):
        n_scanned += 1
        try:
            adata = ad.read_h5ad(path)
        except Exception as e:
            logger.warning("Skip (read failed) %s: %s", path, e)
            continue
        if not needs_enrichment(adata.var):
            continue
        n_need += 1
        logger.info("Enrich %s", path)
        if dry_run:
            continue
        adata = enrich_var_from_reference(adata, ref_var)
        # Same as run_deg.R / aggregating_deg: drop X before save so read_h5ad in
        # aggregating_deg (and older anndata) does not hit IOSpec encoding_type='null'.
        if adata.X is not None:
            adata.X = None
        adata.write_h5ad(path, compression="gzip")
        rewrite_h5ad(path)
        n_written += 1

    return n_scanned, n_need, n_written


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Left-join symbol / is_merged from a reference processed pseudobulk .h5ad "
        "onto all *_de.h5ad files under --deg_dir (recursive)."
    )
    parser.add_argument(
        "--deg_dir",
        type=str,
        required=True,
        help="Directory tree to scan for **/*_de.h5ad (DEG pipeline output root).",
    )
    parser.add_argument(
        "--reference_h5ad",
        type=str,
        required=True,
        help="Processed pseudobulk AnnData; must contain .var columns symbol and is_merged.",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Log actions only; do not write files.",
    )
    args = parser.parse_args()

    deg_out = os.path.abspath(args.deg_dir)
    reference = os.path.abspath(args.reference_h5ad)

    logger.info("DEG directory (scan *_de.h5ad): %s", deg_out)
    logger.info("Reference pseudobulk: %s", reference)

    if not os.path.isdir(deg_out):
        raise SystemExit(f"DEG directory does not exist: {deg_out}")

    scanned, need, written = enrich_tree(deg_out, reference, dry_run=args.dry_run)
    logger.info(
        "Done: scanned=%d, need_enrichment=%d, written=%d (dry_run=%s)",
        scanned,
        need,
        written,
        args.dry_run,
    )


if __name__ == "__main__":
    main()
