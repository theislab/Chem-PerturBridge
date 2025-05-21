#!/usr/bin/env python
"""Sequential pathway scoring for Tahoe DGE tables (tqdm-powered)."""

from __future__ import annotations

import os
import glob
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from gseapy import get_library, prerank

# ───────────────────────────── user-tweakable constants ──────────────────────────
TAHOE_DGE_DIR = Path("../../data/tahoe/dge/mDf2/")
OUTPUT_DIR    = Path("../../data/tahoe/results/pathway_calls/")

GENE_COL    = "gene"
TSTAT_COL   = "t"
COND_COL    = "condition"
FDR_THRESH  = 0.05
GMT_NAME    = "MSigDB_Hallmark_2020"
MIN_OVERLAP = 5
N_PERM      = 10_000
SEED        = 123
PNG_DPI     = 150
# ────────────────────────────────────────────────────────────────────────────────

# Cache gene-set collection once per run
_GMT_CACHE: dict[str, list[str]] | None = None
def _load_gmt() -> dict[str, list[str]]:
    global _GMT_CACHE
    if _GMT_CACHE is None:
        _GMT_CACHE = get_library(
            name=GMT_NAME,
            organism="Human",
            min_size=0,
            max_size=10_000,
        )
    return _GMT_CACHE

# ───────────────────────── single-contrast helper (unchanged) ────────────────────
def pathway_calls_for_one(
    df_subset: pd.DataFrame,
    gene_col: str = GENE_COL,
    t_col: str = TSTAT_COL,
    fdr_thresh: float = FDR_THRESH,
) -> pd.DataFrame:
    gmt = _load_gmt()

    rnk = (
        df_subset[[gene_col, t_col]]
        .dropna()
        .groupby(gene_col)[t_col]
        .mean()
        .sort_values(ascending=False)
        .reset_index()
    )

    if rnk[gene_col].nunique() < MIN_OVERLAP:
        return pd.DataFrame(columns=["pathway", "NES", "FDR", "call"])

    my_genes = set(rnk[gene_col])
    overlap_counts = {p: len(set(genes) & my_genes) for p, genes in gmt.items()}
    if max(overlap_counts.values(), default=0) < MIN_OVERLAP:
        return pd.DataFrame(columns=["pathway", "NES", "FDR", "call"])

    res = prerank(
        rnk=rnk,
        gene_sets=gmt,
        permutation_num=N_PERM,
        outdir=None,
        seed=SEED,
        min_size=MIN_OVERLAP,
        max_size=10_000,
    ).res2d

    res = (
        res.rename(columns={"Term": "pathway", "NES": "NES", "FDR q-val": "FDR"})
        .reset_index(drop=True)
    )

    res["call"] = np.where(
        res["FDR"] < fdr_thresh,
        np.where(res["NES"] > 0, "up", "down"),
        "no_change",
    )
    return res[["pathway", "NES", "FDR", "call"]]

# ───────────────────────── all-conditions helper ────────────────────────────────
def score_all_conditions(limma_df: pd.DataFrame) -> pd.DataFrame:
    results: list[pd.DataFrame] = []
    for cond, sub in limma_df.groupby(COND_COL, sort=False):
        calls = pathway_calls_for_one(sub)
        if not calls.empty:
            calls.insert(0, "condition", cond)
            results.append(calls)
    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()

# ───────────────────────── per-cell-line worker ─────────────────────────────────
def process_cell_line(csv_path: str | os.PathLike) -> str:
    csv_path = Path(csv_path)
    cell_line = csv_path.stem.replace("_differential_expression_results", "")
    out_dir = OUTPUT_DIR / cell_line
    out_dir.mkdir(parents=True, exist_ok=True)

    limma_df = pd.read_csv(csv_path)

    pathway_df = score_all_conditions(limma_df)
    if pathway_df.empty:
        return f"{cell_line}: no pathways -> skipped"

    pathway_csv = out_dir / "pathway_calls.csv"
    pathway_df.to_csv(pathway_csv, index=False)

    sig_calls = (
        pathway_df[pathway_df["call"] != "no_change"]
        .groupby("condition")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    plt.figure(figsize=(12, 6))
    plt.bar(sig_calls["condition"], sig_calls["count"])
    plt.xticks(rotation=90)
    plt.xlabel("Condition")
    plt.ylabel("# Significant Pathway Calls (up | down)")
    plt.title(f"{cell_line}: Significant Pathway Calls per Condition")
    plt.tight_layout()

    png_path = out_dir / "significant_calls.png"
    plt.savefig(png_path, dpi=PNG_DPI)
    plt.close()

    return f"{cell_line}: saved {pathway_csv.name} & {png_path.name}"

# ─────────────────────────────────── main ────────────────────────────────────────
def main():  # noqa: D401
    global TAHOE_DGE_DIR, OUTPUT_DIR

    parser = argparse.ArgumentParser(description="Sequential pathway scoring for Tahoe DGE data.")
    parser.add_argument("--dge_dir", type=Path, default=TAHOE_DGE_DIR,
                        help="Directory with *_differential_expression_results.csv files")
    parser.add_argument("--out_dir", type=Path, default=OUTPUT_DIR,
                        help="Output directory for pathway calls & plots")
    args = parser.parse_args()

    TAHOE_DGE_DIR = args.dge_dir
    OUTPUT_DIR = args.out_dir
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    csv_files = glob.glob(str(TAHOE_DGE_DIR / "*_differential_expression_results.csv"))
    if not csv_files:
        raise SystemExit(f"No CSVs found in {TAHOE_DGE_DIR}")

    print(f"Found {len(csv_files)} cell lines → crunching sequentially…")

    for msg in tqdm(map(process_cell_line, csv_files), total=len(csv_files), unit="cell line"):
        print(msg)

if __name__ == "__main__":
    main()