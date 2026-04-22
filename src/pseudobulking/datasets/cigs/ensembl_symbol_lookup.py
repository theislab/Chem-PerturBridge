"""
Gene symbol ↔ Ensembl ID lookup for the CIGS dataset.

CIGS count matrices ship with HGNC gene symbols as column names, so the main
entry point for this module is :func:`fetch_ensembl_ids_for_symbols`, which
maps symbols → Ensembl gene IDs via the Ensembl REST API. The reverse
mapping (:func:`fetch_symbols_for_ensembl_ids`) is kept for completeness.

Both functions cache results on disk so repeated runs are cheap.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Optional

import requests
from requests.exceptions import HTTPError

from src.utils.parsing_utils import logger


_ENSEMBL_SERVER = "https://rest.ensembl.org"
_BATCH_SIZE     = 1000
_N_RETRIES      = 10
_SLEEP_S        = 10


# ---------------------------------------------------------------------------
# Low-level batch helpers (one POST per call, with retries)
# ---------------------------------------------------------------------------

def _post_batch(
    endpoint: str,
    payload: dict,
    keys: List[str],
    server: str,
    n_retries: int,
    sleep_s: int,
    extract,
) -> Dict[str, Optional[str]]:
    """Generic POST helper with retry/backoff, shared by both endpoints."""
    url     = f"{server}{endpoint}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for attempt in range(n_retries):
        try:
            response = requests.post(url, headers=headers, json=payload, timeout=60)
            response.raise_for_status()
            data = response.json()
            return {key: extract(data.get(key)) for key in keys}
        except HTTPError as e:
            code = e.response.status_code if e.response is not None else None
            retry_after = (
                int(e.response.headers.get("Retry-After", sleep_s))
                if e.response is not None else sleep_s
            )
            logger.warning(
                f"Ensembl {endpoint} HTTP {code} (attempt {attempt + 1}/{n_retries}): {e}"
            )
            time.sleep(retry_after)
        except Exception as e:
            logger.warning(
                f"Ensembl {endpoint} error (attempt {attempt + 1}/{n_retries}): {e}"
            )
            if attempt < n_retries - 1:
                time.sleep(sleep_s)
            else:
                raise

    logger.error(f"Failed to query {endpoint} after {n_retries} attempts for {len(keys)} keys")
    return {key: None for key in keys}


def _lookup_ids_batch(ensembl_ids: List[str], **kw) -> Dict[str, Optional[str]]:
    """POST /lookup/id — Ensembl ID → display_name (gene symbol)."""
    return _post_batch(
        endpoint="/lookup/id",
        payload={"ids": ensembl_ids},
        keys=ensembl_ids,
        extract=lambda entry: entry.get("display_name") if isinstance(entry, dict) else None,
        **kw,
    )


def _lookup_symbols_batch(symbols: List[str], **kw) -> Dict[str, Optional[str]]:
    """POST /lookup/symbol/homo_sapiens — gene symbol → Ensembl gene ID."""
    return _post_batch(
        endpoint="/lookup/symbol/homo_sapiens",
        payload={"symbols": symbols},
        keys=symbols,
        extract=lambda entry: entry.get("id") if isinstance(entry, dict) else None,
        **kw,
    )


# ---------------------------------------------------------------------------
# Cache-aware public API
# ---------------------------------------------------------------------------

def _load_cache(cache_path: Optional[str]) -> Dict[str, Optional[str]]:
    if not cache_path or not Path(cache_path).is_file():
        return {}
    try:
        with open(cache_path) as f:
            cache = json.load(f)
        logger.info(f"Loaded Ensembl lookup cache: {len(cache):,} entries from {cache_path}")
        return cache
    except Exception as e:
        logger.warning(f"Could not load Ensembl cache from {cache_path}: {e}")
        return {}


def _save_cache(cache: Dict[str, Optional[str]], cache_path: Optional[str]) -> None:
    if not cache_path:
        return
    Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(cache, f)
    logger.info(f"Saved Ensembl lookup cache: {len(cache):,} entries to {cache_path}")


def fetch_symbols_for_ensembl_ids(
    ensembl_ids: List[str],
    cache_path: Optional[str] = None,
    batch_size: int = _BATCH_SIZE,
    server: str = _ENSEMBL_SERVER,
) -> Dict[str, Optional[str]]:
    """
    Convert Ensembl gene IDs to gene symbols via ``POST /lookup/id``.

    Version suffixes (``ENSG00000141510.12``) are stripped for the API and
    re-mapped back to the original IDs on return.
    """
    cache = _load_cache(cache_path)

    stripped_to_original: Dict[str, str] = {}
    for eid in ensembl_ids:
        stripped = eid.split(".")[0] if "." in eid else eid
        stripped_to_original.setdefault(stripped, eid)

    to_fetch = [s for s in stripped_to_original if s not in cache]
    logger.info(
        f"Ensembl ID→symbol lookup: {len(stripped_to_original):,} unique, "
        f"{len(stripped_to_original) - len(to_fetch):,} cached, {len(to_fetch):,} to fetch"
    )

    batches = [to_fetch[i:i + batch_size] for i in range(0, len(to_fetch), batch_size)]
    for i, batch in enumerate(batches, 1):
        logger.info(f"  Batch {i}/{len(batches)}: querying {len(batch):,} IDs")
        cache.update(_lookup_ids_batch(
            batch, server=server, n_retries=_N_RETRIES, sleep_s=_SLEEP_S,
        ))

    if to_fetch:
        _save_cache(cache, cache_path)

    return {original: cache.get(stripped) for stripped, original in stripped_to_original.items()}


def fetch_ensembl_ids_for_symbols(
    symbols: List[str],
    cache_path: Optional[str] = None,
    batch_size: int = _BATCH_SIZE,
    server: str = _ENSEMBL_SERVER,
) -> Dict[str, Optional[str]]:
    """
    Convert gene symbols to Ensembl gene IDs via
    ``POST /lookup/symbol/homo_sapiens``.

    Parameters
    ----------
    symbols:
        List of gene symbols (HGNC). Duplicates and case are preserved in the
        return mapping, but the API is queried for unique, stripped symbols.
    cache_path:
        Optional path to a JSON file for persistent caching. Already-cached
        symbols are not re-queried.
    batch_size:
        Number of symbols per API request (Ensembl limit is 1000).
    server:
        Ensembl REST server base URL.

    Returns
    -------
    dict mapping each input symbol to its Ensembl gene ID (or ``None`` when
    the API returns no result for that symbol).
    """
    cache = _load_cache(cache_path)

    unique_symbols = sorted({s for s in symbols if s})
    to_fetch = [s for s in unique_symbols if s not in cache]
    logger.info(
        f"Ensembl symbol→ID lookup: {len(unique_symbols):,} unique, "
        f"{len(unique_symbols) - len(to_fetch):,} cached, {len(to_fetch):,} to fetch"
    )

    batches = [to_fetch[i:i + batch_size] for i in range(0, len(to_fetch), batch_size)]
    for i, batch in enumerate(batches, 1):
        logger.info(f"  Batch {i}/{len(batches)}: querying {len(batch):,} symbols")
        cache.update(_lookup_symbols_batch(
            batch, server=server, n_retries=_N_RETRIES, sleep_s=_SLEEP_S,
        ))

    if to_fetch:
        _save_cache(cache, cache_path)

    return {s: cache.get(s) for s in symbols}
