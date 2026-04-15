"""
Ensembl ID → gene symbol lookup for VCPI Ginkgo datasets.

Uses the Ensembl REST API POST /lookup/id endpoint to convert
Ensembl gene IDs (e.g. ENSG00000141510) to display names (gene symbols,
e.g. TP53). Results are cached to a JSON file to avoid redundant API calls.

Usage
-----
    from src.pseudobulking.datasets.vcpi_ginkgo.ensembl_symbol_lookup import (
        fetch_symbols_for_ensembl_ids,
    )

    symbol_map = fetch_symbols_for_ensembl_ids(
        ensembl_ids=list(adata.var_names),
        cache_path="data/vcpi_ginkgo/raw/ensembl_symbol_cache.json",
    )
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


def _lookup_batch(
    ensembl_ids: List[str],
    server: str = _ENSEMBL_SERVER,
    n_retries: int = _N_RETRIES,
    sleep_s: int = _SLEEP_S,
) -> Dict[str, Optional[str]]:
    """
    POST /lookup/id for a single batch (≤ 1000 IDs).

    Returns a dict mapping each queried Ensembl ID to its ``display_name``
    (gene symbol), or ``None`` when the API returns no entry for that ID.
    """
    url     = f"{server}/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {"ids": ensembl_ids}

    for attempt in range(n_retries):
        try:
            response = requests.post(url, headers=headers, json=payload, timeout=60)
            response.raise_for_status()
            data = response.json()
            return {
                eid: (data[eid].get("display_name") if isinstance(data.get(eid), dict) else None)
                for eid in ensembl_ids
            }
        except HTTPError as e:
            code = e.response.status_code if e.response is not None else None
            retry_after = int(e.response.headers.get("Retry-After", sleep_s)) if e.response is not None else sleep_s
            logger.warning(
                "Ensembl /lookup/id HTTP %s (attempt %d/%d): %s",
                code, attempt + 1, n_retries, e,
            )
            time.sleep(retry_after)
        except Exception as e:
            logger.warning(
                "Ensembl /lookup/id error (attempt %d/%d): %s",
                attempt + 1, n_retries, e,
            )
            if attempt < n_retries - 1:
                time.sleep(sleep_s)
            else:
                raise

    logger.error("Failed to fetch symbols after %d attempts for %d IDs", n_retries, len(ensembl_ids))
    return {eid: None for eid in ensembl_ids}


def fetch_symbols_for_ensembl_ids(
    ensembl_ids: List[str],
    cache_path: Optional[str] = None,
    batch_size: int = _BATCH_SIZE,
    server: str = _ENSEMBL_SERVER,
) -> Dict[str, Optional[str]]:
    """
    Convert Ensembl gene IDs to gene symbols via the Ensembl REST API.

    Queries ``POST /lookup/id`` in batches of ``batch_size`` (max 1000).
    The ``display_name`` field of each response entry is used as the symbol.

    Parameters
    ----------
    ensembl_ids:
        List of Ensembl gene IDs (version suffixes like ``.12`` are stripped
        before querying and re-mapped back to the original IDs).
    cache_path:
        Optional path to a JSON file for persistent caching. Already-cached
        IDs are not re-queried.
    batch_size:
        Number of IDs per API request (Ensembl limit is 1000).
    server:
        Ensembl REST server base URL.

    Returns
    -------
    dict mapping each input Ensembl ID to its symbol string (or ``None``
    when the API returns no result for that ID).
    """
    # Load existing cache
    cache: Dict[str, Optional[str]] = {}
    if cache_path and Path(cache_path).is_file():
        try:
            with open(cache_path) as f:
                cache = json.load(f)
            logger.info("Loaded Ensembl symbol cache: %d entries from %s", len(cache), cache_path)
        except Exception as e:
            logger.warning("Could not load symbol cache from %s: %s", cache_path, e)

    # Strip version suffixes (ENSG00000141510.12 → ENSG00000141510) for the API
    stripped_to_original: Dict[str, str] = {}
    for eid in ensembl_ids:
        stripped = eid.split(".")[0] if "." in eid else eid
        stripped_to_original.setdefault(stripped, eid)

    to_fetch = [s for s in stripped_to_original if s not in cache]
    logger.info(
        "Ensembl symbol lookup: %d unique IDs total, %d already cached, %d to fetch",
        len(stripped_to_original), len(stripped_to_original) - len(to_fetch), len(to_fetch),
    )

    # Fetch in batches
    batches = [to_fetch[i: i + batch_size] for i in range(0, len(to_fetch), batch_size)]
    for i, batch in enumerate(batches, 1):
        logger.info("  Batch %d/%d: querying %d IDs", i, len(batches), len(batch))
        result = _lookup_batch(batch, server=server)
        cache.update(result)

    # Persist updated cache
    if cache_path and to_fetch:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(cache, f)
        logger.info("Saved Ensembl symbol cache: %d entries to %s", len(cache), cache_path)

    # Map back to original IDs (with version suffix)
    return {
        original: cache.get(stripped)
        for stripped, original in stripped_to_original.items()
        for original in [stripped_to_original[stripped]]
    }
