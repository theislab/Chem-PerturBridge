import json
import os
import re
import time
from typing import Any, Dict, Iterable, Optional

import requests

import pandas as pd

from src.utils.parsing_utils import *

def ols_search_clo(query: str, rows: int = 100, start: int = 0):
    url = "https://www.ebi.ac.uk/ols4/api/search"
    params = {
        "q": query,
        "ontology": "clo",
        "rows": rows,
        "start": start,
        "isDefiningOntology": "true",  # matches the UI param you used earlier
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    return r.json().get("response", {}).get("docs", [])

def cellosaurus_search_by_clo(clo_id: str, size: int = 10):
    url = "https://api.cellosaurus.org/search/cell-line"
    # 'dr' = cross-references; searching for "CLO;CLO_...." is a common pattern
    q = f'dr:"CLO;{clo_id}"'
    params = {"q": q, "format": "json", "size": size}
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    return r.json()


def _normalize_cell_label(label: str) -> str:
    if not label:
        return ""
    label = label.lower().removesuffix(" cell").strip(" ").replace("-", "").replace(" ", "")
    return label


def _select_best_clo_doc(docs: list, query: str) -> Optional[Dict[str, Any]]:
    if not docs:
        return None
    query_norm = _normalize_cell_label(query)
    exact = []
    for doc in docs:
        label = doc.get("label", "")
        if _normalize_cell_label(label) == query_norm:
            exact.append(doc)
    return exact[0] if exact else None


def _extract_cellosaurus_accession(data: Dict[str, Any]) -> Optional[str]:
    try:
        cell_line_list = data["Cellosaurus"]["cell-line-list"]
        return cell_line_list[0]["accession-list"][0]["value"]
    except Exception:
        return None


def _call_with_retry(func, *args, n_retries: int = 3, sleep_s: int = 2, **kwargs):
    last_exc = None
    for attempt in range(n_retries):
        try:
            return func(*args, **kwargs)
        except Exception as exc:
            last_exc = exc
            if attempt + 1 < n_retries:
                time.sleep(sleep_s)
    return None


def _lookup_clo_id(cell_line: str, rows: int = 100, n_retries: int = 3) -> Optional[str]:
    """
    Look up CLO ID from OLS API.
    
    Returns:
    --------
    Optional[str]
        CLO ID (e.g., "CLO_0001272") or None if not found
    """
    docs = _call_with_retry(ols_search_clo, cell_line, rows=rows, n_retries=n_retries, sleep_s=2)
    doc = _select_best_clo_doc(docs or [], cell_line)
    clo_id = doc.get("short_form") if doc else None
    return clo_id


def _lookup_cellosaurus_id(clo_id: str, size: int = 10, n_retries: int = 3) -> Optional[str]:
    """
    Look up Cellosaurus ID from CLO ID.
    
    Parameters:
    -----------
    clo_id : str
        CLO identifier (e.g., "CLO_0001272")
    
    Returns:
    --------
    Optional[str]
        Cellosaurus ID (e.g., "CVCL_0132") or None if not found
    """
    cellosaurus_data = _call_with_retry(
        cellosaurus_search_by_clo, clo_id, size=size, n_retries=n_retries, sleep_s=2
    )
    return _extract_cellosaurus_accession(cellosaurus_data) if cellosaurus_data else None


def _search_cellosaurus_by_name(cell_line: str, size: int = 10, n_retries: int = 3) -> Optional[str]:
    """
    Search Cellosaurus directly by cell line name.
    
    Parameters:
    -----------
    cell_line : str
        Cell line name (e.g., "A375")
    
    Returns:
    --------
    Optional[str]
        Cellosaurus ID (e.g., "CVCL_0132") or None if not found or ambiguous
    """
    def _cellosaurus_name_search(name: str, size: int) -> Dict[str, Any]:
        url = "https://api.cellosaurus.org/search/cell-line"
        params = {"q": name, "format": "json", "size": size}
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        return r.json()
    
    data = _call_with_retry(_cellosaurus_name_search, cell_line, size, n_retries=n_retries, sleep_s=2)
    if not data:
        return None
    
    # Only extract if there's exactly one result (unambiguous match)
    try:
        cell_line_list = data.get("Cellosaurus", {}).get("cell-line-list", [])
        if len(cell_line_list) == 1:
            return _extract_cellosaurus_accession(data)
    except Exception:
        pass
    
    return None


def annotate_cell_lines(
    df: pd.DataFrame,
    cell_id_col: str = "cell_id",
    base_cell_id_col: str = "base_cell_id",
    rows: int = 100,
    size: int = 10,
    n_retries: int = 3,
    request_delay_s: float = 0.2,
    manual_map: Optional[Dict[str, str]] = None,
    progress: bool = True,
) -> "pd.DataFrame":
    """
    Annotate cell lines with CLO and Cellosaurus IDs.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Cell line metadata dataframe
    cell_id_col : str, default="cell_id"
        Column containing cell line identifiers for filtering
    base_cell_id_col : str, default="base_cell_id"
        Column containing base cell line names to annotate
    rows : int, default=100
        Number of rows to fetch from OLS API
    size : int, default=10
        Number of results to fetch from Cellosaurus API
    n_retries : int, default=3
        Number of retries for API calls
    request_delay_s : float, default=0.2
        Delay in seconds between OLS/Cellosaurus API requests (fair-use throttling).
    manual_map : Optional[Dict[str, str]], default=None
        Manual mapping of cell line names to Cellosaurus IDs
    progress : bool, default=True
        Whether to log progress information
    
    Strategy:
    1. Filter rows where cell_id == base_cell_id
    2. Look up CLO via OLS API → Look up Cellosaurus via CLO ID
    3. Direct Cellosaurus search by cell line name (fallback)
    4. Use manual mapping if available (last resort)
    
    Returns:
    --------
    pd.DataFrame
        Original dataframe with added columns: clo_id, cellosaurus_id
    """
    if manual_map is None:
        manual_map = globals().get("NAME_TO_CVCL", {})

    # Extract unique cell lines only where cell_id == base_cell_id
    filtered_df = df[df[cell_id_col] == df[base_cell_id_col]]
    cell_lines = filtered_df[base_cell_id_col].dropna().unique()
    n_cell_lines = len(cell_lines)
    
    if progress:
        logger.info(f"Annotating {n_cell_lines} unique cell types")

    results = []
    empty_record = {"clo_id": None, "cellosaurus_id": None}
    iteration_count = 0
    
    for cell_line in cell_lines:
        # Handle missing/empty
        if pd.isna(cell_line) or not str(cell_line).strip():
            results.append({"cell_line": cell_line, **empty_record})
            continue
        
        cell_line = str(cell_line).strip()
        clo_id, cellosaurus_id = None, None
        
        # Strategy 1: Look up via OLS → CLO → Cellosaurus
        clo_id = _lookup_clo_id(cell_line, rows=rows, n_retries=n_retries)
        if clo_id:
            cellosaurus_id = _lookup_cellosaurus_id(clo_id, size=size, n_retries=n_retries)
        
        # Strategy 2: Direct Cellosaurus search by name (fallback)
        if not cellosaurus_id:
            cellosaurus_id = _search_cellosaurus_by_name(cell_line, size=size, n_retries=n_retries)
        
        # Strategy 3: Use manual mapping if still not found
        if not cellosaurus_id and manual_map and cell_line in manual_map:
            cellosaurus_id = manual_map[cell_line]
        
        # Build result
        results.append({"cell_line": cell_line, "cellosaurus_id": cellosaurus_id})
        
        # Progress logging
        iteration_count += 1
        if progress and iteration_count % 10 == 0:
            n_mapped_so_far = sum(1 for r in results if r["cellosaurus_id"] is not None)
            logger.info(f"Processed {iteration_count}/{n_cell_lines} cell types ({n_mapped_so_far} mapped so far)")
        
        if request_delay_s > 0:
            time.sleep(request_delay_s)

    # Create annotations dataframe
    annotations_df = pd.DataFrame(results)
    
    # Final progress summary
    if progress:
        n_mapped = sum(1 for r in results if r["cellosaurus_id"] is not None)
        logger.info(f"Completed: {n_mapped}/{n_cell_lines} cell types successfully mapped to Cellosaurus IDs")
    
    # Merge annotations back into original dataframe
    result_df = df.merge(
        annotations_df,
        left_on=cell_id_col,
        right_on="cell_line",
        how="left"
    )
    
    # Drop the temporary cell_line column from merge
    result_df = result_df.drop(columns=["cell_line"])
    
    return result_df
