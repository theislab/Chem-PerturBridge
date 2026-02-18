import requests
import time
import json
import os
import pandas as pd
from typing import Dict, Optional, List, Union, Callable

from src.utils.parsing_utils import *


class SymbolToEnsemblClient:
    """
    Client for mapping gene symbols to Ensembl gene IDs using HGNC API.
    Supports multiple lookup strategies:
    1. Current symbol lookup
    2. Previous symbol lookup
    
    Includes rate limiting (1 requests/second) and caching.
    """
    
    HGNC_BASE = "https://rest.genenames.org"
    
    def __init__(self, cache_path: Optional[str] = None, requests_per_second: float = 1.0, n_retries: int = 5):
        """
        Initialize gene symbol to Ensembl mapping client.
        
        Parameters
        ----------
        cache_path : Optional[str]
            Path to JSON file for persistent cache storage
        requests_per_second : float
            Maximum number of requests per second (default: 1.0; HGNC limit is 1/s
            and we use up to 2 requests per symbol)
        """
        self.cache: Dict[str, Optional[str]] = {}
        self.cache_path = cache_path
        self.min_interval = 1.0 / requests_per_second
        self.last_request_time = 0.0
        self.n_retries = n_retries
        
        # Load cache if path provided
        if cache_path:
            self._load_cache()
    
    def _load_cache(self) -> None:
        """Load cache from JSON file."""
        if self.cache_path and os.path.exists(self.cache_path):
            try:
                with open(self.cache_path, 'r') as f:
                    self.cache = json.load(f)
                logger.info(f"Loaded symbol-to-Ensembl cache from {self.cache_path} with {len(self.cache)} entries")
            except Exception as e:
                logger.warning(f"Failed to load cache from {self.cache_path}: {e}")
                self.cache = {}
    
    def _save_cache(self) -> None:
        """Save cache to JSON file."""
        if self.cache_path:
            try:
                cache_dir = os.path.dirname(self.cache_path)
                if cache_dir:
                    os.makedirs(cache_dir, exist_ok=True)
                with open(self.cache_path, 'w') as f:
                    json.dump(self.cache, f, indent=2)
                logger.debug(f"Saved symbol-to-Ensembl cache to {self.cache_path} with {len(self.cache)} entries")
            except Exception as e:
                logger.warning(f"Failed to save symbol-to-Ensembl cache to {self.cache_path}: {e}")
    
    def _rate_limit(self) -> None:
        """Enforce rate limiting by sleeping if necessary."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < self.min_interval:
            sleep_time = self.min_interval - time_since_last
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
    
    def symbol_to_ensg(self, symbol: Union[str, int], mode: str = 'symbol') -> Optional[str]:
        """
        Convert gene symbol to Ensembl gene ID using HGNC API.
        
        Parameters
        ----------
        symbol : str | int
            Gene symbol (will be converted to string)
        mode : str, default='symbol'
            HGNC API lookup mode. Must be one of:
            - 'symbol': search by current gene symbol
            - 'prev_symbol': search by previous/alias gene symbol
            
        Returns
        -------
        Optional[str]
            Ensembl gene ID or None if not found
        """
        if not mode in ['symbol', 'prev_symbol']:
            raise ValueError(f"Invalid mode: {mode}. Must be 'symbol' or 'prev_symbol'.")

        symbol_str = str(symbol).strip()
        
        # Check cache (only successful lookups are stored, so a cache hit means we already found the Ensembl ID)
        if symbol_str in self.cache:
            return self.cache[symbol_str]
        
        
        
        # Make HGNC API request, retrying on transient errors (rate limit, timeout, server errors)
        url = f"{self.HGNC_BASE}/fetch/{mode}/{symbol_str}"
        cnt = 0
        result = None
        while cnt < self.n_retries:
            try:
                # Rate limit
                self._rate_limit()
                r = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
                r.raise_for_status()

                docs = r.json()["response"]["docs"]
                if not docs:
                    result = None
                else:
                    # Usually exactly one HGNC record for a gene symbol
                    if len(docs) > 1:
                        logger.warning(f"Multiple HGNC records found for symbol {symbol_str}. Returning None.")
                        result = None
                    else:
                        result = docs[0].get("ensembl_gene_id")
                break

            except requests.exceptions.HTTPError as e:
                if r.status_code in (429, 500, 502, 503, 504):
                    logger.warning(f"HGNC API failed for symbol {symbol_str}: {e}. Retry {cnt}.")
                    cnt += 1
                    time.sleep(5)
                else:
                    logger.warning(f"HGNC API failed for symbol {symbol_str}: {e}")
                    break
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                logger.warning(f"HGNC API failed for symbol {symbol_str}: {e}. Retry {cnt}.")
                cnt += 1
                time.sleep(5)
            except Exception as e:
                logger.warning(f"HGNC API failed for symbol {symbol_str}: {e}")
                break

        # Only cache successful lookups so that a None from symbol mode
        # does not block the subsequent prev_symbol lookup for the same symbol
        if result is not None:
            self.cache[symbol_str] = result
            if len(self.cache) % 100 == 0 and self.cache_path:
                self._save_cache()

        return result
    
    def save_cache(self) -> None:
        """Manually save cache to disk."""
        self._save_cache()


def fetch_ensg_ids_from_symbols(
    symbols: List[Union[str, int]], 
    cache_path: Optional[str] = None,
    requests_per_second: float = 5.0,
    manual_mapping_func: Optional[Callable[[], Dict]] = None,
    dataset_key: Optional[str] = None
) -> List[Optional[str]]:
    """
    Fetch Ensembl gene IDs for a list of gene symbols using HGNC API.
    
    Uses multiple strategies in order:
    1. Manual mapping (if provided)
    2. Current symbol lookup via HGNC
    3. Previous symbol lookup via HGNC
    
    Parameters
    ----------
    symbols : List[str | int]
        List of gene symbols to map
    cache_path : str, optional
        Path to JSON file for persistent cache storage. If None, no persistent
        cache is used. Callers should pass a path relative to project/data root
        (e.g. from paths["ensembl_cache"]) so all runs and scripts share
        the same cache.
    requests_per_second : float
        Maximum number of requests per second (default: 5.0)
    manual_mapping_func : Optional[Callable[[], Dict]], default=None
        Function that returns manual Ensembl ID mappings. If provided, should return
        either a dict directly (e.g., {'symbol': 'ENSG...'}) or a dict with dataset keys
        (e.g., {'dataset_name': {'symbol': 'ENSG...'}}). If None, no manual mapping is used.
    dataset_key : Optional[str], default=None
        Key to extract from the dict returned by manual_mapping_func if it returns
        a nested dict structure. If None and manual_mapping_func returns a nested dict,
        uses the first dataset in that mapping.
        
    Returns
    -------
    List[Optional[str]]
        List of Ensembl gene IDs (or None if not found) in same order as input
    """
    ensg_ids = []
    iteration_count = 0
    n_symbols = len(symbols)
    
    # Get manual mappings
    if manual_mapping_func is not None:
        mapping_result = manual_mapping_func()
        # Handle both flat dicts and nested dicts
        if isinstance(mapping_result, dict):
            if dataset_key is not None:
                # Extract specific key from nested dict
                manual_mapping = mapping_result.get(dataset_key, {})
            elif len(mapping_result) == 1:
                # If single key, use it automatically
                only_value = list(mapping_result.values())[0]
                manual_mapping = only_value if isinstance(only_value, dict) else mapping_result
            else:
                # Check if it's a nested dict (values are dicts) or flat dict (values are strings)
                first_value = list(mapping_result.values())[0] if mapping_result else None
                if isinstance(first_value, dict):
                    # Nested dict but no dataset_key specified - use first value as fallback
                    manual_mapping = first_value
                else:
                    # Flat dict with symbol: ensembl_id mappings
                    manual_mapping = mapping_result
        else:
            manual_mapping = {}
    else:
        manual_mapping = {}
    
    client = SymbolToEnsemblClient(cache_path=cache_path, 
                                   requests_per_second=requests_per_second)
    
    for symbol in symbols:
        iteration_count += 1
        ensg_id = None
        # Strategy 1: Manual mapping
        if manual_mapping and symbol in manual_mapping:
            ensg_id = manual_mapping[symbol]
        
        # Strategy 2: Current symbol lookup
        if not ensg_id:
            ensg_id = client.symbol_to_ensg(symbol, mode='symbol')
        
        # Strategy 3: Previous symbol lookup
        if not ensg_id:
            ensg_id = client.symbol_to_ensg(symbol, mode='prev_symbol')
            
        ensg_ids.append(ensg_id)

        if iteration_count % 50 == 0:
            n_mapped_so_far = sum(1 for i in ensg_ids if i is not None)
            logger.info(f"Processed {iteration_count}/{n_symbols} gene symbols ({n_mapped_so_far} mapped so far)")

    client.save_cache()
    n_mapped_so_far = sum(1 for i in ensg_ids if i is not None)
    logger.info(f"Processed {iteration_count}/{n_symbols} gene symbols ({n_mapped_so_far} mapped so far)")
    return ensg_ids
