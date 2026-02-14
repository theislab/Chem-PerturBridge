import requests
import time
import json
import os
from typing import Dict, Optional, List, Union

from src.utils.parsing_utils import *


class SymbolToEnsemblClient:
    """
    Client for mapping gene symbols to Ensembl gene IDs using HGNC API.
    Includes rate limiting (10 requests/second) and caching.
    """
    
    HGNC_BASE = "https://rest.genenames.org"
    
    def __init__(self, cache_path: Optional[str] = None, requests_per_second: float = 10.0):
        """
        Initialize gene symbol to Ensembl mapping client.
        
        Parameters
        ----------
        cache_path : Optional[str]
            Path to JSON file for persistent cache storage
        requests_per_second : float
            Maximum number of requests per second (default: 10.0)
        """
        self.cache: Dict[str, Optional[str]] = {}
        self.cache_path = cache_path
        self.min_interval = 1.0 / requests_per_second
        self.last_request_time = 0.0
        
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
    
    def symbol_to_ensg(self, symbol: Union[str, int]) -> Optional[str]:
        """
        Convert gene symbol to Ensembl gene ID using HGNC API.
        
        Parameters
        ----------
        symbol : str | int
            Gene symbol (will be converted to string)
            
        Returns
        -------
        Optional[str]
            Ensembl gene ID or None if not found
        """
        symbol_str = str(symbol).strip()
        
        # Check cache first
        if symbol_str in self.cache:
            return self.cache[symbol_str]
        
        # Rate limit
        self._rate_limit()
        
        # Make HGNC API request
        url = f"{self.HGNC_BASE}/fetch/symbol/{symbol_str}"
        try:
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
            r.raise_for_status()
            
            docs = r.json()["response"]["docs"]
            if not docs:
                result = None
            else:
                # Usually exactly one HGNC record for a gene symbol
                if len(docs) > 1:
                    logger.warning(f"Multiple HGNC records found for symbol {symbol_str}. Using the first one.")
                result = docs[0].get("ensembl_gene_id")
            
            # Store in cache
            self.cache[symbol_str] = result
            
            # Save cache periodically (every 100 entries)
            if len(self.cache) % 100 == 0 and self.cache_path:
                self._save_cache()
            
            return result
            
        except Exception as e:
            logger.warning(f"HGNC API failed for symbol {symbol_str}: {e}")
            return None
    
    def save_cache(self) -> None:
        """Manually save cache to disk."""
        self._save_cache()


def fetch_ensg_ids_from_symbols(
    symbols: List[Union[str, int]], 
    cache_path: Optional[str] = None,
    requests_per_second: float = 10.0
) -> List[Optional[str]]:
    """
    Fetch Ensembl gene IDs for a list of gene symbols using HGNC API.
    
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
        Maximum number of requests per second (default: 10.0)
        
    Returns
    -------
    List[Optional[str]]
        List of Ensembl gene IDs (or None if not found) in same order as input
    """
    ensg_ids = []
    iteration_count = 0
    n_symbols = len(symbols)
    
    client = SymbolToEnsemblClient(cache_path=cache_path, 
                                   requests_per_second=requests_per_second)
    
    for symbol in symbols:
        iteration_count += 1
        ensg_id = client.symbol_to_ensg(symbol)
        ensg_ids.append(ensg_id)

        if iteration_count % 50 == 0:
            n_mapped_so_far = sum(1 for i in ensg_ids if i is not None)
            logger.info(f"Processed {iteration_count}/{n_symbols} gene symbols ({n_mapped_so_far} mapped so far)")

    client.save_cache()
    n_mapped_so_far = sum(1 for i in ensg_ids if i is not None)
    logger.info(f"Processed {iteration_count}/{n_symbols} gene symbols ({n_mapped_so_far} mapped so far)")
    return ensg_ids
