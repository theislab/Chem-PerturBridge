import requests
import time
import json
import os
from typing import Dict, Optional, List

from src.utils.parsing_utils import *



class EntrezToEnsemblClient:
    """
    Client for mapping Entrez gene IDs to Ensembl gene IDs.
    Uses HGNC API as primary source, with NCBI Datasets API as fallback.
    Includes rate limiting (10 requests/second) and caching.
    """
    
    HGNC_BASE = "https://rest.genenames.org"
    NCBI_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    
    def __init__(self, cache_path: Optional[str] = None, requests_per_second: float = 10.0):
        """
        Initialize Entrez to Ensembl mapping client.
        
        Parameters:
        -----------
        cache_path : Optional[str]
            Path to JSON file for persistent cache storage
        requests_per_second : float
            Maximum number of requests per second (default: 10.0)
        """
        self.cache: Dict[str, Optional[str]] = {}
        self.cache_path = cache_path
        self.min_interval = 1.0 / requests_per_second  # 0.1 seconds for 10 req/s
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
                logger.info(f"Loaded Entrez-to-Ensembl cache from {self.cache_path} with {len(self.cache)} entries")
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
                logger.debug(f"Saved Entrez-to-Ensembl cache to {self.cache_path} with {len(self.cache)} entries")
            except Exception as e:
                logger.warning(f"Failed to save Entrez-to-Ensembl cache to {self.cache_path}: {e}")
    
    def _rate_limit(self) -> None:
        """Enforce rate limiting by sleeping if necessary."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < self.min_interval:
            sleep_time = self.min_interval - time_since_last
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
    
    def _entrez_to_ensg_hgnc(self, entrez_id_str: str) -> Optional[str]:
        """Convert Entrez gene ID to Ensembl gene ID using HGNC API."""
        # Rate limit
        self._rate_limit()
        
        # Make API request
        url = f"{self.HGNC_BASE}/fetch/entrez_id/{entrez_id_str}"
        try:
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
            r.raise_for_status()
            
            docs = r.json()["response"]["docs"]
            if not docs:
                return None
            
            # Usually exactly one HGNC record for a human Entrez gene ID
            if len(docs) > 1:
                logger.warning(f"Multiple HGNC records found for Entrez ID {entrez_id_str}. Using the first one.")
            return docs[0].get("ensembl_gene_id")
        except Exception as e:
            logger.warning(f"HGNC API failed for Entrez ID {entrez_id_str}: {e}")
            return None
    
    def _entrez_to_ensg_ncbi_datasets(self, entrez_id_str: str) -> Optional[str]:
        """
        Convert Entrez gene ID to Ensembl gene ID using NCBI Datasets API.
        Used as fallback when HGNC returns None.
        """
        # Rate limit
        self._rate_limit()
        
        url = f"{self.NCBI_BASE}/gene/id/{entrez_id_str}"
        headers = {"Accept": "application/json"}
        
        try:
            r = requests.get(url, headers=headers, timeout=60)
            r.raise_for_status()
            data = r.json()
            
            # NCBI Datasets response structure: data['reports'][0]['gene']['ensembl_gene_ids']
            # The response is NOT keyed by gene ID for single gene queries
            reports = data.get("reports", [])
            if not reports:
                return None
            
            # Get first report's gene data
            gene = reports[0].get("gene", {})
            ensg_list = gene.get("ensembl_gene_ids", []) or []
            
            # Filter for Ensembl gene IDs (ENSG...) and return first one
            ensg_list = [e for e in ensg_list if isinstance(e, str) and e.startswith("ENSG")]
            return ensg_list[0] if ensg_list else None
        except Exception as e:
            logger.warning(f"NCBI Datasets API failed for Entrez ID {entrez_id_str}: {e}")
            return None
    
    def entrez_to_ensg(self, entrez_id: str | int) -> Optional[str]:
        """
        Convert Entrez gene ID to Ensembl gene ID.
        Tries HGNC API first, then NCBI Datasets API if HGNC returns None.
        
        Parameters:
        -----------
        entrez_id : str | int
            Entrez gene ID
            
        Returns:
        --------
        Optional[str]
            Ensembl gene ID or None if not found
        """
        entrez_id_str = str(entrez_id)
        
        # Check cache first (only return if not None)
        if entrez_id_str in self.cache and self.cache[entrez_id_str] is not None:
            return self.cache[entrez_id_str]
        
        # Strategy 1: Try HGNC API
        result = self._entrez_to_ensg_hgnc(entrez_id_str)
        
        # Strategy 2: If HGNC returned None, try NCBI Datasets API
        if result is None:
            result = self._entrez_to_ensg_ncbi_datasets(entrez_id_str)
        
        # Store in cache only if result is not None
        if result is not None:
            self.cache[entrez_id_str] = result
        
        # Save cache periodically (every 100 entries)
        if len(self.cache) % 100 == 0 and self.cache_path:
            self._save_cache()
        
        return result
    
    def save_cache(self) -> None:
        """Manually save cache to disk."""
        self._save_cache()


def fetch_ensg_ids(entrez_ids: List[str | int], 
                   cache_path: Optional[str] = None,
                   requests_per_second: float = 10.0) -> List[Optional[str]]:
    """
    Fetch Ensembl gene IDs for a list of Entrez gene IDs.
    
    Parameters:
    -----------
    entrez_ids : List[str | int]
        List of Entrez gene IDs to map
    cache_path : Optional[str]
        Path to JSON file for persistent cache storage (default: None)
    requests_per_second : float
        Maximum number of requests per second (default: 10.0)
        
    Returns:
    --------
    List[Optional[str]]
        List of Ensembl gene IDs (or None if not found) in same order as input
    """
    ensg_ids = []
    iteration_count = 0
    n_ids = len(entrez_ids)
    
    # Use default cache path if not provided
    if cache_path is None:
        cache_path = 'hgnc_cache.json'
    
    client = EntrezToEnsemblClient(cache_path=cache_path, 
                                   requests_per_second=requests_per_second)
    for idx in entrez_ids:
        iteration_count += 1
        eid = client.entrez_to_ensg(idx)
        ensg_ids.append(eid)

        if iteration_count % 50 == 0:
            n_mapped_so_far = sum(1 for i in ensg_ids if i is not None)
            logger.info(f"Processed {iteration_count}/{n_ids} entrez ids ({n_mapped_so_far} mapped so far)")

    client.save_cache()
    n_mapped_so_far = sum(1 for i in ensg_ids if i is not None)
    logger.info(f"Processed {iteration_count}/{n_ids} entrez ids ({n_mapped_so_far} mapped so far)")
    return ensg_ids