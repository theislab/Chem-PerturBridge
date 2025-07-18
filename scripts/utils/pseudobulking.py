import os
import logging
import argparse
import pandas as pd
import scanpy as sc
import decoupler as dc
import pubchempy as pcp
from anndata import AnnData
from typing import Dict, Optional

logger = logging.getLogger(__name__)

class Pseudobulk:
    def __init__(self, file_input: list,
                       file_output: str,
                       bulk_fields: list,
                       dropna: list = [],
                       drop_cols: list = [],
                       select_cols: list = [],
                       min_cells: int = 0):   
        self.file_input = file_input
        dir_output = os.path.dirname(file_output)
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)
        self.file_output = file_output
        self.bulk_fields = bulk_fields
        self.dropna = dropna
        self.drop_cols = drop_cols
        self.select_cols = select_cols
        self.min_cells = min_cells
        self.padata = None
        
    def run_pseudobulking(self, ):
        #for plate, file in enumerate(self.files_input):
        logger.info("Read anndata")
        adata = sc.read_h5ad(self.file_input)
        logger.info("Process anndata")
        adata = self.tiny_preprocessing(adata)
        logger.info("Build pseudobulk")
        self.padata = self.build_pseudoulk(adata)
        logger.info("Filter columns")
        self.filter_cols()
        
            
    def save(self, ):
            if not self.padata is None:
                self.padata.write_h5ad(self.file_output, compression="gzip")
            
    def tiny_preprocessing(self, adata: AnnData):
        adata.layers['counts'] = adata.X.copy()
        if len(self.dropna) != 0:
            adata = adata[~adata.obs[self.dropna].isna().all(axis=1)].copy() 
        return adata
    
    def define_sample_group_cols(self, adata: AnnData):
        from pandas.api.types import is_float_dtype
        for f in self.bulk_fields:
            if f not in adata.obs:
                raise KeyError(f"adata.obs lacks '{f}'")
            else:
                if is_float_dtype(adata.obs[f]):
                    adata.obs[f] = adata.obs[f].astype(int)
        
        adata.obs['sample_id'] = adata.obs[self.bulk_fields] \
            .astype(str) \
            .apply(lambda row: '_'.join(row.values), axis=1)
        adata.obs['pseudo_group'] = 'all'

    def set_sample_idx(self, padata: AnnData):
        padata.obs['sample_id'] = padata.obs[self.bulk_fields] \
            .astype(str) \
            .apply(lambda row: '_'.join(row.values), axis=1)
        padata.obs.set_index('sample_id', inplace=True)
        padata.obs.drop(columns=['pseudo_group'], inplace=True)
        
    def filter_cells(self, padata: AnnData):
        padata = padata[padata.obs.psbulk_cells >= self.min_cells].copy()
        return padata
        
    def build_pseudoulk(self, adata: AnnData):
        self.define_sample_group_cols(adata)
        
        padata = dc.pp.pseudobulk(adata,
                                  sample_col='sample_id',
                                  groups_col='pseudo_group',
                                  layer='counts')
        self.set_sample_idx(padata)
        if self.min_cells > 0:
            padata = self.filter_cells(padata)
        del padata.layers["psbulk_props"]
        return padata

    def filter_cols(self):
        if len(self.drop_cols) > 0:
            self.padata.obs.drop(columns=self.drop_cols, inplace=True)
        if len(self.select_cols) > 0:
            self.padata.obs = self.padata.obs[self.select_cols + ['psbulk_cells', 'psbulk_counts']]


    def add_pubchem_cids_to_padata(self, cache: Dict[str, Optional[int]], drug_col='perturbation') -> None:
        """
        Add a `pubchem_cid` column to adata.obs based on the 'perturbation' field.
        """
        from tqdm import tqdm
        def get_pubchem_cid(drug_name: str, cache: Dict[str, Optional[int]]) -> Optional[int]:
            """
            Fetch PubChem CID for a given drug name, using cache to skip repeat lookups.
            """
            if pd.isna(drug_name) or not drug_name:
                return None
        
            if drug_name in cache:
                return cache[drug_name]
            try:
                compounds = pcp.get_compounds(drug_name, 'name')
                cid = compounds[0].cid if compounds else None
                logger.debug(f"CID for '{drug_name}': {cid}")
            except Exception as e:
                logger.warning(f"PubChem lookup failed for '{drug_name}': {e}")
                cid = None
            cache[drug_name] = cid
            return cid
            
        unique_drugs = self.padata.obs[drug_col].dropna().unique().tolist()
        to_fetch = [d for d in unique_drugs if d not in cache]
        if to_fetch:
            logger.info(f"Looking up {len(to_fetch)} new perturbations on PubChem...")
            for drug in tqdm(to_fetch, desc="PubChem CIDs"):
                get_pubchem_cid(drug, cache)
    
        self.padata.obs['pubchem_cid'] = self.padata.obs[drug_col].map(cache)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='type')
    parser_conf = subparsers.add_parser(name='parse_config')
    parser_conf.add_argument('--config', type=str)
    parser_arg = subparsers.add_parser(name='parse_args')
    parser_arg = parser_arg.add_argument_group()
    parser_arg.add_argument('--input', type=str, required=True)
    parser_arg.add_argument('--output', type=str, required=True)
    parser_arg.add_argument('--bulk_fields', nargs='+', required=True)
    parser_arg.add_argument('--dropna', nargs='+', default=[])
    parser_arg.add_argument('--drop_cols', nargs='+', default=[])
    parser_arg.add_argument('--select_cols', nargs='+', default=[])
    parser_arg.add_argument('--min_cells', type=int, default=0)
    parser_arg.add_argument('--drug_col', type=str, required=True)
    args = parser.parse_args()
    if args.type == 'parse_args':
        pseudo = Pseudobulk(file_input = args.input,
                            file_output = args.output,
                            bulk_fields = args.bulk_fields,
                            dropna = args.dropna,
                            drop_cols = args.drop_cols,
                            select_cols = args.select_cols,
                            min_cells = args.min_cells,
                           )
        pseudo.run_pseudobulking()
        pseudo.add_pubchem_cids_to_padata({}, drug_col=args.drug_col)
        logger.info("Save dataset")
        pseudo.save()

if __name__ == "__main__":
    main()
