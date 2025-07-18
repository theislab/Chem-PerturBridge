import re
import logging
import argparse
import anndata
import pandas as pd
import scanpy as sc
from os import listdir
from os.path import isfile, join
from anndata import AnnData
from tqdm import tqdm

logger = logging.getLogger(__name__)

def get_files(dir_input: str):
    def extract_num(s, p, ret=0):
        search = p.search(s)
        if search:
            return int(search.groups()[0])
        else:
            return ret

    p = re.compile(r'file(\d+)')

    files = [file for file in listdir(dir_input) if isfile(join(dir_input, file))]
    return sorted(files, key=lambda s: extract_num(s, p, float('inf')))

def unite_adatas(dir_input: str):
    files = get_files(dir_input)
    adatas = []
    for i, file in enumerate(tqdm(files)):
        adata = sc.read_h5ad(join(dir_input, file))
        adata.obs['plate'] = i + 1
        adatas.append(adata)
    return anndata.concat(adatas)

def filter_cells(ads: AnnData,
                 groupby,
                 min_cells):
    vc = ads.obs[[groupby, 'psbulk_cells']].groupby(groupby, observed=False).sum()
    c_kept = vc[vc.psbulk_cells>min_cells].index.values
    adss = ads[ads.obs.drugname_drugconc.isin(c_kept)].copy()
    return adss

def save2csv(adss: AnnData, cell_line: str, dir_output: str):
    pd.DataFrame(adss.X).to_csv(join(dir_output, f"{cell_line}_X.csv"), index=False)
    adss.var.to_csv(join(dir_output, f"{cell_line}_var.csv"))
    adss.obs.to_csv(join(dir_output, f"{cell_line}_obs.csv"))

def save_by_cell_lines(adata: AnnData, 
                       dir_output: str, 
                       groupby: str = 'drugname_drugconc', 
                       min_cells: int = 100):
    
    cell_lines = sorted(set(adata.obs.cell_name))
    for cell_line in cell_lines:
        ads = adata[adata.obs['cell_name']==cell_line].copy()
        cell_line = re.sub(r'[^a-zA-Z0-9]', '', cell_line)
        adss = filter_cells(ads, groupby, min_cells)
        save2csv(adss, cell_line, dir_output)

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='type')
    parser_conf = subparsers.add_parser(name='parse_config')
    parser_conf.add_argument('--config', type=str)
    parser_arg = subparsers.add_parser(name='parse_args')
    parser_arg = parser_arg.add_argument_group()
    parser_arg.add_argument('--input', type=str, required=True)
    parser_arg.add_argument('--output', type=str, required=True)
    parser_arg.add_argument('--groupby', type=str)
    parser_arg.add_argument('--min_cells', type=int, default=100)
    args = parser.parse_args()
    if args.type == 'parse_args':
        logger.info("Unite datasets from several plates")
        adata = unite_adatas(args.input)
        logger.info("Save by cell lines")
        save_by_cell_lines(adata, args.output, min_cells=args.min_cells)

if __name__ == "__main__":
    main()
