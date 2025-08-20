from string import ascii_lowercase
import re
import scanpy as sc
from magnitude import mg, new_mag, Magnitude

new_mag('l', Magnitude(0.001, m=3))
new_mag('M', mg(1, 'mol/l'))

STANDARD_ADATA_OBS_COLS = [
                 'plate', 
                 'well', 
                 'ngenes', 
                 'ncounts', 
                 'pcnt_mito',
                 'cell_type', 
                 'perturbagen',
                 'pert_type',
                 'is_control',
                 'pert_dose',
                 'pert_time',
                 'suspension_type',
                 'tissue',
                 'tissue_type', 
                 'disease',
                 'library',
                 'stimulation',
                 'guide',
                 'dataset',
                 'assay',
                 'development_stage',
                 'organism',
                 'sex',
                 'self_reported_ethnicity',
                 ]

CAT_COLS = [
                 'plate', 
                 'well', 
                 'cell_type', 
                 'perturbagen',
                 'pert_type',
                 'suspension_type',
                 'tissue',
                 'tissue_type', 
                 'disease',
                 'library',
                 'stimulation',
                 'guide',
                 'dataset',
                 'development_stage',
                 'organism',
                 'sex',
                 'self_reported_ethnicity'
                 ]


STANDARD_VAR_COLS = ['symbol']

def split_val_units(s):
    r = "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
    try:
        val = float(re.findall(r, s)[0])
    except Exception as e:
        val = None
    try:
        units = re.split(r, s)[-1]
    except Exception as e:
        units = None
    return [val, units]

def standardize_doze(value, unit, default_unit='uM'):
    try:
        value_ = mg(value, unit).toval(default_unit)     
    except Exception as e:
        value_ = None
    return [value_, default_unit]

def standardize_obs_tahoe(adata):
    rename_cols = {'gene_count': 'ngenes', 
                   'tscp_count': 'ncounts',
                   'cell_line_ontology_id': 'cell_type',
                   'pert_compound': 'perturbagen',
                  }
    categories = {c: 'category' for c in CAT_COLS}
    wellnum2id_mapping = wellnum2id()
    adata.obs['well'] = adata.obs.index.to_series().str.split('_').str[0].astype(int).map(wellnum2id_mapping)
    #adata.obs['well'] = adata.obs['plate'] + '_' + adata.obs['well']
    adata.obs.rename(columns=rename_cols, inplace=True)
    adata.obs['is_control'] = adata.obs['perturbagen'] == 'DMSO_TF'
    adata.obs['library'] = None
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'tahoe100'
    adata.obs = adata.obs[STANDARD_ADATA_OBS_COLS]
    adata.obs = adata.obs.astype(categories)
    return adata


def standardize_var_tahoe(adata):
    rename_cols = {'ensembl_gene_id': 'ensembl_id'}
    adata.var.rename(columns=rename_cols, inplace=True)
    adata.var.set_index('ensembl_id', inplace=True)
    return adata

def wellnum2id(nwells=96, ncols=12):
    mapping = {}
    for i in range(nwells):
        mapping[i + 1] = ascii_lowercase[i // ncols].upper() + str(i % ncols + 1)
    return mapping

def standardize_obs_sciplex(adata):
    rename_cols = {'n_genes': 'ngenes', 
                   'pct_counts_mt': 'pcnt_mito', 
                   'cell_type': 'cell_type_', 
                   'cell_line': 'cell_type',
                   'pert_compound': 'perturbagen'
                  }
    categories = {c: 'category' for c in CAT_COLS}

    map_clnames2ids = {'A549': 'CVCL_0023',
                       'K-562': 'CVCL_0004',
                       'MCF7': 'CVCL_0031',
                       }
    
    sc.pp.filter_cells(adata, min_genes=0)
    adata.var["mt"] = adata.var['symbol'].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
    adata.obs['pct_counts_mt'] /= 100
    adata.obs.rename(columns=rename_cols, inplace=True)
    adata.obs['cell_type'] = adata.obs['cell_type'].cat.rename_categories(map_clnames2ids)
    adata.obs['is_control'] = adata.obs['perturbagen'] == 'control'
    adata.obs['perturbagen'] = adata.obs['perturbagen'].cat.add_categories(['DMSO'])
    adata.obs.loc[adata.obs['perturbagen'] == 'control', 'perturbagen'] = 'DMSO'
    adata.obs['perturbagen'] = adata.obs['perturbagen'].cat.remove_categories(['control'])
    adata.obs['library'] = None
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'srivatsan20_sciplex3'
    adata.obs = adata.obs[STANDARD_ADATA_OBS_COLS]
    adata.obs = adata.obs.astype(categories)
    return adata

def standardize_var_sciplex(adata):
    adata.var = adata.var[STANDARD_VAR_COLS]
    return adata


dataset_standardization = {'tahoe': {'standardize_obs': standardize_obs_tahoe,
                                     'standardize_var': standardize_var_tahoe},
                           'sciplex': {'standardize_obs': standardize_obs_sciplex,
                                     'standardize_var': standardize_var_sciplex}
                          }
