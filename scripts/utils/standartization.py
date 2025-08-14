from string import ascii_lowercase
import scanpy as sc


STANDARD_OBS_COLS = ['plate',
                 'well',
                 'pert_time',
                 'ngenes',
                 'ncounts',
                 'pcnt_mito',
                 'cell_type',
                 'perturbagen',
                 'pert_type',
                 'pert_dose',
                 'tissue_type',
                 'tissue',
                 'disease',
                 'is_control',
                 'library',
                 'stimulation',
                 'guide',
                 'dataset',
                ]

STANDARD_VAR_COLS = ['symbol']


def standardize_obs_tahoe(adata):
    rename_cols = {'gene_count': 'ngenes', 
                   'tscp_count': 'ncounts',
                   'cell_line_ontology_id': 'cell_type',
                   'pert_compound': 'perturbagen',
                   'sublibrary': 'library',
                   
                  }
    wellnum2id_mapping = wellnum2id()
    adata.obs['well'] = adata.obs.index.to_series().str.split('_').str[0].astype(int).map(wellnum2id_mapping)
    #adata.obs['well'] = adata.obs['plate'] + '_' + adata.obs['well']
    adata.obs.rename(columns=rename_cols, inplace=True)
    adata.obs['is_control'] = adata.obs['perturbagen'] == 'DMSO_TF'
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'tahoe100'
    adata.obs = adata.obs[STANDARD_OBS_COLS]
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
    adata.obs['is_control'] = adata.obs['perturbagen'].str.lower() == 'control'
    adata.obs['library'] = None
    adata.obs['stimulation'] = None
    adata.obs['guide'] = None
    adata.obs['dataset'] = 'srivatsan20_sciplex3'
    adata.obs = adata.obs[STANDARD_OBS_COLS]
    
    return adata

def standardize_var_sciplex(adata):
    adata.var = adata.var[STANDARD_VAR_COLS]
    return adata

dataset_standardization = {'tahoe': {'standardize_obs': standardize_obs_tahoe,
                                     'standardize_var': standardize_var_tahoe},
                           'sciplex': {'standardize_obs': standardize_obs_sciplex,
                                     'standardize_var': standardize_var_sciplex}
                          }
