import os
import pandas as pd
import scanpy as sc
import anndata as ad

from src.utils.parsing_utils import *

def create_perturbation_label(is_control: bool, 
                              pert: str, 
                              dose: float, 
                              time: float,
                              plate: str,
                              design_param: str) -> str:
    
    if is_control:
        if design_param == 'group_all_replicates':
            return str(pert) \
                + '_' + str(time) + 'h'
        elif design_param == 'separate_replicates':
            return str(pert) \
                + '_' + str(time) + 'h' \
                + '_' + str(plate)
    else:
        return str(pert) \
            + '_' + str(dose) \
            + 'uM_' + str(time) + 'h'

def create_dir_if_not_exists(file_output: str) -> None:
    dir_output = os.path.dirname(file_output)
    if not os.path.exists(dir_output):
        try:
            os.makedirs(dir_output)
        except FileExistsError as e:
            logger.warning('%s', str(e))

def add_perturbation_label_to_padata(file_input: str,
                                     file_output: str,
                                     design_param: str,
                                     ) -> None:
    logger.info('Read pseudobulk file')

    padata = ad.read_h5ad(file_input)
    obs = padata.obs.copy()
    obs['pert_dose_uM'] = obs['pert_dose_uM'].apply(lambda x: format(x, ".15g"))
    obs['pert_time_h'] = obs['pert_time_h'].apply(lambda x: format(x, ".15g"))
    
    logger.info('Add perturbation label')
    padata.obs['perturbation_label'] =  obs.apply(lambda x: create_perturbation_label(x.is_control,
                                                   x.perturbagen,
                                                   x.pert_dose_uM,
                                                   x.pert_time_h,
                                                   x.plate,
                                                   design_param), axis=1).astype("category")
    
    logger.info('Save data')
    create_dir_if_not_exists(file_output)
    padata.write_h5ad(file_output, compression='gzip')


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str)
    parser.add_argument('--input', type=str)
    parser.add_argument('--output', type=str)
    parser.add_argument('--design_param', choices=['group_all_replicates',
                                                   'separate_replicates'])
    args = parser.parse_args()
    d_args = vars(args).copy()
    del d_args['config']

    if not args.config is None:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
        d_args = merge_args(d_args, config)

    for key in ['input', 'output', 'design_param']:
        if not d_args.get(key):
            raise Exception(f'The argument {key} is not set')

    
    add_perturbation_label_to_padata(d_args['input'],
                                     d_args['output'],
                                     d_args['design_param'])



if __name__ == '__main__':
    main()
