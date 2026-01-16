# Guide: Adding a New Dataset to the OP3_v2 Pipeline

This document describes step-by-step integration of the new dataset to the OP3_v2 pipeline.

## General description

The code could be divided into 2 parts: 
* `.py` and `.R` scripts which are targeted to work with the content of datasets, they are located in `op3_v2/src` directory (and related config files)
* `.sh` wrapper scripts located in `op3_v2/pipelines` and `op3_v2/run_pipelines` (and related config files)


`.py` and `.R` scripts are mainly focused on the data downloading, processing, data annotation, aggregation, while .sh are responsible for setting the parameters for running `.py` and `.R` scripts, initiating the run for the pipeline and stacking the scripts in one workflow if it is needed or parallelizing the jobs.

If you want to add the dataset-specific functions <u>to process a dataset</u> and modifies them, you need to use `op3_v2/src` directory: 
* add the scripts to `op3_v2/src/pseudobulking/datasets/mynewdataset` (e.g. `standardization.py`)
* modify `op3_v2/src/configs/datasets.json` pointing the location for the dataset-specific module (e.g. `op3_v2/src/pseudobulking/datasets/mynewdataset/standardization.py`).

Once you add the processing functions unique for your dataset, you need to create/update wrapper scripts responsible for <u>running the pipelines</u> in `op3_v2/pipelines` and `op3_v2/run_pipelines` and add config files as well.

**NB!** Data and Log directories are created automatically when you execute the pipelines.

**NB!** The pseudobulking pipeline is oriented mainly to Laminlabs annotated datasets as the main source for data, therefore before adding the datasets the user needs to check if the dataset has been alredy annotated by Laminlabs and represented on their [site](https://lamin.ai/laminlabs/pertdata). Otherwise the dataset-specific loader needs to be added to `op3_v2/src/pseudobulking/datasets/mynewdataset`.

## Overview

### 1. Update a pseudobulking step:

- Add the dataset folder `mynewdataset` to `src/pseudobulking/datasets/` which stores scripts with the dataset-specific standardization functions

- Update `src/configs/datasets.json` file which stores the paths to the dataset-specific functions mentioned above

- Update shell scripts:
  - Update `run_pipelines/run_pseudobulking.sh`
  - Create the dataset folder `mynewdataset` in `pipelines/` and add `mynewdataset_pseudobulking.sh` script to this folder
  
- Add `pipelines/mynewdataset/configs/` folder

### 2. Update a DEG step:

- Update shell script:
  - Update `run_pipelines/run_deg.sh`

- If you are planning to process additionally the obtained pseudobulk, then you need to add post-processing functions to `src/pseudobulking/datasets/mynewdataset`

- Create `pipelines/mynewdataset/configs/deg/` directory

- Add `config.json` to the `pipelines/mynewdataset/configs/deg/` which keeps the parameters for DEG analysis

