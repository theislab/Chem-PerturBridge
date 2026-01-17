# Guide: Adding a New Dataset to the OP3_v2 Pipeline

This document describes step-by-step integration of a new dataset with an identifier `mynewdataset` into the OP3_v2 pipeline.

## General description

The whole pipeline could be divided into two main steps:
* Pseudobulking
* DEG analysis

The executed code could be divided into two parts:

* `.py` and `.R` scripts which are targeted to work with the content of datasets, they are located in `op3_v2/src` directory (and related config files)
* `.sh` wrapper scripts located in `op3_v2/pipelines` and `op3_v2/run_pipelines` (and related config files)

`.py` and `.R` scripts are mainly focused on data downloading, processing, and annotation, while `.sh` scripts are responsible for setting parameters for running `.py` and `.R` scripts, initiating the pipeline, and stacking scripts into a workflow when needed or parallelizing the jobs.

**NB!** Data and Log directories are created automatically when you execute the pipelines.

## Overview

### 1 Update a pseudobulking step:

**1.1** **`.py` scripts and config**

> **1.1.1** **directories**
> * Create a directory `op3_v2/src/pseudobulking/datasets/mynewdataset` which will store modules with the dataset-specific functions.
>
> **1.1.2** **scripts**
> * Add python modules to the directory created in **1.1.1** containing specific helping functions to prepare the dataset which you want to add.
> * **NB!** The pseudobulking pipeline is oriented mainly to Laminlabs annotated datasets as the main source for data, therefore before adding a dataset you need to check if the dataset has been already annotated by Laminlabs and represented on their [site](https://lamin.ai/laminlabs/pertdata). Otherwise the script which downloads and assembles your dataset, and saves in `.h5ad` format needs to be added to `op3_v2/src/pseudobulking/datasets/mynewdataset`.
> * **NB!** If you already has the bulked/pseudobulked version of the data, there is no need to compute pseudobulks. But you need to add the assembling modules to `op3_v2/src/pseudobulking/datasets/mynewdataset`, and execute them later with `op3_v2/run_pipelines/run_pseudobulking.sh`. Look at **L1000** case for the detailes.
>
> **1.1.3** **configs**
> * Update a file `op3_v2/src/configs/datasets.json` which stores locations to the modules from **1.1.2**.

**1.2** **`.sh` scripts and configs**
>
> **1.2.1** **directories**
> * Create a directory `op3_v2/pipelines/mynewdataset` which will store dataset-specific wrapping script to start the execution of the python modules on the cluster.
>
> **1.2.2** **scripts**
> * Add `mynewdataset_pseudobulking.sh` to the directory from **1.2.1**.
> * Update `op3_v2/run_pipelines/run_pseudobulking.sh`:
>   - Add the name of the dataset (`mynewdataset` in our case) to `VALID_CHOICES` variable.
>   - Add the condition specifying the new dataset and command to initiate the execution of `op3_v2/pipelines/mynewdataset/mynewdataset_pseudobulking.sh`.
>
> **1.2.3** **configs**
> * No need to add configs on this step for the current version of the repository. In the future there might be updates, but not now.

### **2** Update a DEG anaysis step:

**2.1** **`.py` scripts and config**
> **2.1.1** **directories**
> * If you followed steps from **1.1** you already have the dataset-specific directory. Otherwise, optionally, follow the step **1.1.1** in case you are planning to add the unique module to post-process your bulked/pseudobulked dataset.
>
> **2.1.2** **scripts**
> * Optionally, add python modules to the directory from **2.1.1** in case you need to post-process bulked/pseudobulked data. The dataset-specific post-processing will be the part of `op3_v2/src/deg/run_processing_pseudobulk.py`.
>    
> **2.1.3** **configs**
> * Optionally, follow the step **1.1.3** if you added the module on the step **2.1.2**.
>
**2.2** **`.sh` scripts and configs**
> **2.2.1** **directories**
> * Create a directory `op3_v2/pipelines/mynewdataset/configs`.
> * Create a directory `op3_v2/pipelines/mynewdataset/configs/deg`.
>
> **2.2.2** **scripts**
> * Update `run_pipelines/run_deg.sh`:
>   - Add the name of the dataset (`mynewdataset` in our case) to `VALID_CHOICES` variable.
>
> **2.2.3** **configs**
> * Create a file `op3_v2/pipelines/mynewdataset/configs/deg/config.json` which keeps the parameters for DEG analysis.