# Adding a New Dataset to the Chem-PerturBridge Pipeline

This document describes step-by-step integration of a new dataset with an identifier `mynewdataset` into the **Chem-PerturBridge** pipeline.

## General description

The whole **pipeline** could be divided into two main steps:
* Pseudobulking
* DEG analysis

The **executed code** could be divided into two parts:

* `.py` and `.R` scripts which are targeted to work with the content of datasets, they are located in `Chem-PerturBridge/src` directory (and related config files)
* `.sh` wrapper scripts located in `Chem-PerturBridge/pipelines` and `Chem-PerturBridge/run_pipelines` (and related config files)

`.py` and `.R` scripts are mainly focused on data downloading, processing, aggregation, and annotation, while `.sh` scripts are responsible for setting parameters for running `.py` and `.R` scripts, initiating the pipeline, and stacking scripts into a workflow when needed or parallelizing the jobs.

**NB!** Data and log directories are created automatically when you execute the pipelines.

## Quick start guide

**NB!** Use `sciplex` or `tahoe` as reference examples by replacing `mynewdataset` in the scripts and directories.

### 1 Update a pseudobulking step:

**1.1** **`.py` scripts and config**

> **1.1.1** **directories**
> * Create a directory `Chem-PerturBridge/src/pseudobulking/datasets/mynewdataset` which will store modules containing dataset-specific functions.
>
> **1.1.2** **scripts**
> * Add python modules to the directory created in **1.1.1** containing specific helping functions to prepare the dataset which you want to add.
> * **NB!** Please, check out [this description](https://github.com/theislab/Chem-PerturBridge/blob/main/docs/Format_Pseudobulk.ipynb) of a core standardized schema for bulked/pseudobulked data which you need to receive at the end of pseudobulking pipeline. This schema is not entirely strict, so you can add some other columns (e.g. a `split` column including train/test/validation label) which you think contains useful information.
> * **NB!** The pseudobulking pipeline is oriented mainly to Laminlabs annotated datasets as the main source for data, therefore before adding a dataset you need to check if the dataset has been already annotated by Laminlabs and represented on their [site](https://lamin.ai/laminlabs/pertdata). Otherwise the script that downloads and assembles your dataset and saves it in `.h5ad` format needs to be added to `Chem-PerturBridge/src/pseudobulking/datasets/mynewdataset`. Look at [**L1000** case](https://github.com/theislab/Chem-PerturBridge/tree/main/src/pseudobulking/datasets/l1000) for the details.
> * **NB!** If you already have the bulked/pseudobulked version of the data, there is no need to compute pseudobulks. But you need to add the assembling modules to `Chem-PerturBridge/src/pseudobulking/datasets/mynewdataset`, and execute them later with `Chem-PerturBridge/run_pipelines/run_pseudobulking.sh`. Look at [**L1000** case](https://github.com/theislab/Chem-PerturBridge/tree/main/src/pseudobulking/datasets/l1000) for the details.
>
> **1.1.3** **configs**
> * Update a file `Chem-PerturBridge/src/configs/datasets.json` which stores paths to the modules from **1.1.2**.

**1.2** **`.sh` scripts and configs**
>
> **1.2.1** **directories**
> * Create a directory `Chem-PerturBridge/pipelines/mynewdataset` which will store dataset-specific wrapping script to start the execution of the python modules on the cluster.
>
> **1.2.2** **scripts**
> * Add `mynewdataset_pseudobulking.sh` to the directory from **1.2.1**.
> * Update `Chem-PerturBridge/run_pipelines/run_pseudobulking.sh`:
>   - Add the name of the dataset (`mynewdataset` in our case) to `VALID_CHOICES` variable.
>   - Add the condition specifying the new dataset and command to initiate the execution of `Chem-PerturBridge/pipelines/mynewdataset/mynewdataset_pseudobulking.sh`.
>
> **1.2.3** **configs**
> * No need to add configs on this step for the current version of the repository. In the future there might be updates, but not now.

### **2** Update a DEG analysis step:

**2.1** **`.py` scripts and config**
> **2.1.1** **directories**
> * If you followed steps from **1.1** you already have the dataset-specific directory. Otherwise, optionally, follow the step **1.1.1** in case you are planning to add the unique module to post-process your bulked/pseudobulked dataset.
>
> **2.1.2** **scripts**
> * Optionally, add python modules to the directory from **2.1.1** in case you need to post-process bulked/pseudobulked data. The dataset-specific post-processing will be the part of `Chem-PerturBridge/src/deg/run_processing_pseudobulk.py`.
>    
> **2.1.3** **configs**
> * Optionally, follow the step **1.1.3** if you added the module on the step **2.1.2**.
>
**2.2** **`.sh` scripts and configs**
> **2.2.1** **directories**
> * Create a directory `Chem-PerturBridge/pipelines/mynewdataset/configs`.
> * Create a directory `Chem-PerturBridge/pipelines/mynewdataset/configs/deg`.
>
> **2.2.2** **scripts**
> * Update `run_pipelines/run_deg.sh`:
>   - Add the name of the dataset (`mynewdataset` in our case) to `VALID_CHOICES` variable.
> * **NB!** If you already have the bulked/pseudobulked version of the data, you need to figure out if your data is already normalized or it has raw counts. Use `-n` option while executing `run_deg.sh` if your data is already nomalized, otherwise for raw counts skip setting `-n` option. Look at the example of [**L1000** case](https://github.com/theislab/Chem-PerturBridge?tab=readme-ov-file#3-run-scripts) in README.md for the details.
>
> **2.2.3** **configs**
> * Create a file `Chem-PerturBridge/pipelines/mynewdataset/configs/deg/config.json` which keeps the parameters for DEG analysis.


### **3** Validate your integration

**3.1** **Test with subsample `-s` first**

**NB!** Always test with a subsample before running on the full dataset.

> **3.1.1** **run the script**
>
>```bash
># E.g. Test pseudobulking
>cd /path/to/Chem-PerturBridge
>./run_pipelines/run_pseudobulking.sh -s -d mynewdataset
>```

> **3.1.2** **check logs**
>
>```bash
># E.g. Monitor logs for pseudobulking
>tail -f logs/mynewdataset/subsample/pseudobulk_mynewdataset.*.out
>tail -f logs/mynewdataset/subsample/pseudobulk_mynewdataset.*.err
>```

> **3.1.3** **verify pseudobulk output**
>
>```bash
># E.g. Check that the pseudobulk data was created correctly
>ls -lh data/mynewdataset/pseudobulk/subsample/
>```

**3.2** **Run on the full dataset...**

For full command usage and examples, refer to [`../README.md`](https://github.com/theislab/Chem-PerturBridge/blob/main/README.md) in the `main` branch.

### **4** **Optionally...**

If you are a contributor, it would be great to update `README.md` file by adding information on the dataset which you processed (url link to dataset paper, a license of a processed dataset, and a license for additional sources used in dataset processing). Just follow our sample dataset descriptions in main [`../README.md`](https://github.com/theislab/Chem-PerturBridge/blob/main/README.md).

## Additional Resources

- Main README: [`../README.md`](https://github.com/theislab/Chem-PerturBridge/blob/main/README.md)
- Check existing files and directories for the following datasets: **sciplex**, **tahoe**, **L1000**, **op3** in the `main` branch as reference examples.
- [Core pseudobulk schema](https://github.com/theislab/Chem-PerturBridge/blob/main/docs/Format_Pseudobulk.ipynb)
- `Chem-PerturBridge/docs/STANDARDIZATION.md` - a more detailed description of how to standardize data.
