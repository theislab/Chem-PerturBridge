# Description

**OP3_v2** is a set of pipelines for analyzing single-cell RNA sequencing data from perturbation experiments. 

The current version of OP3_v2 includes the scripts for processing and analyzing [**Sci-Plex**](https://www.science.org/doi/10.1126/science.aax6234) and [**Tahoe**](https://www.biorxiv.org/content/10.1101/2025.02.20.639398v1.full) datasets.

It consists of the following steps:

* **Pseudobulking**: 
Aggregates raw single-cell RNA sequencing data to pseudobulk samples for downstream analysis

* **DEG Analysis**: 
Identifies differentially expressed genes between treatment and control conditions in pseudobulk samples

# License

This repository is a collection; each component retains its original license; our processed integration is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode.txt).

The datasets used in this project are provided under the following licenses:

- **Sciplex3** – Creative Commons Attribution 4.0 International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode.txt)) - [scPerturb Single-Cell Perturbation Data](https://zenodo.org/records/13350497) additionally annotated by [Laminlabs](https://lamin.ai/laminlabs/pertdata/artifacts?filter%5Band%5D%5B0%5D%5Bor%5D%5B0%5D%5Bbranch.name%5D%5Beq%5D=main&filter%5Band%5D%5B1%5D%5Bor%5D%5B0%5D%5Bis_latest%5D%5Beq%5D=true&filter%5Band%5D%5B2%5D%5Bor%5D%5B0%5D%5Bprojects.name%5D%5Beq%5D=scPerturb)
- **Tahoe** - Creative Commons Zero v1.0 Universal ([CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/legalcode.txt)) - [Arc Virtual Cell Atlas](https://arcinstitute.org/tools/virtualcellatlas) additionally annotated by [Laminlabs](https://lamin.ai/laminlabs/pertdata/artifacts?filter%5Band%5D%5B0%5D%5Bor%5D%5B0%5D%5Bbranch.name%5D%5Beq%5D=main&filter%5Band%5D%5B1%5D%5Bor%5D%5B0%5D%5Bis_latest%5D%5Beq%5D=true&filter%5Band%5D%5B2%5D%5Bor%5D%5B0%5D%5Bprojects.name%5D%5Beq%5D=Tahoe-100M)
- **L1000** - Contains data from GEO accessions [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) and [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138). Upstream terms apply; we do not assert CC BY for upstream L1000. See provenance + preprocessing notes. Datasets were additionally annotated by [Laminlabs](https://lamin.ai/laminlabs/pertdata/artifacts?filter%5Band%5D%5B0%5D%5Bor%5D%5B0%5D%5Bbranch.name%5D%5Beq%5D=main&filter%5Band%5D%5B1%5D%5Bor%5D%5B0%5D%5Bis_latest%5D%5Beq%5D=true&filter%5Band%5D%5B2%5D%5Bor%5D%5B0%5D%5Bprojects.name%5D%5Beq%5D=LINCS).
- **OP3** – Creative Commons Attribution 4.0 International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode.txt)) - [Open Problems Perturbation Prediction benchmark](https://openreview.net/forum?id=WTI4RJYSVm)

# Quick start guide
You can run scripts on the HPC cluster with the workload manager (SLURM)

**NB!** Before running some bash scripts, you may need to run `chmod u+x ./path/script.sh` before their execution.

**NB!** After the first trial of fetching the dataset from lamindb with the pipeline, the full version of data has been downloaded into your `~/.cache` directory as well as the version from your request (full/subsample) has been downloaded into the specified directory. Then, during the next re-runs the data is retrieved from `~/.cache` synchronization to the external database.

### 1. Set up prerequisites
Before executing scripts we need to install packages into the certain environment. Here we provide the illustrative scripts for using a package manager called **Mamba**. Also you can adapt our scripts to use with other managers such as Micromamba or Conda.

**1.1** **Install and prepare Mamba**

**1.1.1** To install Mamba you can follow instructions described in the [Miniforge](https://github.com/conda-forge/miniforge) repository, and provide `yes` answers to installer's questions:

```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

**1.1.2** To stop running base environment automatically after running the terminal, you need to execute:
```
conda deactivate
conda config --set auto_activate_base false
```

**1.1.3** Additionally, once you install Mamba and Conda with Miniforge, add the path of their location to `~/.bashrc` file:
```
export PATH="/you_path/miniforge3/bin:$PATH"
```
And then run:
```
source ~/.bashrc
```
**NB!** If you installed Mamba and have a problem with its execution at the next log in, you might have a problem with the automatical execution of `~/.bashrc`.

To handle it, you need to create/open `~/.bash_profile` and add following lines:
```
if [ -f ~/.bashrc ]; then
        . ~/.bashrc
fi
```
**1.2** **Create environment**

**1.2.1** If you need to create the folder with the Python environment (in case you do not have it) and install the packages, please, run:
```
cd /op3_v2
mamba create -f env.yaml -p path_to_env/venv --yes
```
**1.3** **Create symlinks**

In order to organize directories for scripts' outputs, you need to create folders in the desired locations for **data**; and then link **data** and **venv** folders to a script. 

To do that, please, execute:

**1.3.1** **For environment**
```
cd /op3_v2
ln -s path_to_env/venv venv
```
**1.3.2** **For data**
```
mkdir path_to_data/data
cd /op3_v2
ln -s path_to_data/data data
```

### 2. Set up the SLURM parameters for the pipeline script. 
If you want to run the dataset pseudobulking pipeline, e.g. for the Sci-plex dataset: `./pipelines/sciplex/sciplex_pseudobulking.sh`, you need to determine the parameters such as the number of cpus or the requested memory for the node under the commented lines by editing `#SBATCH` lines inside the script.

However, for Tahoe you need to change the arguments under `sbatch` command in the script `./pipelines/tahoe/tahoe_pseudobulking_parallel.sh`. 

For example:

```
sbatch -W -J pseudobulk_tahoe \
        -t 10:00:00 \
        -n 1 \
        --array=1-14 \
        --qos=cpu_normal \
        --mem=250G \
        --partition=cpu_p \
        --cpus-per-task=2 \
						...
```
as `sbatch` commands are executed on the different nodes for different plates.

Currently, scripts have pre-defined parameters.

### 3. Run scripts 

**3.1** To start a SLURM-based pseudobulking pipeline, you need to execute a wrapper script (`./run_pipelines/run_pseudobulking.sh`), which allow to run dataset-specific pipelines:

You can execute it with following arguments:
```
-s if this flag is included, then the script is executed for a subsample of a dataset, default=false
-h if this flag is included, then the help is printed, default=false
-d it is a required flag speifying which dataset should be processed: sciplex | tahoe ...
```

For example, for subsample of Sci-Plex dataset, you need to run:
```
cd /op3_v2
./run_pipelines/run_pseudobulking.sh -s -d sciplex
```

**3.2** To start a SLURM-based DGE (Differential Gene Expression) pipeline, you need to execute a wrapper script (`./run_pipelines/run_deg.sh`), which allows to run dataset-specific DGE pipelines:

You can execute it with following arguments:
```
-s if this flag is included, then the script is executed for a subsample of a dataset, default=false
-j if this flag is included, then the script is executed in parallel mode (array jobs per cell type), default=false
-h if this flag is included, then the help is printed, default=false
-f (value <int>) Min number of cells in pseudobulk to filter samples with the lower number, default=0 (no filtering)
-q if this flag is included, then samples that did not pass quality control are filtered out, default=false
-d it is a required flag specifying which dataset should be processed: sciplex | tahoe
-p it is a required flag specifying the parameter for DEG pipeline: group_all_replicates | separate_replicates
```
For example, for the subsample of Sci-Plex dataset with group_all_replicates parameter in parallel mode you need to run:
```
cd /op3_v2
./run_pipelines/run_deg.sh -s -j -d sciplex -p group_all_replicates
```

# Project structure:
The structure of the repo:
```
 tree .
.
├── data -> /lustre/groups/ml01/workspace/olga.novitskaia/data
├── env.yaml
├── logs
│   └── dataset_i
│       ├── full (or subsample)
│       │   └── deg
│       │       └── parameter_name (group_all_replicates or separate_replicates)
│       │           └── qc_true (or qc_false)
│       │               └── filter_min_cells_f
│       │                   └── deg_*.out, deg_*.err
│       └── ...
├── pipelines
│   ├── common
│   │   └── deg.sh
│   ├── sciplex
│   │   └─── configs/
│   │       ├── deg/
│   │       │   └── config.json
│   │       └── sciplex_pseudobulking.sh
│   └── tahoe
│       └─── configs/
│           ├── deg/
│           │   └── config.json
│           └── tahoe_pseudobulking_parallel.sh
├── README.md
├── requirements.txt
├── run_pipelines
│   ├── run_pseudobulking.sh
│   └── run_deg.sh
├── src
│   ├── configs
│   │   └── datasets.json
│   ├── deg
│   │   ├── run_deg.R
│   │   ├── run_processing_pseudobulk.py
│   │   └── subsampling.R
│   ├── downloading
│   │   └── run_downloading_datasets.py
│   ├── pseudobulking
│   │   ├── common
│   │   │   ├── pseudobulk.py
│   │   │   ├── run_combining_datasets.py
│   │   │   └── run_pseudobulking.py
│   │   └── datasets
│   │       ├── sciplex
│   │       │   ├── pubchem_imputation.py
│   │       │   └── standardization.py
│   │       └── tahoe
│   │           ├── pubchem_imputation.py
│   │           └── standardization.py
│   └── utils
│       └─── parsing_utils.py
└── venv -> /home/icb/olga.novitskaia/venv/
```

# Data
The structure of the `data` folder:
```
├──dataset_i
│   ├──deg_data
│   │   └──group_rep (or sep_rep)
│   │       ├──full
│   │       │   ├──qc_true
│   │       │   │   ├──filter_min_cells_0
│   │       │   │   │   ├──celltype1_de.h5ad
│   │       │   │   │   ├──celltype2_de.h5ad
│   │       │   │   │   └──...
│   │       │   │   ├──filter_min_cells_50
│   │       │   │   │   └──...
│   │       │   │   └──...
│   │       │   └──qc_false
│   │       │       └──...
│   │       └──subsampling
│   │           ├──qc_true
│   │           │   └──...
│   │           └──qc_false
│   │               └──filter_min_cells_0
│   │                   ├──celltype1_de.h5ad
│   │                   └──...
│   ├──raw
│   │   ├──dataset_i.h5ad
│   │   ├──dataset_i_subsample.h5ad (optionally)
│   │   ├──dataset_i_obs.parquet
│   │   └──dataset_i_var.parquet
│   ├──pseudobulk
│   │   ├──full
│   │   │   └──dataset_i.h5ad
│   │   └──subsample
│   │       └──dataset_i_subsample.h5ad
│   ├──pseudobulk_to_merge (optionally for tahoe)
│   │   ├──full
│   │   │   └──dataset_i.h5ad
│   │   └──subsample
│   │       └──dataset_i_subsample.h5ad
│   └──pseudobulk_processed
│       ├──group_rep (or sep_rep)
│       │   ├──dataset_i_processed.h5ad
│       │   └──by_celltype
│       │       ├──celltype1_processed.h5ad
│       │       ├──celltype2_processed.h5ad
│       │       └──...
│       └──...
└──...        

```
		
