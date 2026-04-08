# Description

**OP3_v2** is a set of pipelines for analyzing single-cell/bulk RNA sequencing data from perturbation experiments. 

The current version of OP3_v2 includes the scripts for processing and analyzing following datasets: 
- [**Sci-Plex3**](https://www.science.org/doi/10.1126/science.aax6234)
- [**Tahoe**](https://www.biorxiv.org/content/10.1101/2025.02.20.639398v1.full)
- [**L1000**](https://www.cell.com/cell/fulltext/S0092-8674(17)31309-0)
- [**OP3**](https://openreview.net/forum?id=WTI4RJYSVm)
- [**Novartis**](https://www.nature.com/articles/s41467-018-06500-x)

It consists of the following steps:

* **Pseudobulking**: 
Aggregates raw single-cell RNA-seq counts into pseudobulk samples, or standardizes bulk RNA-seq count data and enriches them with metadata for downstream analysis.

* **DEG Analysis**: 
Identifies differentially expressed genes between treatment and control conditions in pseudobulk/bulk samples

* **Enrich DEG `.var` (optional)**: 
After DEG, some `*_de.h5ad` objects may lose gene metadata (`symbol`, `is_merged`) when concatenated or re-saved. This step left-joins those columns from the processed pseudobulk reference and rewrites files in a form compatible with **aggregation** (`aggregating_deg.py`). Run only for the DEG output tree you care about (e.g. `results/` and matching `intermediate/`).

# License

This repository is a collection; each component retains its original license; our processed integration is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode.txt).

The licences of the datasets used in this project are provided by their source:

- **Sci-Plex** – [scPerturb Single-Cell Perturbation Data](https://zenodo.org/records/13350497) additionally annotated by [Laminlabs](https://lamin.ai/laminlabs/pertdata/artifacts?filter%5Band%5D%5B0%5D%5Bor%5D%5B0%5D%5Bbranch.name%5D%5Beq%5D=main&filter%5Band%5D%5B1%5D%5Bor%5D%5B0%5D%5Bis_latest%5D%5Beq%5D=true&filter%5Band%5D%5B2%5D%5Bor%5D%5B0%5D%5Bprojects.name%5D%5Beq%5D=scPerturb)
- **Tahoe** - [Arc Virtual Cell Atlas](https://arcinstitute.org/tools/virtualcellatlas) additionally annotated by [Laminlabs](https://lamin.ai/laminlabs/pertdata/artifacts?filter%5Band%5D%5B0%5D%5Bor%5D%5B0%5D%5Bbranch.name%5D%5Beq%5D=main&filter%5Band%5D%5B1%5D%5Bor%5D%5B0%5D%5Bis_latest%5D%5Beq%5D=true&filter%5Band%5D%5B2%5D%5Bor%5D%5B0%5D%5Bprojects.name%5D%5Beq%5D=Tahoe-100M)
- **L1000** - Contains data from GEO accessions [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) and [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138). Upstream terms apply; we do not assert CC BY for upstream L1000. See provenance + preprocessing notes. Datasets were additionally annotated by [Laminlabs](https://lamin.ai/laminlabs/pertdata/artifacts?filter%5Band%5D%5B0%5D%5Bor%5D%5B0%5D%5Bbranch.name%5D%5Beq%5D=main&filter%5Band%5D%5B1%5D%5Bor%5D%5B0%5D%5Bis_latest%5D%5Beq%5D=true&filter%5Band%5D%5B2%5D%5Bor%5D%5B0%5D%5Bprojects.name%5D%5Beq%5D=LINCS); small-molecule annotations from the [LINCS Data Portal](https://lincsportal.ccs.miami.edu/dcic-portal/#/terms) were used.
- **OP3** - [Open Problems Perturbation Prediction dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945)
- **Novartis** - [Novartis DRUG-seq MoABox Dataset](https://zenodo.org/records/14291446)

The terms of use/license for each database used for annotations are defined by its original data provider:

- **PubChem** - [NCBI](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- **NCBI datasets** - [NCBI](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- **Cellosaurus** - [Cellosaurus database](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- **OLS4** - [EMBL-EBI Data Resources and Tools](https://www.ebi.ac.uk/about/terms-of-use/)
- **HGNC** - [HGNC resources](https://www.genenames.org/about/license/)
- **Ensembl** - [Ensembl data](https://www.ensembl.org/info/about/legal/disclaimer.html)


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

**1.2.2** If you need to add environment variable (e.g. API_TOKEN), you can enter environment and then run:

```
mamba activate path_to_env/venv # Enter environment
conda env config vars set MY_API_KEY="SECRET_KEY" 
```
Or (which is better to hide key information)

```
mamba activate path_to_env/venv # Enter environment
conda env config vars set MY_API_KEY=$(cat secret_key.txt) #secret_key.txt contains your SECRET_KEY
```

To make your changes take effect please reactivate your environment:

```
mamba deactivate
mamba activate path_to_env/venv
echo $MY_API_KEY #If you want to print your MY_API_KEY
mamba deactivate
```

**NB!** Make sure that a path to your environment (as well as a path to your secret_key.txt file) is included in .gitignore, .cursorignore, .cursorindexingignore, .codexignore and .codexindexingignor to try not to allow to read/index API KEYS by Cursor/Codex. But still, there might be a leak of KEYS to AI-agents.

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
-d it is a required flag speifying which dataset should be processed: e.g. sciplex | tahoe | ...
```

For example, for subsample of Sci-Plex dataset, you need to run:
```
cd /op3_v2
./run_pipelines/run_pseudobulking.sh -s -d sciplex
```

**3.2** To start a SLURM-based DGE (Differential Gene Expression) pipeline, you need to execute a wrapper script (`./run_pipelines/run_deg.sh`), which allows to run dataset-specific DGE pipelines:

The DEG pipeline consists of three main steps:
1. **Preprocessing**: Processes pseudobulk data and adds perturbation labels (`run_processing_pseudobulk.py`)
2. **DGE Analysis**: Performs differential gene expression analysis using limma/edgeR (`run_deg.R`)
3. **Aggregation**: Aggregates batched results into final output files (`aggregating_deg.py`)

You can execute it with following arguments:
```
-s if this flag is included, then the script is executed for a subsample of a dataset, default=false
-j if this flag is included, then the script is executed in parallel mode (array jobs per cell type), default=false
-g Use GPU QOS/partition (gpu_normal/gpu_p), default=CPU
-h if this flag is included, then the help is printed, default=false
-f (value <int>) Min number of cells in pseudobulk to filter samples with the lower number, default=0 (no filtering)
-q if this flag is included, then samples that did not pass quality control are filtered out, default=false
-n if this flag is included, then the dataset is already normalized and normalization steps are skipped, default=false
-d it is a required flag specifying which dataset should be processed: sciplex | tahoe | ...
-p it is a required flag specifying the parameter for DEG pipeline: group_all_replicates | separate_replicates
```
**NB!** The -g option was added to make it possible to distribute DEG calculations across a larger number of nodes.

For example, for the subsample of Sci-Plex dataset with group_all_replicates parameter in parallel mode on the CPU nodes you need to run:
```
cd /op3_v2
./run_pipelines/run_deg.sh -s -j -d sciplex -p group_all_replicates
```

For a **normalized dataset** (e.g. L1000) with separate_replicates parameter on the CPU nodes:
```
cd /op3_v2
./run_pipelines/run_deg.sh -d l1000_phase1 -p separate_replicates -j -n
```

For a **normalized dataset** (e.g. L1000) with separate_replicates parameter on the GPU nodes:
```
cd /op3_v2
./run_pipelines/run_deg.sh -d l1000_phase1 -p separate_replicates -j -n -g
```

**3.3** **Enrich pseudobulk `.var` on DEG outputs (optional)**

Use `./run_pipelines/run_enrich_pseudobulk_var.sh` when you need `symbol` / `is_merged` back on `*_de.h5ad` before aggregation or downstream Python. It submits `./pipelines/common/enrich_pseudobulk_var.sh`, which runs **`sbatch -W`** (Slurm) and scans all `*_de.h5ad` under `-i` recursively.

Arguments:
```
-d  required: dataset name (sciplex | tahoe | l1000_phase1 | l1000_phase2 | op3 | novartis)
-i  required: root directory containing DEG *_de.h5ad files (e.g. .../deg_data/.../results or a parent tree)
-r  required: processed pseudobulk .h5ad used as reference (must contain .var columns symbol and is_merged)
--dry-run  log only; no file writes
-h  help
```

Example (paths depend on your layout and `group_rep` / `sep_rep`):
```
cd /op3_v2
./run_pipelines/run_enrich_pseudobulk_var.sh \
  -d novartis \
  -i ./data/novartis/deg_data/group_rep/full/qc_false/filter_min_cells_0/results \
  -r ./data/novartis/pseudobulk_processed/group_rep/novartis_standardized_processed.h5ad
```

Logs: `./logs/<dataset>/full/enrich_pseudobulk_var/` — launcher `enrich_<dataset>.PID*.out|.err` and Slurm `enrich_<dataset>.<jobid>.out|.err`.

### 4. Aggregate logs

After running pipelines, aggregate scattered log files into combined logs for easier review.

**Recursive mode** (all datasets):
```
cd /op3_v2
./run_pipelines/run_combining_logs.sh -r -l logs
```

**Recursive mode** (single dataset):
```
cd /op3_v2
./run_pipelines/run_combining_logs.sh -r -l logs/tahoe/
```

**Targeted mode** (specific directory):

Pseudobulk logs:
```
cd /op3_v2
./run_pipelines/run_combining_logs.sh -l logs/tahoe/full -t pseudobulk
```

DEG logs:
```
cd /op3_v2
./run_pipelines/run_combining_logs.sh -l logs/sciplex/full/deg/separate_replicates/qc_false/filter_min_cells_10 -t deg
```

Enrich (pseudobulk `.var`) logs:
```
cd /op3_v2
./run_pipelines/run_combining_logs.sh -l logs/novartis/full/enrich_pseudobulk_var -t enrich
```

### 5. Add datasets
To add a new dataset to the pipeline, follow the instructions in [`./docs/README.md`](https://github.com/theislab/op3_v2/blob/readme_new_dataset/docs/README.md)

# Project structure:
The structure of the repo:
```
 tree .
.
├── data -> /lustre/groups/ml01/workspace/olga.novitskaia/data
├── docs
├── env.yaml
├── logs
│   └── dataset_i
│       ├── full (or subsample)
│       │   ├── enrich_pseudobulk_var/
│       │   │   └── enrich_<dataset>.*.out, *.err
│       │   └── deg
│       │       └── parameter_name (group_all_replicates or separate_replicates)
│       │           └── qc_true (or qc_false)
│       │               └── filter_min_cells_f
│       │                   ├── preprocessing/
│       │                   │   └── deg_processing_pseudobulk.*.out, *.err
│       │                   ├── deg/
│       │                   │   ├── celltype1/
│       │                   │   │   └── deg_analysis.*.out, *.err
│       │                   │   ├── celltype2/
│       │                   │   │   └── deg_analysis.*.out, *.err
│       │                   │   └── ...
│       │                   └── aggregation/
│       │                       └── deg_aggregation.*.out, *.err
│       └── ...
├── pipelines
│   ├── common
|   |   ├── combining_logs.sh
│   │   ├── deg.sh
│   │   └── enrich_pseudobulk_var.sh
│   ├── dataset_i
│   │   └─── configs/
│   │       ├── deg/
│   │       │   └── config.json
│   │       └── dataset_i_pseudobulking.sh
│   └── ...
├── README.md
├── run_pipelines
│   ├── run_combining_logs.sh
│   ├── run_deg.sh
│   ├── run_enrich_pseudobulk_var.sh
│   └── run_pseudobulking.sh
├── LICENSE
├── src
│   ├── configs
│   │   └── datasets.json
│   ├── deg
│   │   ├── run_deg.R
│   │   ├── run_processing_pseudobulk.py
│   │   ├── aggregating_deg.py
│   │   ├── enrich_pseudobulk_var_metadata.py
│   │   ├── ensembl_mapping.py
│   │   └── subsampling.R
│   ├── downloading
│   │   └── run_downloading_datasets.py
│   ├── pseudobulking
│   │   ├── common
│   │   │   ├── pseudobulk.py
│   │   │   ├── pubchem.py
│   │   │   ├── run_combining_datasets.py
│   │   │   └── run_pseudobulking.py
│   │   └── datasets
│   │       ├── dataset_i
│   │       │   ├── pubchem_imputation.py
│   │       │   └── ...
│   │       └── ...
│   └── utils
│       ├──  aggregate_logs.py
│       └─── parsing_utils.py
└── venv -> /lustre/groups/ml01/workspace/olga.novitskaia/venv
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
│   │       │   │   │   ├──intermediate/
│   │       │   │   │   │   ├──celltype1/
│   │       │   │   │   │   │   ├──batch_1.h5ad
│   │       │   │   │   │   │   ├──batch_2.h5ad
│   │       │   │   │   │   │   └──...
│   │       │   │   │   │   ├──celltype2/
│   │       │   │   │   │   │   └──...
│   │       │   │   │   │   └──...
│   │       │   │   │   └──results/
│   │       │   │   │       ├──celltype1_de.h5ad
│   │       │   │   │       ├──celltype2_de.h5ad
│   │       │   │   │       └──...
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
│   │                   ├──intermediate/
│   │                   │   └──...
│   │                   └──results/
│   │                       ├──celltype1_de.h5ad
│   │                       └──...
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
│   ├──pseudobulk_to_merge (optionally)
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
