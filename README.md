You can run scripts on the HPC cluster with the workload manager (SLURM)

**NB!** Before running some bash scripts, you may need to run `chmod u+x ./path/script.sh` before their execution.
### 1. Set up prerequisites
Before executing scripts we need to install packages into the certain environment. Here we provide the illustrative scripts for using a package manager called **Mamba**. Also you can adapt our scripts to use other managers such as Micromamba or Conda.

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
mkdir path_to_env/.venv
cd /op3_v2
mamba create -f env.yaml -p path_to_env/.venv --yes
```
**1.3** **Create symlinks**

In order to organize directories for scripts' outputs, you need to create folders in the desired locations for **logs** and **data**; and then link **logs**, **data** and **.venv** folders to a script. 

To do that, please, execute:

**1.3.1** **For environment**
```
cd /op3_v2
ln -s path_to_env/.venv .venv
```
**1.3.2** **For logs**
```
mkdir path_to_logs/logs
cd /op3_v2
ln -s path_to_logs/logs logs
```
**1.3.3** **For data**
```
mkdir path_to_data/data
cd /op3_v2
ln -s path_to_data/data data
```

### 2. Set up the SLURM parameters for the pipeline script. 
If you want to run the dataset pseudobulking pipeline, e.g. the Sci-plex subsample pseudobulking: `./pipelines/sciplex/sciplex_pseudobulking.sh`, you need to determine the parameters such as the number of cpus or the requested memory for the node under the commented lines by editing `#SBATCH` lines inside the script.

However, for Tahoe you need to change the arguments under `sbatch` command in the script `./pipelines/tahoe/tahoe_pseudobulking_parallel.sh`:

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

### 3. Set up a Python environment and run scripts 

**3.1** To start a SLURM-based pipeline, you need to execute

**For Sci-plex**
```
cd /op3_v2
sbatch -e ./logs/pseudobulk_sciplex.%j.err -o ./logs/pseudobulk_sciplex.%j.out ./pipelines/sciplex/sciplex_pseudobulking.sh
```

**For Tahoe**
```
cd /op3_v2
./pipelines/tahoe/tahoe_pseudobulking_parallel.sh > ./logs/pseudobulk_tahoe.out 2>./logs/pseudobulk_tahoe.err &
```
**NB!** If you want to run on the subsample you need to add `subsampling` right after the bash script. For example:
```
./pipelines/tahoe/tahoe_pseudobulking_parallel.sh subsampling > ~/logs/pseudobulk_tahoe.out 2>~/logs/pseudobulk_tahoe.err &
```

**3.2** The script `run_with_slurm.sh` includes two versions of pipelines for datasets, full and subsampling for code debugging. You can edit/comment/uncomment lines to choose the dataset for pseudobulking and execute it as:
```
cd /op3_v2
./run_pipelines/run_with_slurm.sh
```

## Project structure:
The structure of the repo:
```
tree .
.
├── pipelines
│   ├── sciplex
│   │   └── sciplex_pseudobulking.sh
│   └── tahoe
│       └── tahoe_pseudobulking_parallel.sh
├── README.md
├── requirements.txt
├── run_pipelines
│   ├── run_with_slurm.sh
│   └── set_env.sh
└── src
    ├── configs
    │   └── datasets.json
    ├── downloading
    │   └── run_downloading_datasets.py
    ├── pseudobulking
    │   ├── common
    │   │   ├── pseudobulk.py
    │   │   ├── run_combining_datasets.py
    │   │   └── run_pseudobulking.py
    │   └── datasets
    │       ├── sciplex
    │       │   ├── pubchem_imputation.py
    │       │   └── standardization.py
    │       └── tahoe
    │           ├── pubchem_imputation.py
    │           └── standardization.py
    └── utils
        ├── parsing_utils.py
```

## Data
The structure of the data folder:
```
├──dataset_i
    ├──raw
        ├──dataset_i.h5ad
        ├──dataset_i_subsample.h5ad (optionally)
        ├──dataset_i_obs.parquet
        ├──dataset_i_var.parquet
    ├──pseudobulk
        ├──full
            ├──dataset_i.h5ad
        ├──subsample
            ├──dataset_i_subsample.h5ad
    ├──pseudobulk_to_merge (optionally for tahoe)
        ├──full
            ├──dataset_i.h5ad
        ├──subsample
            ├──dataset_i_subsample.h5ad
```
		
