You can run scripts on the HPC cluster with the workload manager (SLURM)

**NB!** Before running some bash scripts, you may need to run `chmod u+x ./path/script.sh` before their execution.

### 1.1 Set up the SLURM parameters for the pipeline script. 
If you want to run the dataset pseudobulking pipeline, e.g. the Sci-plex subsample pseudobulking: `./pipelines/sciplex/sciplex_pseudobulking_subsample.sh`, you need to determine the parameters under the commented lines, for example, by editing `#SBATCH` lines inside the script.

However, for Tahoe you need to change the arguments under `sbatch` command in the script `./pipelines/tahoe/tahoe_pseudobulking_subsample_parallel.sh`:

```
sbatch -W -J pseudobulk_tahoe \
                        -t 10:00:00 \
                        -n 1 \
                        --qos=cpu_normal \
                        --mem=250G \
                        --partition=cpu_p \
                        --cpus-per-task=2 \
						...
```
as `sbatch` commands are executed on the different nodes for different plates.

### 1.2. Set up a Python environment and run scripts 

**1.2.1** Firstly, you need to create the folder with the Python environment (in case you do not have it) and install the packages by running `set_env.sh`:
```
cd /op3_v2
./run_pipelines/set_env.sh
```
**1.2.2** To start a SLURM-based pipeline, you need to execute

**For Sci-plex**
```
mkdir -p ~/logs
cd /op3_v2
sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./pipelines/sciplex/sciplex_pseudobulking.sh
```

**For Tahoe**
```
mkdir -p ~/logs
cd /op3_v2
./pipelines/tahoe/tahoe_pseudobulking_parallel.sh > ~/logs/pseudobulk_tahoe.out 2>~/logs/pseudobulk_tahoe.err &
```
**1.2.3** The script `run_with_slurm.sh` unites steps from **1.2.1** and **1.2.2** and includes two versions of pipelines for datasets, full and subsample for code debugging. You can edit/comment/uncomment lines to choose the dataset for pseudobulking and execute it as:
```
cd /op3_v2
./run_pipelines/run_with_slurm.sh
```

## Project structure:
The structure of the script folder:
```
tree .
├── pipelines
│   ├── sciplex
│   │   ├── configs
│   │   │   ├── sciplex_downloading.json
│   │   │   └── sciplex_pseudobulking.json
│   │   ├── sciplex_pseudobulking.sh
│   │   └── sciplex_pseudobulking_subsample.sh
│   └── tahoe
│       ├── configs
│       │   ├── tahoe_combining_datasets.json
│       │   ├── tahoe_downloading.json
│       │   └── tahoe_pseudobulking.json
│       ├── tahoe_pseudobulking_parallel.sh
│       └── tahoe_pseudobulking_subsample_parallel.sh
├── README.md
├── requirements.txt
├── run_pipelines
│   ├── run_with_slurm.sh
│   └── set_env.sh
└── src
    ├── downloading
    │   └── run_downloading_datasets.py
    ├── pseudobulking
    │   ├── pseudobulk.py
    │   ├── pubchem_imputation.py
    │   ├── run_combining_datasets.py
    │   ├── run_pseudobulking.py
    │   └── standardization.py
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
		
