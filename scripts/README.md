There are two options to run scripts: 
* inside the Docker container, for example, on the Amazon server
* on the HPC cluster with the workload manager (SLURM)

## 1. Docker

**NB!** It is recommended to run everything below under `tmux`, especially processing pipelines.

**NB!** Before running some bash scripts, you may need to run `chmod u+x ./path/script.sh` before their execution.

### 1.1 Setup the container:
**1.1.1** Go to the `docker` directory: 
```
cd /op3_v2/scripts/docker 
```
**1.1.2** Look at the `Dockerfile` and set the working directory, e.g.: 
```
ENV HOME=/home/ubuntu
```
**1.1.3** Run: 
```
docker-compose up -d to build a container
```
**1.1.4** To get into the container shell, run: 
```
docker exec -it pseudobulking bash
```

**1.1.5** To stop the container: 
```
docker stop pseudobulking
```

**1.1.6** To remove the container: 
```
docker rm pseudobulking
```

### 1.2 Set up Python environment and run scripts
**1.2.1** Once you have gotten into a container shell, you need to create the folder with the Python environment (in case you do not have it) and install the packages by running `set_env.sh`:
```
cd /op3_v2/scripts
./run_pipelines/set_env.sh
```
**1.2.2** To start a pipeline, you need to execute the bash script, which is located in the dataset-specific folder that you want to process. For example, running
```
cd /op3_v2/scripts
./pipelines/sciplex/sciplex_pseudobulking_subsample.sh 
```
you will execute a pseudobulking procedure for the subsample of the Sci-plex dataset.

**1.2.3** The script `run_inside_docker.sh` unites steps from **1.2.1** and **1.2.2**. You can edit/comment/uncomment lines to choose the dataset for pseudobulking and execute it as:
```
cd /op3_v2/scripts
./run_pipelines/run_inside_docker.sh
```
## 2. SLURM
### 2.1 Set up the SLURM parameters for the pipeline script. 
For instance, if you want to run `./pipelines/sciplex/sciplex_pseudobulking_subsample.sh`, you need to determine the parameters under the commented lines, for example, by editing `#SBATCH` lines inside it.
### 2.2. Set up a Python environment and run scripts. 
The set of steps is essentially the same as for **1.2**:

**2.2.1** Create the environment as it was described in **1.2.1**.

**2.2.2** To start a SLURM-based pipeline, you need to execute
```
cd /op3_v2/scripts
sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./pipelines/sciplex/sciplex_pseudobulking_subsample.sh
```
**2.2.3** The script `run_with_slurm.sh` unites steps from 2.2.1 and 2.2.2. You can edit/comment/uncomment lines to choose the dataset for pseudobulking and execute it as:
```
cd /op3_v2/scripts
./run_pipelines/run_with_slurm.sh
```

## Project structure:
The structure of the script folder:
```
tree scripts/
scripts/
├── docker
│   ├── Dockerfile
│   └── docker-compose.yml
├── pipelines
│   ├── sciplex
│   │   ├── configs
│   │   │   ├── sciplex_downloading.json
│   │   │   └── sciplex_pseudobulking.json
│   │   ├── sciplex_pseudobulking.sh
│   │   └── sciplex_pseudobulking_subsample.sh
│   └── tahoe
│       ├── configs
│       │   ├── tahoe_combining_datasets.json
│       │   ├── tahoe_downloading.json
│       │   └── tahoe_pseudobulking.json
│       ├── tahoe_pseudobulking.sh
│       └── tahoe_pseudobulking_subsample.sh
├── run_pipelines
│   ├── run_inside_docker.sh
│   ├── run_with_slurm.sh
│   └── set_env.sh
└── src
    ├── downloading
    │   └── run_downloading_datasets.py
    ├── pseudobulking
    │   ├── pseudobulk.py
    │   ├── run_combining_datasets.py
    │   ├── run_pseudobulking.py
    │   └── standartization.py
    └── utils
        └── parsing_utils.py
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
        ├──dataset_i.h5ad
        ├──dataset_i_subsample.h5ad
    ├──pseudobulk_to_merge (optionally for tahoe)
        ├──full
            ├──dataset_i.h5ad
        ├──subsample
            ├──dataset_i_subsample.h5ad
```
		
