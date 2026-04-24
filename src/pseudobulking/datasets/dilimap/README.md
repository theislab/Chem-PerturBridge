# DILImap Training Dataset

This module handles two variants of the DILImap dataset for the pseudobulking pipeline:

- **dilimap_train**: Training data only
- **dilimap_train_val**: Training + validation data combined

Keeping them separate prevents test-set contamination when computing DGE signatures.

## Data provenance

### Training data

The training dataset is not publicly available. It was obtained by request from the
DILImap authors (see [paper](https://www.nature.com/articles/s41467-025-65690-3),
section "Data availability").

The file is shared via Google Drive. Place it at:

```
Chem-PerturBridge/data/adata_training_counts.h5ad
```

### Validation data

The validation dataset is publicly available via the DILImap S3 bucket:

```python
import dilimap as dmap
adata = dmap.s3.read('validation_data_counts.h5ad', package_name='public/data')
```

Place it at:

```
Chem-PerturBridge/data/adata_validation_counts.h5ad
```

## Running the pipeline

```bash
# Training data only
./run_pipelines/run_pseudobulking.sh -d dilimap_train

# Training + validation combined
./run_pipelines/run_pseudobulking.sh -d dilimap_train_val
```
