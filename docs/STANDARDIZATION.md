## How to standardize?

You might have some questions on how to standardize diverse datasets.
Therefore, we added some comments on potential issues to help you avoid them.

To standardize you just need to convert dataset column content or column names to the specific schema described [here](https://github.com/theislab/op3_v2/blob/main/docs/Format_Pseudobulk.ipynb).

* For example, your dataset probably has columns `PLATE_ID` and `PLATE_NAME`, both could be converted to `plate` column from our schema. How to choose? Probably, take a look if authors of the dataset have some preferences by using content of `PLATE_ID` and `PLATE_NAME` in other columns or looking at content readability (could numbers in `PLATE_ID` be confused with other IDs, e.g. `pubchem_cid`?).

* Additionally, pay attention to differences in units between the original version and the standardized one, just to convert them correctly.

## Which columns are important?
The most important columns to keep an eye on are `plate`, `well`, `cell_type`, `perturbagen`, `pert_type`, `is_control`, `pert_dose_uM`, `pert_time_h`, `pubchem_cid`. Description of them is located [here](https://github.com/theislab/op3_v2/blob/main/docs/Format_Pseudobulk.ipynb).

**NB!** Some cell types are not represented as cell lines. It is ok to leave them with the original names/ids, but also if you have a chance it would be great to map them to **CL:** IDs via [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl). In our columns we keep it in the format: `CL_` instead of `CL:`.

**NB!** Additionally, check that the `perturbagen` column has no duplicates in compound names with slightly different chemical structures. If so, this case should be analysed individually, but it is similar to the further preprocessing run for the L1000 dataset before DEG (we added [`combine_perturbagen_pubchem_cid` option](https://github.com/theislab/op3_v2/blob/44ad730e8a8341d0be6d4838922e5529ca4be403/src/deg/run_processing_pseudobulk.py#L1171) to distinguish compounds with the same names but different structures in DEG analysis).

## Which columns are not important?

Some columns, especially related to the information regarding **donors**, might not be so highly important for future analysis and therefore not included in the original dataset, but we included this information into our schema to store details as well. **Donor** information connected to the specific cell line could be found on `www.cellosaurus.org` website, or in the publication associated with the dataset release.

## Data
### Single-Cell data
Examples: `sciplex3`, `tahoe`

For Single-Cell data pseudobulk aggregation is needed to run. Dataset-specific standardization is executed during [`pseudobulk.py`](https://github.com/theislab/op3_v2/blob/main/src/pseudobulking/common/pseudobulk.py) script running. To include standardization scripts into `pseudobulk.py` pipeline, and make them imported you need to:

* include `standardization.py` scripts in dataset-specific directory `op3_v2/src/pseudobulking/datasets/mynewdataset`
* add the path to `standardization` module in the `config` file: `/op3_v2/src/configs/datasets.json`
* Additionally, `pseudobulk.py` includes `pubchem_cid` enrichment function to map drug names to PubChem CIDs. Note that automatic enrichment of CIDs does not always find CIDs for drug names, therefore it might be worth adding a dictionary of manual mappings for unmapped drug names (see [`pubchem_imputation.py`](https://github.com/theislab/op3_v2/blob/main/src/pseudobulking/datasets/sciplex/pubchem_imputation.py) for **Sci-Plex**) and specify it in the mentioned `config` file: `/op3_v2/src/configs/datasets.json`.

### Bulk data
Bulk data also might include both normalized values and raw counts.
Examples: `l1000`, `op3`*

*originally `op3` is a single cell dataset, but we took a prepared pseudobulk version for the current pipeline.

For bulk data you probably have more space for creating different scripts with different namings inside the dataset-specific directory, e.g. [L1000](https://github.com/theislab/op3_v2/tree/main/src/pseudobulking/datasets/l1000). But the resulting dataset schema should be similar to the one described [here](https://github.com/theislab/op3_v2/blob/main/docs/Format_Pseudobulk.ipynb), slight differences might be possible.

Note, to look up `pubchem_cid` for bulk data you can include `lookup_pubchem_cids` from [`pubchem.py`](https://github.com/theislab/op3_v2/blob/main/src/pseudobulking/common/pubchem.py) module into your processing/assembling pipeline.

**NB!** Some columns that might be represented in pseudobulk data are outputs of the `decoupler` library, such as `psbulk_cells` and `psbulk_counts`. If we do not know information on the number of cells or counts in the bulk sample - it is better to fill these columns with -666. For intensities or normalized data we cannot determine the number of counts in the sample, therefore `psbulk_counts` should probably be filled with -666; for bulk data with raw counts we can calculate the values for this column as `.X.sum(1)`.

**NB!** For single-cell data we also included some filtration to get rid of noisy signal. This filtration is already included into our `pseudobulk.py` pipeline. However, bulked data might have different recommendations for filtration, and these recommendations might be found in the papers associated with released datasets. So, it is better to check the publication if you have some concerns about what to do with the filtration step.
