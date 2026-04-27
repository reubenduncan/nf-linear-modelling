# linear-modelling

A Nextflow pipeline for linear modelling of microbiome data, combining CODA-GLMNET penalized regression and best-subset regression.

## Introduction

This pipeline runs two complementary modelling approaches:

- **CODA-GLMNET:** Fits penalized (LASSO/ridge/elastic net) regression models of environmental covariates on log-ratio transformed microbiome compositions. A separate model is fitted per group × covariate combination, identifying which taxa are most predictive of each environmental variable.
- **Subset regression:** Selects the best linear model of a dependent variable (e.g. an alpha diversity metric or any metadata column) from a set of candidate environmental predictors using exhaustive, forward, backward, or sequential replacement search (via `leaps`), with model selection by cross-validation, AIC, or BIC.

Both analyses run in parallel. An optional `--merge_parquet` flag consolidates all outputs into a single Parquet file.

## Quick start

```bash
nextflow run main.nf \
  --feature_table            /path/to/feature_table.biom \
  --meta_table               /path/to/meta_table.csv \
  --groups_column            Treatment \
  --environmental_covariates pH,Temperature \
  --explanatory_variables    pH,Temperature,Conductivity \
  --dependent_variable       Shannon \
  --dependent_source         csv \
  --dependent_csv            /path/to/alpha_diversity.csv \
  --label                    my_analysis
```

## Parameters

### Input / Output

| Parameter | Default | Description |
|---|---|---|
| `--feature_table` | *(required)* | Path to feature table (BIOM, TSV, or GTDB format) |
| `--meta_table` | *(required)* | Path to sample metadata CSV (first column = sample IDs) |
| `--taxonomy_table` | `""` | Taxonomy TSV (required for `tsv`/`gtdb` input formats) |
| `--input_format` | `biom` | `biom` \| `tsv` \| `gtdb` |
| `--output_dir` | `results/` | Directory for output files |
| `--label` | `analysis` | Label prepended to all output file names |

### Filtering

| Parameter | Default | Description |
|---|---|---|
| `--min_library_size` | `5000` | Minimum per-sample read depth; samples below this are dropped |
| `--exclude_column` | `""` | Metadata column used to identify samples for exclusion |
| `--exclude_values` | `""` | Comma-separated values in `exclude_column` to remove |

### Grouping

| Parameter | Default | Description |
|---|---|---|
| `--groups_column` | `""` | Metadata column defining groups for per-group CODA-GLMNET models |
| `--groups_paste_columns` | `""` | Comma-separated columns pasted together to form groups |

### CODA-GLMNET

| Parameter | Default | Description |
|---|---|---|
| `--environmental_covariates` | *(required)* | Comma-separated metadata columns to use as response variables |
| `--coda_taxon_rank` | `Feature` | Taxonomic level for feature aggregation |
| `--coda_normalisation` | `TSS+CLR` | Normalisation: `TSS+CLR` \| `TSS+ILR` \| `logrelative` |
| `--coda_occupancy` | `1` | Minimum number of samples a feature must appear in |
| `--coda_top_n` | `100` | Retain only the top N most abundant features |
| `--coda_lambda` | `lambda.min` | Regularisation selection: `lambda.min` \| `lambda.1se` |
| `--coda_alpha` | `1.0` | GLMNET mixing: `0` = ridge, `1` = lasso, `0.5` = elastic net |
| `--coda_binary_outcome` | `false` | Auto-detect binary covariates and fit logistic CODA-GLMNET |
| `--coda_groups_paste_columns` | `""` | CODA-specific override for groups pasting |
| `--covariates_global` | `false` | Fit a single global model across all samples rather than per group |

### Subset regression

| Parameter | Default | Description |
|---|---|---|
| `--explanatory_variables` | *(required)* | Comma-separated metadata columns used as candidate predictors |
| `--dependent_variable` | `Shannon` | Column name of the response variable |
| `--dependent_source` | `csv` | Source of the dependent variable: `csv` \| `metadata` |
| `--dependent_csv` | `""` | Path to CSV containing the dependent variable (e.g. alpha diversity output). Required when `dependent_source=csv` |
| `--regression_method` | `forward` | Subset search algorithm: `exhaustive` \| `backward` \| `forward` \| `seqrep` |
| `--test_method` | `cv` | Model selection criterion: `cv` \| `aic` \| `bic` |
| `--cv_folds` | `5` | Number of cross-validation folds (used when `test_method=cv`) |
| `--scale_predictors` | `true` | Standardise predictors before regression |
| `--really_big` | `false` | Pass `really.big=TRUE` to `regsubsets()` for large predictor sets with exhaustive search |

### Output options

| Parameter | Default | Description |
|---|---|---|
| `--merge_parquet` | `false` | Merge all output CSVs into a single Parquet file |

## Outputs

Files are written to subdirectories of `--output_dir`.

### CODA-GLMNET (`coda_glmnet/`)

| File | Description |
|---|---|
| `CODA_coefficients_{label}_{group}_{covariate}.csv` | Non-zero GLMNET coefficients per group × covariate model |
| `CODA_predictions_{label}_{group}_{covariate}.csv` | Predicted vs observed values per group × covariate model |
| `CODA_model_summary_{label}.csv` | Cross-validation performance summary for all models |

### Subset regression (`subset_regression/`)

| File | Description |
|---|---|
| `SubsetReg_best_model_{label}.csv` | Coefficients of the selected best model |
| `SubsetReg_cv_errors_{label}.csv` | Cross-validation error by model size |
| `SubsetReg_all_models_{label}.csv` | All candidate models with their selection criteria |

### Merged output

| File | Description |
|---|---|
| `linear_modelling_{label}.parquet` | All CSVs merged with `analysis` and `table` metadata columns (`--merge_parquet` only) |

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (default executor — environment built automatically from `environment.yml`)
- **or** Docker with `-profile docker`
- **or** Singularity with `-profile singularity`
- **or** a local R installation with: `optparse`, `leaps`, `dplyr`, `caret`, `purrr`, `stringr`, `coda4microbiome`, `phyloseq`, `mixOmics`, `arrow`

## Running with a local R installation

Add `-profile` to select your execution environment (conda is used by default if no profile is specified):

```bash
nextflow run main.nf \
  -c nextflow.config \
  --feature_table            /path/to/table.biom \
  --meta_table               /path/to/meta.csv \
  --groups_column            Treatment \
  --environmental_covariates pH,Temperature \
  --explanatory_variables    pH,Temperature,Conductivity \
  --dependent_variable       pH \
  --dependent_source         metadata \
  --label                    my_analysis
```

Available profiles: `conda` (default), `docker`, `singularity`.