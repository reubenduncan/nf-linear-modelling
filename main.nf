#!/usr/bin/env nextflow
// main.nf — Linear Modelling of Microbiome Data
// Processes: CODA_GLMNET, SUBSET_REGRESSION (from CSV or metadata)

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// PROCESS: CODA-GLMNET penalized regression
// ---------------------------------------------------------------------------
process CODA_GLMNET {
    tag "${params.label}"
    publishDir "${params.output_dir}/coda_glmnet", mode: 'copy'

    input:
    path feature_table
    path meta_table
    path taxonomy_table   // may be an empty sentinel file

    output:
    path "*.csv", emit: coda_csvs

    script:
    def tax_arg     = (params.taxonomy_table && params.taxonomy_table != '')
                      ? "--taxonomy_table ${taxonomy_table}"
                      : ''
    def excl_col    = params.exclude_column    ? "--exclude_column '${params.exclude_column}'" : ''
    def excl_vals   = params.exclude_values    ? "--exclude_values '${params.exclude_values}'" : ''
    def grp_col     = params.groups_column     ? "--groups_column '${params.groups_column}'" : ''
    def grp_paste   = params.coda_groups_paste_columns
                      ? "--groups_paste_columns '${params.coda_groups_paste_columns}'" : ''
    def global_flag = params.covariates_global ? '--covariates_global' : ''
    def binary_flag = params.coda_binary_outcome ? '--coda_binary_outcome' : ''

    """
    Rscript ${params.scripts_dir}/src/R/lm_coda_glmnet.R \\
        --feature_table       "${feature_table}" \\
        --meta_table          "${meta_table}" \\
        --input_format        "${params.input_format}" \\
        --output_dir          "." \\
        --which_level         "${params.coda_which_level}" \\
        --label               "${params.label}" \\
        --min_library_size    ${params.min_library_size} \\
        --normalisation_method "${params.coda_normalisation}" \\
        --occupancy_threshold  ${params.coda_occupancy} \\
        --top_n               ${params.coda_top_n} \\
        --coda_lambda         "${params.coda_lambda}" \\
        --coda_alpha          ${params.coda_alpha} \\
        --environmental_covariates "${params.environmental_covariates}" \\
        ${tax_arg} \\
        ${excl_col} \\
        ${excl_vals} \\
        ${grp_col} \\
        ${grp_paste} \\
        ${global_flag} \\
        ${binary_flag}
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Subset regression — dependent variable from a CSV file
// ---------------------------------------------------------------------------
process SUBSET_REGRESSION {
    tag "${params.label}"
    publishDir "${params.output_dir}/subset_regression", mode: 'copy'

    input:
    path meta_table
    path dependent_csv

    output:
    path "*.csv", emit: subset_csvs

    script:
    def scale_pred = params.scale_predictors ? 'TRUE' : 'FALSE'
    def big_flag   = params.really_big       ? '--really_big' : ''

    """
    Rscript ${params.scripts_dir}/src/R/subset_regression.R \\
        --meta_table          "${meta_table}" \\
        --dependent_csv       "${dependent_csv}" \\
        --dependent_source    "csv" \\
        --output_dir          "." \\
        --label               "${params.label}" \\
        --dependent_variable  "${params.dependent_variable}" \\
        --explanatory_variables "${params.explanatory_variables}" \\
        --regression_method   "${params.regression_method}" \\
        --test_method         "${params.test_method}" \\
        --cv_folds            ${params.cv_folds} \\
        --scale_predictors    ${scale_pred} \\
        ${big_flag}
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Subset regression — dependent variable from metadata column
// ---------------------------------------------------------------------------
process SUBSET_REGRESSION_FROM_META {
    tag "${params.label}"
    publishDir "${params.output_dir}/subset_regression", mode: 'copy'

    input:
    path meta_table

    output:
    path "*.csv", emit: subset_csvs

    script:
    def scale_pred = params.scale_predictors ? 'TRUE' : 'FALSE'
    def big_flag   = params.really_big       ? '--really_big' : ''

    """
    Rscript ${params.scripts_dir}/src/R/subset_regression.R \\
        --meta_table          "${meta_table}" \\
        --dependent_source    "metadata" \\
        --output_dir          "." \\
        --label               "${params.label}" \\
        --dependent_variable  "${params.dependent_variable}" \\
        --explanatory_variables "${params.explanatory_variables}" \\
        --regression_method   "${params.regression_method}" \\
        --test_method         "${params.test_method}" \\
        --cv_folds            ${params.cv_folds} \\
        --scale_predictors    ${scale_pred} \\
        ${big_flag}
    """
}

// ---------------------------------------------------------------------------
// PROCESS: Merge all CSVs into a single Parquet file
// ---------------------------------------------------------------------------
process MERGE_PARQUET {
    tag "${params.label}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path csvs

    output:
    path "linear_modelling_${params.label}.parquet"

    script:
    """
    Rscript ${params.scripts_dir}/src/R/merge_parquet.R \\
        --label      '${params.label}' \\
        --output_dir '.'
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    // Validate required parameters
    if (!params.environmental_covariates) {
        error "ERROR: params.environmental_covariates must be set (comma-sep metadata columns)"
    }
    if (!params.explanatory_variables) {
        error "ERROR: params.explanatory_variables must be set (comma-sep metadata columns)"
    }

    feat_ch = Channel.fromPath(params.feature_table, checkIfExists: true)
    meta_ch = Channel.fromPath(params.meta_table,    checkIfExists: true)

    // Taxonomy table: emit a sentinel empty file if not provided
    if (params.taxonomy_table && params.taxonomy_table != '') {
        tax_ch = Channel.fromPath(params.taxonomy_table, checkIfExists: true)
    } else {
        tax_ch = Channel.fromPath("${projectDir}/NO_TAXONOMY_FILE", checkIfExists: false)
                        .ifEmpty(file("${projectDir}/.no_taxonomy"))
    }

    // CODA-GLMNET always runs
    coda_out = CODA_GLMNET(feat_ch, meta_ch, tax_ch)

    // Subset regression: depends on dependent_source param
    if (params.dependent_csv && params.dependent_csv != '') {
        dep_ch     = Channel.fromPath(params.dependent_csv, checkIfExists: true)
        subreg_out = SUBSET_REGRESSION(meta_ch, dep_ch)
    } else {
        subreg_out = SUBSET_REGRESSION_FROM_META(meta_ch)
    }

    if (params.merge_parquet) {
        all_csvs = coda_out.coda_csvs.flatten()
            .mix(subreg_out.subset_csvs.flatten())
            .collect()
        MERGE_PARQUET(all_csvs)
    }
}
