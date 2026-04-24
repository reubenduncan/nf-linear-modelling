#!/usr/bin/env Rscript
# lm_coda_glmnet.R
# CODA-GLMNET penalized regression on microbiome feature tables.
# Outputs CSV files only; no plotting.

suppressPackageStartupMessages({
  library(optparse)
  library(coda4microbiome)
  library(phyloseq)
  library(mixOmics)
  library(stringr)
})

# ---------------------------------------------------------------------------
# CLI options
# ---------------------------------------------------------------------------
option_list <- list(
  # Input
  make_option(c("--feature_table"),          type = "character", default = NULL,
              help = "Path to feature table (BIOM, TSV, or GTDB). Alias: --biom_file"),
  make_option(c("--biom_file"),              type = "character", default = NULL,
              help = "Alias for --feature_table (kept for backwards compatibility)"),
  make_option(c("--meta_table"),             type = "character", default = NULL,
              help = "Path to metadata CSV (samples x variables)"),
  make_option(c("--taxonomy_table"),         type = "character", default = "",
              help = "Path to taxonomy TSV (required for tsv/gtdb formats)"),
  make_option(c("--input_format"),           type = "character", default = "biom",
              help = "Feature table format: biom | tsv | gtdb [default: biom]"),
  make_option(c("--output_dir"),             type = "character", default = ".",
              help = "Directory for output CSVs [default: .]"),
  # Filtering / grouping
  make_option(c("--which_level"),            type = "character", default = "Otus",
              help = "Taxonomy level to aggregate to [default: Otus]"),
  make_option(c("--label"),                  type = "character", default = "Hypothesis1",
              help = "Label prefix for output files [default: Hypothesis1]"),
  make_option(c("--min_library_size"),       type = "integer",   default = 5000L,
              help = "Minimum read depth per sample [default: 5000]"),
  make_option(c("--exclude_column"),         type = "character", default = "",
              help = "Metadata column to filter samples on"),
  make_option(c("--exclude_values"),         type = "character", default = "",
              help = "Comma-separated values to exclude from exclude_column"),
  make_option(c("--groups_column"),          type = "character", default = "",
              help = "Metadata column defining groups for per-group models"),
  make_option(c("--groups_paste_columns"),   type = "character", default = "",
              help = "Comma-sep columns to paste together to form a Groups column"),
  make_option(c("--normalisation_method"),   type = "character", default = "TSS+CLR",
              help = "Normalisation: logrelative | TSS+ILR | TSS+CLR [default: TSS+CLR]"),
  make_option(c("--occupancy_threshold"),    type = "integer",   default = 1L,
              help = "Min samples a feature must be present in [default: 1]"),
  make_option(c("--top_n"),                  type = "integer",   default = 100L,
              help = "Retain top N most abundant features [default: 100]"),
  make_option(c("--environmental_covariates"), type = "character", default = "",
              help = "Comma-sep metadata columns to model as response variables (required)"),
  # CODA-GLMNET tuning
  make_option(c("--coda_lambda"),            type = "character", default = "lambda.min",
              help = "GLMNET lambda: lambda.min | lambda.1se [default: lambda.min]"),
  make_option(c("--coda_alpha"),             type = "double",    default = 1.0,
              help = "GLMNET alpha: 0=ridge, 1=lasso, 0.5=elastic net [default: 1]"),
  make_option(c("--coda_binary_outcome"),    action = "store_true", default = FALSE,
              help = "If set, auto-detect binary covariates (else all treated as continuous)"),
  make_option(c("--covariates_global"),      action = "store_true", default = FALSE,
              help = "If set, run ONE model across all groups combined instead of per-group")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --biom_file alias
if (is.null(opt$feature_table) && !is.null(opt$biom_file)) {
  opt$feature_table <- opt$biom_file
}

# ---------------------------------------------------------------------------
# Validate required inputs
# ---------------------------------------------------------------------------
if (is.null(opt$feature_table) || !nzchar(opt$feature_table)) {
  stop("[lm_coda_glmnet] --feature_table (or --biom_file) is required")
}
if (is.null(opt$meta_table) || !nzchar(opt$meta_table)) {
  stop("[lm_coda_glmnet] --meta_table is required")
}
if (!file.exists(opt$feature_table)) {
  stop("[lm_coda_glmnet] Feature table not found: ", opt$feature_table)
}
if (!file.exists(opt$meta_table)) {
  stop("[lm_coda_glmnet] Metadata table not found: ", opt$meta_table)
}
if (!opt$input_format %in% c("biom", "tsv", "gtdb")) {
  stop("[lm_coda_glmnet] --input_format must be one of: biom, tsv, gtdb")
}
if (!opt$coda_lambda %in% c("lambda.min", "lambda.1se")) {
  stop("[lm_coda_glmnet] --coda_lambda must be: lambda.min or lambda.1se")
}
if (opt$coda_alpha < 0 || opt$coda_alpha > 1) {
  stop("[lm_coda_glmnet] --coda_alpha must be between 0 and 1")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Source loader
# ---------------------------------------------------------------------------
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(script_dir) || script_dir == ".") script_dir <- "."
loader_path <- file.path(script_dir, "load_feature_table.R")
if (!file.exists(loader_path)) {
  # fallback: look relative to opt$output_dir or /opt/ecology-scripts
  candidates <- c(
    file.path(dirname(opt$output_dir), "src", "R", "load_feature_table.R"),
    "/opt/ecology-scripts/src/R/load_feature_table.R",
    "src/R/load_feature_table.R"
  )
  loader_path <- Filter(file.exists, candidates)[1]
  if (is.na(loader_path)) stop("[lm_coda_glmnet] Cannot find load_feature_table.R")
}
source(loader_path)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
tax_tbl <- if (nzchar(opt$taxonomy_table) && file.exists(opt$taxonomy_table)) opt$taxonomy_table else NULL
ft <- load_feature_table(opt$feature_table, opt$input_format, tax_tbl)
abund_table_full <- ft$abund_table   # samples x features
OTU_taxonomy     <- ft$OTU_taxonomy

message(sprintf("[lm_coda_glmnet] Feature table: %d samples x %d features",
                nrow(abund_table_full), ncol(abund_table_full)))

meta_table <- read.csv(opt$meta_table, header = TRUE, row.names = 1, check.names = FALSE)

# ---------------------------------------------------------------------------
# Aggregate to requested taxonomy level
# ---------------------------------------------------------------------------
which_level <- opt$which_level  # e.g. "Otus", "Genus", "Family", ...

if (which_level != "Otus" && which_level %in% colnames(OTU_taxonomy)) {
  message("[lm_coda_glmnet] Aggregating to level: ", which_level)
  level_labels <- OTU_taxonomy[[which_level]]
  level_labels[level_labels == "" | is.na(level_labels)] <- "Unknown"
  # Sum counts per level label
  unique_labels <- unique(level_labels)
  agg_mat <- vapply(unique_labels, function(lbl) {
    cols <- which(level_labels == lbl)
    if (length(cols) == 1) abund_table_full[, cols]
    else rowSums(abund_table_full[, cols, drop = FALSE])
  }, numeric(nrow(abund_table_full)))
  rownames(agg_mat) <- rownames(abund_table_full)
  abund_table_full <- agg_mat
} else if (which_level != "Otus") {
  message("[lm_coda_glmnet] Warning: which_level '", which_level,
          "' not found in taxonomy; using Otus (no aggregation)")
}

# ---------------------------------------------------------------------------
# Align samples between abundance table and metadata
# ---------------------------------------------------------------------------
common_samples <- intersect(rownames(abund_table_full), rownames(meta_table))
if (length(common_samples) == 0) {
  stop("[lm_coda_glmnet] No overlapping sample IDs between feature table and metadata")
}
abund_table_full <- abund_table_full[common_samples, , drop = FALSE]
meta_table       <- meta_table[common_samples, , drop = FALSE]
message(sprintf("[lm_coda_glmnet] %d samples after aligning with metadata", length(common_samples)))

# ---------------------------------------------------------------------------
# Library size filter
# ---------------------------------------------------------------------------
lib_sizes  <- rowSums(abund_table_full)
keep_samps <- lib_sizes >= opt$min_library_size
n_removed  <- sum(!keep_samps)
if (n_removed > 0) {
  message(sprintf("[lm_coda_glmnet] Removing %d samples below min_library_size (%d)",
                  n_removed, opt$min_library_size))
}
abund_table_full <- abund_table_full[keep_samps, , drop = FALSE]
meta_table       <- meta_table[keep_samps, , drop = FALSE]

# ---------------------------------------------------------------------------
# Exclude samples by column/value
# ---------------------------------------------------------------------------
if (nzchar(opt$exclude_column) && opt$exclude_column %in% colnames(meta_table)) {
  excl_vals <- trimws(strsplit(opt$exclude_values, ",")[[1]])
  keep_excl <- !as.character(meta_table[[opt$exclude_column]]) %in% excl_vals
  n_excl    <- sum(!keep_excl)
  if (n_excl > 0) {
    message(sprintf("[lm_coda_glmnet] Excluding %d samples by %s ∈ {%s}",
                    n_excl, opt$exclude_column, paste(excl_vals, collapse = ", ")))
  }
  abund_table_full <- abund_table_full[keep_excl, , drop = FALSE]
  meta_table       <- meta_table[keep_excl, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Build Groups column
# ---------------------------------------------------------------------------
if (nzchar(opt$groups_paste_columns)) {
  paste_cols <- trimws(strsplit(opt$groups_paste_columns, ",")[[1]])
  missing_pc <- setdiff(paste_cols, colnames(meta_table))
  if (length(missing_pc) > 0) {
    warning("[lm_coda_glmnet] groups_paste_columns not found in metadata: ",
            paste(missing_pc, collapse = ", "))
    paste_cols <- intersect(paste_cols, colnames(meta_table))
  }
  if (length(paste_cols) > 0) {
    meta_table$Groups <- apply(meta_table[, paste_cols, drop = FALSE], 1,
                               function(r) paste(r, collapse = "_"))
  }
} else if (nzchar(opt$groups_column) && opt$groups_column %in% colnames(meta_table)) {
  meta_table$Groups <- meta_table[[opt$groups_column]]
} else {
  # No grouping: treat all samples as one group
  meta_table$Groups <- "All"
}
meta_table$Groups <- as.factor(as.character(meta_table$Groups))

# ---------------------------------------------------------------------------
# Validate environmental covariates
# ---------------------------------------------------------------------------
if (!nzchar(opt$environmental_covariates)) {
  stop("[lm_coda_glmnet] --environmental_covariates is required (comma-sep metadata columns)")
}
environmental_covariates <- trimws(strsplit(opt$environmental_covariates, ",")[[1]])
environmental_covariates <- environmental_covariates[nzchar(environmental_covariates)]

missing_covs <- setdiff(environmental_covariates, colnames(meta_table))
if (length(missing_covs) > 0) {
  stop("[lm_coda_glmnet] The following environmental_covariates are not in metadata: ",
       paste(missing_covs, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Occupancy filter + Top-N selection
# ---------------------------------------------------------------------------
occ_threshold <- function(m, threshold, max_absent = 0) {
  occs <- colSums(m > max_absent)
  m[, occs >= threshold, drop = FALSE]
}

abund_table_full <- occ_threshold(abund_table_full, opt$occupancy_threshold)
if (ncol(abund_table_full) == 0) {
  stop("[lm_coda_glmnet] No features remain after occupancy filtering")
}
# Top-N by total abundance
top_n_keep <- order(colSums(abund_table_full), decreasing = TRUE)[
  seq_len(min(opt$top_n, ncol(abund_table_full)))]
abund_table_full <- abund_table_full[, top_n_keep, drop = FALSE]
message(sprintf("[lm_coda_glmnet] After occupancy filter + top-N: %d features",
                ncol(abund_table_full)))

# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------
normalise_table <- function(m, method) {
  TSS.divide <- function(x) x / sum(x)
  if (method == "logrelative") {
    log((m + 1) / (rowSums(m) + ncol(m)))
  } else if (method == "TSS+ILR") {
    as(logratio.transfo(t(apply(m + 1, 1, TSS.divide)), logratio = "ILR"), "matrix")
  } else if (method == "TSS+CLR") {
    as(logratio.transfo(t(apply(m + 1, 1, TSS.divide)), logratio = "CLR"), "matrix")
  } else {
    warning("[lm_coda_glmnet] Unknown normalisation_method '", method, "'; using raw counts")
    m
  }
}

normalised_table_full <- normalise_table(abund_table_full, opt$normalisation_method)

# ---------------------------------------------------------------------------
# Determine group levels for modelling
# ---------------------------------------------------------------------------
if (opt$covariates_global) {
  group_levels <- "All_Combined"
} else {
  group_levels <- levels(meta_table$Groups)
}

# ---------------------------------------------------------------------------
# Summary accumulator
# ---------------------------------------------------------------------------
model_summary_rows <- list()

# ---------------------------------------------------------------------------
# Per-group (or global) CODA-GLMNET
# ---------------------------------------------------------------------------
for (grp in group_levels) {

  if (opt$covariates_global) {
    mt <- meta_table
    at <- abund_table_full[rownames(mt), , drop = FALSE]
  } else {
    mt <- meta_table[meta_table$Groups == grp, , drop = FALSE]
    at <- abund_table_full[rownames(mt), , drop = FALSE]
  }

  # Drop all-zero features within this group
  at <- at[, colSums(at) > 0, drop = FALSE]

  for (cov in environmental_covariates) {

    # Skip if too few non-NA observations
    mt2 <- mt[!is.na(mt[[cov]]), , drop = FALSE]
    at2 <- at[rownames(mt2), , drop = FALSE]

    if (nrow(mt2) < 10) {
      message(sprintf(
        "[lm_coda_glmnet] SKIP group='%s' covariate='%s': only %d samples with non-NA values (min 10 required)",
        grp, cov, nrow(mt2)
      ))
      model_summary_rows[[length(model_summary_rows) + 1]] <- data.frame(
        group              = grp,
        covariate          = cov,
        n_samples          = nrow(mt2),
        n_taxa_in_signature = NA_integer_,
        lambda_used        = NA_character_,
        r_squared          = NA_real_,
        converged          = FALSE,
        skip_reason        = sprintf("n_samples=%d < 10", nrow(mt2)),
        stringsAsFactors   = FALSE
      )
      next
    }

    if (ncol(at2) < 3) {
      message(sprintf(
        "[lm_coda_glmnet] SKIP group='%s' covariate='%s': only %d features (min 3 required for GLMNET)",
        grp, cov, ncol(at2)
      ))
      model_summary_rows[[length(model_summary_rows) + 1]] <- data.frame(
        group              = grp,
        covariate          = cov,
        n_samples          = nrow(mt2),
        n_taxa_in_signature = NA_integer_,
        lambda_used        = NA_character_,
        r_squared          = NA_real_,
        converged          = FALSE,
        skip_reason        = sprintf("n_features=%d < 3", ncol(at2)),
        stringsAsFactors   = FALSE
      )
      next
    }

    if (length(unique(mt2[[cov]])) < 2) {
      message(sprintf(
        "[lm_coda_glmnet] SKIP group='%s' covariate='%s': zero variance (only one unique value)",
        grp, cov
      ))
      model_summary_rows[[length(model_summary_rows) + 1]] <- data.frame(
        group              = grp,
        covariate          = cov,
        n_samples          = nrow(mt2),
        n_taxa_in_signature = NA_integer_,
        lambda_used        = NA_character_,
        r_squared          = NA_real_,
        converged          = FALSE,
        skip_reason        = "zero variance in covariate",
        stringsAsFactors   = FALSE
      )
      next
    }

    tryCatch({
      y_val <- mt2[[cov]]

      # Determine if binary (and coda_binary_outcome is set)
      is_binary <- opt$coda_binary_outcome &&
                   (length(unique(y_val)) == 2 ||
                    (is.factor(y_val) && nlevels(y_val) == 2))

      res <- NULL
      if (is_binary) {
        if (!is.factor(y_val)) y_val <- as.factor(as.character(y_val))
        message(sprintf("[lm_coda_glmnet] Fitting binary CODA-GLMNET: group='%s' covariate='%s'",
                        grp, cov))
        res <- coda_glmnet(
          x         = at2,
          y         = y_val,
          lambda    = opt$coda_lambda,
          alpha     = opt$coda_alpha,
          showPlots = FALSE
        )
      } else {
        if (!is.numeric(y_val)) y_val <- as.numeric(as.character(y_val))
        message(sprintf("[lm_coda_glmnet] Fitting continuous CODA-GLMNET: group='%s' covariate='%s'",
                        grp, cov))
        res <- coda_glmnet(
          x         = at2,
          y         = y_val,
          lambda    = opt$coda_lambda,
          alpha     = opt$coda_alpha,
          showPlots = FALSE
        )
      }

      # ----- Extract signature -----
      sig_data <- res$`signature plot`$data
      # sig_data columns: first = taxa name, second = coefficient
      coef_df <- data.frame(
        taxa       = as.character(sig_data[[1]]),
        coefficient = as.numeric(sig_data[[2]]),
        sign       = ifelse(as.numeric(sig_data[[2]]) > 0, "positive", "negative"),
        group      = grp,
        covariate  = cov,
        stringsAsFactors = FALSE
      )

      coef_file <- file.path(opt$output_dir,
        paste0("CODA_coefficients_", opt$label, "_", grp, "_", cov, ".csv"))
      write.csv(coef_df, coef_file, row.names = FALSE)
      message("[lm_coda_glmnet] Wrote: ", coef_file)

      # ----- Extract predictions -----
      pred_data <- res$`predictions plot`$data
      pred_df <- data.frame(
        sample    = as.character(rownames(pred_data)),
        predicted = as.numeric(pred_data[[1]]),
        actual    = as.numeric(pred_data[[2]]),
        Groups    = grp,
        covariate = cov,
        stringsAsFactors = FALSE
      )

      pred_file <- file.path(opt$output_dir,
        paste0("CODA_predictions_", opt$label, "_", grp, "_", cov, ".csv"))
      write.csv(pred_df, pred_file, row.names = FALSE)
      message("[lm_coda_glmnet] Wrote: ", pred_file)

      # ----- R-squared (correlation^2 between predicted and actual) -----
      r_sq <- tryCatch(
        cor(pred_df$predicted, pred_df$actual, use = "complete.obs")^2,
        error = function(e) NA_real_
      )

      # ----- Lambda used -----
      lambda_val <- tryCatch(as.character(res$lambda), error = function(e) opt$coda_lambda)

      model_summary_rows[[length(model_summary_rows) + 1]] <- data.frame(
        group               = grp,
        covariate           = cov,
        n_samples           = nrow(mt2),
        n_taxa_in_signature = nrow(coef_df),
        lambda_used         = lambda_val,
        r_squared           = r_sq,
        converged           = TRUE,
        skip_reason         = "",
        stringsAsFactors    = FALSE
      )

    }, error = function(e) {
      message(sprintf(
        "[lm_coda_glmnet] ERROR group='%s' covariate='%s': %s",
        grp, cov, conditionMessage(e)
      ))
      model_summary_rows[[length(model_summary_rows) + 1]] <<- data.frame(
        group               = grp,
        covariate           = cov,
        n_samples           = nrow(mt2),
        n_taxa_in_signature = NA_integer_,
        lambda_used         = NA_character_,
        r_squared           = NA_real_,
        converged           = FALSE,
        skip_reason         = conditionMessage(e),
        stringsAsFactors    = FALSE
      )
    })
  }
}

# ---------------------------------------------------------------------------
# Write model summary CSV
# ---------------------------------------------------------------------------
summary_file <- file.path(opt$output_dir, paste0("CODA_model_summary_", opt$label, ".csv"))
if (length(model_summary_rows) > 0) {
  summary_df <- do.call(rbind, model_summary_rows)
} else {
  summary_df <- data.frame(
    group               = character(0),
    covariate           = character(0),
    n_samples           = integer(0),
    n_taxa_in_signature = integer(0),
    lambda_used         = character(0),
    r_squared           = numeric(0),
    converged           = logical(0),
    skip_reason         = character(0)
  )
}
write.csv(summary_df, summary_file, row.names = FALSE)
message("[lm_coda_glmnet] Wrote model summary: ", summary_file)

# ---------------------------------------------------------------------------
# Placeholder if no outputs were written
# ---------------------------------------------------------------------------
coef_files <- list.files(opt$output_dir, pattern = paste0("CODA_coefficients_", opt$label),
                         full.names = FALSE)
if (length(coef_files) == 0) {
  placeholder_file <- file.path(opt$output_dir,
    paste0("CODA_no_output_", opt$label, ".csv"))
  placeholder_df <- data.frame(
    status = "No CODA-GLMNET models were successfully fitted. See CODA_model_summary for details.",
    label  = opt$label,
    stringsAsFactors = FALSE
  )
  write.csv(placeholder_df, placeholder_file, row.names = FALSE)
  message("[lm_coda_glmnet] No successful models; wrote placeholder: ", placeholder_file)
}

message("[lm_coda_glmnet] Done.")
