#!/usr/bin/env Rscript
# subset_regression.R
# Best-subset regression (via leaps) with cross-validation or information criterion
# model selection. Outputs CSV files only; no plotting.

suppressPackageStartupMessages({
  library(optparse)
  library(leaps)
  library(dplyr)
  library(caret)
  library(purrr)
})

# ---------------------------------------------------------------------------
# CLI options
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--meta_table"),            type = "character", default = NULL,
              help = "Path to metadata CSV (samples x variables)"),
  make_option(c("--dependent_csv"),         type = "character", default = "",
              help = "Path to CSV with dependent variable (e.g. alpha diversity output)"),
  make_option(c("--output_dir"),            type = "character", default = ".",
              help = "Directory for output CSVs [default: .]"),
  make_option(c("--label"),                 type = "character", default = "Hypothesis1_Shannon",
              help = "Label prefix for output files [default: Hypothesis1_Shannon]"),
  make_option(c("--dependent_variable"),    type = "character", default = "Shannon",
              help = "Column name of the dependent variable [default: Shannon]"),
  make_option(c("--dependent_source"),      type = "character", default = "csv",
              help = "Where to read dependent variable: csv | metadata [default: csv]"),
  make_option(c("--explanatory_variables"), type = "character", default = "",
              help = "Comma-sep metadata columns to use as predictors (required)"),
  make_option(c("--regression_method"),     type = "character", default = "forward",
              help = "Subset selection method: exhaustive | backward | forward | seqrep [default: forward]"),
  make_option(c("--test_method"),           type = "character", default = "cv",
              help = "Model selection criterion: cv | aic | bic [default: cv]"),
  make_option(c("--really_big"),            action = "store_true", default = FALSE,
              help = "Pass really_big=TRUE to regsubsets (for large predictor sets)"),
  make_option(c("--cv_folds"),             type = "integer",   default = 5L,
              help = "Number of cross-validation folds (used when test_method=cv) [default: 5]"),
  make_option(c("--scale_predictors"),      type = "logical",   default = TRUE,
              help = "Standardize predictors before regression [default: TRUE]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if (is.null(opt$meta_table) || !nzchar(opt$meta_table)) {
  stop("[subset_regression] --meta_table is required")
}
if (!file.exists(opt$meta_table)) {
  stop("[subset_regression] meta_table not found: ", opt$meta_table)
}
if (!opt$dependent_source %in% c("csv", "metadata")) {
  stop("[subset_regression] --dependent_source must be: csv or metadata")
}
if (!opt$test_method %in% c("cv", "aic", "bic")) {
  stop("[subset_regression] --test_method must be: cv, aic, or bic")
}
if (!opt$regression_method %in% c("exhaustive", "backward", "forward", "seqrep")) {
  stop("[subset_regression] --regression_method must be: exhaustive, backward, forward, or seqrep")
}
if (!nzchar(opt$explanatory_variables)) {
  stop("[subset_regression] --explanatory_variables is required. Provide a comma-separated list of metadata column names.")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load metadata
# ---------------------------------------------------------------------------
meta_table <- read.csv(opt$meta_table, header = TRUE, row.names = 1, check.names = FALSE)

# ---------------------------------------------------------------------------
# Load dependent variable
# ---------------------------------------------------------------------------
dep_var <- opt$dependent_variable

if (opt$dependent_source == "csv") {
  if (!nzchar(opt$dependent_csv)) {
    stop("[subset_regression] --dependent_csv is required when --dependent_source=csv")
  }
  if (!file.exists(opt$dependent_csv)) {
    stop("[subset_regression] dependent_csv not found: ", opt$dependent_csv)
  }
  message("[subset_regression] Reading dependent variable from CSV: ", opt$dependent_csv)
  dep_table <- read.csv(opt$dependent_csv, header = TRUE, row.names = 1, check.names = FALSE)

  if (!dep_var %in% colnames(dep_table)) {
    stop(sprintf(
      "[subset_regression] Dependent variable '%s' not found in dependent_csv. Available columns: %s",
      dep_var, paste(colnames(dep_table), collapse = ", ")
    ))
  }

  # Align samples
  common_samps <- intersect(rownames(meta_table), rownames(dep_table))
  if (length(common_samps) == 0) {
    stop("[subset_regression] No overlapping sample IDs between meta_table and dependent_csv")
  }
  meta_table   <- meta_table[common_samps, , drop = FALSE]
  dep_values   <- dep_table[common_samps, dep_var, drop = TRUE]
  names(dep_values) <- common_samps

} else {
  # dependent_source == "metadata"
  if (!dep_var %in% colnames(meta_table)) {
    stop(sprintf(
      "[subset_regression] Dependent variable '%s' not found in meta_table. Available columns: %s",
      dep_var, paste(colnames(meta_table), collapse = ", ")
    ))
  }
  message("[subset_regression] Reading dependent variable from metadata column: ", dep_var)
  dep_values <- meta_table[[dep_var]]
  names(dep_values) <- rownames(meta_table)
}

dep_values <- as.numeric(dep_values)

# ---------------------------------------------------------------------------
# Parse and validate explanatory variables
# ---------------------------------------------------------------------------
expl_vars <- trimws(strsplit(opt$explanatory_variables, ",")[[1]])
expl_vars <- expl_vars[nzchar(expl_vars)]

missing_expl <- setdiff(expl_vars, colnames(meta_table))
if (length(missing_expl) > 0) {
  stop(sprintf(
    "[subset_regression] The following explanatory_variables are not in meta_table: %s",
    paste(missing_expl, collapse = ", ")
  ))
}

# ---------------------------------------------------------------------------
# Build modelling dataset
# ---------------------------------------------------------------------------
lm_dat <- meta_table[, expl_vars, drop = FALSE]
lm_dat[[dep_var]] <- dep_values

# Remove rows with NA in dependent variable
na_dep <- is.na(lm_dat[[dep_var]])
if (any(na_dep)) {
  message(sprintf("[subset_regression] Removing %d rows with NA in dependent variable '%s'",
                  sum(na_dep), dep_var))
  lm_dat <- lm_dat[!na_dep, , drop = FALSE]
}

if (nrow(lm_dat) < 5) {
  stop("[subset_regression] Fewer than 5 complete observations — cannot fit models")
}

# ---------------------------------------------------------------------------
# Drop zero-variance or all-NA predictors
# ---------------------------------------------------------------------------
dropped_vars <- character(0)
for (v in expl_vars) {
  col <- lm_dat[[v]]
  # Convert to numeric if possible
  if (!is.numeric(col)) {
    col_num <- suppressWarnings(as.numeric(as.character(col)))
    if (sum(!is.na(col_num)) > 0) {
      lm_dat[[v]] <- col_num
      col <- col_num
    }
  }
  if (all(is.na(col))) {
    message(sprintf("[subset_regression] Dropping predictor '%s': all values are NA", v))
    dropped_vars <- c(dropped_vars, v)
  } else if (is.numeric(col) && var(col, na.rm = TRUE) == 0) {
    message(sprintf("[subset_regression] Dropping predictor '%s': zero variance", v))
    dropped_vars <- c(dropped_vars, v)
  }
}

if (length(dropped_vars) > 0) {
  warning(sprintf(
    "[subset_regression] Dropped %d predictor(s) due to all-NA or zero variance: %s",
    length(dropped_vars), paste(dropped_vars, collapse = ", ")
  ))
  expl_vars <- setdiff(expl_vars, dropped_vars)
  lm_dat    <- lm_dat[, c(expl_vars, dep_var), drop = FALSE]
}

if (length(expl_vars) == 0) {
  stop("[subset_regression] No explanatory variables remain after dropping zero-variance / all-NA columns")
}

# ---------------------------------------------------------------------------
# Remove rows with NA in any predictor
# ---------------------------------------------------------------------------
na_rows <- apply(lm_dat[, expl_vars, drop = FALSE], 1, function(r) any(is.na(r)))
if (any(na_rows)) {
  message(sprintf("[subset_regression] Removing %d rows with NA in predictors", sum(na_rows)))
  lm_dat <- lm_dat[!na_rows, , drop = FALSE]
}

if (nrow(lm_dat) < 5) {
  stop("[subset_regression] Fewer than 5 complete observations after removing NA rows")
}

message(sprintf("[subset_regression] Modelling dataset: %d observations x %d predictors",
                nrow(lm_dat), length(expl_vars)))

# ---------------------------------------------------------------------------
# Scale predictors
# ---------------------------------------------------------------------------
if (isTRUE(opt$scale_predictors)) {
  message("[subset_regression] Standardizing predictors (scale_predictors=TRUE)")
  for (v in expl_vars) {
    col <- lm_dat[[v]]
    if (is.numeric(col)) {
      sd_col <- sd(col, na.rm = TRUE)
      if (!is.na(sd_col) && sd_col > 0) {
        lm_dat[[v]] <- as.numeric(scale(col))
      }
    }
  }
}

# ---------------------------------------------------------------------------
# Best-subset regression via leaps::regsubsets
# ---------------------------------------------------------------------------
n_vars  <- length(expl_vars)
formula_str <- paste(dep_var, "~", paste(expl_vars, collapse = " + "))

message(sprintf("[subset_regression] Running regsubsets: method=%s, nvmax=%d",
                opt$regression_method, n_vars))

models <- regsubsets(
  as.formula(formula_str),
  data      = lm_dat,
  nvmax     = n_vars,
  method    = opt$regression_method,
  really.big = opt$really_big
)

model_summary_obj <- summary(models)
n_models          <- nrow(model_summary_obj$which)

# ---------------------------------------------------------------------------
# Helper: build lm formula for a given model size from regsubsets
# ---------------------------------------------------------------------------
get_model_formula <- function(n_vars_in_model, models_obj, dep_var) {
  which_mat  <- summary(models_obj)$which
  which_row  <- which_mat[n_vars_in_model, ]
  sel_vars   <- names(which_row)[which_row]
  sel_vars   <- sel_vars[sel_vars != "(Intercept)"]
  as.formula(paste(dep_var, "~", paste(sel_vars, collapse = " + ")))
}

# ---------------------------------------------------------------------------
# Model selection
# ---------------------------------------------------------------------------
best_n <- NULL

if (opt$test_method == "cv") {
  message(sprintf("[subset_regression] Running %d-fold CV for model selection", opt$cv_folds))
  set.seed(42)
  folds     <- createFolds(lm_dat[[dep_var]], k = opt$cv_folds, list = TRUE)
  cv_errors <- numeric(n_models)

  for (m in seq_len(n_models)) {
    fold_errors <- numeric(opt$cv_folds)
    fml <- get_model_formula(m, models, dep_var)

    for (f in seq_along(folds)) {
      test_idx   <- folds[[f]]
      train_dat  <- lm_dat[-test_idx, , drop = FALSE]
      test_dat   <- lm_dat[test_idx, , drop = FALSE]

      fit  <- tryCatch(lm(fml, data = train_dat), error = function(e) NULL)
      if (is.null(fit)) { fold_errors[f] <- NA_real_; next }
      pred <- tryCatch(predict(fit, newdata = test_dat), error = function(e) NULL)
      if (is.null(pred)) { fold_errors[f] <- NA_real_; next }

      actual         <- test_dat[[dep_var]]
      fold_errors[f] <- mean((actual - pred)^2, na.rm = TRUE)
    }
    cv_errors[m] <- mean(fold_errors, na.rm = TRUE)
  }

  cv_se  <- vapply(seq_len(n_models), function(m) {
    fml <- get_model_formula(m, models, dep_var)
    fold_errors <- numeric(opt$cv_folds)
    for (f in seq_along(folds)) {
      test_idx  <- folds[[f]]
      train_dat <- lm_dat[-test_idx, , drop = FALSE]
      test_dat  <- lm_dat[test_idx, , drop = FALSE]
      fit  <- tryCatch(lm(fml, data = train_dat), error = function(e) NULL)
      if (is.null(fit)) { fold_errors[f] <- NA_real_; next }
      pred <- tryCatch(predict(fit, newdata = test_dat), error = function(e) NULL)
      if (is.null(pred)) { fold_errors[f] <- NA_real_; next }
      actual <- test_dat[[dep_var]]
      fold_errors[f] <- mean((actual - pred)^2, na.rm = TRUE)
    }
    sd(fold_errors, na.rm = TRUE) / sqrt(opt$cv_folds)
  }, numeric(1))

  cv_df <- data.frame(
    nvars            = seq_len(n_models),
    cv_error         = cv_errors,
    se               = cv_se,
    selection_method = "cv"
  )

  best_n <- which.min(cv_errors)

} else if (opt$test_method == "aic") {
  message("[subset_regression] Using AIC for model selection")
  aic_vals <- vapply(seq_len(n_models), function(m) {
    fml <- get_model_formula(m, models, dep_var)
    fit <- tryCatch(lm(fml, data = lm_dat), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    AIC(fit)
  }, numeric(1))

  cv_df <- data.frame(
    nvars            = seq_len(n_models),
    cv_error         = aic_vals,
    se               = NA_real_,
    selection_method = "aic"
  )
  best_n <- which.min(aic_vals)

} else {
  # bic
  message("[subset_regression] Using BIC for model selection")
  bic_vals <- -model_summary_obj$bic   # regsubsets reports negative BIC; lower is better in std BIC
  # Use actual BIC from lm() for consistency
  bic_lm <- vapply(seq_len(n_models), function(m) {
    fml <- get_model_formula(m, models, dep_var)
    fit <- tryCatch(lm(fml, data = lm_dat), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    BIC(fit)
  }, numeric(1))

  cv_df <- data.frame(
    nvars            = seq_len(n_models),
    cv_error         = bic_lm,
    se               = NA_real_,
    selection_method = "bic"
  )
  best_n <- which.min(bic_lm)
}

# Fallback
if (is.null(best_n) || is.na(best_n)) best_n <- 1L

message(sprintf("[subset_regression] Best model: %d predictor(s)", best_n))

# ---------------------------------------------------------------------------
# OUTPUT 1: CV/AIC/BIC errors
# ---------------------------------------------------------------------------
cv_file <- file.path(opt$output_dir, paste0("SubsetReg_cv_errors_", opt$label, ".csv"))
write.csv(cv_df, cv_file, row.names = FALSE)
message("[subset_regression] Wrote: ", cv_file)

# ---------------------------------------------------------------------------
# OUTPUT 2: Best model coefficients
# ---------------------------------------------------------------------------
best_fml <- get_model_formula(best_n, models, dep_var)
best_fit  <- lm(best_fml, data = lm_dat)
best_coef <- summary(best_fit)$coefficients

best_coef_df <- data.frame(
  variable    = rownames(best_coef),
  coefficient = as.numeric(best_coef[, "Estimate"]),
  std_error   = as.numeric(best_coef[, "Std. Error"]),
  t_statistic = as.numeric(best_coef[, "t value"]),
  pvalue      = as.numeric(best_coef[, "Pr(>|t|)"]),
  stringsAsFactors = FALSE
)

best_file <- file.path(opt$output_dir, paste0("SubsetReg_best_model_", opt$label, ".csv"))
write.csv(best_coef_df, best_file, row.names = FALSE)
message("[subset_regression] Wrote: ", best_file)

# ---------------------------------------------------------------------------
# OUTPUT 3: All candidate models summary
# ---------------------------------------------------------------------------
which_mat <- model_summary_obj$which

# Variables included (excluding intercept) as a semicolon-joined string
vars_included <- vapply(seq_len(n_models), function(m) {
  row    <- which_mat[m, ]
  sel    <- names(row)[row]
  sel    <- sel[sel != "(Intercept)"]
  paste(sel, collapse = ";")
}, character(1))

# AIC and BIC for each model size
aic_all <- vapply(seq_len(n_models), function(m) {
  fml <- get_model_formula(m, models, dep_var)
  fit <- tryCatch(lm(fml, data = lm_dat), error = function(e) NULL)
  if (is.null(fit)) NA_real_ else AIC(fit)
}, numeric(1))

bic_all <- vapply(seq_len(n_models), function(m) {
  fml <- get_model_formula(m, models, dep_var)
  fit <- tryCatch(lm(fml, data = lm_dat), error = function(e) NULL)
  if (is.null(fit)) NA_real_ else BIC(fit)
}, numeric(1))

all_models_df <- data.frame(
  nvars             = seq_len(n_models),
  variables_included = vars_included,
  adjr2             = model_summary_obj$adjr2,
  bic               = bic_all,
  cp                = model_summary_obj$cp,
  aic               = aic_all,
  stringsAsFactors  = FALSE
)

all_file <- file.path(opt$output_dir, paste0("SubsetReg_all_models_", opt$label, ".csv"))
write.csv(all_models_df, all_file, row.names = FALSE)
message("[subset_regression] Wrote: ", all_file)

message("[subset_regression] Done.")
