# load_feature_table.R
# Loads a feature table from BIOM, TSV, or GTDB format.
# Returns a list with:
#   $abund_table   - samples x features matrix (numeric)
#   $feature_taxonomy  - data.frame with columns Kingdom,Phylum,Class,Order,Family,Genus,Feature

suppressPackageStartupMessages({
  library(phyloseq)
  library(stringr)
})

# ---------------------------------------------------------------------------
# Helper: strip rank prefixes from taxonomy strings
# ---------------------------------------------------------------------------
.strip_qiime2_prefixes <- function(x) {
  # QIIME2 style: D_0__, D_1__, etc.  or  d__, p__, c__, o__, f__, g__, s__
  x <- gsub("D_[0-9]+__", "", x)
  x <- gsub("^[dpcofgs]__", "", x)
  x <- trimws(x)
  x
}

.strip_gtdb_prefixes <- function(x) {
  # GTDB style: d__, p__, c__, o__, f__, g__, s__
  x <- gsub("^[dpcofgs]__", "", x)
  x <- trimws(x)
  x
}

# ---------------------------------------------------------------------------
# Helper: parse a semicolon-separated taxonomy string into 7 ranks
# ---------------------------------------------------------------------------
.parse_taxonomy_string <- function(tax_str, strip_fn = identity) {
  parts <- strsplit(as.character(tax_str), ";")[[1]]
  parts <- strip_fn(trimws(parts))
  # pad to 7 ranks
  length(parts) <- 7
  parts[is.na(parts)] <- ""
  names(parts) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Feature")
  parts
}

# ---------------------------------------------------------------------------
# Helper: apply common taxonomy filters
#   - Remove features with Unassigned / empty Kingdom
#   - Remove features with blank Phylum
#   - Remove Chloroplast (Order) and Mitochondria (Family)
# ---------------------------------------------------------------------------
.filter_taxonomy <- function(abund_table, tax_df) {
  keep <- rep(TRUE, nrow(tax_df))

  # Unassigned or missing Kingdom
  unassigned_k <- tax_df$Kingdom == "" |
                  is.na(tax_df$Kingdom) |
                  grepl("^[Uu]nassigned", tax_df$Kingdom)
  keep[unassigned_k] <- FALSE

  # Blank Phylum
  blank_p <- tax_df$Phylum == "" | is.na(tax_df$Phylum)
  keep[blank_p] <- FALSE

  # Chloroplast Order
  chloro <- grepl("[Cc]hloroplast", tax_df$Order)
  keep[chloro] <- FALSE

  # Mitochondria Family
  mito <- grepl("[Mm]itochondri", tax_df$Family)
  keep[mito] <- FALSE

  n_removed <- sum(!keep)
  if (n_removed > 0) {
    message(sprintf(
      "[load_feature_table] Pruned %d features (Unassigned/blank Kingdom, blank Phylum, Chloroplast, Mitochondria)",
      n_removed
    ))
  }

  list(
    abund_table  = abund_table[, keep, drop = FALSE],
    feature_taxonomy = tax_df[keep, , drop = FALSE]
  )
}

# ---------------------------------------------------------------------------
# BIOM loader
# ---------------------------------------------------------------------------
.load_biom <- function(feature_table) {
  message("[load_feature_table] Reading BIOM file: ", feature_table)
  ps <- import_biom(feature_table)

  # Abundance matrix: phyloseq stores features x samples → transpose
  abund_mat <- as.matrix(otu_table(ps))
  if (taxa_are_rows(ps)) {
    abund_mat <- t(abund_mat)   # now samples x features
  }
  storage.mode(abund_mat) <- "numeric"

  # Taxonomy
  tax_raw <- tax_table(ps)
  tax_df  <- as.data.frame(tax_raw, stringsAsFactors = FALSE)

  # Rename columns to standard names and strip QIIME2 prefixes
  col_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Feature")
  if (ncol(tax_df) >= 7) {
    colnames(tax_df)[1:7] <- col_names[1:7]
  } else {
    # pad
    for (cn in col_names[seq(ncol(tax_df) + 1, 7)]) {
      tax_df[[cn]] <- ""
    }
    colnames(tax_df)[seq_len(min(ncol(tax_df), 7))] <- col_names[seq_len(min(ncol(tax_df), 7))]
  }
  tax_df <- tax_df[, col_names, drop = FALSE]

  for (cn in col_names) {
    tax_df[[cn]] <- .strip_qiime2_prefixes(tax_df[[cn]])
  }

  list(abund_table = abund_mat, feature_taxonomy = tax_df)
}

# ---------------------------------------------------------------------------
# TSV / GTDB loader
# ---------------------------------------------------------------------------
.load_tsv <- function(feature_table, taxonomy_table = NULL, strip_fn = identity) {
  message("[load_feature_table] Reading TSV file: ", feature_table)
  raw <- read.table(
    feature_table,
    sep       = "\t",
    header    = TRUE,
    row.names = 1,
    check.names = FALSE,
    comment.char = "#"
  )

  # Auto-detect orientation: if nrow > ncol → likely features-first → transpose
  if (nrow(raw) > ncol(raw)) {
    message("[load_feature_table] Detected features-first orientation; transposing to samples x features")
    raw <- t(raw)
  }

  abund_mat <- as.matrix(raw)
  storage.mode(abund_mat) <- "numeric"

  col_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Feature")

  # Parse taxonomy table if provided
  if (!is.null(taxonomy_table) && nchar(taxonomy_table) > 0 && file.exists(taxonomy_table)) {
    message("[load_feature_table] Reading taxonomy table: ", taxonomy_table)
    tax_raw <- read.table(
      taxonomy_table,
      sep           = "\t",
      header        = TRUE,
      row.names     = 1,
      check.names   = FALSE,
      stringsAsFactors = FALSE
    )
    # Expect: col1 = featureID (rowname), col2 (or named "Taxon"/"taxonomy") = semicolon-sep string
    tax_col_idx <- grep("taxon|taxonomy|classification", colnames(tax_raw), ignore.case = TRUE)
    if (length(tax_col_idx) == 0) tax_col_idx <- 1L

    tax_strings <- tax_raw[, tax_col_idx[1]]
    names(tax_strings) <- rownames(tax_raw)

    parsed <- lapply(tax_strings, .parse_taxonomy_string, strip_fn = strip_fn)
    tax_df  <- as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)

    # Align to feature IDs in abund_mat
    feat_ids   <- colnames(abund_mat)
    common_ids <- intersect(feat_ids, rownames(tax_df))
    if (length(common_ids) == 0) {
      warning("[load_feature_table] No overlapping IDs between feature table and taxonomy table. Returning empty taxonomy.")
      tax_df <- data.frame(
        Kingdom = character(length(feat_ids)),
        Phylum  = character(length(feat_ids)),
        Class   = character(length(feat_ids)),
        Order   = character(length(feat_ids)),
        Family  = character(length(feat_ids)),
        Genus   = character(length(feat_ids)),
        Feature = feat_ids,
        stringsAsFactors = FALSE
      )
      rownames(tax_df) <- feat_ids
    } else {
      # Subset to common IDs
      abund_mat <- abund_mat[, common_ids, drop = FALSE]
      tax_df    <- tax_df[common_ids, col_names, drop = FALSE]
    }

  } else {
    # No taxonomy table: use feature IDs as Feature, leave other ranks blank
    message("[load_feature_table] No taxonomy table provided; using feature IDs as Feature column")
    feat_ids <- colnames(abund_mat)
    tax_df <- data.frame(
      Kingdom = rep("", length(feat_ids)),
      Phylum  = rep("", length(feat_ids)),
      Class   = rep("", length(feat_ids)),
      Order   = rep("", length(feat_ids)),
      Family  = rep("", length(feat_ids)),
      Genus   = rep("", length(feat_ids)),
      Feature = feat_ids,
      stringsAsFactors = FALSE
    )
    rownames(tax_df) <- feat_ids
  }

  list(abund_table = abund_mat, feature_taxonomy = tax_df)
}

# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------
load_feature_table <- function(feature_table,
                               input_format,
                               taxonomy_table = NULL) {

  input_format <- tolower(trimws(input_format))

  result <- switch(
    input_format,
    "biom" = .load_biom(feature_table),
    "tsv"  = .load_tsv(feature_table, taxonomy_table, strip_fn = .strip_qiime2_prefixes),
    "gtdb" = .load_tsv(feature_table, taxonomy_table, strip_fn = .strip_gtdb_prefixes),
    stop("[load_feature_table] Unknown input_format '", input_format,
         "'. Must be one of: biom, tsv, gtdb")
  )

  # Apply common taxonomy filters (only when we actually have taxonomy)
  has_taxonomy <- any(result$feature_taxonomy$Kingdom != "" & !is.na(result$feature_taxonomy$Kingdom))
  if (has_taxonomy) {
    result <- .filter_taxonomy(result$abund_table, result$feature_taxonomy)
  }

  message(sprintf(
    "[load_feature_table] Loaded %d samples x %d features",
    nrow(result$abund_table), ncol(result$abund_table)
  ))

  result
}
