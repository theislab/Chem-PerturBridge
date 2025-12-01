#!/usr/bin/env Rscript

# Step 1: Filter - Filter genes and normalize counts

# Configure reticulate to use conda Python
suppressPackageStartupMessages(library(reticulate))
conda_env <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_env)) use_python(file.path(conda_env, "bin/python"), required = TRUE)

requireNamespace("anndata", quietly=TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(argparse)
  library(jsonlite)
  library(rhdf5)
})

source("./src/deg/subsampling.R")
# helper to sanitize names for model.matrix/contrasts
clean <- function(x) gsub("[^[:alnum:]_]", "_", x)

#' Filter samples with low cell counts
filter_cells <- function(adata, min_cells) {
  stopifnot(!missing(adata), !missing(min_cells), is.numeric(min_cells), min_cells >= 0,
            "psbulk_cells" %in% names(adata$obs))
  
  n_obs_prev <- adata$n_obs
  cat("Filter samples with the low number of cells\n")
  cat("    Filter threshold (min_cells): ", min_cells, "\n")
  
  keep_mask <- (adata$obs$psbulk_cells == -666) | (adata$obs$psbulk_cells >= min_cells)
  n_sentinel <- sum(adata$obs$psbulk_cells == -666, na.rm = TRUE)
  
  if (n_sentinel > 0) {
    cat("    Found ", n_sentinel, " samples with sentinel value (-666), keeping all\n")
  }
  
  adata <- adata[keep_mask]
  cat("    n_obs: ", n_obs_prev, "--> ", adata$n_obs, "\n")
  return(adata)
}

#' Filter samples based on quality control
filter_qc <- function(adata, qc_enabled) {
  stopifnot(!missing(adata), !missing(qc_enabled), is.logical(qc_enabled))
  
  if (!qc_enabled) {
    return(adata)
  }
  
  qc_column <- "qc_pass"
  
  if (qc_column %in% names(adata$obs)) {
    n_obs_prev <- adata$n_obs
    cat("Filter samples which did not pass qc\n")
    cat("    QC column: ", qc_column, "\n")
    adata <- adata[adata$obs[[qc_column]] == TRUE, ]
    cat("    n_obs: ", n_obs_prev, "--> ", adata$n_obs, "\n")
  } else {
    warning("QC filtering requested but no QC column found in data. ",
            "Checked for: ", qc_column, ". ",
            "Proceeding without QC filtering.")
  }
  return(adata)
}

#' Filter redundant controls from dataset
filter_controls <- function(ad, col) {
  stopifnot(!missing(ad), !missing(col), is.character(col), col %in% names(ad$obs))
  
  cat("    Filter controls by ", col, " column\n")
  n_obs_prev <- ad$n_obs
    
  control <- ad$obs %>%
    filter(is_control==TRUE) %>%
    pull(col) %>%
    unique()

  treated <- ad$obs %>%
    filter(is_control!=TRUE) %>%
    pull(col) %>%
    unique()
  
  redundant_controls <- setdiff(control, treated)
  if (length(redundant_controls) > 0) {
    ad <- ad[!(ad$obs[, col] %in% redundant_controls)]
  }
  
  cat("        n_obs: ", n_obs_prev, "--> ",ad$n_obs, "\n\n")
  return(ad)
}

#' Check for missing controls in dataset
check_for_controls <- function(obs, col) {
  stopifnot(!missing(obs), !missing(col), is.character(col),
            col %in% names(obs), "is_control" %in% names(obs))
  
  control <- obs %>%
    filter(is_control==TRUE) %>%
    pull(col) %>%
    unique()

  treated <- obs %>%
    filter(is_control!=TRUE) %>%
    pull(col) %>%
    unique()
    
  missing_controls <- setdiff(treated, control)
  if (length(missing_controls) > 0) {
    warning("Missing controls for '", col, "' columns: ", paste(missing_controls, collapse=", "))
  }
}

#' Filter genes and normalize counts
#' 
#' @param adata AnnData object
#' @param args Arguments list containing parameters
#' @param input_basename Basename for output file naming
#' @return List containing filtered DGEList, design matrix, and metadata
filter_and_prepare_dge <- function(adata,
                                   args,
                                   input_basename) {
  # Parameter validation
  stopifnot(!missing(adata),
            !missing(args),
            !missing(input_basename),
            is.character(input_basename))
  
  # Get cell type (should be only one since input files are split by cell type)
  cell_types <- unique(as.character(adata$obs$cell_type))
  if (length(cell_types) != 1) {
    stop("Expected exactly one cell type in input file, found: ", paste(cell_types, collapse = ", "))
  }
  cl <- cell_types[1]
  
  cat("\n▶︎ Processing cell_type: ", cl, "\n")
  
  # Filter controls
  ad <- filter_controls(adata, "pert_time_h")
  
  obs <- ad$obs %>%
    mutate(
      raw_cond = perturbation_label,
      cond = factor(clean(perturbation_label)),
    )
  
  # Check for required controls
  check_for_controls(obs, "pert_time_h")

  # Build DGEList and design matrix
  cat("    Constructing design matrix...\n")
  start_time <- Sys.time()
  counts <- Matrix::t(ad$X)
  dge    <- DGEList(counts=counts)
  design <- model.matrix(
    ~ 0 + cond + plate,
    data = obs
  )
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("    Design matrix is constructed in %.1f seconds\n", diff_time))
  cat("    Design matrix dimensions:", nrow(design), "samples ×", ncol(design), "coefficients\n")
  
  # Filter genes and normalize
  cat("    Filtering genes and normalizing counts...\n")
  start_time <- Sys.time()
  keep <- filterByExpr(dge, design)
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  dge  <- dge[keep, , keep.lib.sizes=FALSE] %>% calcNormFactors()
  cat(sprintf("    Genes are filtered and normalized in %.1f seconds\n", diff_time))
  
  # Prepare intermediate data
  intermediate_data <- list(
    dge = dge,
    design = design,
    obs = obs,
    ad_var = ad$var,
    ad_uns = if (!is.null(adata$uns) && length(adata$uns) > 0) as.list(adata$uns) else list(),
    parameters = list(
      min_cells = args$min_cells,
      qc = args$qc,
      subsampling = args$subsampling
    )
  )
  
  return(intermediate_data)
}

#' Merge command line arguments with configuration
merge_config <- function(args, config) {
  stopifnot(!missing(args), !missing(config), is.list(args), is.list(config))
  
  config_merged <- args
  for (name in names(config)) {
    if (is.null(args[[name]])) {
      config_merged[[name]] <- config[[name]]
    }
  }
  return(config_merged)
}

#' Main function for Step 1: Filter genes and normalize
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--input", type="character", required=TRUE,
                     help="Path to input h5ad file")
  parser$add_argument("--output_dir", type="character", required=TRUE,
                     help="Directory to save intermediate filtered DGEList")
  parser$add_argument("--min_cells", type="integer", default=0)
  parser$add_argument("--qc", action="store_true")
  parser$add_argument("--subsampling", action="store_true")
  parser$add_argument("--max_cell_types", type="integer")
  parser$add_argument("--max_perturbations", type="integer")
  parser$add_argument("--max_genes", type="integer")
  parser$add_argument("--specific_times", nargs = "+", type="integer", default=c(24))
  parser$add_argument("--specific_perturbagens", type="character", default=NULL)
  parser$add_argument("--config", default="{}")

  args <- parser$parse_args()
  config <- jsonlite::fromJSON(args$config)
  args$config = NULL
  args <- merge_config(args, config)  
  
  if (is.null(args$min_cells)) {
    args$min_cells <- 0
  }
  stopifnot(is.numeric(args$min_cells), args$min_cells >= 0)
  
  if (is.null(args$qc)) {
    args$qc <- FALSE
  }
  stopifnot(is.logical(args$qc))
  
  deg_start <- Sys.time()
  cat("\nDEG Step 1: Filter - Started...\n")
  
  warning_messages <- list()
  
  input_basename <- tools::file_path_sans_ext(basename(args$input))
  input_basename <- gsub("_processed_", "_", input_basename)
  input_basename <- gsub("_processed$", "", input_basename)
  
  adata <- anndata::read_h5ad(args$input)
  cat("Loaded input file: ", args$input, "\n")
  cat("Initial dimensions: ", adata$n_obs, " samples × ", adata$n_vars, " genes\n")
  
  adata <- filter_cells(adata, args$min_cells)
  adata <- filter_qc(adata, args$qc)
  adata <- subsampling(adata, args)
  
  # Process with warning capture
  withCallingHandlers({
    intermediate_data <- filter_and_prepare_dge(adata, args, input_basename)
    
    if (!dir.exists(args$output_dir)) {
      dir.create(args$output_dir, recursive = TRUE)
    }
    
    intermediate_file <- file.path(args$output_dir, paste0(input_basename, "_filtered.rds"))
    cat("\n    Saving filtered DGEList: ", intermediate_file, "\n")
    
    saveRDS(intermediate_data, intermediate_file, compress = "xz")
    
  }, warning = function(w) {
    warning_messages[[length(warning_messages) + 1]] <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  
  deg_time <- difftime(Sys.time(), deg_start, units = "secs")
  cat(sprintf("\nStep 1 Total runtime: %.1f seconds\n", deg_time))
  
  if (length(warning_messages) > 0) {
    message(sprintf("\n=== %d Warnings encountered ===\n", length(warning_messages)))
    unique_warns <- unique(unlist(warning_messages))
    n_show <- min(20, length(unique_warns))
    for (i in 1:n_show) {
      message(sprintf("%d. %s", i, unique_warns[i]))
    }
    if (length(unique_warns) > n_show) {
      message(sprintf("... and %d more unique warning types (total %d warnings)", 
                  length(unique_warns) - n_show, length(warning_messages)))
    }
  }
  
  cat("\nDEG Step 1: Filter - Completed!\n")
}

main()

