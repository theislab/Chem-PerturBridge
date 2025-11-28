#!/usr/bin/env Rscript

# Step 1: Prepare - Run lmFit and save intermediate results for batch processing

# Configure reticulate to use conda Python
suppressPackageStartupMessages(library(reticulate))
conda_env <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_env)) use_python(file.path(conda_env, "bin/python"), required = TRUE)


requireNamespace("anndata", quietly=TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(argparse)
  library(jsonlite)
  library(rhdf5)
  library(parallel)
})

source("./src/deg/subsampling.R")
# helper to sanitize names for model.matrix/contrasts
clean <- function(x) gsub("[^[:alnum:]_]", "_", x)

#' Filter samples with low cell counts
#' 
#' @param adata AnnData object containing pseudobulk data
#' @param min_cells Minimum number of cells required per sample
#' @return Filtered AnnData object
filter_cells <- function(adata,
                         min_cells) {
  # Parameter validation
  stopifnot(!missing(adata),
            !missing(min_cells),
            is.numeric(min_cells), min_cells >= 0,
            "psbulk_cells" %in% names(adata$obs))
  
  n_obs_prev <- adata$n_obs
  cat("Filter samples with the low number of cells\n")
  cat("    Filter threshold (min_cells): ", min_cells, "\n")
  
  # Handle -666 as a sentinel value meaning "keep this sample regardless of threshold"
  # Samples with -666 are kept, others are filtered by min_cells threshold
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
#' 
#' @param adata AnnData object containing pseudobulk data
#' @param qc_enabled Whether to filter by quality control
#' @return Filtered AnnData object (or original if QC column not found)
filter_qc <- function(adata,
                      qc_enabled) {
  # Parameter validation
  stopifnot(!missing(adata),
            !missing(qc_enabled),
            is.logical(qc_enabled))
  
  # If QC filtering is not enabled, return original data
  if (!qc_enabled) {
    return(adata)
  }
  
  # Check for QC column
  qc_column <- "qc_pass"
  
  if (qc_column %in% names(adata$obs)) {
    n_obs_prev <- adata$n_obs
    cat("Filter samples which did not pass qc\n")
    cat("    QC column: ", qc_column, "\n")
    # Use [[ to access column by variable name
    adata <- adata[adata$obs[[qc_column]] == TRUE, ]
    cat("    n_obs: ", n_obs_prev, "--> ", adata$n_obs, "\n")
  } else {
    # If no QC column found, warn and return original data
    warning("QC filtering requested but no QC column found in data. ",
            "Checked for: ", qc_column, ". ",
            "Proceeding without QC filtering.")
  }
  return(adata)
}

#' Check for missing controls in dataset
#' 
#' @param obs Observation dataframe
#' @param col Column name to check for controls
#' @return NULL (prints warnings for missing controls)
check_for_controls <- function(obs, 
                               col) {
  # Parameter validation
  stopifnot(!missing(obs),
            !missing(col),
            is.character(col),
            col %in% names(obs),
            "is_control" %in% names(obs))
  
  # Check for required controls
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

#' Filter redundant controls from dataset
#' 
#' @param ad AnnData object
#' @param col Column name to filter controls by
#' @return Filtered AnnData object
filter_controls <- function(ad,
                            col) {
  # Parameter validation
  stopifnot(!missing(ad),
            !missing(col),
            is.character(col),
            col %in% names(ad$obs))
  
  # Filter controls
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

#' Run differential gene expression analysis (voom + lmFit)
#' 
#' @param ad AnnData object
#' @param obs Observations dataframe
#' @return Linear model fit object
run_dge <- function(ad,
                    obs) {
  # Parameter validation
  stopifnot(!missing(ad),
            !missing(obs),
            "cond" %in% names(obs),
            "plate" %in% names(obs),
            !is.null(ad$X))
  
  # build DGEList + design
  counts <- Matrix::t(ad$X)
  dge    <- DGEList(counts=counts)

  cat("    Constructing design matrix...\n")
  start_time <- Sys.time()
  design <- model.matrix(
    ~ 0 + cond + plate,
    data = obs
    )
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("    Design matrix is constructed in %.1f seconds\n", diff_time))
  cat("    Design matrix dimensions:", nrow(design), "samples ×", ncol(design), "coefficients\n")

  # filter genes and normalize
  cat("    Filtering genes and normalizing counts...\n")
  start_time <- Sys.time()
  keep <- filterByExpr(dge, design)
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  dge  <- dge[keep, , keep.lib.sizes=FALSE] %>% calcNormFactors()
  cat(sprintf("    Genes are filtered and normalized in %.1f seconds\n", diff_time))
  
  
  # voom + lmFit
  cat("    Running voom...\n")
  start_time <- Sys.time()
  v   <- voom(dge, design, plot=FALSE)
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("    Voom is completed in %.1f seconds\n", diff_time))

  cat("    Running lmFit...\n")
  start_time <- Sys.time()
  fit <- lmFit(v, design)
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("    LmFit is completed in %.1f seconds\n", diff_time))
  return(fit)
}

#' Merge command line arguments with configuration
#' 
#' @param args Command line arguments list
#' @param config Configuration list
#' @return Merged parameters list
merge_config <- function(args, config) {
  # Parameter validation
  stopifnot(!missing(args),
            !missing(config),
            is.list(args), is.list(config))
  
  config_merged <- args
  for (name in names(config)) {
    if (is.null(args[[name]])) {
      config_merged[[name]] <- config[[name]]
    }
  }
  return(config_merged)
}

#' Main function for Step 1: Prepare intermediate data for batch processing
#' 
#' Processes input data, runs voom+lmFit, and saves intermediate results
#' @return NULL (saves results to files)
main <- function() {
  

  parser <- ArgumentParser()
  parser$add_argument("--input", type="character")
  parser$add_argument("--output_dir", type="character", required=TRUE)
  parser$add_argument("--subsampling", action="store_true")
  parser$add_argument("--max_cell_types", type="integer")
  parser$add_argument("--max_perturbations", type="integer")
  parser$add_argument("--max_genes", type="integer")
  parser$add_argument("--specific_times", nargs = "+", type="integer", default=c(24))
  parser$add_argument("--specific_perturbagens", type="character", default=NULL)
  parser$add_argument("--min_cells", type="integer", default=0)
  parser$add_argument("--qc", action="store_true")
  parser$add_argument("--config", default="{}")

  args <- parser$parse_args()
  config <- jsonlite::fromJSON(args$config)
  args$config = NULL
  args <- merge_config(args, config)  
  
  # Set default min_cells if not provided
  if (is.null(args$min_cells)) {
    args$min_cells <- 0
  }
  stopifnot(is.numeric(args$min_cells), args$min_cells >= 0)
  
  # Set default qc_enabled if not provided
  if (is.null(args$qc)) {
    args$qc <- FALSE
  }
  stopifnot(is.logical(args$qc))
  
  # Start timer
  deg_start <- Sys.time()
  cat("\nDEG Step 1: Prepare - Started...\n")
  
  # Capture warnings during processing
  warning_messages <- list()
  
  # Get basename of input file (without extension) for intermediate file naming
  # Remove "_processed" from the basename if present (handles both "_processed" at end or followed by underscore)
  input_basename <- tools::file_path_sans_ext(basename(args$input))
  input_basename <- gsub("_processed_", "_", input_basename)  # Replace "_processed_" with "_"
  input_basename <- gsub("_processed$", "", input_basename)    # Remove "_processed" at the end
  
  # load the full pseudobulk AnnData
  adata <- anndata::read_h5ad(args$input)
  cat("Loaded input file: ", args$input, "\n")
  cat("Initial dimensions: ", adata$n_obs, " samples × ", adata$n_vars, " genes\n")
  
  adata <- filter_cells(adata, args$min_cells)
  adata <- filter_qc(adata, args$qc)
  adata <- subsampling(adata, args)
  
  # Get cell type (should be only one since input files are split by cell type)
  cell_types <- unique(as.character(adata$obs$cell_type))
  if (length(cell_types) != 1) {
    stop("Expected exactly one cell type in input file, found: ", paste(cell_types, collapse = ", "))
  }
  cl <- cell_types[1]
  
  
  # Process with warning capture
  withCallingHandlers({
    cat("\n▶︎ Processing cell_type: ", cl, "\n")
    
    # Filter controls
    ad <- filter_controls(adata, "pert_time_h")
    
    obs <- ad$obs %>%
      mutate(
        raw_cond = perturbation_label,
        cond = factor(clean(perturbation_label)),  # sanitize values here!
      )
    
    # Check for required controls
    check_for_controls(obs, "pert_time_h")
  
    cond <- obs$cond
    start_time <- Sys.time()
    fit <- run_dge(ad, obs)
    delta_time <- difftime(Sys.time(), start_time, units = "secs")
    cat(sprintf("DGE fit is completed in %.1f seconds\n", delta_time))

    # separate the control and treated conds (we won't DE on control-vs-control)
    control_obs <- obs %>%
        distinct(cond, .keep_all = TRUE) %>%
        filter(is_control==TRUE)
    
    treated_obs <- obs %>%
        distinct(cond, .keep_all = TRUE) %>%
        filter(is_control!=TRUE)
    
    # Use output_dir as the intermediate directory
    intermediate_dir <- args$output_dir
    
    # Create intermediate directory if it doesn't exist
    if (!dir.exists(intermediate_dir)) {
      dir.create(intermediate_dir, recursive = TRUE)
    }
    
    # Save intermediate data for batch processing (use basename of input file)
    intermediate_file <- file.path(intermediate_dir, paste0(input_basename, ".rds"))
    cat("\n    Saving intermediate data: ", intermediate_file, "\n")
    
    # Save fit object and metadata
    intermediate_data <- list(
      fit = fit,
      control_obs = control_obs,
      treated_obs = treated_obs,
      obs = obs,
      ad_var = ad$var,  # Save var for later use
      ad_uns = if (!is.null(adata$uns) && length(adata$uns) > 0) as.list(adata$uns) else list(),
      parameters = list(
        min_cells = args$min_cells,
        qc = args$qc,
        subsampling = args$subsampling
      ),
      genes = rownames(fit$coefficients)  # Save gene names
    )
    
    saveRDS(intermediate_data, intermediate_file, compress = "xz")
    
  }, warning = function(w) {
    warning_messages[[length(warning_messages) + 1]] <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  
  # Show runtime
  deg_time <- difftime(Sys.time(), deg_start, units = "secs")
  cat(sprintf("\nStep 1 Total runtime: %.1f seconds\n", deg_time))
  
  # Print any warnings that occurred (to stderr)
  if (length(warning_messages) > 0) {
    message(sprintf("\n=== %d Warnings encountered ===\n", length(warning_messages)))
    # Print first 20 unique warnings
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
  
  cat("\nDEG Step 1: Prepare - Completed!\n")
}

main()

