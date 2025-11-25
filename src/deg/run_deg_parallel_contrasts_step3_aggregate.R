#!/usr/bin/env Rscript

# Step 3: Aggregate - Combine batch results and create final h5ad output

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
  library(argparse)
  library(rhdf5)
})

# helper to sanitize names for model.matrix/contrasts
clean <- function(x) gsub("[^[:alnum:]_]", "_", x)

# Define which DE results to store in layers
res_cols <- c("logFC", 
              "stdev.unscaled", 
              "stdev.scaled", 
              "CI.L", 
              "CI.R", 
              "AveExpr", 
              "t", 
              "P.Value", 
              "adj.P.Value.within_one_contrast", 
              "adj.P.Value.across_all_contrasts", 
              "B"
             )

#' Rewrite h5ad file by completely removing X matrix
#' 
#' Completely removes the X matrix from an h5ad file using direct HDF5 access.
#' This is equivalent to the Python rewrite_h5ad function and ensures older
#' Python anndata versions can read the file.
#' 
#' @param path2file Path to the h5ad file
#' @return NULL (modifies file in place)
rewrite_h5ad <- function(path2file) {
  # Parameter validation
  stopifnot(!missing(path2file),
            is.character(path2file),
            file.exists(path2file))
  
  cat("\n    Rewriting h5ad file:", path2file, "\n")
  
  # Check if file is valid HDF5
  if (!rhdf5::H5Fis_hdf5(path2file)) {
    stop("File is not a valid HDF5 file")
  }
  
  # Open file in read-write mode
  fid <- rhdf5::H5Fopen(path2file, "H5F_ACC_RDWR")
  
  # Check if X group exists and remove it completely
  if (rhdf5::H5Lexists(fid, "X")) {
    rhdf5::H5Ldelete(fid, "X")
    cat("    Removed X matrix\n")
  } else {
    cat("    X matrix not found or already removed\n")
  }
  
  # Close file
  rhdf5::H5Fclose(fid)
  
  cat("    Rewriting h5ad file is done!\n")
}

#' Extract observations for differential expression results
#' 
#' @param de_df Differential expression results dataframe
#' @param treated_obs Treated observations dataframe
#' @return Processed observations dataframe
get_obs <- function(de_df,
                    treated_obs) {
  # Parameter validation
  stopifnot(!missing(de_df),
            !missing(treated_obs),
            "control" %in% names(de_df), "cond" %in% names(de_df),
            "cond" %in% names(treated_obs))
  
  # FIX 1: Ensure dose_value and time are preserved as actual values, not factors
  de_unique <- de_df %>% distinct(control, cond)
  obs_out <- treated_obs %>%
    rownames_to_column("id") %>%
    right_join(de_unique, by = "cond") %>%                     # only contrasts with results; keep same order as de_df
    mutate(
      contrast=factor(paste0(cond,' - ', control)),
      pert_dose_uM = as.numeric(as.character(pert_dose_uM)),       # ensure numeric
      pert_time_h       = as.numeric(as.character(pert_time_h))
    ) %>%
    remove_rownames() %>%
    column_to_rownames("contrast") %>%
    select(-control) %>%
    select(-raw_cond) %>%
    select(-cond)
  return(obs_out)
}

#' Extract gene information for differential expression results
#' 
#' @param de_df Differential expression results dataframe
#' @param ad_var Variable information dataframe from original data
#' @return Gene information dataframe
get_var <- function(de_df,
                    ad_var) {
  # Parameter validation
  stopifnot(!missing(de_df),
            !missing(ad_var),
            "gene" %in% names(de_df),
            !is.null(ad_var), is.data.frame(ad_var))
  
  genes   <- unique(de_df$gene)
  # FIX 2: Carry over all columns from the original var dataframe
  # Get the original var data for these genes
  var_out <- ad_var[genes, , drop = FALSE]
  
  # If for some reason the genes aren't in the original var, create a basic var
  if (nrow(var_out) == 0 || !all(genes %in% rownames(var_out))) {
    message("  Warning: Some genes not found in original var, creating basic var dataframe")
    var_out <- data.frame(gene = genes, row.names = genes)
    # Try to merge with available var data
    if (nrow(ad_var) > 0) {
      available_genes <- intersect(genes, rownames(ad_var))
      if (length(available_genes) > 0) {
        var_out[available_genes, ] <- ad_var[available_genes, ]
      }
    }
  }
  return(var_out)
}

#' Create layers for differential expression statistics
#' 
#' @param de_df Differential expression results dataframe
#' @param obs_out Processed observations dataframe
#' @return List of matrices for each DE statistic
get_layers <- function(de_df,
                       obs_out) {
  # Parameter validation
  stopifnot(!missing(de_df),
            !missing(obs_out),
            "gene" %in% names(de_df), "cond" %in% names(de_df))
  
  # Create layers for each DE statistic
  layers <- map(res_cols, function(m) {
    if (!m %in% names(de_df)) {
      warning("Column ", m, " not found in DE results")
      return(NULL)
    }
    
    de_df %>%
      select(gene, cond, !!sym(m)) %>%
      pivot_wider(names_from = gene, values_from = !!sym(m)) %>%
      arrange(match(cond, rownames(obs_out))) %>%
      select(-cond) %>%
      as.matrix()
  }) %>% 
    set_names(res_cols) %>%
    compact()  # Remove NULL entries
  return(layers)
}

#' Main function for Step 3: Aggregate batch results
#' 
#' Combines all batch results, applies global p-value adjustment, and creates final h5ad
#' @return NULL (saves results to files)
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--input_file", type="character", required=TRUE,
                     help="Path to intermediate RDS file from step 1")
  parser$add_argument("--input_dir", type="character", required=TRUE,
                     help="Directory containing batch result files")
  parser$add_argument("--output_file", type="character", required=TRUE,
                     help="Path to output h5ad file")

  args <- parser$parse_args()
  
  # Start timer
  start_time <- Sys.time()
  cat("\nDEG Step 3: Aggregate - Started...\n")
  
  # Capture warnings during processing
  warning_messages <- list()
  
  # Process with warning capture
  withCallingHandlers({
    # Load intermediate data for metadata
    cat("Loading intermediate data: ", args$input_file, "\n")
    intermediate_data <- readRDS(args$input_file)
    
    treated_obs <- intermediate_data$treated_obs
    ad_var <- intermediate_data$ad_var
    ad_uns <- intermediate_data$ad_uns
    parameters <- intermediate_data$parameters
    
    # Find all batch result files
    batch_files <- list.files(args$input_dir, pattern = "^batch_[0-9]+\\.rds$", full.names = TRUE)
    n_batches <- length(batch_files)
    
    if (n_batches == 0) {
      stop("No batch result files found in: ", args$input_dir)
    }
    
    cat(sprintf("Found %d batch result files\n", n_batches))
    
    # Load and combine all batch results
    cat("Loading and combining batch results...\n")
    batch_results <- list()
    
    for (batch_file in batch_files) {
      cat("  Loading: ", basename(batch_file), "\n")
      batch_data <- readRDS(batch_file)
      if (!is.null(batch_data) && nrow(batch_data) > 0) {
        batch_results[[length(batch_results) + 1]] <- batch_data
      }
    }
    
    if (length(batch_results) == 0) {
      stop("No valid batch results found")
    }
    
    # Combine all batch results
    cat("Combining all batches...\n")
    de_res <- bind_rows(batch_results)
    
    cat(sprintf("Combined results: %d rows\n", nrow(de_res)))
    
    # Apply global p-value adjustment across all contrasts
    cat("Applying global p-value adjustment...\n")
    de_df <- de_res %>%
      mutate(
        adj.P.Value.across_all_contrasts = p.adjust(P.Value, method="BH")
      ) %>%
      rename(adj.P.Value.within_one_contrast = adj.P.Val)
    
    cat(sprintf("Creating output structures...\n"))
    obs_out <- get_obs(de_df, treated_obs)
    var_out <- get_var(de_df, ad_var)
    layers <- get_layers(de_df, obs_out)
    
    # Prepare uns metadata
    new_uns <- ad_uns
    new_uns$threshold_filter_cells <- parameters$min_cells
    new_uns$qc_filtering_enabled <- parameters$qc
    
    # Assemble and write AnnData
    cat("Creating AnnData object...\n")
    out_adata <- anndata::AnnData(
      obs    = obs_out,
      var    = var_out,
      layers = layers,
      uns    = new_uns
    )
    
    # Use output_file as provided
    outfile <- args$output_file
    
    # Create output directory if it doesn't exist
    output_dir <- dirname(outfile)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    cat("\n    Writing: ", outfile, "\n")
    out_adata$write_h5ad(outfile, compression = "gzip")
    
    # Remove X matrix for compatibility with older Python anndata versions
    rewrite_h5ad(outfile)
  }, warning = function(w) {
    warning_messages[[length(warning_messages) + 1]] <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  
  # Show runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "secs")
  cat(sprintf("\nStep 3 completed in %.1f seconds\n", runtime))
  
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
  
  cat("\nDEG Step 3: Aggregate - Completed!\n")
}

main()

