#!/usr/bin/env Rscript

# Step 2: Batch - Process contrasts in parallel batches on SLURM

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(limma)
  library(argparse)
})

# helper to sanitize names for model.matrix/contrasts
clean <- function(x) gsub("[^[:alnum:]_]", "_", x)

#' Find control conditions for a specific time point
#' 
#' @param control_obs Control observations dataframe
#' @param time Time point in hours
#' @return Vector of control condition names or NULL if none found
find_controls <- function(control_obs,
                          time) {
  # Parameter validation
  stopifnot(!missing(control_obs),
            !missing(time),
            is.numeric(time),
            "pert_time_h" %in% names(control_obs))
    
    select_ctrs <- (control_obs$pert_time_h == time)

    if (!any(select_ctrs)) {
        warning("No control found for time ", time, "h")
        return(NULL)
        }
    
    raw_controls <- unique(control_obs[select_ctrs, ]$cond)
    return(raw_controls)
    
}

#' Derive differential expression results table (single contrast)
#' 
#' @param fit Linear model fit object
#' @param raw_control Control condition name
#' @param raw_cond Treatment condition name
#' @param perturbagen Perturbagen name
#' @param pert_dose_uM Dose concentration
#' @param pert_time_h Time point
#' @param subsampling Whether subsampling mode is enabled
#' @return Dataframe with differential expression results or NULL if error
derive_top_table <- function(fit,
                             raw_control,
                             raw_cond,
                             perturbagen,
                             pert_dose_uM,
                             pert_time_h,
                             subsampling) {
  # Parameter validation
  stopifnot(!missing(fit),
            !missing(raw_control),
            !missing(raw_cond),
            !missing(perturbagen),
            !missing(pert_dose_uM),
            !missing(pert_time_h))
  
    contrast <- paste0("cond", clean(raw_cond), " - cond", clean(raw_control))
    tryCatch({
        ctr <- makeContrasts(contrasts = contrast, levels = colnames(coef(fit)))
        
        fit2 <- contrasts.fit(fit, ctr) %>% 
          eBayes(robust = !isTRUE(subsampling))  # Skip robust for speed in subsampling mode
        
        top_table <- topTable(fit2, number = Inf, sort = "none", adjust.method="BH", confint=TRUE) %>%
          rownames_to_column("gene") %>%
          mutate(
            control = clean(raw_control),
            cond = clean(raw_cond),
            perturbagen = perturbagen,
            pert_dose_uM = pert_dose_uM,
            pert_time_h = pert_time_h
          )

        
        stdev_unscaled <- data.frame(fit2$stdev.unscaled, check.names = FALSE) %>%
            rownames_to_column("gene") %>%
            rename(stdev.unscaled = !!sym(contrast))

        stdev_scaled <- data.frame(fit2$stdev.unscaled * sqrt(fit2$s2.post), check.names = FALSE) %>%
            rownames_to_column("gene") %>%
            rename(stdev.scaled = !!sym(contrast))

        top_table <- left_join(top_table, stdev_unscaled, by = "gene")
        top_table <- left_join(top_table, stdev_scaled, by = "gene")
        return(top_table)
        
      }, error = function(e) {
        warning("Error in contrast for ", raw_cond, ": ", e$message)
        return(NULL)
      })
}

#' Process a batch of contrasts
#' 
#' @param fit Linear model fit object
#' @param control_obs Control observations dataframe
#' @param treated_obs Treated observations dataframe (all rows)
#' @param batch_indices Vector of indices to process in this batch
#' @param subsampling Whether subsampling mode is enabled
#' @return Dataframe with differential expression results for this batch
process_batch <- function(fit,
                         control_obs,
                         treated_obs,
                         batch_indices,
                         subsampling) {
  # Parameter validation
  stopifnot(!missing(fit),
            !missing(control_obs),
            !missing(treated_obs),
            !missing(batch_indices))
  
  cat(sprintf("    Processing batch with %d contrasts (indices: %d to %d)\n", 
              length(batch_indices), min(batch_indices), max(batch_indices)))
  
  # Pre-compute control lookup for efficiency
  control_lookup <- split(control_obs, control_obs$pert_time_h)
  
  batch_results <- list()
  n_contrasts <- length(batch_indices)
  contrast_num <- 0
  
  for (i in batch_indices) {
    contrast_num <- contrast_num + 1
    
    # Print progress every 50 contrasts or at the end
    if (contrast_num %% 50 == 1 || contrast_num == n_contrasts) {
      cat(sprintf("    Processing contrast %d/%d (index %d)\n", contrast_num, n_contrasts, i))
    }
    
    row <- treated_obs[i, ]
    cond <- row$cond
    perturbagen <- row$perturbagen
    pert_dose_uM <- row$pert_dose_uM
    pert_time_h <- row$pert_time_h
    plate <- row$plate
    raw_cond <- row$raw_cond
    
    # Get controls for this time point (using pre-computed lookup)
    raw_controls <- control_lookup[[as.character(pert_time_h)]]
    if (is.null(raw_controls) || nrow(raw_controls) == 0) {
      next
    }
    
    # Extract unique control condition names
    raw_control_names <- unique(raw_controls$cond)
    if (length(raw_control_names) == 0) {
      next
    }
    
    # Process all controls for this contrast
    top_tables <- lapply(raw_control_names, function(raw_control) {
      derive_top_table(
        fit,
        raw_control,
        raw_cond,
        perturbagen,
        pert_dose_uM,
        pert_time_h,
        subsampling
      )
    })
    
    top_tables <- compact(top_tables)
    if (length(top_tables) > 0) {
      result <- tryCatch({
        bind_rows(top_tables)
      }, error = function(e) {
        warning(sprintf("Error combining results for contrast %d (%s): %s", 
                       i, as.character(cond), e$message))
        return(NULL)
      })
      
      if (!is.null(result)) {
        batch_results[[length(batch_results) + 1]] <- result
      }
    }
  }
  
  # Combine all results for this batch
  if (length(batch_results) > 0) {
    return(bind_rows(batch_results))
  } else {
    return(NULL)
  }
}

#' Main function for Step 2: Process a single batch of contrasts
#' 
#' Loads intermediate data, processes assigned batch, and saves results
#' @return NULL (saves results to files)
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--intermediate_file", type="character", required=TRUE,
                     help="Path to intermediate RDS file from step 1")
  parser$add_argument("--output_dir", type="character", required=TRUE,
                     help="Directory to save batch results")
  parser$add_argument("--batch_id", type="integer", required=TRUE,
                     help="Batch ID (0-based index)")
  parser$add_argument("--n_batches", type="integer", required=TRUE,
                     help="Total number of batches")

  args <- parser$parse_args()
  
  # Start timer
  start_time <- Sys.time()
  cat(sprintf("\nDEG Step 2: Batch %d/%d - Started...\n", args$batch_id + 1, args$n_batches))
  
  # Load intermediate data
  cat("Loading intermediate data: ", args$intermediate_file, "\n")
  intermediate_data <- readRDS(args$intermediate_file)
  
  fit <- intermediate_data$fit
  control_obs <- intermediate_data$control_obs
  treated_obs <- intermediate_data$treated_obs
  subsampling <- intermediate_data$parameters$subsampling
  
  n_contrasts <- nrow(treated_obs)
  cat(sprintf("Total contrasts to process: %d\n", n_contrasts))
  
  # Calculate batch indices
  batch_size <- ceiling(n_contrasts / args$n_batches)
  start_idx <- args$batch_id * batch_size + 1
  end_idx <- min(start_idx + batch_size - 1, n_contrasts)
  
  if (start_idx > n_contrasts) {
    cat(sprintf("Batch %d has no contrasts to process (start_idx %d > n_contrasts %d)\n", 
                args$batch_id, start_idx, n_contrasts))
    # Create empty output file to indicate completion
    output_file <- file.path(args$output_dir, sprintf("batch_%04d.rds", args$batch_id))
    saveRDS(NULL, output_file, compress = "xz")
    return()
  }
  
  batch_indices <- start_idx:end_idx
  cat(sprintf("Processing indices %d to %d (%d contrasts)\n", 
              start_idx, end_idx, length(batch_indices)))
  
  # Capture warnings during processing
  warning_messages <- list()
  
  # Process this batch with warning capture
  de_res <- withCallingHandlers(
    process_batch(
      fit = fit,
      control_obs = control_obs,
      treated_obs = treated_obs,
      batch_indices = batch_indices,
      subsampling = subsampling
    ),
    warning = function(w) {
      warning_messages[[length(warning_messages) + 1]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  
  # Save batch results
  output_file <- file.path(args$output_dir, sprintf("batch_%04d.rds", args$batch_id))
  cat("Saving batch results: ", output_file, "\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  
  saveRDS(de_res, output_file, compress = "xz")
  
  # Show runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "secs")
  cat(sprintf("\nBatch %d completed in %.1f seconds\n", args$batch_id, runtime))
  
  # Print any warnings that occurred
  if (length(warning_messages) > 0) {
    cat(sprintf("\n=== %d Warnings encountered ===\n", length(warning_messages)))
    # Print first 20 unique warnings
    unique_warns <- unique(unlist(warning_messages))
    n_show <- min(20, length(unique_warns))
    for (i in 1:n_show) {
      cat(sprintf("%d. %s\n", i, unique_warns[i]))
    }
    if (length(unique_warns) > n_show) {
      cat(sprintf("... and %d more unique warning types (total %d warnings)\n", 
                  length(unique_warns) - n_show, length(warning_messages)))
    }
  }
  
  cat(sprintf("\nDEG Step 2: Batch %d/%d - Completed!\n", args$batch_id + 1, args$n_batches))
}

main()

