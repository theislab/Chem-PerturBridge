#!/usr/bin/env Rscript

# Step 3: lmFit - Fit linear model

suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
  library(argparse)
})

# helper to sanitize names for model.matrix/contrasts
clean <- function(x) gsub("[^[:alnum:]_]", "_", x)

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

#' Fit linear model
#' 
#' @param v Voom object
#' @param design Design matrix
#' @param obs Observations dataframe
#' @return List containing fit object, control_obs, and treated_obs
run_lmfit <- function(v,
                      design,
                      obs) {
  # Parameter validation
  stopifnot(!missing(v),
            !missing(design),
            !missing(obs),
            is.matrix(design),
            "is_control" %in% names(obs))
  
  # Check for required controls
  check_for_controls(obs, "pert_time_h")
  
  # Fit linear model
  cat("    Running lmFit...\n")
  start_time <- Sys.time()
  fit <- lmFit(v, design)
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("    LmFit is completed in %.1f seconds\n", diff_time))
  
  # Separate the control and treated conds
  control_obs <- obs %>%
    distinct(cond, .keep_all = TRUE) %>%
    filter(is_control==TRUE)
  
  treated_obs <- obs %>%
    distinct(cond, .keep_all = TRUE) %>%
    filter(is_control!=TRUE)
  
  return(list(
    fit = fit,
    control_obs = control_obs,
    treated_obs = treated_obs
  ))
}

#' Main function for Step 3: Fit linear model
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--input_file", type="character", required=TRUE,
                     help="Path to voom RDS file from step 2")
  parser$add_argument("--output_dir", type="character", required=TRUE,
                     help="Directory to save fit object")

  args <- parser$parse_args()
  
  deg_start <- Sys.time()
  cat("\nDEG Step 3: lmFit - Started...\n")
  
  cat("Loading voom object: ", args$input_file, "\n")
  voom_data <- readRDS(args$input_file)
  
  v <- voom_data$voom
  design <- voom_data$design
  obs <- voom_data$obs
  
  # Fit linear model
  fit_results <- run_lmfit(v, design, obs)
  fit <- fit_results$fit
  control_obs <- fit_results$control_obs
  treated_obs <- fit_results$treated_obs
  
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  
  input_basename <- tools::file_path_sans_ext(basename(args$input_file))
  input_basename <- gsub("_voom$", "", input_basename)
  output_file <- file.path(args$output_dir, paste0(input_basename, ".rds"))
  cat("\n    Saving fit object: ", output_file, "\n")
  
  intermediate_data <- list(
    fit = fit,
    control_obs = control_obs,
    treated_obs = treated_obs,
    obs = obs,
    ad_var = voom_data$ad_var,
    ad_uns = voom_data$ad_uns,
    parameters = voom_data$parameters,
    genes = rownames(fit$coefficients)
  )
  
  saveRDS(intermediate_data, output_file, compress = "xz")
  
  deg_time <- difftime(Sys.time(), deg_start, units = "secs")
  cat(sprintf("\nStep 3 Total runtime: %.1f seconds\n", deg_time))
  
  cat("\nDEG Step 3: lmFit - Completed!\n")
}

main()

