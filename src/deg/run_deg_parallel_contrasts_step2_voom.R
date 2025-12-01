#!/usr/bin/env Rscript

# Step 2: Voom - Transform counts using voom

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(argparse)
})

#' Apply voom transformation
#' 
#' @param dge DGEList object
#' @param design Design matrix
#' @return Voom object
run_voom <- function(dge,
                     design) {
  # Parameter validation
  stopifnot(!missing(dge),
            !missing(design),
            class(dge) == "DGEList",
            is.matrix(design))
  
  cat("    Running voom...\n")
  start_time <- Sys.time()
  v <- voom(dge, design, plot=FALSE)
  diff_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("    Voom is completed in %.1f seconds\n", diff_time))
  
  return(v)
}

#' Main function for Step 2: Apply voom transformation
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--input_file", type="character", required=TRUE,
                     help="Path to filtered DGEList RDS file from step 1")
  parser$add_argument("--output_dir", type="character", required=TRUE,
                     help="Directory to save voom object")

  args <- parser$parse_args()
  
  deg_start <- Sys.time()
  cat("\nDEG Step 2: Voom - Started...\n")
  
  cat("Loading filtered DGEList: ", args$input_file, "\n")
  intermediate_data <- readRDS(args$input_file)
  
  dge <- intermediate_data$dge
  design <- intermediate_data$design
  
  # Apply voom transformation
  v <- run_voom(dge, design)
  
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  
  input_basename <- tools::file_path_sans_ext(basename(args$input_file))
  input_basename <- gsub("_filtered$", "", input_basename)
  output_file <- file.path(args$output_dir, paste0(input_basename, "_voom.rds"))
  cat("\n    Saving voom object: ", output_file, "\n")
  
  voom_data <- list(
    voom = v,
    design = design,
    obs = intermediate_data$obs,
    ad_var = intermediate_data$ad_var,
    ad_uns = intermediate_data$ad_uns,
    parameters = intermediate_data$parameters
  )
  
  saveRDS(voom_data, output_file, compress = "xz")
  
  deg_time <- difftime(Sys.time(), deg_start, units = "secs")
  cat(sprintf("\nStep 2 Total runtime: %.1f seconds\n", deg_time))
  
  cat("\nDEG Step 2: Voom - Completed!\n")
}

main()

