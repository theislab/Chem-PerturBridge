#!/usr/bin/env Rscript

# Configure reticulate to use conda Python
suppressPackageStartupMessages(library(reticulate))
conda_env <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_env)) use_python(file.path(conda_env, "bin/python"), required = TRUE)

requireNamespace("anndata", quietly = TRUE)
suppressPackageStartupMessages({
  library(Matrix)
  library(anndata)
  library(argparse)
})


#' Convert a single plate entry to a sparse count matrix and obs data.frame
#'
#' @param plate_data List with elements UMI.counts (genes x samples) and
#'   Annotation (samples x metadata)
#' @return List with:
#'   \item{X}{sparse Matrix of shape (n_samples x n_genes)}
#'   \item{obs}{data.frame of shape (n_samples x 18), rownames = sample barcodes}
plate_to_parts <- function(plate_data) {
  counts <- plate_data$UMI.counts   # genes x samples dense matrix
  obs    <- plate_data$Annotation   # samples x metadata

  # Reorder Annotation rows to match the column order of UMI.counts
  obs <- obs[colnames(counts), , drop = FALSE]

  # Transpose to samples x genes and convert to sparse format
  X_sparse <- Matrix(t(counts), sparse = TRUE)

  list(X = X_sparse, obs = obs)
}


#' Build a combined AnnData from the full Exp list
#'
#' Exp structure:
#'   Exp[[batch_id]][[plate_barcode]]$UMI.counts  # genes x samples
#'   Exp[[batch_id]][[plate_barcode]]$Annotation  # samples x 18 metadata cols
#'
#' @param Exp Named list of batches, each containing named plates
#' @param output_path Optional path to write the .h5ad file
#' @return AnnData object (samples x genes)
exp_to_anndata <- function(Exp, output_path = NULL) {

  X_list    <- list()
  obs_list  <- list()
  var_names <- NULL

  batch_ids <- names(Exp)
  n_batches <- length(batch_ids)

  for (i in seq_along(batch_ids)) {
    batch_id <- batch_ids[[i]]
    batch    <- Exp[[batch_id]]
    message(sprintf("[%d/%d] batch %s", i, n_batches, batch_id))

    for (plate_id in names(batch)) {
      key         <- paste(batch_id, plate_id, sep = "_")
      parts       <- plate_to_parts(batch[[plate_id]])
      plate_genes <- colnames(parts$X)

      # Fix gene order: use first plate as reference, reorder subsequent plates
      if (is.null(var_names)) {
        var_names <- plate_genes
      } else {
        parts$X <- parts$X[, var_names, drop = FALSE]
      }

      X_list[[key]]   <- parts$X
      obs_list[[key]] <- parts$obs
    }
  }

  # Stack all samples
  X_combined   <- do.call(rbind, X_list)    # total_samples x genes
  obs_combined <- do.call(rbind, obs_list)  # total_samples x 18

  var_df <- data.frame(row.names = var_names)

  adata <- AnnData(
    X   = X_combined,
    obs = obs_combined,
    var = var_df
  )

  if (!is.null(output_path)) {
    adata$write_h5ad(output_path, compression = "gzip")
    message("Saved AnnData to ", output_path)
  }

  adata
}


#' Main entry point
#'
#' Parses command-line arguments and runs the RData -> AnnData conversion
#' and gene annotation export.
#' @return NULL (saves results to files)
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--input_file",       type = "character", required = TRUE,
                      help = "Path to input .RData file containing the Exp object")
  parser$add_argument("--output_file",      type = "character", required = TRUE,
                      help = "Path for the output .h5ad file")
  parser$add_argument("--annotation_file",  type = "character", required = TRUE,
                      help = "Path to drugseq_ensembl_v98_annotation_and_entrez_mapping.RData")
  parser$add_argument("--genes_csv",        type = "character", required = TRUE,
                      help = "Path for the output drugseq_ensg_v98 CSV file")

  args <- parser$parse_args()

  start_time <- Sys.time()
  cat("\nRData -> AnnData conversion started...\n")

  tryCatch({
    message("Loading ", args$input_file)
    load(args$input_file)   # loads Exp into the current environment

    adata <- exp_to_anndata(Exp, output_path = args$output_file)

    cat(sprintf("obs : %d samples x %d metadata columns\n",
                nrow(adata$obs), ncol(adata$obs)))
    cat(sprintf("var : %d genes\n", nrow(adata$var)))
    cat("\nConversion completed!\n")
  }, error = function(e) {
    message("Error in conversion: ", e$message)
    stop("Conversion failed")
  })

  tryCatch({
    message("Loading gene annotation from ", args$annotation_file)
    load(args$annotation_file)   # loads drugseq_ensg_v98 into the current environment

    write.csv(drugseq_ensg_v98, args$genes_csv, row.names = FALSE)
    message("Saved gene annotation CSV to ", args$genes_csv)
  }, error = function(e) {
    message("Error saving gene annotation: ", e$message)
    stop("Gene annotation export failed")
  })

  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("\nTotal runtime: %.1f seconds\n", elapsed))
}

main()
