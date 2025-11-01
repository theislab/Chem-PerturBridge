#!/usr/bin/env Rscript

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
})

source("./src/deg/subsampling.R")
# helper to sanitize names for model.matrix/contrasts
clean <- function(x) gsub("[^[:alnum:]_]", "_", x)

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
  
  cat("Rewriting h5ad file:", path2file, "\n")
  
  # Check if file is valid HDF5
  if (!rhdf5::H5Fis_hdf5(path2file)) {
    stop("File is not a valid HDF5 file")
  }
  
  # Open file in read-write mode
  fid <- rhdf5::H5Fopen(path2file, "H5F_ACC_RDWR")
  
  # Check if X group exists and remove it completely
  if (rhdf5::H5Lexists(fid, "X")) {
    rhdf5::H5Ldelete(fid, "X")
    cat("Removed X matrix\n")
  } else {
    cat("X matrix not found or already removed\n")
  }
  
  # Close file
  rhdf5::H5Fclose(fid)
  
  cat("Rewriting h5ad file is done!\n")
}

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
  adata <- adata[adata$obs$psbulk_cells >= min_cells]
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

#' Derive differential expression results table
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

#' Run differential expression contrasts
#' 
#' @param fit Linear model fit object
#' @param control_obs Control observations dataframe
#' @param treated_obs Treated observations dataframe
#' @param par Parameters list
#' @return Dataframe with differential expression results
run_contrasts <- function(fit,
                          control_obs,
                          treated_obs,
                          par) {
  # Parameter validation
  stopifnot(!missing(fit),
            !missing(control_obs),
            !missing(treated_obs),
            !missing(par),
            is.list(par))
  
  # run one contrast per treated cond vs. the matching control@same time
  n_contrasts <- nrow(treated_obs)
    contrast_num <- 0
    
    de_res <- pmap_dfr(treated_obs, function(cond, perturbagen, pert_dose_uM, pert_time_h, plate, raw_cond, ...) {
      contrast_num <<- contrast_num + 1
      if (contrast_num %% 5 == 1 || contrast_num == n_contrasts) {
        cat(sprintf("    Running contrast %d/%d\n", contrast_num, n_contrasts))
      }

      raw_controls <- find_controls(control_obs,
                    pert_time_h)
      # Check if control exists for this time point
      if (is.null(raw_controls)) {
        return(NULL)
      }
      top_tables <- list()
        
      for (raw_control in raw_controls) {
          l <- length(top_tables)
          top_tables[[l + 1]] <- derive_top_table(
                                  fit,
                                  raw_control,
                                  raw_cond,
                                  perturbagen,
                                  pert_dose_uM,
                                  pert_time_h,
                                  par$subsampling)
        }
      top_tables <- compact(top_tables)
      if (length(top_tables) > 0) {
        tryCatch({
          return(bind_rows(top_tables))
        }, error = function(e) {
          warning("Error binding top_tables: ", e$message)
          return(NULL)
        })
      } else {
        return(NULL)
      }
    })
    
    
    return(de_res)
    
}

#' Run differential gene expression analysis
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
  design <- model.matrix(
    ~ 0 + cond + plate,
    data = obs
    )
  # filter genes and normalize
  keep <- filterByExpr(dge, design)
  dge  <- dge[keep, , keep.lib.sizes=FALSE] %>% calcNormFactors()
  
  # voom + lmFit
  v   <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
  return(fit)
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
#' @param ad AnnData object
#' @return Gene information dataframe
get_var <- function(de_df,
                    ad) {
  # Parameter validation
  stopifnot(!missing(de_df),
            !missing(ad),
            "gene" %in% names(de_df),
            !is.null(ad$var), is.data.frame(ad$var))
  
  genes   <- unique(de_df$gene)
  # FIX 2: Carry over all columns from the original var dataframe
  # Get the original var data for these genes
  original_var <- ad$var
  var_out <- original_var[genes, , drop = FALSE]
  
  # If for some reason the genes aren't in the original var, create a basic var
  if (nrow(var_out) == 0 || !all(genes %in% rownames(var_out))) {
    message("  Warning: Some genes not found in original var, creating basic var dataframe")
    var_out <- data.frame(gene = genes, row.names = genes)
    # Try to merge with available var data
    if (nrow(original_var) > 0) {
      available_genes <- intersect(genes, rownames(original_var))
      if (length(available_genes) > 0) {
        var_out[available_genes, ] <- original_var[available_genes, ]
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



#' Main differential gene expression pipeline
#' 
#' @param par Parameters list containing input/output paths and analysis settings
#' @return NULL (saves results to files)
run_dge_pipeline <- function(par) {
  # Parameter validation
  stopifnot(!missing(par),
            is.list(par),
            "input" %in% names(par), !is.null(par$input),
            "output_dir" %in% names(par), !is.null(par$output_dir),
            is.character(par$input), file.exists(par$input),
            is.character(par$output_dir))
  
  # Set default min_cells if not provided
  if (is.null(par$min_cells)) {
    par$min_cells <- 0
  }
  stopifnot(is.numeric(par$min_cells), par$min_cells >= 0)
  
  # Set default qc_enabled if not provided
  if (is.null(par$qc)) {
    par$qc <- FALSE
  }
  stopifnot(is.logical(par$qc))
  
  # load the full pseudobulk AnnData
  adata <- anndata::read_h5ad(par$input)
  cat("Loaded input file: ", par$input, "\n")
  cat("Initial dimensions: ", adata$n_obs, " samples × ", adata$n_vars, " genes\n")
  
  adata <- filter_cells(adata, par$min_cells)
  adata <- filter_qc(adata, par$qc)
  adata <- subsampling(adata, par)
  cell_types_to_process <- unique(adata$obs$cell_type)
  n_cell_types <- length(cell_types_to_process)
  cl_num <- 0

  for (cl in cell_types_to_process) {
    cl_num <- cl_num + 1
    cl_start <- Sys.time()
    cat("\n▶︎ Processing cell_type ", cl_num, "/", n_cell_types, ": ", cl, "\n")
  
    # subset to this cell_type
    ad  <- adata[adata$obs$cell_type == cl, ]
    ad <- filter_controls(ad, "pert_time_h")
    
    obs <- ad$obs %>%
      mutate(
      raw_cond = perturbation_label,
      cond = factor(clean(perturbation_label)),  # sanitize values here!
      )
    
    # Check for required controls
    check_for_controls(obs, "pert_time_h")
  
    cond <- obs$cond
    fit <- run_dge(ad,
            obs)
  
    # separate the control and treated conds (we won't DE on control-vs-control)
    control_obs <- obs %>%
        distinct(cond, .keep_all = TRUE) %>%
        filter(is_control==TRUE)
    
    treated_obs <- obs %>%
        distinct(cond, .keep_all = TRUE) %>%
        filter(is_control!=TRUE)

    de_res <- run_contrasts(fit,
                          control_obs,
                          treated_obs,
                          par)

    # Skip if no valid DE results
    if (is.null(de_res) || nrow(de_res) == 0) {
      warning("No valid DE results for cell line ", cl)
      next
    }
  
    # adjust p-values globally across all contrasts
    de_df <- de_res %>%
      mutate(
        adj.P.Value.across_all_contrasts = p.adjust(P.Value, method="BH")
      )  %>%
      rename(adj.P.Value.within_one_contrast = adj.P.Val)
  
    obs_out <- get_obs(de_df, treated_obs)
    var_out <- get_var(de_df, ad)
    layers <- get_layers(de_df, obs_out)
    
    # Get uns metadata before subsetting to avoid ImplicitModificationWarning
    # Access uns early while adata is still a proper object (not a view)
    if (!is.null(adata$uns) && length(adata$uns) > 0) {
      new_uns <- as.list(adata$uns)  # Create a copy
    } else {
      new_uns <- list()
    }
  
    new_uns$threshold_filter_cells <- par$min_cells
    new_uns$qc_filtering_enabled <- par$qc
    

    
    # assemble and write
    out_adata <- anndata::AnnData(
      obs    = obs_out,
      var    = var_out,
      layers = layers,
      uns    = new_uns
    )
  
    # Create output directory if it doesn't exist
    if (!is.null(par$subsampling) && par$subsampling) {
      output_dir <- file.path(par$output_dir, "subsampling")
    }
    else {
      output_dir <- file.path(par$output_dir, "full")
    }
    # Add QC folder based on whether QC filtering is enabled
    qc_folder <- if (par$qc) "qc_true" else "qc_false"
    output_dir <- file.path(output_dir, qc_folder)
    
    # Add filter folder based on min_cells threshold
    filter_folder <- paste0("filter_min_cells_", par$min_cells)
    output_dir <- file.path(output_dir, filter_folder)
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  
    outfile <- file.path(output_dir, paste0(cl, "_de.h5ad"))
    cat("\n    Writing: ", outfile, "\n")
    out_adata$write_h5ad(outfile, compression = "gzip")
    
    # Remove X matrix for compatibility with older Python anndata versions
    rewrite_h5ad(outfile)
  
    # Show time for this cell line
    cl_time <- difftime(Sys.time(), cl_start, units = "secs")
    cat(sprintf("  ✓ Cell line completed in %.1f seconds\n", cl_time))
  }
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

#' Main function for differential gene expression analysis
#' 
#' Processes command line arguments and runs the complete DGE pipeline
#' @return NULL (saves results to files)
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--input", type="character")
  parser$add_argument("--output_dir", type="character")
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
  # Start timer
  start_time <- Sys.time()
  cat("\nDE analysis started...\n")
	
  tryCatch({
    	run_dge_pipeline(args)
      cat("DE analysis completed!\n")
  }, error = function(e) {
    	message("Error in running DGE: ", e$message)
  	message("Stack trace:")
	traceback()
	stop("DGE pipeline failed")
  })
  
	# Show runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "secs")
  cat(sprintf("\nTotal runtime: %.1f seconds\n", runtime))

  if (!is.null(args$subsampling) && args$subsampling) {
  	cat("\n📝 Note: This was a TEST RUN with reduced data.\n")
  	cat("Set subsampling = FALSE for full analysis.\n\n")
  }

}

main()
