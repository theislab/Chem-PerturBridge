#!/usr/bin/env Rscript
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
})

source("/src/deg/subsampling.R")
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
	
filter_cells <- function(adata,
                         min_cells) {
  n_obs_prev <- adata$n_obs
  message("Filter samples with the low number of cells")
  adata <- adata[adata$obs$psbulk_cells >= min_cells]
  message("    n_obs: ", n_obs_prev, "--> ", adata$n_obs)
  return(adata)
}

check_for_controls <- function(obs, 
                               col) {
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

filter_controls <- function(ad,
                            col) {
  # Filter controls
  message("Filter controls by ", col, " column")
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
  
  message("    n_obs: ", n_obs_prev, "--> ",ad$n_obs)
  return(ad)
}

find_controls <- function(control_obs,
                          time) {
    
    select_ctrs <- (control_obs$pert_time_h == time)

    if (!any(select_ctrs)) {
        warning("No control found for time ", pert_time_h, "h, skipping ", raw_cond)
        return(NULL)
        }
    
    raw_controls <- unique(control_obs[select_ctrs, ]$cond)
    return(raw_controls)
    
}

derive_top_table <- function(fit,
                             raw_control,
                             raw_cond,
                             perturbagen,
                             pert_dose_uM,
                             pert_time_h) {
    contrast <- paste0("cond", clean(raw_cond), " - cond", clean(raw_control))
    tryCatch({
        ctr <- makeContrasts(contrasts = contrast, levels = colnames(coef(fit)))
        
        fit2 <- contrasts.fit(fit, ctr) %>% 
          eBayes(robust = !isTRUE(par$subsampling))  # Skip robust for speed in subsampling mode
        
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

run_contrasts <- function(fit,
                          control_obs,
                          treated_obs) {
    # run one contrast per treated cond vs. the matching control@same time
    n_contrasts <- nrow(treated_obs)
    contrast_num <- 0
    
    de_res <- pmap_dfr(treated_obs, function(cond, perturbagen, pert_dose_uM, pert_time_h, plate, raw_cond, ...) {
      contrast_num <<- contrast_num + 1
      if (contrast_num %% 5 == 1 || contrast_num == n_contrasts) {
        message(sprintf("    Running contrast %d/%d", contrast_num, n_contrasts))
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
                                  pert_time_h)
        }
      return(bind_rows(top_tables))
    })
    
    
    return(de_res)
    
}

run_dge <- function(ad,
                    obs) {
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

get_obs <- function(de_df,
                    treated_obs) {
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

get_var <- function(de_df,
                    ad) {
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

get_layers <- function(de_df) {
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



run_dge_pipeline <- function(par) {
  
  # load the full pseudobulk AnnData
  adata <- anndata::read_h5ad(par$input)
  adata <- filter_cells(adata, par$min_cells)
  adata <- subsampling(adata, par)
  cell_types_to_process <- unique(adata$obs$cell_type)
  n_cell_types <- length(cell_types_to_process)
  cl_num <- 0

  for (cl in cell_types_to_process) {
    cl_num <- cl_num + 1
    cl_start <- Sys.time()
    message("\n▶︎ Processing cell_type ", cl_num, "/", n_cell_types, ": ", cl)
  
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
                          treated_obs)
    
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
  
    layers <- get_layers(de_df)
  
    # carry over global uns if you like
    
    new_uns <- adata$uns
    print(par$output_dir) 
    
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
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  
    outfile <- file.path(output_dir, paste0(cl, "_de.h5ad"))
      
    message("  Writing: ", outfile)
    out_adata$write_h5ad(outfile, compression = "gzip")
  
    # Show time for this cell line
    cl_time <- difftime(Sys.time(), cl_start, units = "secs")
    message(sprintf("  ✓ Cell line completed in %.1f seconds", cl_time))
  }
}

merge_config <- function(args, config) {
  config_merged <- args
  for (name in names(config)) {
    if (is.null(args[[name]])) {
      config_merged[[name]] <- config[[name]]
    }
  }
  return(config_merged)
}

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
	parser$add_argument("--config", default="{}")

	args <- parser$parse_args()
	config <- jsonlite::fromJSON(args$config)
	args$config = NULL
	args <- merge_config(args, config)

	# Start timer
	start_time <- Sys.time()
	message("DE analysis started...")
	run_dge_pipeline(args)

	message("DE analysis completed!")

	# Show runtime
	end_time <- Sys.time()
	runtime <- difftime(end_time, start_time, units = "secs")
	message(sprintf("\nTotal runtime: %.1f seconds", runtime))

	if (!is.null(args$subsampling) && args$subsampling) {
  		message("\n📝 Note: This was a TEST RUN with reduced data.")
  		message("Set subsampling = FALSE for full analysis.")
	}

}

main()
