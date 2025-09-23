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
})

subsampling <- function(adata, par) {
  	# Subset cell lines
  	cell_types <- unique(adata$obs$cell_type)
	if (!is.null(par$max_cell_types) && length(cell_types) > par$max_cell_types) {
    	cell_types <- cell_types[1:par$max_cell_types]
    	message("Subsetting to ", par$max_cell_types, " cell lines: ", paste(cell_types, collapse=", "))
    	adata <- adata[adata$obs$cell_type %in% cell_types, ]
  	}

  
 	# Subset time points
  	if (!is.null(par$specific_times)) {
    	message("Filtering to time points: ", paste(par$specific_times, collapse=", "), " hours")
    	adata <- adata[adata$obs$pert_time_h %in% par$specific_times, ]
  	}
  
  	# Subset perturbations 
  	if (!is.null(par$specific_perturbagens)) {
    	# Use specific perturbagens if provided
    	message("Using specific perturbagens: ", paste(par$specific_perturbagens, collapse=", "))
   		adata[(adata$obs$is_control == TRUE)|(adata$obs$perturbagen %in% par$specific_perturbagens)]
  	} else if (!is.null(par$max_perturbations)) {
    	# Otherwise limit to max_perturbations
    	perturbagens <- unique(adata$obs$perturbagen)
    	control_perturbagens <- unique(adata[(adata$obs$is_control == TRUE)]$obs$perturbagen)
    	treated_perts <- setdiff(perturbagens, control_perturbagens)
    	if (length(treated_perts) > par$max_perturbations) {
      		keep_perts <- c(levels(control_perturbagens), treated_perts[1:par$max_perturbations])
      		message("Subsetting to ", par$max_perturbations, " perturbations (+ control): ", 
              		paste(keep_perts, collapse=", "))
      		adata <- adata[adata$obs$perturbagen %in% keep_perts, ]
    	}
  	}
  
  	# Subset genes - keep most variable
  	if (!is.null(par$max_genes) && ncol(adata) > par$max_genes) {
		message("Calculating gene variance for subsetting...")
		gene_vars <- apply(as.matrix(adata$X), 2, var)
		top_genes <- order(gene_vars, decreasing = TRUE)[1:par$max_genes]
		adata <- adata[, top_genes]
		message("Subset to top ", par$max_genes, " most variable genes")
  	}
  
	message("Test mode dataset: ", nrow(adata$obs), " samples x ", ncol(adata), " genes")
  
  	# Show analysis overview
  	analysis_summary <- adata$obs %>%
		group_by(cell_type, perturbagen, pert_time_h) %>%
		summarise(n_samples = n(), .groups = "drop") %>%
		arrange(cell_type, pert_time_h, perturbagen)
  
  	message("\nConditions to analyze:")
  	print(analysis_summary)
  	message("")
}
	


start_deg <- function(par) {
	# Example test configurations:
	# 1. Quick test (current settings) - ~1-2 minutes
	# 2. Single drug test:
	#    specific_perturbagens = c("Belinostat"), max_cell_types = 1
	# 3. Time course test:
	#    specific_times = c(8, 24, 72), max_perturbations = 2
	# 4. Full analysis:
	#    subsampling = FALSE (or set all limits to NULL)

	# helper to sanitize names for model.matrix/contrasts
	clean <- function(x) gsub("[^[:alnum:]_]", "_", x)

	# Define which DE results to store in layers
	res_cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Value", "B")

	# load the full pseudobulk AnnData
	adata <- anndata::read_h5ad(par$input)
	if (!is.null(par$subsampling) && par$subsampling) {
  		message("\n⚡ TEST MODE ENABLED ⚡")
  		message("Original dataset dimensions: ", nrow(adata$obs), " samples x ", ncol(adata$obs), " genes")
		adata <- subsampling(adata, par)
	}

	cell_types_to_process <- unique(adata$obs$cell_type)
	n_cell_types <- length(cell_types_to_process)
	cl_num <- 0

	for (cl in cell_types_to_process) {

		cl_num <- cl_num + 1
		cl_start <- Sys.time()
		message("\n▶︎ Processing cell_type ", cl_num, "/", n_cell_types, ": ", cl)
		
		# subset to this cell_type
		ad  <- adata[adata$obs$cell_type == cl, ]
		obs <- ad$obs %>%
		mutate(
			# collapsed + cleaned conditions
			raw_cond = if_else(
			is_control == TRUE,
			paste0(perturbagen, "_", as.character(pert_time_h), "h"),
			paste0(perturbagen, "_", as.character(pert_dose_uM), "uM_", as.character(pert_time_h), "h")
			),
			cond = factor(clean(raw_cond)),  # sanitize values here!
		)

		# Check for required controls
		control_times <- obs %>%
		filter(is_control==TRUE) %>%
		pull(pert_time_h) %>%
		unique()
		
		treated_times <- obs %>%
		filter(is_control!=TRUE) %>%
		pull(pert_time_h) %>%
		unique()
		
		missing_controls <- setdiff(treated_times, control_times)
		if (length(missing_controls) > 0) {
			warning("Missing controls for time points: ", paste(missing_controls, collapse=", "))
		}

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
		
		# pick out only the treated conds (we won't DE on control-vs-control)
		new_obs <- obs %>%
		distinct(cond, .keep_all = TRUE) %>%
		filter(is_control!=TRUE)
		
		# run one contrast per treated cond vs. the matching control@same time
		n_contrasts <- nrow(new_obs)
		contrast_num <- 0
		de_res <- pmap_dfr(new_obs, function(cond, perturbagen, pert_dose_uM, pert_time_h, raw_cond, ...) {
		contrast_num <<- contrast_num + 1
		if (contrast_num %% 5 == 1 || contrast_num == n_contrasts) {
			message(sprintf("    Running contrast %d/%d", contrast_num, n_contrasts))
		}
		
		#raw_control <- paste0(par$control, "_", pert_time_h, "h")
		
		# Check if control exists for this time point
		if (!any(obs[obs$is_control==TRUE, ]$pert_time_h == pert_time_h)) {
			warning("No control found for time ", pert_time_h, "h, skipping ", raw_cond)
			return(NULL)
		}
		
		raw_control <- paste0(unique(obs[(obs$is_control==TRUE) & (obs$pert_time_h == 24), ]$perturbagen)[1], "_", pert_time_h, "h")
			
		contrast <- paste0("cond", clean(raw_cond), " - cond", clean(raw_control))
		
		tryCatch({
			ctr <- makeContrasts(contrasts = contrast, levels = colnames(coef(fit)))
			fit2 <- contrasts.fit(fit, ctr) %>% 
			eBayes(robust = !isTRUE(par$subsampling))  # Skip robust for speed in test mode
		
			topTable(fit2, number = Inf, sort = "none", adjust.method="BH") %>%
			rownames_to_column("gene") %>%
			mutate(
				cond = clean(raw_cond),
				perturbagen = perturbagen,
				pert_dose_uM = pert_dose_uM,
				pert_time_h = pert_time_h
			) 
		}, error = function(e) {
			warning("Error in contrast for ", raw_cond, ": ", e$message)
			return(NULL)
		})
		})

		# Skip if no valid DE results
		if (is.null(de_res) || nrow(de_res) == 0) {
			warning("No valid DE results for cell line ", cl)
		next
		}
		
		# adjust p-values globally across all contrasts
		de_df <- de_res %>%
			mutate(
				adj.P.Value = p.adjust(P.Value, method="BH")
			)

		# FIX 1: Ensure dose_value and time are preserved as actual values, not factors
		obs_out <- new_obs %>%
			filter(cond %in% de_df$cond) %>%                           # only contrasts with results
				arrange(factor(cond, levels = unique(de_df$cond))) %>%     # keep same order as de_df
					mutate(
						pert_dose_uM = as.numeric(as.character(pert_dose_uM)),       # ensure numeric
						pert_time_h       = as.numeric(as.character(pert_time_h))
						) %>%
							remove_rownames() %>%
								column_to_rownames("cond")
		
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
		
		# carry over global uns if you like
		new_uns <- adata$uns
		
		# assemble and write
		out_adata <- anndata::AnnData(
			obs    = obs_out,
			var    = var_out,
			layers = layers,
			uns    = new_uns
		)
		
		# Create output directory if it doesn't exist
		if (!dir.exists(par$output_dir)) {
			dir.create(par$output_dir, recursive = TRUE)
		}
		
		outfile <- file.path(par$output_dir, paste0(cl, "_de.h5ad"))
		message("  Writing: ", outfile)
		out_adata$write_h5ad(outfile, compression = "gzip")
		
		# Show time for this cell line
		cl_time <- difftime(Sys.time(), cl_start, units = "secs")
		message(sprintf("  ✓ Cell line completed in %.1f seconds", cl_time))
		}
}

main <- function() {
	parser <- ArgumentParser()
	# by default ArgumentParser will add an help option 
	parser$add_argument("--input", type="string")
	parser$add_argument("--output_dir", type="string")
	parser$add_argument("--subsampling", )
	parser$add_argument("--max_cell_types", type="integer", default=3)
	parser$add_argument("--max_perturbations", type="integer", default=3)
	parser$add_argument("--max_genes", type="integer", default=1000)
	parser$add_argument("--specific_times", nargs = "+", type="integer", default=c(24))
	parser$add_argument("--specific_perturbagens", type="string", default=NULL)

	args <- parser$parse_args()


	# Start timer
	start_time <- Sys.time()
	start_deg(args)

	message("DE analysis complete!")

	# Show runtime
	end_time <- Sys.time()
	runtime <- difftime(end_time, start_time, units = "secs")
	message(sprintf("\nTotal runtime: %.1f seconds", runtime))

	if (!is.null(par$subsampling) && par$subsampling) {
  		message("\n📝 Note: This was a TEST RUN with reduced data.")
  		message("Set subsampling = FALSE for full analysis.")
	}

}

main()
