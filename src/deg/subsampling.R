subsampling <- function(adata, par) {
    stopifnot(!missing(adata),
              !missing(par))
    # Test mode subsetting
    if (!is.null(par$subsampling) && par$subsampling) {
      cat("\n⚡ TEST MODE ENABLED ⚡\n")
      cat("Original dataset dimensions: ", nrow(adata$obs), " samples x ", ncol(adata$obs), " genes\n")
      
      # Subset cell lines
      cell_types <- unique(adata$obs$cell_type)
      if (!is.null(par$max_cell_types) && length(cell_types) > par$max_cell_types) {
        cell_types <- cell_types[1:par$max_cell_types]
        cat("Subsetting to ", par$max_cell_types, " cell lines: ", paste(cell_types, collapse=", "), "\n")
        adata <- adata[adata$obs$cell_type %in% cell_types, ]
      }
      
      # Subset time points
      if (!is.null(par$specific_times)) {
        cat("Filtering to time points: ", paste(par$specific_times, collapse=", "), " hours\n")
        adata <- adata[adata$obs$pert_time_h %in% par$specific_times, ]
      }
      
      # Subset perturbations 
      if (!is.null(par$specific_perturbagens)) {
        # Use specific perturbagens if provided
        control_perts <- as.character(unique(adata[(adata$obs$is_control == TRUE)]$obs$perturbagen))
        keep_perts <- c(control_perts, par$specific_perturbagens)
        cat("Using specific perturbagens: ", paste(keep_perts, collapse=", "), "\n")
        adata <- adata[adata$obs$perturbation %in% keep_perts, ]
      } else if (!is.null(par$max_perturbations)) {
        # Otherwise limit to max_perturbations
        perturbagens <- as.character(unique(adata$obs$perturbagen))
        control_perts <- as.character(unique(adata[(adata$obs$is_control == TRUE)]$obs$perturbagen))
        treated_perts <- setdiff(perturbagens, control_perts)
        if (length(treated_perts) > par$max_perturbations) {
          keep_perts <- c(control_perts, treated_perts[1:par$max_perturbations])
          cat("Subsetting to ", par$max_perturbations, " perturbations (+ control): ", 
                  paste(keep_perts, collapse=", "), "\n")
          adata <- adata[adata$obs$perturbagen %in% keep_perts, ]
        }
      }
      
      # Subset genes - keep most variable
      if (!is.null(par$max_genes) && ncol(adata) > par$max_genes) {
        cat("Calculating gene variance for subsetting...\n")
        gene_vars <- apply(as.matrix(adata$X), 2, var)
        top_genes <- order(gene_vars, decreasing = TRUE)[1:par$max_genes]
        adata <- adata[, top_genes]
        cat("Subset to top ", par$max_genes, " most variable genes\n")
      }
      
      cat("Test mode dataset: ", nrow(adata$obs), " samples x ", ncol(adata), " genes\n")
      
      # Show analysis overview
      analysis_summary <- adata$obs %>%
        group_by(cell_type, perturbagen, pert_time_h) %>%
        summarise(n_samples = n(), .groups = "drop") %>%
        arrange(cell_type, pert_time_h, perturbagen)
      
      cat("\nConditions to analyze:\n")
      print(analysis_summary)
      message("")
    }
    return(adata)
}
