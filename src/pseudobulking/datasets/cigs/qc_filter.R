#!/usr/bin/env Rscript

# Configure reticulate to use conda Python (kept for consistency with the
# other R helpers in this repo, e.g. novartis/rdata_to_adata.R).
suppressPackageStartupMessages(library(reticulate))
conda_env <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_env)) use_python(file.path(conda_env, "bin/python"), required = TRUE)

suppressPackageStartupMessages({
  library(argparse)
})


# Paper-level plate exclusions hard-coded in
# Wang-lab302/CIGS/Code/00._DataPre.Rmd for the MCE HEK293T subset.  The
# authors drop these three plates without explanation; they likely
# represent a library-construction batch (MCE4) with catastrophic
# sequencing QC issues that the per-sample IQR rule cannot detect (see the
# paper's Extended Data Fig. 1c/f discussing HTS2 mapping-rate variance
# and inter-batch reproducibility).
PAPER_EXCLUDED_PLATES <- c(
  "MCE4_293T_24H_X1",
  "MCE4_293T_24H_X2",
  "MCE4_293T_24H_X3"
)


#' Flag BLANK/RNA control wells, paper-excluded plates, and per-subset IQR
#' low-QC outliers.
#'
#' Mirrors the filtering block of the official CIGS paper R code
#' (Wang-lab302/CIGS, Code/00._DataPre.Rmd):
#'
#'   deltaQ    <- quantile(log2(reads), 0.75) - quantile(log2(reads), 0.25)
#'   threshold <- quantile(log2(reads), 0.25) - 1.5 * deltaQ
#'   status    <- ifelse(log2(reads) < threshold, "lowqc", "pass")
#'
#' applied independently to each of the 6 (library x cell line x dose)
#' subsets.  BLANK/RNA controls are dropped first so they do not distort
#' the IQR, matching the paper's `meta[meta$Treat != "BLANK", ]` /
#' `meta[meta$treatment != "RNA" & meta$treatment != "Blank", ]` step.
#' The three MCE4 HEK293T plates that the paper excludes by name
#' (see PAPER_EXCLUDED_PLATES) are tagged as "excluded_plate" so they are
#' also dropped and so they do not influence the IQR.
#'
#' @param df data.frame with columns: sample_id, subset_key, total_reads,
#'   treatment, sample_plate
#' @return df augmented with a `status` column: "pass", "blank",
#'   "excluded_plate", or "lowqc".
flag_samples <- function(df) {
  df$status <- "pass"

  treat_up <- toupper(trimws(as.character(df$treatment)))
  df$status[treat_up %in% c("BLANK", "RNA")] <- "blank"

  plates <- trimws(as.character(df$sample_plate))
  df$status[plates %in% PAPER_EXCLUDED_PLATES] <- "excluded_plate"

  for (key in unique(df$subset_key)) {
    idx <- which(df$subset_key == key & df$status == "pass")
    if (length(idx) < 4) next

    vals      <- log2(df$total_reads[idx])
    q1        <- as.numeric(quantile(vals, 0.25))
    q3        <- as.numeric(quantile(vals, 0.75))
    deltaQ    <- q3 - q1
    threshold <- q1 - 1.5 * deltaQ
    low_idx   <- idx[vals < threshold]
    df$status[low_idx] <- "lowqc"

    message(sprintf(
      "  %s: IQR threshold log2(lib)=%.4f, lowqc=%d/%d",
      key, threshold, length(low_idx), length(idx)
    ))
  }
  df
}


#' Main entry point.
#'
#' Reads a TSV with columns (sample_id, subset_key, total_reads, treatment),
#' applies flag_samples(), and writes a TSV with columns (sample_id, status).
main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--input_tsv",  type = "character", required = TRUE,
                      help = "Input TSV: sample_id, subset_key, total_reads, treatment, sample_plate")
  parser$add_argument("--output_tsv", type = "character", required = TRUE,
                      help = "Output TSV: sample_id, status (pass/blank/lowqc)")
  args <- parser$parse_args()

  df <- read.table(
    args$input_tsv,
    header           = TRUE,
    sep              = "\t",
    stringsAsFactors = FALSE,
    check.names      = FALSE,
    quote            = "",
    comment.char     = "",
    na.strings       = c("", "NA")
  )
  if (!"sample_plate" %in% colnames(df)) df$sample_plate <- ""

  n_blank_ctrl <- sum(toupper(trimws(as.character(df$treatment))) %in% c("BLANK", "RNA"), na.rm = TRUE)
  n_excluded   <- sum(trimws(as.character(df$sample_plate)) %in% PAPER_EXCLUDED_PLATES, na.rm = TRUE)
  message(sprintf(
    "Loaded %d samples across %d subsets (%d BLANK/RNA wells, %d paper-excluded-plate wells)",
    nrow(df), length(unique(df$subset_key)), n_blank_ctrl, n_excluded
  ))

  df <- flag_samples(df)

  tab <- table(df$status)
  message("Status counts:")
  for (s in names(tab)) message(sprintf("  %s: %d", s, tab[[s]]))

  write.table(
    df[, c("sample_id", "status")],
    file      = args$output_tsv,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  message("Wrote ", args$output_tsv)
}

main()
