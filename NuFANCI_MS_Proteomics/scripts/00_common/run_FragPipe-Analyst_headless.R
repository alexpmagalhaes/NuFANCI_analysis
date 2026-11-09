#!/usr/bin/env Rscript
# =============================================================================
# FragPipe-Analyst LFQ analysis without the Shiny app (command line)
#
# Runs the same steps as the FragPipe-Analyst web app ("LFQ" data type) by
# calling the app's own R functions (R/customized.R, R/functions.R, R/filter.R),
# so results are identical to the GUI with the same settings:
#
#   read combined_protein(_f).tsv -> remove contaminants -> make SE (log2)
#   -> missing-value filter -> normalisation -> imputation
#   -> limma (all pairwise contrasts, BH) -> significance (p.adj, |log2FC|)
#   -> tables (DE results, matrices), dual comparisons, run log
#
# FragPipe-Analyst: https://github.com/MonashProteomics/FragPipe-Analyst
# Used version: commit 3e1af92
#   git clone https://github.com/MonashProteomics/FragPipe-Analyst.git
#   cd FragPipe-Analyst && git checkout 3e1af92
#
# Usage:
#   Rscript run_FragPipe-Analyst_headless.R \
#     [--fpa_dir=~/FragPipe-Analyst] [--use_renv=TRUE] \
#     --input=combined_protein_f.tsv \
#     --annotation=experiment_annotation.tsv \
#     --outdir=. --suffix=_minImputation \
#     [--lfq_type=MaxLFQ] [--normalization=vsn] [--imputation=min] \
#     [--fdr=BH] [--p=0.05] [--lfc=1] [--de_type=all] [--control=COND] \
#     [--min_global=0] [--min_condition=0] [--remove_contaminants=TRUE] \
#     [--dual=CONTRAST1,CONTRAST2 --dual=...] [--dual_alpha=0.05] \
#     [--dual_lfc=1] [--dual_p_adj=TRUE]
#
# Defaults are the app defaults used for the NuFANCI analysis
# (MaxLFQ, VSN, min imputation, BH, all pairs, p.adj <= 0.05, |log2FC| >= 1).
#
# Outputs (in --outdir; <sfx> = --suffix):
#   DE_results<sfx>.csv         same as "DE results" download
#   Original_matrix<sfx>.csv, Filtered_matrix<sfx>.csv,
#   Normalized_matrix<sfx>.csv, Imputed_matrix<sfx>.csv
#   dual_comparison_<c1>_and_<c2>.tsv / dual_comparison_scatter_<c1>_and_<c2>.pdf
#   FragPipe-Analyst_log<sfx>.txt
#
# Requires a local FragPipe-Analyst installation (default ~/FragPipe-Analyst).
# If it has an renv library (renv/activate.R), that library is used, so the packages
# are exactly those of the app. Otherwise these R/Bioconductor packages must be
# installed: SummarizedExperiment, limma, MSnbase, vsn, dplyr, tidyr, tibble, purrr,
# stringr, ggplot2, ggrepel, assertthat, matrixStats, conflicted.
# =============================================================================

# ---- command-line arguments (--key=value; --dual may be repeated) -----------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hits <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hits) == 0) return(default)
  sub(paste0("^--", key, "="), "", hits)
}
as_bool <- function(x) toupper(x) %in% c("TRUE", "T", "1", "YES")

opt <- list(
  fpa_dir        = get_arg("fpa_dir", "~/FragPipe-Analyst"),
  use_renv       = as_bool(get_arg("use_renv", "TRUE")),
  input          = get_arg("input", "combined_protein_f.tsv"),
  annotation     = get_arg("annotation", "experiment_annotation.tsv"),
  outdir         = get_arg("outdir", "."),
  suffix         = get_arg("suffix", ""),
  lfq_type       = get_arg("lfq_type", "MaxLFQ"),          # Intensity | MaxLFQ | Spectral Count
  normalization  = get_arg("normalization", "vsn"),        # none | vsn | MD
  imputation     = get_arg("imputation", "min"),           # none | man | knn | MLE | min | zero
  fdr            = get_arg("fdr", "BH"),                   # BH | fdrtool
  p              = as.numeric(get_arg("p", "0.05")),
  lfc            = as.numeric(get_arg("lfc", "1")),
  de_type        = get_arg("de_type", "all"),              # all | control | others
  control        = get_arg("control", ""),
  min_global     = as.numeric(get_arg("min_global", "0")),
  min_condition  = as.numeric(get_arg("min_condition", "0")),
  remove_contam  = as_bool(get_arg("remove_contaminants", "TRUE")),
  dual           = get_arg("dual", character(0)),
  dual_alpha     = as.numeric(get_arg("dual_alpha", "0.05")),
  dual_lfc       = as.numeric(get_arg("dual_lfc", "1")),
  dual_p_adj     = as_bool(get_arg("dual_p_adj", "TRUE"))
)
opt$fpa_dir <- normalizePath(path.expand(opt$fpa_dir), mustWork = FALSE)
if (!file.exists(file.path(opt$fpa_dir, "R", "customized.R"))) {
  stop("FragPipe-Analyst not found in '", opt$fpa_dir, "'. Set --fpa_dir to your FragPipe-Analyst folder.")
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Use the same R package library as the local FragPipe-Analyst app (renv), if present,
# so the analysis runs with exactly the package versions the app uses.
renv_activate <- file.path(opt$fpa_dir, "renv", "activate.R")
if (opt$use_renv && file.exists(renv_activate)) {
  owd <- setwd(opt$fpa_dir)
  source(renv_activate)
  setwd(owd)
  message("Using the FragPipe-Analyst renv library: ", .libPaths()[1])
}

# Sort condition names like macOS does (case-insensitive: "HepG2" before "HSLiver"),
# so contrast/column order matches the app. Only affects ordering, not values.
for (loc in c("en_US.UTF-8", "en_US.utf8")) {
  if (nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc)))) break
}
out <- function(name) file.path(opt$outdir, name)

# ---- logging (same message format as the app's run log) ---------------------
log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(...))
  log_lines <<- c(log_lines, msg)
  message(msg)
}

# ---- packages and FragPipe-Analyst functions --------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(stringr)
  library(ggplot2); library(ggrepel)
  library(limma); library(MSnbase); library(vsn)
  library(assertthat); library(matrixStats)
  library(conflicted)
})
# same conflict preferences as the app (global.R)
conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)
conflict_prefer("desc", "dplyr", quiet = TRUE)
conflict_prefer("setdiff", "base", quiet = TRUE)
conflict_prefer("unite", "tidyr", quiet = TRUE)
conflict_prefer("intersect", "base", quiet = TRUE)
conflict_prefer("colMedians", "matrixStats", quiet = TRUE)
conflicts_prefer(base::as.factor, .quiet = TRUE)

for (f in c("R/functions.R", "R/customized.R", "R/filter.R")) {
  source(file.path(opt$fpa_dir, f))
}

log_msg("Analysis started. Data type: LFQ")

# ---- 1. read FragPipe protein table (as in server.R, LFQ) -------------------
fp <- read.table(opt$input, header = TRUE, fill = TRUE, sep = "\t", quote = "",
                 comment.char = "", blank.lines.skip = FALSE, check.names = FALSE)
colnames(fp) <- make.unique.2(colnames(fp), "_")
colnames(fp) <- gsub("-", ".", colnames(fp))          # dashes in experiment names
if (opt$remove_contam) fp <- fp[!grepl("contam", fp$Protein), ]

# ---- 2. read experiment annotation (as in server.R, LFQ) --------------------
ed <- read.table(opt$annotation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(ed) <- tolower(colnames(ed))
ed$sample_name <- as.character(ed$sample_name)
ed$condition <- make.names(ed$condition)
ed$sample <- gsub("-", ".", ed$sample)
ed$label <- ed$sample
if (opt$lfq_type == "Intensity")      ed$label <- paste(ed$label, "Intensity", sep = " ")
if (opt$lfq_type == "MaxLFQ")         ed$label <- paste(ed$label, "MaxLFQ.Intensity", sep = " ")
if (opt$lfq_type == "Spectral Count") ed$label <- paste(ed$label, "Spectral.Count", sep = " ")

# ---- 3. SummarizedExperiment (log2 intensities) ------------------------------
data_unique <- make_unique(fp, "Gene", "Protein ID")
if (opt$lfq_type == "Intensity") {
  lfq_columns <- setdiff(grep("Intensity", colnames(data_unique)), grep("MaxLFQ", colnames(data_unique)))
  lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data_unique)))
  lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data_unique)))
} else if (opt$lfq_type == "MaxLFQ") {
  lfq_columns <- grep("MaxLFQ", colnames(data_unique))
} else {
  lfq_columns <- grep("Spectral", colnames(data_unique))
  lfq_columns <- setdiff(lfq_columns, grep("Total Spectral Count", colnames(data_unique)))
  lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(data_unique)))
}
spectral <- opt$lfq_type == "Spectral Count"
se <- make_se_customized(data_unique, lfq_columns, ed, log2transform = !spectral,
                         exp = "LFQ", lfq_type = opt$lfq_type, level = "protein")
log_msg(sprintf("Data loaded: %d features x %d samples", nrow(se), ncol(se)))
log_msg(sprintf("Conditions: %s", paste(unique(colData(se)$condition), collapse = ", ")))

# ---- 4. missing-value filtering ---------------------------------------------
filtered <- se
if (opt$min_global != 0)    filtered <- global_filter(se, 100 - opt$min_global)
if (opt$min_condition != 0) filtered <- filter_by_condition(filtered, opt$min_condition)
log_msg(sprintf("Filtering: %d features retained, %d removed (global >= %g%%, per-condition >= %g%%)",
                nrow(filtered), nrow(se) - nrow(filtered), opt$min_global, opt$min_condition))

# ---- 5. normalisation --------------------------------------------------------
normalized <- filtered
if (opt$normalization == "vsn" && !spectral) {
  normalized <- normalize_vsn(filtered)
} else if (opt$normalization == "MD") {
  assay(normalized) <- sweep(assay(filtered), 2, colMedians(assay(filtered), na.rm = TRUE))
}
log_msg(sprintf("Normalization: %s",
                switch(opt$normalization, none = "none", vsn = "VSN", MD = "Median-centered", opt$normalization)))

# ---- 6. imputation -----------------------------------------------------------
if (opt$imputation == "none") {
  imputed <- normalized
  rowData(imputed)$imputed <- FALSE
  rowData(imputed)$num_NAs <- rowSums(is.na(assay(normalized)))
} else {
  imputed <- impute_customized(normalized, opt$imputation)
}
log_msg(sprintf("Imputation: %s", opt$imputation))

# ---- 7. differential abundance -----------------------------------------------
data <- imputed
if (spectral) assay(data) <- log2(assay(data))
control <- if (grepl(";", opt$control)) unlist(str_split(opt$control, ";")) else opt$control
if (opt$fdr == "BH") {
  diff_all <- switch(opt$de_type,
    all     = test_limma_customized(data, type = "all", paired = FALSE),
    control = test_limma_customized(data, type = "control", control = control, paired = FALSE),
    others  = test_limma_customized(data, type = "others", paired = FALSE))
} else {
  diff_all <- switch(opt$de_type,
    all     = test_diff_customized(data, type = "all"),
    control = test_diff_customized(data, type = "control", control = control),
    others  = test_diff_customized(data, type = "others"))
}
dep <- add_rejections_customized(diff_all, alpha = opt$p, lfc = opt$lfc)

sig_cols <- grep("_significant$", colnames(rowData(dep)), value = TRUE)
n_sig <- sum(rowSums(as.data.frame(rowData(dep)[, sig_cols, drop = FALSE]), na.rm = TRUE) > 0)
log_msg(sprintf("DE analysis complete (FDR method: %s, p <= %g, |LFC| >= %g): %d significant feature(s) across all contrasts",
                opt$fdr, opt$p, opt$lfc, n_sig))
log_msg(sprintf("Contrasts tested: %s", paste(sub("_significant$", "", sig_cols), collapse = ", ")))

# ---- 8. tables (same format as the app's downloads) --------------------------
matrix_table <- function(x, tag) {
  m <- assay(x)
  colnames(m) <- paste(colnames(m), paste0(tag, if (spectral) "_spectral_count" else "_intensity"), sep = "_")
  as.data.frame(cbind(ProteinID = rownames(m), m))
}
write_csv_app <- function(df, name) write.table(df, out(name), col.names = TRUE, row.names = FALSE, sep = ",")
sfx <- opt$suffix
write_csv_app(get_results_proteins_customized(dep), paste0("DE_results", sfx, ".csv"))
write_csv_app(matrix_table(se, "original"),          paste0("Original_matrix", sfx, ".csv"))
write_csv_app(matrix_table(filtered, "filtered"),    paste0("Filtered_matrix", sfx, ".csv"))
write_csv_app(matrix_table(normalized, "normalized"), paste0("Normalized_matrix", sfx, ".csv"))
write_csv_app(matrix_table(imputed, "imputed"),      paste0("Imputed_matrix", sfx, ".csv"))

# ---- 9. dual comparisons (as the app's "Dual Comparison" tab) ----------------
dual_comparison <- function(c1, c2) {
  rd <- as.data.frame(rowData(dep))
  p1 <- paste0(c1, if (opt$dual_p_adj) "_p.adj" else "_p.val")
  p2 <- paste0(c2, if (opt$dual_p_adj) "_p.adj" else "_p.val")
  l1 <- paste0(c1, "_diff"); l2 <- paste0(c2, "_diff")
  stopifnot(all(c(l1, l2, p1, p2) %in% colnames(rd)))
  sig1 <- !is.na(rd[[p1]]) & abs(rd[[l1]]) >= opt$dual_lfc & rd[[p1]] <= opt$dual_alpha
  sig2 <- !is.na(rd[[p2]]) & abs(rd[[l2]]) >= opt$dual_lfc & rd[[p2]] <= opt$dual_alpha
  status <- dplyr::case_when(sig1 & sig2 ~ "Both significant",
                             sig1 & !sig2 ~ paste0("Only ", c1),
                             !sig1 & sig2 ~ paste0("Only ", c2),
                             TRUE ~ "Not significant")
  label_col <- if ("Gene" %in% colnames(rd)) "Gene" else "name"
  df <- data.frame(lfc1 = rd[[l1]], lfc2 = rd[[l2]], status = status,
                   name = rd[[label_col]], stringsAsFactors = FALSE)
  colnames(df)[1:2] <- c(paste0("log2FC_", c1), paste0("log2FC_", c2))
  for (cc in c(c1, c2)) df[[paste0("p.val_", cc)]] <- rd[[paste0(cc, "_p.val")]]
  for (cc in c(c1, c2)) df[[paste0("p.adj_", cc)]] <- rd[[paste0(cc, "_p.adj")]]
  df
}
dual_plot <- function(df, c1, c2) {
  x <- paste0("log2FC_", c1); y <- paste0("log2FC_", c2)
  lv <- c("Both significant", paste0("Only ", c1), paste0("Only ", c2), "Not significant")
  lv <- lv[lv %in% unique(df$status)]
  cols <- c("Both significant" = "#C0392B", "Not significant" = "grey70")
  cols[paste0("Only ", c1)] <- "#E67E22"; cols[paste0("Only ", c2)] <- "#2980B9"
  df$status <- factor(df$status, levels = lv)
  d <- df[!is.na(df[[x]]) & !is.na(df[[y]]), ]
  lim <- max(abs(c(d[[x]], d[[y]])), na.rm = TRUE) * 1.05
  lab <- d[d$status == "Both significant", ]
  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]], color = status)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = cols, name = "Significance") +
    geom_hline(yintercept = c(-opt$dual_lfc, opt$dual_lfc), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-opt$dual_lfc, opt$dual_lfc), linetype = "dashed", color = "grey40") +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = paste0("log2FC (", c1, ")"), y = paste0("log2FC (", c2, ")"),
         title = "Marker Nomination: Dual Comparison") +
    theme_bw() + theme(legend.position = "right")
  if (nrow(lab) > 0 && nrow(lab) <= 40) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(x = .data[[x]], y = .data[[y]], label = name),
                                      color = "#C0392B", size = 3,
                                      box.padding = unit(0.2, "lines"), max.overlaps = 20)
  }
  p
}
for (pair in opt$dual) {
  cc <- strsplit(pair, ",")[[1]]
  base <- paste0(cc[1], "_and_", cc[2])
  df <- dual_comparison(cc[1], cc[2])
  write.table(df, out(paste0("dual_comparison_", base, ".tsv")),
              col.names = TRUE, row.names = FALSE, sep = "\t")
  pdf(out(paste0("dual_comparison_scatter_", base, ".pdf")))
  print(dual_plot(df, cc[1], cc[2]))
  dev.off()
}

# ---- 10. log -----------------------------------------------------------------
writeLines(log_lines, out(paste0("FragPipe-Analyst_log", sfx, ".txt")))
message("Done. Results in: ", normalizePath(opt$outdir))
