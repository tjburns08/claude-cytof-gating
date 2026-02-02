#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ExperimentHub)
  library(HDCytoData)
  library(SummarizedExperiment)
})

cache_dir <- file.path(getwd(), "ehub_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
setExperimentHubOption("CACHE", cache_dir)
setExperimentHubOption("ASK", FALSE)

se <- Samusik_01_SE()

assay_name <- if ("exprs" %in% assayNames(se)) "exprs" else assayNames(se)[1]
if (is.null(assay_name) || is.na(assay_name) || length(assay_name) == 0) {
  stop("No assays found in Samusik_01_SE.")
}

exprs <- SummarizedExperiment::assay(se, assay_name)
rd <- SummarizedExperiment::rowData(se)
rd_df <- as.data.frame(rd)
if (nrow(exprs) == nrow(rd_df)) {
  exprs_df <- as.data.frame(exprs, check.names = FALSE)
} else if (ncol(exprs) == nrow(rd_df)) {
  exprs_df <- as.data.frame(t(exprs), check.names = FALSE)
} else {
  stop("Could not align expression data with rowData.")
}

preferred_cols <- c(
  "population_id", "population_name", "population_label",
  "population", "Population", "cell_type", "celltype", "CellType", "Cell_Type",
  "cell_label", "cell_labels", "label", "Label", "gate", "Gate", "gating", "Gating",
  "annotation", "Annotation", "manual_gate", "ManualGate", "manual_gating", "ManualGating"
)

select_pop_col <- function(df) {
  candidates <- intersect(preferred_cols, colnames(df))
  if (length(candidates) == 0) {
    target_patterns <- c("CD4", "CD8", "B", "NK", "Mono", "Monocyte", "DC", "Dendritic", "Treg")
    candidates <- names(df)[vapply(df, function(x) {
      if (!is.character(x) && !is.factor(x)) return(FALSE)
      vals <- unique(as.character(x))
      any(grepl(paste(target_patterns, collapse = "|"), vals, ignore.case = TRUE))
    }, logical(1))]
  }

  if (length(candidates) == 0) return(NA_character_)
  if (length(candidates) == 1) return(candidates)

  score_col <- function(col) {
    vals <- unique(as.character(df[[col]]))
    sum(grepl("CD4|CD8|B|NK|Mono|Monocyte|DC|Dendritic|Treg", vals, ignore.case = TRUE))
  }

  scores <- vapply(candidates, score_col, numeric(1))
  candidates[which.max(scores)]
}

pop_col <- select_pop_col(rd_df)
if (is.na(pop_col)) {
  stop(
    "No population metadata column found in rowData. Available columns: ",
    paste(colnames(rd_df), collapse = ", ")
  )
}

spitzer_gate <- as.character(rd_df[[pop_col]])
if (!any(grepl("CD4", spitzer_gate, ignore.case = TRUE))) {
  stop(
    "Population labels do not appear to include CD4 entries. ",
    "Selected column: ", pop_col, ". Unique labels (first 20): ",
    paste(head(unique(spitzer_gate), 20), collapse = ", ")
  )
}

exprs_df$spitzer_gate <- spitzer_gate

saveRDS(exprs_df, file = "samusik_01.rds")

message("Saved samusik_01.rds using assay '", assay_name, "' and population column '", pop_col, "'.")
