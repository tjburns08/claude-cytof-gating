#!/usr/bin/env Rscript
# Calculate per-population precision, recall, and F1 scores
# Compares LLM gate labels against Spitzer reference labels

proj_root <- here::here()

llm_labels <- readRDS(file.path(proj_root, "llm_gate.rds"))
ref_labels <- readRDS(file.path(proj_root, "data/spitzer_gate.rds"))

populations <- sort(unique(ref_labels))

results <- do.call(rbind, lapply(populations, function(pop) {
  llm_pos <- llm_labels == pop
  ref_pos <- ref_labels == pop
  tp <- sum(llm_pos & ref_pos)
  fp <- sum(llm_pos & !ref_pos)
  fn <- sum(!llm_pos & ref_pos)
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall    <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1        <- ifelse(precision + recall == 0, 0,
                      2 * precision * recall / (precision + recall))
  data.frame(
    Population = pop,
    Precision  = round(precision * 100, 1),
    Recall     = round(recall * 100, 1),
    F1         = round(f1 * 100, 1),
    LLM_Count  = sum(llm_pos),
    Ref_Count  = sum(ref_pos),
    TP         = tp,
    stringsAsFactors = FALSE
  )
}))

results <- results[order(-results$F1), ]
rownames(results) <- NULL

cat("Per-population precision, recall, and F1 scores:\n\n")
print(results, row.names = FALSE)
