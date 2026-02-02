# Can an LLM Replace a Human at Manual Gating?

An experiment in fully autonomous CyTOF manual gating using Claude Code (Opus 4.5). No human touched the data. The LLM wrote and iteratively debugged an R script implementing a 24-population hierarchical gating pipeline — from QC gates through terminal cell type assignment.

## The Setup

- **Dataset**: Samusik 01 bone marrow CyTOF (86,864 events, 39 phenotypic markers)
- **Reference**: Expert-curated labels from Spitzer et al.
- **Agent**: Claude Code (Opus 4.5), fully autonomous (zero human intervention)
- **Pipeline**: R script with asinh(x/5) transformation, density-based threshold detection, hierarchical gating
- **Iteration**: 11 autonomous runs with self-debugging and QA validation

## Results

| Metric | Value |
|--------|-------|
| Overall exact-match accuracy | 48.4% |
| Accuracy on assigned cells | 60.5% |
| Lineage-level accuracy | 82.4% |
| Weighted F1 score | 40.2% |
| Populations gated | 24 terminal types |
| Human intervention | Zero |

Top performers: IgD+IgM+ B cells (86.5% F1), CD8 T cells (74.6% F1), IgD-IgM+ B cells (72.5% F1).

Bottom performers: CLP (0% F1), HSC (0.3% F1), Macrophages (0.5% F1) — rare or poorly resolved populations.

## What Went Right

- **82.4% lineage-level accuracy**: The LLM rarely confuses major lineages (B vs T vs myeloid vs progenitor)
- **Self-correcting**: Autonomously diagnosed 6 failure modes across 11 iterations, including identifying that CD115 was uninformative in this panel and that its own QA metric was mathematically flawed
- **B cell gating**: Near-expert performance on IgD/IgM B cell subsets
- **Complete artifact suite**: 31 gate plots, UMAP visualizations, confusion matrices, precision/recall tables, and full reasoning trace

## What Went Wrong

- **Monocyte undercount**: 13k vs 27k reference (CD11b threshold too strict)
- **Rare progenitor confusion**: CLP, HSC, and other rare populations (<1% of cells) were poorly resolved by density-based thresholding
- **NK cell overcount**: 846 vs 189 reference (NKp46 is dim in this panel)
- **47.5% unassigned rate** vs 38.8% reference

## Repository Structure

```
├── data/                          # Raw CyTOF data (.rds)
├── scripts/
│   ├── run_manual_gating.R        # Main gating pipeline (881 lines, written by Claude)
│   └── calculate_precision_recall.R
├── figures/                       # 48 PNG files (gate plots, UMAP, QA)
├── objects/                       # Intermediate data (confusion matrix, QA tables)
├── prompt.md                      # Agent instructions given to Claude
├── qa_checks.md                   # QA requirements specification
├── reasoning_trace_report.md      # Full iteration log (11 runs)
├── gate_summary_report.md/.html   # Main results report
├── gate_comparison_report.md/.html # LLM vs reference comparison
├── precision_recall.md/.html/.pdf # Per-population metrics
├── executive_summary.md/.html/.pdf # High-level summary
├── llm_gate.rds                   # Final cell labels (86,864 events)
└── llm_gate_metadata.csv          # Event-by-event label table
```

## Key Takeaway

Claude can autonomously implement a complex bioinformatics pipeline with full QA and iteration, achieving strong lineage-level accuracy (82.4%). It is not yet a replacement for expert gating at the subtype level (48.4% overall), particularly for rare populations and markers with poor bimodal separation. The autonomous debugging capability — independently diagnosing failure modes and adapting — is the most promising signal for future work.

## License

GPL-3.0
