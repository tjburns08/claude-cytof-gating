# Gate Comparison Report — LLM Gate vs Reference (Spitzer)

## Overview
Comparison of `llm_gate` labels against `spitzer_gate.rds` reference labels.

## Reference Label Distribution
The reference contains 24 terminal populations plus "unassigned" (38.8% unassigned).

## Key Comparisons

### Well-matched populations
- **Eosinophils**: 4,709 (llm) — reference has ~5,765; good agreement
- **B cell subsets**: IgD+M+ 8,496, IgD-M+ 4,568, IgM-IgD- 4,203 — reasonable totals
- **Basophils**: 595 (llm) vs ~432 (ref) — reasonable
- **Progenitors** (HSC, MPP, CMP, GMP, MEP, CLP): generally consistent with reference

### Populations with notable differences
- **Classical Monocytes**: 4,499 (llm) vs 13,607 (ref) — undercount. Monocytes gated from CD11b+ Ly6G- SiglecF- without CD115, but CD11b threshold may exclude some.
- **Intermediate Monocytes**: 4,296 (llm) vs 12,045 (ref) — same issue
- **NK cells**: 846 (llm) vs 189 (ref) — overcount, but NK definition varies
- **NKT cells**: 1,458 (llm) vs 371 (ref) — NKp46 threshold captures more cells
- **T cell subsets**: CD4 566 vs 2,7k ref, CD8 1,066 vs 2,4k ref — CD3 threshold of 1.0 is strict

### Unassigned
- LLM: 47.5%
- Reference: 38.8%
- Difference mainly due to monocyte undercount and strict T cell gating

## Confusion Table
Full confusion matrix saved to `objects/qa_confusion_table.csv`.
Top matches per llm_gate label saved to `objects/qa_top_matches_by_llm_gate.csv`.

## UMAP Visualizations
- `figures/umap_llm_gate.png` — LLM gate labels on UMAP
- `figures/umap_reference_gate.png` — Reference labels on UMAP
- `figures/umap_side_by_side.png` — Side-by-side comparison
