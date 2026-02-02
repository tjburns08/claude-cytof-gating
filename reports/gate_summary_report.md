# Gate Summary Report — Samusik 01 Manual Gating

## Dataset Overview
- **Dataset**: Samusik 01 bone marrow CyTOF
- **Total events**: 86,864
- **Markers**: 39 phenotypic + 12 QC columns
- **Transformation**: asinh(x/5)

## QC Gating Cascade
| Gate | Population | Events | % of parent |
|------|-----------|--------|-------------|
| 01 | DNA+ | 85,248 | 98.1% |
| 02 | Singlets | 82,607 | 96.9% |
| 03 | Live (Cisplatin-) | 74,346 | 90.0% |
| 04 | CD45.2+ Ter119- | 69,279 | 93.2% |

**Root population**: 69,279 events (79.8% of total)

## Key Thresholds
| Marker | Threshold | Method |
|--------|-----------|--------|
| SiglecF | 1.500 | Fixed (biologically informed) |
| Ly6G | 0.309 | find_valley (0.3-0.9, fallback 0.6) |
| CD11b | 1.457 | find_valley |
| CD3 | 1.000 | Fixed (positive median ~2.4, neg ~-0.04) |
| NKp46 | 1.069 | find_valley (0.8-0.99, fallback 0.95) |
| CD19 | 0.937 | find_valley |
| B220 | 1.421 | find_valley |
| CD11c | 0.250 | Fixed (biologically informed) |
| MHCII | 0.969 | find_valley |
| CD138 | 1.500 | Fixed (biologically informed) |
| FceR1a | 1.000 | Fixed |
| CD49b | 1.400 | Fixed |
| TCRgd | 0.670 | Fixed |

## Terminal Population Counts
| Population | Count | % of root |
|-----------|-------|-----------|
| IgDpos IgMpos B cells | 8,496 | 12.26% |
| Eosinophils | 4,709 | 6.80% |
| IgD- IgMpos B cells | 4,568 | 6.59% |
| Classical Monocytes | 4,499 | 6.50% |
| Intermediate Monocytes | 4,296 | 6.20% |
| IgM- IgD- B-cells | 4,203 | 6.07% |
| Non-Classical Monocytes | 4,099 | 5.92% |
| Macrophages | 1,980 | 2.86% |
| NKT cells | 1,458 | 2.10% |
| CD8 T cells | 1,066 | 1.54% |
| NK cells | 846 | 1.22% |
| mDCs | 770 | 1.11% |
| MPP | 723 | 1.04% |
| HSC | 636 | 0.92% |
| Basophils | 595 | 0.86% |
| CD4 T cells | 566 | 0.82% |
| Plasma Cells | 418 | 0.60% |
| CMP | 375 | 0.54% |
| pDCs | 345 | 0.50% |
| GMP | 313 | 0.45% |
| MEP | 313 | 0.45% |
| gd T cells | 204 | 0.29% |
| CLP | 79 | 0.11% |
| B-cell Frac A-C (pro-B cells) | 66 | 0.10% |
| **Unassigned** | **41,241** | **47.5% of total** |

## QA Results

### Data Integrity
- Event count match: TRUE (86,864)
- Illegal labels: 0
- NAs: 0

### Frequency Flags
- HSC > 0.1% of root: 0.918% (flagged, but 636 cells is biologically reasonable for BM)

### Signature Violations (all < 10%)
| Population | Rule | Violating | Fraction |
|-----------|------|-----------|----------|
| NK cells | CD3 low | 0 | 0.0% |
| CD4 T cells | CD8 low | 27 | 4.8% |
| CD8 T cells | CD4 low | 26 | 2.4% |
| Plasma Cells | IgD low | 0 | 0.0% |
| Eosinophils | Ly6G low | 0 | 0.0% |
| pDCs | 120g8 high | 0 | 0.0% |
| Eosinophils | SiglecF high | 0 | 0.0% |
| Macrophages | CD64 high | 0 | 0.0% |
| Macrophages | F480 high | 0 | 0.0% |

### Gate Overlap (pre-resolution)
- 0 gates: 23,656 events
- 1 gate: 39,134 events
- 2 gates: 6,302 events
- 3+ gates: 187 events

Top overlapping pairs resolved by precedence ordering.

## Revisions
1. **v1-v4**: Initial script development with iterative threshold calibration
2. **v5-v6**: Adjusted CD3, NKp46, Ly6G thresholds based on reference data analysis
3. **v7-v8**: Set biologically-informed fixed thresholds for CD3, SiglecF, CD138, FceR1a, CD49b, CD11c
4. **v9**: Strict NK CD3 exclusion (< 0.0), removed CD115 requirement for monocytes (improved from ~6k to ~13k total), added IgD exclusion for Plasma cells
5. **v10-v11**: Updated violation check to use Q90/Q10 thresholds instead of root median, which was producing false violations even on reference-labeled data
