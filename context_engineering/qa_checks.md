# qa_checks.md — Manual Gating Quality Assurance (CyTOF)

This file defines **required QA checks** for the manual gating run in `output/claude/`.  
You must execute these checks and report results in:

- `output/claude/gate_summary_report.md`
- and (if reference labels exist) `output/claude/gate_comparison_report.md`

Where possible, also save machine-readable outputs into:
- `output/claude/objects/qa_results.rds`
- `output/claude/objects/qa_tables/`

---

## 0) Definitions

- **Transformed space**: asinh(x / 5) for all marker channels used in gating/plots.
- **Root population**: `CD45.2+ Ter119−` after QC (time/beads/DNA/singlets/live).
- **Terminal labels**: the final subset names defined in `gating_strategy_for_agents.md`.
- **Unassigned**: events not placed in any terminal subset.

---

## 1) File / artifact completeness checks (hard requirements)

### 1.1 Required directories exist
Confirm these directories exist (create if missing):
- `output/claude/figures/`
- `output/claude/scripts/`
- `output/claude/objects/`

### 1.2 Required core outputs exist
Must exist at end of run:
- `output/claude/llm_gate.rds`
- `output/claude/llm_gate_metadata.csv`
- `output/claude/gate_summary_report.md`
- `output/claude/gate_summary_report.html`
- `output/claude/figures/umap_llm_gate.png`

If reference labels exist (e.g., `data/spitzer_gate.rds`), must also exist:
- `output/claude/gate_comparison_report.md`
- `output/claude/gate_comparison_report.html`
- `output/claude/figures/umap_reference_gate.png`
- `output/claude/figures/umap_side_by_side.png`

### 1.3 Gate plot completeness
For **every gate node** in `gating_strategy_for_agents.md`, there must be a corresponding plot:
- `output/claude/figures/gate_<NN>_<gate_name>.png`

Also must exist:
- `output/claude/figures/gating_tree_overview.png`

**QA action:** produce a manifest table:
- `output/claude/objects/qa_gate_plot_manifest.csv`
with columns:
- `gate_name`, `expected_file`, `exists`, `notes`

---

## 2) Data integrity checks (hard requirements)

### 2.1 Event ordering consistency
- `llm_gate` length must equal number of events in the loaded dataset.
- `llm_gate` must be in the same event order as the dataset object used downstream.

**QA action:** write these scalars to QA results:
- `n_events_dataset`, `n_events_llm_gate`, `identical_lengths (TRUE/FALSE)`

### 2.2 Allowed label set
- Every entry of `llm_gate` must be either:
  - one of the terminal subset labels, or
  - exactly `unassigned`

**QA action:** save:
- `output/claude/objects/qa_illegal_labels.csv`
listing any illegal labels and their counts (should be empty).

### 2.3 No missing labels
- No `NA` in `llm_gate`.

---

## 3) Coverage + frequency sanity checks

### 3.1 Basic label counts
Produce a table of counts and proportions:
- `table(llm_gate)`
- `% = count / n_events`

Save:
- `output/claude/objects/qa_llm_gate_counts.csv`

### 3.2 Unassigned rate
Compute:
- `unassigned_fraction = mean(llm_gate == "unassigned")`

**Guideline (soft):**
- If `unassigned_fraction > 0.40`, you must explain why and identify which branch is absorbing the cells.
- If `unassigned_fraction > 0.60`, treat as a likely gating failure and revisit thresholds.

### 3.3 Rare population sanity
These should be **rare** in most BM-like contexts:

- `HSC` should be extremely small.
- `CLP_PROXY` should be small.
- `MPP` small-to-moderate (depends).
- DCs generally smaller than monocytes.

**Guideline (soft bounds, not strict truth):**
- If `HSC > 0.1%` of root, flag.
- If `HSC == 0`, flag (may be real but usually indicates too strict gates).
- If `pDCs + mDCs > 10%` of root, flag.

**QA action:** add a short section “Frequency flags” listing any triggers.

---

## 4) Overlap / precedence audit (critical for correctness)

Before final precedence resolution, compute overlaps among terminal gates:
- For each event, count how many terminal gates it falls into.

**QA action:** report:
- number of events in 0, 1, 2, 3+ terminal gates pre-resolution
- list the most common overlapping label pairs (top 10)

Save:
- `output/claude/objects/qa_terminal_overlap_summary.csv`

**Expectation:**
- Vast majority should be in **0 or 1** terminal gate.
- If many are in 2+, gates are not mutually exclusive; refine.

Also document the **precedence order** used to resolve overlaps:
- `output/claude/objects/qa_precedence_order.txt`

---

## 5) Marker-signature validation (the most important QA)

For each terminal population, compute:
- median (and optionally IQR) of key markers in transformed space.

Save:
- `output/claude/objects/qa_marker_medians_by_label.csv`

### 5.1 Expected signatures (must generally hold)

#### Eosinophils
- **SiglecF**: high
- **Ly6G**: low
- Often **CD11b**: positive

**Flag if:** SiglecF median not clearly elevated vs root.

#### Basophils
- **FceR1a**: high
- **CD49b**: high
- **cKit**: low (because you excluded cKit+)

**Flag if:** cKit median is high, or FceR1a not elevated.

#### Macrophages
- **F480**: high
- **CD64**: high
- Often **CD11b**: positive
- **Ly6C**: low (typically)

**Flag if:** CD64 or F480 not elevated.

#### Classical / Intermediate / Non-classical Monocytes
All monocyte subsets should generally be:
- **CD11b** positive
- **CD115** positive
- **Ly6G** low
- **SiglecF** low

And differ mainly by **Ly6C**:
- Classical: Ly6C highest
- Intermediate: Ly6C intermediate
- Non-classical: Ly6C lowest

**Flag if:** Ly6C medians are not ordered: Classical > Intermediate > Non-classical.

#### pDCs
- **120g8**: high
- Usually **MHCII**: positive
- Often **B220**: positive-ish
- **CD11c**: intermediate

**Flag if:** 120g8 median not elevated.

#### mDCs
- **CD11c**: high
- **MHCII**: high
- **120g8**: low

**Flag if:** 120g8 median is high or CD11c not elevated.

#### NK cells
- **NKp46**: high
- **CD3**: low
- Often **CD49b**: positive

**Flag if:** CD3 median is high.

#### NKT cells
- **NKp46**: positive
- **CD3**: positive
- Often **TCRb**: positive

**Flag if:** either NKp46 or CD3 not elevated.

#### CD8 T cells
- **CD3**: high
- **TCRb**: high
- **CD8**: high
- **CD4**: low

**Flag if:** CD4 median is high or CD8 not elevated.

#### CD4 T cells
- **CD3**: high
- **TCRb**: high
- **CD4**: high
- **CD8**: low

**Flag if:** CD8 median is high or CD4 not elevated.

#### γδ T cells
- **CD3**: positive
- **TCRgd**: high
- **TCRb**: low-ish

**Flag if:** TCRgd median not elevated.

#### B cells (IgD/IgM subsets)
All B subsets should generally be:
- **CD19** and/or **B220** positive
- **CD3** low
- **NKp46** low

Then:
- IgDpos IgMpos: IgD high, IgM high
- IgD− IgMpos: IgD low, IgM high
- IgM− IgD−: both low

**Flag if:** IgD/IgM medians don’t match the subset name.

#### Plasma cells
- **CD138**: high
- Often **B220** low, **CD19** low
- **IgD** low

**Flag if:** CD138 not elevated.

#### Pro-B (Frac A–C proxy)
- **B220** positive
- **CD43** positive
- **IgM** low
- **IgD** low

**Flag if:** IgM median is high.

#### CMP/GMP/MEP
From LK:
- CMP: CD34 high, CD16_32 low
- GMP: CD34 high, CD16_32 high
- MEP: CD34 low, CD16_32 low

**Flag if:** these medians do not match their quadrant identity.

#### MPP / HSC (LSK-based)
- Both should be Lin−, Sca1+, cKit+
- HSC: CD150 high, CD34 low
- MPP: CD34 high (often CD150 low)

**Flag if:** HSC has high CD34 or low CD150 median.

---

## 6) “Impossible co-expression” spot checks (quick fails)

Run these checks on final labels:

- NK cells should not be **CD3 high**
- CD4 T should not be **CD8 high**
- CD8 T should not be **CD4 high**
- Plasma should not be **IgD high**
- Eosinophils should not be **Ly6G high**
- pDCs should not be **120g8 low**
- Macrophages should not be **CD64 low** AND **F480 low**

**QA action:** produce a table of “violations”:
- `output/claude/objects/qa_signature_violations.csv`
with rows:
- `label`, `rule`, `n_violating`, `fraction_violating`

**Guideline:** if any rule has >10% violations within a label, flag for gate revision.

---

## 7) Visual QA (required figures)

Create and save these additional plots for QA:

### 7.1 Marker heatmap (median by label)
- Heatmap of median transformed expression per label for key markers:
  - CD3, TCRb, TCRgd, CD4, CD8
  - CD19, B220, IgD, IgM, CD138, CD43
  - CD11b, CD115, Ly6C, Ly6G, SiglecF, F480, CD64
  - CD11c, MHCII, 120g8
  - Sca1, cKit, CD150, CD34, CD16_32

Save:
- `output/claude/figures/qa_marker_heatmap.png`

### 7.2 Density overlays for key binary decisions
For each of these, overlay density/histograms by label to confirm separation:
- `Cisplatin` (live vs dead decision)
- `CD45.2` and `Ter119` (leukocyte vs RBC decision)
- `Ly6G` (neutrophil exclusion)
- `SiglecF` (eosinophils)
- `120g8` (pDC)
- `CD138` (plasma)
- `Sca1` and `cKit` (LSK/LK)
- `CD34` vs `CD16_32` (CMP/GMP/MEP)

Save as:
- `output/claude/figures/qa_density_<marker_or_pair>.png`

---

## 8) Reference comparison QA (only if reference exists)

If `spitzer_gate` (or other reference) exists:

1. Tabulate:
   - counts per reference label
   - counts per `llm_gate`
2. Confusion matrix:
   - `table(llm_gate, spitzer_gate)`
3. For each `llm_gate` label, compute:
   - top 3 matching reference labels by proportion

Save:
- `output/claude/objects/qa_confusion_table.csv`
- `output/claude/objects/qa_top_matches_by_llm_gate.csv`

**Guideline:** you are allowed label-name differences; focus on lineage-consistent matches.

---

## 9) What to do if QA fails (required remediation loop)

If any of these conditions occur:
- unassigned_fraction > 0.60
- many overlaps (2+ terminal gates) > 5% of events
- major signature violations (>10% within a label)
- impossible frequency flags (e.g., HSC > 0.1% of root)

Then you must:
1. Identify the gate(s) responsible.
2. Adjust polygon thresholds (still polygon/rect only).
3. Re-run labeling and re-run QA.
4. Document the change in `gate_summary_report.md` under “Revisions”.

Keep all versions of altered gate plots (append `_v2`, `_v3`).

---


