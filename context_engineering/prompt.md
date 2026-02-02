# prompt.md — Agent Instructions for Manual Gating Experiment (CyTOF)

You are operating inside a project folder that contains an `output/claude/` directory.  
**All results you create must be written to `output/claude/`** (scripts, figures, tables, reports, serialized objects).

This project may contain outputs from previous runs. **Do not bias yourself by anything already in `output/claude/`.**  
Treat this run as a fresh experiment.

---

## 1) Inputs and context engineering files (read these first)

### 1.1 Gating strategy (source of truth)
In `output/claude/` you will find:

- `output/claude/xshift_gating_strategy_from_pdf.png`

This file contains an image of the exact hierarchical gating plan (e.g. QC → CD45+ root → immune and progenitors → subsets). Use it to craft an explicit gating strategy and intended cell subsets that you follow from here on out.  
**You must follow it.** Do not invent an alternative gating tree unless a marker is missing from the dataset (in which case, document the deviation explicitly).

### 1.2 QA requirements (also source of truth)
In `output/claude/` you will find:

- `output/claude/qa_checks.md`

This file defines the required QA checks, signature validations, overlap audits, and remediation loop.  
**You must run QA and report it.** If QA fails, revise gates and re-run QA until acceptable (or explain clearly why not).

### 1.3 Dataset location
Assume there is a `data/` directory at the project root containing the CyTOF dataset for this run.

Your first job is to:
1. Locate the input data file(s) under `data/` (FCS, CSV, RDS, etc.).
2. Load the dataset into R.
3. Identify the feature/marker columns present.

If multiple datasets exist, pick the one referenced by the run context (if a clear dataset subfolder exists), otherwise pick the most “primary” dataset (largest, most complete), and document what you picked.

---

## 2) Global requirements / constraints

1. **Transformation:** Apply an **asinh transform with cofactor 5** to all marker channels used for gating and visualization.
   - Keep raw values if needed, but all gating decisions should be made in transformed space.

2. **Labeling policy:** Avoid a giant “other” bucket.
   - Use `unassigned` only when you genuinely cannot place a cell into any defined subset.
   - If a cell is clearly within a known branch (e.g., myeloid, lymphoid), prefer the closest matching label rather than dumping into `unassigned`.

3. **Reproducibility:** Every step must be reproducible from scripts saved into `output/claude/`.
   - Save **all scripts you create** (e.g., `output/claude/scripts/run_manual_gating.R`).
   - Save all intermediate gating objects (e.g., gate definitions, indices).

---

## 3) Deliverables (must produce all of these)

### 3.1 Plots for every gate
For **each gate node** in the gating tree:
- Create a plot of the parent population in the relevant marker space
- Overlay the polygon/rectangle gate
- Save as `output/claude/figures/gate_<NN>_<gate_name>.png`

Minimum expectation:
- Time gate
- Bead removal
- DNA+ gate
- Singlet gate
- Live (cisplatin) gate
- CD45.2+ Ter119− gate
- Every branch and terminal subset gate from `gating_strategy_for_agents.md`

Also produce an **index figure** showing the high-level tree structure (can be a simple diagram or a montage of key gates):
- `output/claude/figures/gating_tree_overview.png`

### 3.2 Create `llm_gate`
Add a new vector/column labeling each event:
- Column name: **`llm_gate`**
- Values: subset names exactly matching the terminal labels from `gating_strategy_for_agents.md`
- Any unlabeled events: `unassigned`

Save the label vector:
- `output/claude/llm_gate.rds` (an RDS containing the character vector in event order)

Also save a slim metadata table:
- `output/claude/llm_gate_metadata.csv` with at least:
  - `event_id` (or row index)
  - `llm_gate`

### 3.3 UMAP visualization and comparison
Subsample to **10,000 cells** (stratified if possible, so rare populations aren’t lost).

Create a UMAP on transformed marker space (exclude QC channels like Time, beadDist, DNA1/2 if appropriate; document which channels you use).

Save:
- `output/claude/figures/umap_llm_gate.png` (UMAP colored by `llm_gate`)

If a reference annotation exists:
- Look for `data/spitzer_gate.rds` first.
- If not present, search `data/` for any `*_gate.rds` or obvious annotation file.
- If none exists, still produce `umap_llm_gate.png` and skip the comparison section (but explain clearly in the report).

If a reference gate vector exists (e.g., `spitzer_gate`):
- `output/claude/figures/umap_reference_gate.png` (UMAP colored by reference gate)
- `output/claude/figures/umap_side_by_side.png` (reference vs llm_gate)

### 3.4 QA artifacts (required)
You must run the QA described in `output/claude/qa_checks.md` and save the resulting artifacts, including at minimum:
- `output/claude/objects/qa_gate_plot_manifest.csv`
- `output/claude/objects/qa_llm_gate_counts.csv`
- `output/claude/objects/qa_marker_medians_by_label.csv`
- `output/claude/objects/qa_signature_violations.csv`
- `output/claude/figures/qa_marker_heatmap.png`
- `output/claude/figures/qa_density_<marker_or_pair>.png` (multiple files)

If reference labels exist:
- `output/claude/objects/qa_confusion_table.csv`
- `output/claude/objects/qa_top_matches_by_llm_gate.csv`

### 3.5 Summary + comparison report (Markdown + knitted HTML)
Create:
1. `output/claude/gate_summary_report.md`
   - Dataset description (file used, number of events)
   - Marker list used for gating
   - Brief note on transformation
   - Table: counts per `llm_gate`
   - QA summary (pass/fail + key flags) referencing `qa_checks.md`
   - Notes: ambiguous gates, tricky thresholds, deviations from strategy
   - Revisions section if QA required adjustments

2. If reference exists: `output/claude/gate_comparison_report.md`
   - Counts per reference gate
   - Counts per `llm_gate`
   - Confusion table: `table(llm_gate, reference_gate)`
   - Short interpretation: where it matches, where it diverges, likely reasons
   - Include a few “spot-check” gate plots if useful

Finally, **knit markdown → HTML**:
- `output/claude/gate_summary_report.html`
- `output/claude/gate_comparison_report.html` (if applicable)

---

## 4) Execution plan (do this sequentially)

### Step 1 — Setup
1. Create output subfolders:
   - `output/claude/figures/`
   - `output/claude/scripts/`
   - `output/claude/objects/`
   - `output/claude/objects/qa_tables/` (optional but recommended)
2. Start an R script: `output/claude/scripts/run_manual_gating.R`
3. Load required libraries (choose sensible ones; examples: `flowCore`, `ggplot2`, `dplyr`, `uwot`, `Matrix`, `rmarkdown`, etc.).
4. Set a fixed seed for reproducibility.

### Step 2 — Load data and inspect markers
1. Load the dataset from `data/`.
2. Print and save:
   - number of events
   - column names
   - identify marker columns vs metadata/QC columns
3. Save this info to `output/claude/objects/dataset_overview.rds` and summarize in `gate_summary_report.md`.

### Step 3 — Transform
1. Apply `asinh(x / 5)` (cofactor 5) to gating-relevant channels.
2. Store transformed values (either overwrite or create a transformed matrix; document what you do).
3. From here onward, gate using transformed values.

### Step 4 — Manual gating (follow `gating_strategy_for_agents.md`)
Implement gates in the exact order:

1. Start from all events (after debarcoding if needed).
2. Apply QC gates: time → beads → DNA+ → singlets → live → CD45.2+ Ter119−.
3. Branch into mature immune vs progenitors (Lin−).
4. Gate terminal subsets exactly as described.

**For each gate:**
- Decide polygon vertices based on density/bi-modality in the transformed space.
- Save gate definitions (vertices, parent population name, axes used) into:
  - `output/claude/objects/gates_<gate_name>.rds` (or a single `gates_all.rds`)
- Save the gate plot figure to `output/claude/figures/`.

### Step 5 — Create `llm_gate`
1. Initialize all events as `unassigned`.
2. Assign terminal subset labels in a deterministic order:
   - Prefer assigning terminal gates first, and ensure no event belongs to two terminal subsets.
   - If overlap occurs, resolve by:
     - gating hierarchy precedence (child overrides parent)
     - document precedence decisions
3. Save:
   - `output/claude/llm_gate.rds`
   - `output/claude/llm_gate_metadata.csv`

### Step 6 — Run QA (mandatory; follow `qa_checks.md`)
1. Execute every check in `output/claude/qa_checks.md`.
2. Save all required QA tables and figures to the locations specified in `qa_checks.md`.
3. Determine pass/fail and list any flags in `gate_summary_report.md`.

**If QA fails:** you must iterate:
- adjust the relevant polygon gates (still polygon/rect only),
- regenerate `llm_gate`,
- re-run QA,
- document revisions (what changed and why) in the report,
- keep prior gate plots with version suffixes (`_v2`, `_v3`, etc.).

### Step 7 — UMAP (after QA is acceptable)
1. Subsample to 10k events (try stratified sampling using `llm_gate` so rare populations appear).
2. Run UMAP (document parameters).
3. Save UMAP plots:
   - `output/claude/figures/umap_llm_gate.png`
   - reference + side-by-side if applicable

### Step 8 — Reports and knitting
1. Write/refresh `output/claude/gate_summary_report.md` (include QA results).
2. If reference exists, write/refresh `output/claude/gate_comparison_report.md`.
3. Knit markdown → html (use `rmarkdown::render`), output HTML into `output/claude/`.

---

## 5) Notes on “proxy” gates (important)
The strategy file may label some gates as **proxy** (e.g., CLP proxy, Hardy A–C proxy).  
If markers needed for canonical definitions are missing, still perform the best-effort proxy gate as written, and clearly label it as proxy in reports and in the `llm_gate` values exactly as specified in the gating strategy file.

---

## 6) What “done” looks like
You are done when:
- All gate figures exist in `output/claude/figures/`
- `output/claude/llm_gate.rds` exists and matches event order
- QA artifacts exist and QA results are summarized in `gate_summary_report.md`
- UMAP plot(s) exist
- Reports exist as both `.md` and `.html`
- All code required to reproduce is in `output/claude/scripts/`

Proceed now.


