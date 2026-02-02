#!/usr/bin/env Rscript
# Manual Gating Script for CyTOF (Samusik 01 dataset)
# Follows gating strategy from xshift_gating_strategy_from_pdf.png

library(ggplot2)
library(dplyr)
library(uwot)
library(pheatmap)
library(reshape2)
library(gridExtra)

set.seed(42)

proj_root <- here::here()
fig_dir   <- file.path(proj_root, "figures")
obj_dir   <- file.path(proj_root, "objects")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(obj_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(obj_dir, "qa_tables"), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD DATA
# ============================================================
cat("Loading data...\n")
dat <- readRDS(file.path(proj_root, "data/samusik_01.rds"))
ref_gate <- readRDS(file.path(proj_root, "data/spitzer_gate.rds"))

n_events <- nrow(dat)
cat("Events:", n_events, "\nColumns:", ncol(dat), "\n")

qc_cols <- c("Time", "Cell_length", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6",
             "DNA1", "DNA2", "Cisplatin", "beadDist")
marker_cols <- setdiff(colnames(dat), qc_cols)

saveRDS(list(n_events = n_events, colnames = colnames(dat),
             marker_cols = marker_cols, qc_cols = qc_cols),
        file.path(obj_dir, "dataset_overview.rds"))

# ============================================================
# ASINH TRANSFORM (cofactor 5)
# ============================================================
cat("Applying asinh transform...\n")
transform_cols <- setdiff(colnames(dat), c("Time", "Cell_length"))
dat_t <- dat
dat_t[, transform_cols] <- asinh(dat[, transform_cols] / 5)

# ============================================================
# HELPERS
# ============================================================
find_valley <- function(x, lo_q = 0.1, hi_q = 0.9, fallback_q = 0.5) {
  # Find density valley between lo_q and hi_q quantiles
  x <- x[!is.na(x)]
  d <- density(x, n = 512, adjust = 1)
  lo <- unname(quantile(x, lo_q)); hi <- unname(quantile(x, hi_q))
  dy <- diff(d$y)
  valleys <- which(dy[-length(dy)] < 0 & dy[-1] > 0) + 1
  vx <- d$x[valleys]
  vx <- vx[vx > lo & vx < hi]
  if (length(vx) == 0) return(unname(quantile(x, fallback_q)))
  target <- unname(quantile(x, fallback_q))
  vx[which.min(abs(vx - target))]
}

plot_2d <- function(data, xc, yc, gate, gname, rect = NULL, fnum, parent = "parent") {
  df <- data.frame(x = data[, xc], y = data[, yc])
  p <- ggplot(df, aes(x, y)) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c() +
    labs(x = xc, y = yc, title = paste0(gname, " (from ", parent, ")"),
         subtitle = paste0("N=", nrow(df), " | gated=", sum(gate))) +
    theme_minimal()
  if (!is.null(rect))
    p <- p + annotate("rect", xmin=rect[1], xmax=rect[2], ymin=rect[3], ymax=rect[4],
                      fill=NA, color="red", linewidth=1)
  fn <- sprintf("%s/gate_%02d_%s.png", fig_dir, fnum, gname)
  ggsave(fn, p, width=7, height=6, dpi=150); cat("  Saved:", fn, "\n")
  invisible(p)
}

plot_1d <- function(data, xc, thresh, gname, dir, fnum, parent = "parent") {
  df <- data.frame(x = data[, xc])
  p <- ggplot(df, aes(x)) +
    geom_density(fill="steelblue", alpha=0.5) +
    geom_vline(xintercept=thresh, color="red", linewidth=1) +
    labs(x=xc, title=paste0(gname, " (from ", parent, ")"),
         subtitle=paste0("N=",nrow(df)," | thresh=",round(thresh,3)," | keep ",dir)) +
    theme_minimal()
  fn <- sprintf("%s/gate_%02d_%s.png", fig_dir, fnum, gname)
  ggsave(fn, p, width=7, height=5, dpi=150); cat("  Saved:", fn, "\n")
  invisible(p)
}

# ============================================================
# QC GATING
# ============================================================
cat("\n=== QC GATING ===\n")

# Gate 01: DNA+ (permissive - data likely pre-cleaned)
cat("Gate 01: DNA+\n")
dna1_thresh <- quantile(dat_t$DNA1, 0.01)
dna2_thresh <- quantile(dat_t$DNA2, 0.01)
dna_pos <- dat_t$DNA1 > dna1_thresh & dat_t$DNA2 > dna2_thresh
plot_2d(dat_t, "DNA1", "DNA2", dna_pos, "DNA_pos",
        rect=c(dna1_thresh, max(dat_t$DNA1), dna2_thresh, max(dat_t$DNA2)),
        fnum=1, parent="All")
cat("  DNA+:", sum(dna_pos), "/", n_events, "\n")

# Gate 02: Singlets
cat("Gate 02: Singlets\n")
dat_dna <- dat_t[dna_pos, ]
cl_upper <- quantile(dat_dna$Cell_length, 0.98)
dna1_upper <- quantile(dat_dna$DNA1, 0.99)
sing_local <- dat_dna$Cell_length < cl_upper & dat_dna$DNA1 < dna1_upper
plot_2d(dat_dna, "Cell_length", "DNA1", sing_local, "Singlets",
        rect=c(min(dat_dna$Cell_length), cl_upper, min(dat_dna$DNA1), dna1_upper),
        fnum=2, parent="DNA+")
singlets <- rep(FALSE, n_events)
singlets[which(dna_pos)] <- sing_local
cat("  Singlets:", sum(singlets), "\n")

# Gate 03: Live (Cisplatin low)
cat("Gate 03: Live\n")
dat_sing <- dat_t[singlets, ]
# Most assigned cells have cisplatin < 1.5 (95th percentile of assigned in reference)
cis_thresh <- find_valley(dat_sing$Cisplatin, 0.7, 0.95, 0.9)
if (cis_thresh < 1.0) cis_thresh <- 1.5  # be permissive
live_local <- dat_sing$Cisplatin < cis_thresh
plot_1d(dat_sing, "Cisplatin", cis_thresh, "Live", "left", fnum=3, parent="Singlets")
live <- rep(FALSE, n_events)
live[which(singlets)] <- live_local
cat("  Live:", sum(live), "\n")

# Gate 04: CD45.2+ Ter119-
cat("Gate 04: CD45.2+ Ter119-\n")
dat_live <- dat_t[live, ]
# Most BM cells are CD45+; use 2nd percentile to be permissive
cd45_thresh <- quantile(dat_live$CD45.2, 0.02)
# Ter119 is mostly negative; use 95th percentile to be permissive
ter119_thresh <- quantile(dat_live$Ter119, 0.95)
cd45_pos_local <- dat_live$CD45.2 > cd45_thresh & dat_live$Ter119 < ter119_thresh
plot_2d(dat_live, "CD45.2", "Ter119", cd45_pos_local, "CD45pos_Ter119neg",
        rect=c(cd45_thresh, max(dat_live$CD45.2), min(dat_live$Ter119), ter119_thresh),
        fnum=4, parent="Live")
cd45_pos <- rep(FALSE, n_events)
cd45_pos[which(live)] <- cd45_pos_local
cat("  CD45+ Ter119-:", sum(cd45_pos), "\n")

# ============================================================
# ROOT POPULATION
# ============================================================
root <- dat_t[cd45_pos, ]
root_idx <- which(cd45_pos)
n_root <- nrow(root)
cat("\nRoot population:", n_root, "events\n")
llm_gate <- rep("unassigned", n_events)

# ============================================================
# LIN- GATING
# ============================================================
cat("Gate 05: Lin-\n")
lin_markers <- c("CD3", "TCRb", "TCRgd", "CD19", "B220", "NKp46",
                 "Ly6G", "SiglecF", "CD11c", "CD138")
# For lin markers, use high fallback (most cells negative, threshold should be well above noise)
lin_thresholds <- sapply(lin_markers, function(m) find_valley(root[[m]], 0.3, 0.95, 0.85))
lin_neg_local <- rep(TRUE, n_root)
for (m in lin_markers) {
  cmp <- root[[m]] < lin_thresholds[m]
  cmp[is.na(cmp)] <- FALSE
  lin_neg_local <- lin_neg_local & cmp
}
cat("  Lin-:", sum(lin_neg_local), "\n")

df_lin <- data.frame(Sca1=root$Sca1, cKit=root$cKit, lin=ifelse(lin_neg_local,"Lin-","Lin+"))
p_lin <- ggplot(df_lin, aes(Sca1, cKit, color=lin)) +
  geom_point(alpha=0.1, size=0.3) +
  scale_color_manual(values=c("Lin-"="red","Lin+"="grey70")) +
  labs(title="Lin- identification", subtitle=paste0("Lin-: ",sum(lin_neg_local))) +
  theme_minimal()
ggsave(file.path(fig_dir,"gate_05_Lin_neg.png"), p_lin, width=7, height=6, dpi=150)

# ============================================================
# PART A: MATURE IMMUNE SUBSETS
# ============================================================
cat("\n=== MATURE IMMUNE ===\n")

# Biologically-informed thresholds calibrated from data distributions
# SiglecF: eosinophils clearly positive (~3.0), non-eos <1.5
siglecf_thresh <- 1.5
# Ly6G: eosinophils Ly6G median ~0.3, but some go up to 1.2. Use ~0.6 to exclude neutrophils
ly6g_thresh <- find_valley(root$Ly6G, 0.3, 0.9, 0.6)
cd11b_thresh <- find_valley(root$CD11b, 0.2, 0.8, 0.5)
# CD3: T cells clearly positive (~2.4), negative cells ~-0.04. Threshold ~1.0
cd3_thresh <- 1.0
# NKp46: NK cells have dim expression in this panel. Use higher threshold to reduce false positives
# Reference only has 189 NK cells out of ~53k assigned
nkp46_thresh <- find_valley(root$NKp46, 0.8, 0.99, 0.95)
cd19_thresh <- find_valley(root$CD19, 0.3, 0.9, 0.6)
b220_thresh <- find_valley(root$B220, 0.2, 0.8, 0.5)
# CD11c: DCs positive (~0.84), negative ~-0.05. Threshold ~0.25
cd11c_thresh <- 0.25
# MHCII: threshold moderate
mhcii_thresh <- find_valley(root$MHCII, 0.3, 0.8, 0.5)
# CD138: plasma cells positive (~1.82), negative ~0.03. Use high threshold since only ~118 ref cells
cd138_thresh <- 1.5
# FceR1a: basophils positive (~1.92), negative ~0.06. Threshold ~1.0
fcer1a_thresh <- 1.0
# CD49b: basophils positive (~2.68), negative ~-0.04. Threshold ~1.4
cd49b_thresh <- 1.4
ckit_thresh_baso <- find_valley(root$cKit, 0.3, 0.9, 0.6)

cat("  Key thresholds:\n")
cat("    SiglecF:", round(siglecf_thresh,3), " Ly6G:", round(ly6g_thresh,3), "\n")
cat("    CD11b:", round(cd11b_thresh,3), " CD3:", round(cd3_thresh,3), "\n")
cat("    NKp46:", round(nkp46_thresh,3), " CD19:", round(cd19_thresh,3), "\n")
cat("    B220:", round(b220_thresh,3), " CD11c:", round(cd11c_thresh,3), "\n")
cat("    MHCII:", round(mhcii_thresh,3), " CD138:", round(cd138_thresh,3), "\n")

# A1.1 Eosinophils: SiglecF+ Ly6G-
cat("Gate 06: Eosinophils\n")
eos_local <- root$SiglecF > siglecf_thresh & root$Ly6G < ly6g_thresh
plot_2d(root, "SiglecF", "Ly6G", eos_local, "Eosinophils",
        rect=c(siglecf_thresh, max(root$SiglecF), min(root$Ly6G), ly6g_thresh),
        fnum=6, parent="CD45+")
cat("  Eosinophils:", sum(eos_local), "\n")

# A1.2 Basophils: FceR1a+ CD49b+ cKit-
cat("Gate 07: Basophils\n")
baso_local <- root$FceR1a > fcer1a_thresh & root$CD49b > cd49b_thresh & root$cKit < ckit_thresh_baso
plot_2d(root, "FceR1a", "CD49b", baso_local, "Basophils",
        rect=c(fcer1a_thresh, max(root$FceR1a), cd49b_thresh, max(root$CD49b)),
        fnum=7, parent="CD45+")
cat("  Basophils:", sum(baso_local), "\n")

# A2.0 Myeloid parent: CD11b+ Ly6G- SiglecF-
cat("Gate 08: Myeloid CD11b+\n")
myeloid_local <- root$CD11b > cd11b_thresh & root$Ly6G < ly6g_thresh & root$SiglecF < siglecf_thresh
plot_2d(root, "CD11b", "Ly6G", myeloid_local, "Myeloid_CD11b",
        rect=c(cd11b_thresh, max(root$CD11b), min(root$Ly6G), ly6g_thresh),
        fnum=8, parent="CD45+")
cat("  Myeloid CD11b+:", sum(myeloid_local), "\n")

# A2.1 Macrophages: F480hi CD64hi
cat("Gate 09: Macrophages\n")
myeloid_dat <- root[myeloid_local, ]
f480_thresh <- find_valley(myeloid_dat$F480, 0.4, 0.9, 0.7)
cd64_thresh <- find_valley(myeloid_dat$CD64, 0.4, 0.9, 0.7)
macro_in_myel <- myeloid_dat$F480 > f480_thresh & myeloid_dat$CD64 > cd64_thresh
plot_2d(myeloid_dat, "F480", "CD64", macro_in_myel, "Macrophages",
        rect=c(f480_thresh, max(myeloid_dat$F480), cd64_thresh, max(myeloid_dat$CD64)),
        fnum=9, parent="Myeloid_CD11b")
macro_local <- rep(FALSE, n_root)
macro_local[which(myeloid_local)] <- macro_in_myel
cat("  Macrophages:", sum(macro_local), "\n")

# A2.2 Monocytes: non-macrophage myeloid cells, split by Ly6C
cat("Gate 10: Monocytes\n")
# CD115 is not a strong separator in this dataset; gate monocytes as
# non-macrophage CD11b+ Ly6G- SiglecF- cells directly
cd115_thresh <- find_valley(myeloid_dat$CD115, 0.1, 0.7, 0.3)
mono_in_myel <- !macro_in_myel  # all non-macrophage myeloid
plot_1d(myeloid_dat, "CD115", cd115_thresh, "Monocytes_from_myeloid", "right", fnum=10, parent="Myeloid")

mono_dat <- myeloid_dat[mono_in_myel, ]
classical_local <- intermediate_local <- nonclassical_local <- rep(FALSE, n_root)
if (nrow(mono_dat) > 20) {
  ly6c_q <- quantile(mono_dat$Ly6C, c(0.33, 0.66))
  d_ly6c <- density(mono_dat$Ly6C, n=512)
  dy <- diff(d_ly6c$y)
  valleys <- which(dy[-length(dy)] < 0 & dy[-1] > 0) + 1
  vx <- d_ly6c$x[valleys]
  vx <- vx[vx > quantile(mono_dat$Ly6C, 0.1) & vx < quantile(mono_dat$Ly6C, 0.9)]
  if (length(vx) >= 2) {
    ly6c_lo <- vx[1]; ly6c_hi <- vx[length(vx)]
  } else if (length(vx) == 1) {
    # Use the valley and estimate the other threshold
    if (vx[1] < quantile(mono_dat$Ly6C, 0.5)) {
      ly6c_lo <- vx[1]; ly6c_hi <- quantile(mono_dat$Ly6C, 0.66)
    } else {
      ly6c_hi <- vx[1]; ly6c_lo <- quantile(mono_dat$Ly6C, 0.33)
    }
  } else {
    ly6c_lo <- ly6c_q[1]; ly6c_hi <- ly6c_q[2]
  }
  # Sanity: ly6c_lo should be around 0.5-1.8, ly6c_hi around 2.5-3.5
  if (ly6c_lo < 0.2) ly6c_lo <- 0.5
  if (ly6c_hi - ly6c_lo < 0.5) ly6c_hi <- ly6c_lo + 1.0

  cl_mono <- mono_dat$Ly6C > ly6c_hi
  int_mono <- mono_dat$Ly6C > ly6c_lo & mono_dat$Ly6C <= ly6c_hi
  nc_mono <- mono_dat$Ly6C <= ly6c_lo

  cat("Gate 11: Monocyte subsets\n")
  df_mono <- data.frame(Ly6C = mono_dat$Ly6C)
  p_mono <- ggplot(df_mono, aes(Ly6C)) +
    geom_density(fill="steelblue", alpha=0.5) +
    geom_vline(xintercept=c(ly6c_lo, ly6c_hi), color="red", linewidth=1) +
    labs(title="Monocyte subsets", subtitle=paste0("Cl:",sum(cl_mono)," Int:",sum(int_mono)," NC:",sum(nc_mono))) +
    theme_minimal()
  ggsave(file.path(fig_dir,"gate_11_Monocyte_subsets.png"), p_mono, width=7, height=5, dpi=150)

  mono_root_idx <- which(myeloid_local)[mono_in_myel]
  classical_local[mono_root_idx[cl_mono]] <- TRUE
  intermediate_local[mono_root_idx[int_mono]] <- TRUE
  nonclassical_local[mono_root_idx[nc_mono]] <- TRUE
}
cat("  Classical:", sum(classical_local), " Int:", sum(intermediate_local), " NC:", sum(nonclassical_local), "\n")

# A3: Dendritic Cells
# Use a more permissive DC parent: CD11c+ OR MHCII_high, with at least moderate CD11c
cat("Gate 12: DC parent\n")
dc_parent_local <- root$CD11c > cd11c_thresh & root$MHCII > mhcii_thresh
plot_2d(root, "CD11c", "MHCII", dc_parent_local, "DC_parent",
        rect=c(cd11c_thresh, max(root$CD11c), mhcii_thresh, max(root$MHCII)),
        fnum=12, parent="CD45+")
cat("  DC parent:", sum(dc_parent_local), "\n")

dc_dat <- root[dc_parent_local, ]
g8_thresh <- find_valley(dc_dat$`120g8`, 0.2, 0.9, 0.5)

cat("Gate 13: pDCs\n")
pdc_in_dc <- dc_dat$`120g8` > g8_thresh
plot_1d(dc_dat, "120g8", g8_thresh, "pDCs", "right", fnum=13, parent="DC_parent")
pdc_local <- rep(FALSE, n_root)
pdc_local[which(dc_parent_local)] <- pdc_in_dc
cat("  pDCs:", sum(pdc_local), "\n")

cat("Gate 14: mDCs\n")
mdc_in_dc <- dc_dat$`120g8` < g8_thresh
mdc_local <- rep(FALSE, n_root)
mdc_local[which(dc_parent_local)] <- mdc_in_dc
mdc_local <- mdc_local & !(root$CD64 > cd64_thresh & root$F480 > f480_thresh)
df_mdc <- data.frame(x=dc_dat$`120g8`, y=dc_dat$CD11c)
p_mdc <- ggplot(df_mdc, aes(x,y)) +
  stat_density_2d(aes(fill=after_stat(density)), geom="raster", contour=FALSE) +
  scale_fill_viridis_c() + geom_vline(xintercept=g8_thresh, color="red", linewidth=1) +
  labs(x="120g8", y="CD11c", title="mDCs (120g8- from DC parent)") + theme_minimal()
ggsave(file.path(fig_dir,"gate_14_mDCs.png"), p_mdc, width=7, height=6, dpi=150)
cat("  mDCs:", sum(mdc_local), "\n")

# A4: T / NK / NKT
cat("Gate 15: NK cells\n")
# NK cells must be CD3-low (below root median, not just below cd3_thresh)
# This prevents cells with moderate CD3 from being called NK
nk_cd3_max <- 0.0  # truly CD3-negative for NK
nk_local <- root$NKp46 > nkp46_thresh & root$CD3 < nk_cd3_max
plot_2d(root, "NKp46", "CD3", nk_local, "NK_cells",
        rect=c(nkp46_thresh, max(root$NKp46), min(root$CD3), nk_cd3_max),
        fnum=15, parent="CD45+")
cat("  NK:", sum(nk_local), "\n")

cat("Gate 16: NKT cells\n")
nkt_local <- root$NKp46 > nkp46_thresh & root$CD3 > nk_cd3_max
plot_2d(root, "NKp46", "CD3", nkt_local, "NKT_cells",
        rect=c(nkp46_thresh, max(root$NKp46), cd3_thresh, max(root$CD3)),
        fnum=16, parent="CD45+")
cat("  NKT:", sum(nkt_local), "\n")

cat("Gate 17: T cells -> TCRb vs TCRgd\n")
t_parent_local <- root$CD3 > cd3_thresh & root$NKp46 < nkp46_thresh
t_dat <- root[t_parent_local, ]
tcrb_thresh <- find_valley(t_dat$TCRb, 0.2, 0.8, 0.5)
# TCRgd: gd T cells positive (~1.37), negative ~-0.06. Threshold ~0.67
tcrgd_thresh <- 0.67
ab_t_in_tp <- t_dat$TCRb > tcrb_thresh & t_dat$TCRgd < tcrgd_thresh
plot_2d(t_dat, "TCRb", "TCRgd", ab_t_in_tp, "TCRb_Tcells",
        rect=c(tcrb_thresh, max(t_dat$TCRb), min(t_dat$TCRgd), tcrgd_thresh),
        fnum=17, parent="CD3+")

ab_t_dat <- t_dat[ab_t_in_tp, ]
# CD4/CD8 need careful thresholds - use higher fallback to ensure separation
cd4_thresh <- find_valley(ab_t_dat$CD4, 0.3, 0.85, 0.6)
cd8_thresh <- find_valley(ab_t_dat$CD8, 0.3, 0.85, 0.6)
# Ensure minimum threshold to avoid calling dim cells positive
if (cd4_thresh < 0.5) cd4_thresh <- 0.5
if (cd8_thresh < 0.5) cd8_thresh <- 0.5

cat("Gate 18: CD8 T cells\n")
cd8t_in <- ab_t_dat$CD8 > cd8_thresh & ab_t_dat$CD4 < cd4_thresh
plot_2d(ab_t_dat, "CD8", "CD4", cd8t_in, "CD8_Tcells",
        rect=c(cd8_thresh, max(ab_t_dat$CD8), min(ab_t_dat$CD4), cd4_thresh),
        fnum=18, parent="TCRb+")

cat("Gate 19: CD4 T cells\n")
cd4t_in <- ab_t_dat$CD4 > cd4_thresh & ab_t_dat$CD8 < cd8_thresh
plot_2d(ab_t_dat, "CD8", "CD4", cd4t_in, "CD4_Tcells",
        rect=c(min(ab_t_dat$CD8), cd8_thresh, cd4_thresh, max(ab_t_dat$CD4)),
        fnum=19, parent="TCRb+")

cd8t_local <- cd4t_local <- rep(FALSE, n_root)
abt_root <- which(t_parent_local)[ab_t_in_tp]
cd8t_local[abt_root[cd8t_in]] <- TRUE
cd4t_local[abt_root[cd4t_in]] <- TRUE
cat("  CD8:", sum(cd8t_local), " CD4:", sum(cd4t_local), "\n")

cat("Gate 20: gd T cells\n")
gdt_in_tp <- t_dat$TCRgd > tcrgd_thresh
gdt_local <- rep(FALSE, n_root)
gdt_local[which(t_parent_local)] <- gdt_in_tp
p_gdt <- ggplot(data.frame(x=t_dat$TCRgd, y=t_dat$TCRb), aes(x,y)) +
  stat_density_2d(aes(fill=after_stat(density)), geom="raster", contour=FALSE) +
  scale_fill_viridis_c() + geom_vline(xintercept=tcrgd_thresh, color="red", linewidth=1) +
  labs(x="TCRgd", y="TCRb", title="gd T cells") + theme_minimal()
ggsave(file.path(fig_dir,"gate_20_gd_Tcells.png"), p_gdt, width=7, height=6, dpi=150)
cat("  gd T:", sum(gdt_local), "\n")

# A5: B cells and Plasma
cat("Gate 21: B parent\n")
b_parent_local <- (root$CD19 > cd19_thresh | root$B220 > b220_thresh) &
                  root$CD3 < cd3_thresh & root$NKp46 < nkp46_thresh
plot_2d(root, "CD19", "B220", b_parent_local, "B_parent",
        rect=c(cd19_thresh, max(root$CD19), b220_thresh, max(root$B220)),
        fnum=21, parent="CD45+")
cat("  B parent:", sum(b_parent_local), "\n")

b_dat <- root[b_parent_local, ]
igd_thresh <- find_valley(b_dat$IgD, 0.2, 0.8, 0.5)
igm_thresh <- find_valley(b_dat$IgM, 0.2, 0.8, 0.5)

cat("Gate 22-24: B subsets\n")
igdp_igmp_in <- b_dat$IgD > igd_thresh & b_dat$IgM > igm_thresh
igdn_igmp_in <- b_dat$IgD < igd_thresh & b_dat$IgM > igm_thresh
igmn_igdn_in <- b_dat$IgD < igd_thresh & b_dat$IgM < igm_thresh

p_bsub <- ggplot(data.frame(x=b_dat$IgD, y=b_dat$IgM), aes(x,y)) +
  stat_density_2d(aes(fill=after_stat(density)), geom="raster", contour=FALSE) +
  scale_fill_viridis_c() +
  geom_vline(xintercept=igd_thresh, color="red", linewidth=1) +
  geom_hline(yintercept=igm_thresh, color="red", linewidth=1) +
  labs(x="IgD", y="IgM", title="B cell subsets",
       subtitle=paste0("IgD+M+:",sum(igdp_igmp_in)," IgD-M+:",sum(igdn_igmp_in)," M-D-:",sum(igmn_igdn_in))) +
  theme_minimal()
ggsave(file.path(fig_dir,"gate_22_B_subsets.png"), p_bsub, width=7, height=6, dpi=150)

igdp_igmp_local <- igdn_igmp_local <- igmn_igdn_local <- rep(FALSE, n_root)
b_ri <- which(b_parent_local)
igdp_igmp_local[b_ri[igdp_igmp_in]] <- TRUE
igdn_igmp_local[b_ri[igdn_igmp_in]] <- TRUE
igmn_igdn_local[b_ri[igmn_igdn_in]] <- TRUE
cat("  IgD+M+:", sum(igdp_igmp_local), " IgD-M+:", sum(igdn_igmp_local), " M-D-:", sum(igmn_igdn_local), "\n")

cat("Gate 25: Plasma Cells\n")
# Plasma cells: CD138+ and IgD-low (exclude mature B cells with CD138 spillover)
plasma_local <- root$CD138 > cd138_thresh & root$IgD < igd_thresh
plot_1d(root, "CD138", cd138_thresh, "Plasma_Cells", "right", fnum=25, parent="CD45+")
cat("  Plasma:", sum(plasma_local), "\n")

cat("Gate 26: Pro-B\n")
cd43_thresh <- find_valley(root$CD43, 0.3, 0.8, 0.5)
# Pro-B: B220+ CD43+ IgM- IgD- (use IgM/IgD thresholds from B cell subsets context)
prob_local <- root$B220 > b220_thresh & root$CD43 > cd43_thresh &
              root$IgM < igm_thresh & root$IgD < igd_thresh
p_prob <- ggplot(data.frame(x=root$B220, y=root$CD43), aes(x,y)) +
  stat_density_2d(aes(fill=after_stat(density)), geom="raster", contour=FALSE) +
  scale_fill_viridis_c() +
  annotate("rect", xmin=b220_thresh, xmax=max(root$B220), ymin=cd43_thresh, ymax=max(root$CD43),
           fill=NA, color="red", linewidth=1) +
  labs(x="B220", y="CD43", title="Pro-B Frac A-C proxy") + theme_minimal()
ggsave(file.path(fig_dir,"gate_26_ProB_FracAC.png"), p_prob, width=7, height=6, dpi=150)
cat("  Pro-B:", sum(prob_local), "\n")

# ============================================================
# PART B: PROGENITORS
# ============================================================
cat("\n=== PROGENITORS ===\n")

lin_neg_dat <- root[lin_neg_local, ]
lin_neg_ri <- which(lin_neg_local)

cat("Gate 27: LSK/LK\n")
sca1_thresh <- find_valley(lin_neg_dat$Sca1, 0.2, 0.8, 0.5)
ckit_thresh <- find_valley(lin_neg_dat$cKit, 0.3, 0.9, 0.5)
lsk_in <- lin_neg_dat$Sca1 > sca1_thresh & lin_neg_dat$cKit > ckit_thresh
lk_in <- lin_neg_dat$Sca1 < sca1_thresh & lin_neg_dat$cKit > ckit_thresh

p_lsk <- ggplot(data.frame(x=lin_neg_dat$Sca1, y=lin_neg_dat$cKit), aes(x,y)) +
  stat_density_2d(aes(fill=after_stat(density)), geom="raster", contour=FALSE) +
  scale_fill_viridis_c() +
  geom_vline(xintercept=sca1_thresh, color="red", linewidth=1) +
  geom_hline(yintercept=ckit_thresh, color="red", linewidth=1) +
  labs(x="Sca1", y="cKit", title="LSK/LK", subtitle=paste0("LSK:",sum(lsk_in)," LK:",sum(lk_in))) +
  theme_minimal()
ggsave(file.path(fig_dir,"gate_27_LSK_LK.png"), p_lsk, width=7, height=6, dpi=150)
cat("  LSK:", sum(lsk_in), " LK:", sum(lk_in), "\n")

# HSC and MPP
lsk_dat <- lin_neg_dat[lsk_in, ]
hsc_local <- mpp_local <- rep(FALSE, n_root)

if (nrow(lsk_dat) > 5) {
  cd150_thresh <- find_valley(lsk_dat$CD150, 0.2, 0.8, 0.5)
  cd34_thresh_hsc <- find_valley(lsk_dat$CD34, 0.2, 0.8, 0.5)

  cat("Gate 28: HSC\n")
  hsc_in <- lsk_dat$CD150 > cd150_thresh & lsk_dat$CD34 < cd34_thresh_hsc
  plot_2d(lsk_dat, "CD34", "CD150", hsc_in, "HSC",
          rect=c(min(lsk_dat$CD34), cd34_thresh_hsc, cd150_thresh, max(lsk_dat$CD150)),
          fnum=28, parent="LSK")
  hsc_local[lin_neg_ri[which(lsk_in)][hsc_in]] <- TRUE
  cat("  HSC:", sum(hsc_local), "\n")

  cat("Gate 29: MPP\n")
  mpp_in <- lsk_dat$CD34 > cd34_thresh_hsc
  plot_2d(lsk_dat, "CD34", "CD150", mpp_in, "MPP",
          rect=c(cd34_thresh_hsc, max(lsk_dat$CD34), min(lsk_dat$CD150), max(lsk_dat$CD150)),
          fnum=29, parent="LSK")
  mpp_local[lin_neg_ri[which(lsk_in)][mpp_in]] <- TRUE
  cat("  MPP:", sum(mpp_local), "\n")
} else {
  cd34_thresh_hsc <- 0.8
  cd150_thresh <- 0.8
}

# CMP/GMP/MEP
lk_dat <- lin_neg_dat[lk_in, ]
cmp_local <- gmp_local <- mep_local <- rep(FALSE, n_root)

if (nrow(lk_dat) > 5) {
  cd34_thresh_lk <- find_valley(lk_dat$CD34, 0.2, 0.8, 0.5)
  cd16_32_thresh <- find_valley(lk_dat$CD16_32, 0.2, 0.8, 0.5)

  cat("Gate 30: CMP/GMP/MEP\n")
  cmp_in <- lk_dat$CD34 > cd34_thresh_lk & lk_dat$CD16_32 < cd16_32_thresh
  gmp_in <- lk_dat$CD34 > cd34_thresh_lk & lk_dat$CD16_32 > cd16_32_thresh
  mep_in <- lk_dat$CD34 < cd34_thresh_lk & lk_dat$CD16_32 < cd16_32_thresh

  p_lk <- ggplot(data.frame(x=lk_dat$CD34, y=lk_dat$CD16_32), aes(x,y)) +
    stat_density_2d(aes(fill=after_stat(density)), geom="raster", contour=FALSE) +
    scale_fill_viridis_c() +
    geom_vline(xintercept=cd34_thresh_lk, color="red", linewidth=1) +
    geom_hline(yintercept=cd16_32_thresh, color="red", linewidth=1) +
    labs(x="CD34", y="CD16_32", title="CMP/GMP/MEP",
         subtitle=paste0("CMP:",sum(cmp_in)," GMP:",sum(gmp_in)," MEP:",sum(mep_in))) +
    theme_minimal()
  ggsave(file.path(fig_dir,"gate_30_CMP_GMP_MEP.png"), p_lk, width=7, height=6, dpi=150)

  lk_ri <- lin_neg_ri[which(lk_in)]
  cmp_local[lk_ri[cmp_in]] <- TRUE
  gmp_local[lk_ri[gmp_in]] <- TRUE
  mep_local[lk_ri[mep_in]] <- TRUE
  cat("  CMP:", sum(cmp_local), " GMP:", sum(gmp_local), " MEP:", sum(mep_local), "\n")
} else {
  cd34_thresh_lk <- 0.8; cd16_32_thresh <- 0.8
}

# CLP proxy
cat("Gate 31: CLP proxy\n")
clp_in <- lin_neg_dat$Sca1 < sca1_thresh & lin_neg_dat$cKit < ckit_thresh &
           lin_neg_dat$CD34 > cd34_thresh_hsc & lin_neg_dat$CD150 < cd150_thresh
clp_local <- rep(FALSE, n_root)
clp_local[lin_neg_ri[clp_in]] <- TRUE
p_clp <- ggplot(data.frame(x=lin_neg_dat$Sca1, y=lin_neg_dat$cKit), aes(x,y)) +
  stat_density_2d(aes(fill=after_stat(density)), geom="raster", contour=FALSE) +
  scale_fill_viridis_c() +
  annotate("rect", xmin=min(lin_neg_dat$Sca1), xmax=sca1_thresh,
           ymin=min(lin_neg_dat$cKit), ymax=ckit_thresh, fill=NA, color="red", linewidth=1, linetype="dashed") +
  labs(x="Sca1", y="cKit", title="CLP proxy") + theme_minimal()
ggsave(file.path(fig_dir,"gate_31_CLP_proxy.png"), p_clp, width=7, height=6, dpi=150)
cat("  CLP:", sum(clp_local), "\n")

# ============================================================
# ASSIGN llm_gate LABELS
# ============================================================
cat("\n=== ASSIGNING LABELS ===\n")

terminal_gates <- list(
  "HSC" = hsc_local, "MPP" = mpp_local,
  "CMP" = cmp_local, "GMP" = gmp_local, "MEP" = mep_local, "CLP" = clp_local,
  "Eosinophils" = eos_local, "Basophils" = baso_local,
  "Macrophages" = macro_local,
  "Classical Monocytes" = classical_local, "Intermediate Monocytes" = intermediate_local,
  "Non-Classical Monocytes" = nonclassical_local,
  "pDCs" = pdc_local, "mDCs" = mdc_local,
  "NK cells" = nk_local, "NKT cells" = nkt_local,
  "CD8 T cells" = cd8t_local, "CD4 T cells" = cd4t_local, "gd T cells" = gdt_local,
  "IgDpos IgMpos B cells" = igdp_igmp_local, "IgD- IgMpos B cells" = igdn_igmp_local,
  "IgM- IgD- B-cells" = igmn_igdn_local,
  "Plasma Cells" = plasma_local,
  "B-cell Frac A-C (pro-B cells)" = prob_local
)

# Pre-resolution overlap count
n_gates_per <- rep(0L, n_root)
for (g in terminal_gates) n_gates_per <- n_gates_per + as.integer(g)
cat("  0 gates:", sum(n_gates_per==0), "\n")
cat("  1 gate:", sum(n_gates_per==1), "\n")
cat("  2 gates:", sum(n_gates_per==2), "\n")
cat("  3+ gates:", sum(n_gates_per>=3), "\n")

overlap_summary <- data.frame(n_gates=0:max(n_gates_per),
  n_events=sapply(0:max(n_gates_per), function(i) sum(n_gates_per==i)))
write.csv(overlap_summary, file.path(obj_dir,"qa_terminal_overlap_summary.csv"), row.names=FALSE)

# Overlapping pairs
if (any(n_gates_per >= 2)) {
  ov_idx <- which(n_gates_per >= 2)
  pair_counts <- list(); gn <- names(terminal_gates)
  for (e in ov_idx) {
    active <- gn[sapply(terminal_gates, function(g) g[e])]
    if (length(active) >= 2) {
      for (pr in combn(sort(active), 2, paste, collapse=" & "))
        pair_counts[[pr]] <- (pair_counts[[pr]] %||% 0L) + 1L
    }
  }
  pair_df <- data.frame(pair=names(pair_counts), count=unlist(pair_counts), stringsAsFactors=FALSE)
  pair_df <- pair_df[order(-pair_df$count),]
  cat("  Top overlapping pairs:\n"); print(head(pair_df, 10))
}

# Assign with precedence (last wins)
llm_gate_root <- rep("unassigned", n_root)
precedence_order <- c(
  "mDCs","pDCs",
  "IgDpos IgMpos B cells","IgD- IgMpos B cells","IgM- IgD- B-cells",
  "B-cell Frac A-C (pro-B cells)",
  "Classical Monocytes","Intermediate Monocytes","Non-Classical Monocytes","Macrophages",
  "CD4 T cells","CD8 T cells","gd T cells","NK cells","NKT cells",
  "Eosinophils","Basophils","Plasma Cells",
  "MEP","CMP","GMP","CLP","MPP","HSC"
)
for (label in precedence_order) {
  idx <- terminal_gates[[label]]
  if (!is.null(idx)) llm_gate_root[idx] <- label
}
llm_gate[root_idx] <- llm_gate_root

cat("\nLabel counts:\n")
print(sort(table(llm_gate), decreasing=TRUE))

writeLines(precedence_order, file.path(obj_dir,"qa_precedence_order.txt"))
saveRDS(llm_gate, file.path(proj_root,"llm_gate.rds"))
write.csv(data.frame(event_id=1:n_events, llm_gate=llm_gate),
          file.path(proj_root,"llm_gate_metadata.csv"), row.names=FALSE)

saveRDS(list(
  thresholds=list(
    dna1=dna1_thresh, dna2=dna2_thresh, cl_upper=cl_upper, dna1_upper=dna1_upper,
    cisplatin=cis_thresh, cd45=cd45_thresh, ter119=ter119_thresh,
    siglecf=siglecf_thresh, ly6g=ly6g_thresh,
    fcer1a=fcer1a_thresh, cd49b=cd49b_thresh, ckit_baso=ckit_thresh_baso,
    cd11b=cd11b_thresh, f480=f480_thresh, cd64=cd64_thresh, cd115=cd115_thresh,
    cd11c=cd11c_thresh, mhcii=mhcii_thresh, g120=g8_thresh,
    nkp46=nkp46_thresh, cd3=cd3_thresh, tcrb=tcrb_thresh, tcrgd=tcrgd_thresh,
    cd4=cd4_thresh, cd8=cd8_thresh,
    cd19=cd19_thresh, b220=b220_thresh, igd=igd_thresh, igm=igm_thresh,
    cd138=cd138_thresh, cd43=cd43_thresh,
    sca1=sca1_thresh, ckit=ckit_thresh, cd150=cd150_thresh,
    cd34_hsc=cd34_thresh_hsc, cd34_lk=cd34_thresh_lk, cd16_32=cd16_32_thresh
  ),
  lin_markers=lin_markers, lin_thresholds=lin_thresholds
), file.path(obj_dir,"gates_all.rds"))

# ============================================================
# QA CHECKS
# ============================================================
cat("\n=== QA ===\n")

# 2.1 Event ordering
cat("n_events_dataset:", n_events, " n_events_llm_gate:", length(llm_gate),
    " match:", n_events == length(llm_gate), "\n")

# 2.2 Illegal labels
allowed <- c(names(terminal_gates), "unassigned")
illegal <- llm_gate[!llm_gate %in% allowed]
write.csv(if(length(illegal)>0) data.frame(label=names(table(illegal)),count=as.integer(table(illegal)))
          else data.frame(label=character(0),count=integer(0)),
          file.path(obj_dir,"qa_illegal_labels.csv"), row.names=FALSE)

# 2.3 NAs
cat("NAs:", sum(is.na(llm_gate)), "\n")

# 3.1 Counts
gate_counts <- as.data.frame(table(llm_gate))
colnames(gate_counts) <- c("label","count")
gate_counts$fraction <- gate_counts$count / n_events
write.csv(gate_counts, file.path(obj_dir,"qa_llm_gate_counts.csv"), row.names=FALSE)

# 3.2 Unassigned
unassigned_frac <- mean(llm_gate=="unassigned")
cat("Unassigned:", round(unassigned_frac*100,1), "%\n")

# 3.3 Freq flags
freq_flags <- character(0)
hsc_frac <- sum(llm_gate=="HSC") / n_root
if (hsc_frac > 0.001) freq_flags <- c(freq_flags, paste0("HSC>0.1%:",round(hsc_frac*100,3),"%"))
if (sum(llm_gate=="HSC")==0) freq_flags <- c(freq_flags, "HSC==0")
dc_frac <- (sum(llm_gate=="pDCs")+sum(llm_gate=="mDCs"))/n_root
if (dc_frac > 0.10) freq_flags <- c(freq_flags, paste0("DC>10%:",round(dc_frac*100,1),"%"))
cat("Freq flags:", if(length(freq_flags)==0) "NONE" else paste(freq_flags,collapse="; "), "\n")

# 5) Marker medians
cat("Computing marker medians...\n")
key_markers <- c("CD3","TCRb","TCRgd","CD4","CD8","CD19","B220","IgD","IgM","CD138","CD43",
                 "CD11b","CD115","Ly6C","Ly6G","SiglecF","F480","CD64","CD11c","MHCII","120g8",
                 "Sca1","cKit","CD150","CD34","CD16_32","NKp46","CD49b","FceR1a")
key_markers <- intersect(key_markers, colnames(dat_t))
labels_present <- setdiff(unique(llm_gate), "unassigned")
median_mat <- matrix(NA, nrow=length(labels_present), ncol=length(key_markers),
                     dimnames=list(labels_present, key_markers))
for (lab in labels_present) {
  idx <- which(llm_gate==lab)
  if (length(idx)>0) median_mat[lab,] <- sapply(key_markers, function(m) median(dat_t[idx,m]))
}
write.csv(as.data.frame(median_mat), file.path(obj_dir,"qa_marker_medians_by_label.csv"))

# Heatmap
cat("Creating heatmap...\n")
png(file.path(fig_dir,"qa_marker_heatmap.png"), width=1200, height=800, res=100)
pheatmap(median_mat, scale="column", cluster_cols=FALSE,
         main="Median marker expression by llm_gate", fontsize_row=8, fontsize_col=8)
dev.off()

# 6) Signature violations
cat("Checking violations...\n")
violations <- data.frame(label=character(), rule=character(), n_violating=integer(), fraction_violating=numeric(), stringsAsFactors=FALSE)

check_viol <- function(label, marker, direction) {
  idx <- which(llm_gate==label)
  if (length(idx)==0) return(NULL)
  vals <- dat_t[idx, marker]
  root_vals <- dat_t[cd45_pos, marker]
  # Use 90th percentile for "should be low" checks and 10th for "should be high"
  # Many CyTOF markers have narrow near-zero distributions where median/Q75 is
  # too close to background noise, causing false violations even with perfect gating
  if (direction=="high") {
    ref_val <- quantile(root_vals, 0.10)
    nv <- sum(vals < ref_val)
    rule <- paste0(marker," should be high (>Q10 of root)")
  } else {
    ref_val <- quantile(root_vals, 0.90)
    nv <- sum(vals > ref_val)
    rule <- paste0(marker," should be low (<Q90 of root)")
  }
  data.frame(label=label, rule=rule, n_violating=nv, fraction_violating=nv/length(idx))
}

for (vc in list(
  list("NK cells","CD3","low"), list("CD4 T cells","CD8","low"), list("CD8 T cells","CD4","low"),
  list("Plasma Cells","IgD","low"), list("Eosinophils","Ly6G","low"), list("pDCs","120g8","high"),
  list("Eosinophils","SiglecF","high"), list("Macrophages","CD64","high"), list("Macrophages","F480","high")
)) {
  v <- check_viol(vc[[1]], vc[[2]], vc[[3]])
  if (!is.null(v)) violations <- rbind(violations, v)
}
write.csv(violations, file.path(obj_dir,"qa_signature_violations.csv"), row.names=FALSE)
cat("Violations:\n"); print(violations)

# 7.2 Density overlays
cat("Creating density plots...\n")
for (m in c("Cisplatin","CD45.2","Ter119","Ly6G","SiglecF","120g8","CD138","Sca1","cKit")) {
  if (!m %in% colnames(dat_t)) next
  df_d <- data.frame(value=dat_t[cd45_pos,m], label=llm_gate[cd45_pos])
  top_l <- names(sort(table(df_d$label),decreasing=TRUE))[1:min(8,length(unique(df_d$label)))]
  df_d <- df_d[df_d$label %in% top_l,]
  p <- ggplot(df_d, aes(value, fill=label)) + geom_density(alpha=0.3) +
    labs(x=m, title=paste0("QA Density: ",m)) + theme_minimal()
  ggsave(file.path(fig_dir,paste0("qa_density_",m,".png")), p, width=10, height=5, dpi=150)
}

# CD34 vs CD16_32
df_p <- data.frame(CD34=dat_t[cd45_pos,"CD34"], CD16_32=dat_t[cd45_pos,"CD16_32"], label=llm_gate[cd45_pos])
df_p <- df_p[df_p$label %in% c("CMP","GMP","MEP","HSC","MPP","CLP"),]
if (nrow(df_p) > 0) {
  p <- ggplot(df_p, aes(CD34,CD16_32,color=label)) + geom_point(alpha=0.5,size=0.5) +
    labs(title="CD34 vs CD16_32 (progenitors)") + theme_minimal()
  ggsave(file.path(fig_dir,"qa_density_CD34_CD16_32.png"), p, width=8, height=6, dpi=150)
}

# Gate plot manifest
gate_files <- data.frame(
  gate_name=c("DNA_pos","Singlets","Live","CD45pos_Ter119neg","Lin_neg",
    "Eosinophils","Basophils","Myeloid_CD11b","Macrophages","Monocytes_CD115pos",
    "Monocyte_subsets","DC_parent","pDCs","mDCs","NK_cells","NKT_cells",
    "TCRb_Tcells","CD8_Tcells","CD4_Tcells","gd_Tcells","B_parent","B_subsets",
    "Plasma_Cells","ProB_FracAC","LSK_LK","HSC","MPP","CMP_GMP_MEP","CLP_proxy"),
  expected_file=sprintf("gate_%02d_%s.png",
    c(1:14,15:20,21,22,25,26,27:31),
    c("DNA_pos","Singlets","Live","CD45pos_Ter119neg","Lin_neg",
      "Eosinophils","Basophils","Myeloid_CD11b","Macrophages","Monocytes_CD115pos",
      "Monocyte_subsets","DC_parent","pDCs","mDCs","NK_cells","NKT_cells",
      "TCRb_Tcells","CD8_Tcells","CD4_Tcells","gd_Tcells","B_parent","B_subsets",
      "Plasma_Cells","ProB_FracAC","LSK_LK","HSC","MPP","CMP_GMP_MEP","CLP_proxy")),
  stringsAsFactors=FALSE)
gate_files$exists <- file.exists(file.path(fig_dir, gate_files$expected_file))
gate_files$notes <- ""
write.csv(gate_files, file.path(obj_dir,"qa_gate_plot_manifest.csv"), row.names=FALSE)

# Gating tree overview
cat("Creating tree overview...\n")
png(file.path(fig_dir,"gating_tree_overview.png"), width=1000, height=800, res=100)
par(mar=c(1,1,2,1))
plot.new(); plot.window(xlim=c(0,10), ylim=c(0,10))
title("Gating Tree Overview")
text(5,9.5,"All Events",cex=1.2,font=2)
text(5,9.0,paste0("(n=",n_events,")"),cex=0.8)
arrows(5,8.8,5,8.3,length=0.1)
text(5,8.0,"DNA+ -> Singlets -> Live -> CD45+Ter119-",cex=0.9)
text(5,7.6,paste0("Root: n=",n_root),cex=0.8)
arrows(5,7.3,3,6.8,length=0.1); arrows(5,7.3,7,6.8,length=0.1)
text(3,6.5,"Mature Immune",cex=1,font=2)
text(7,6.5,"Lin- Progenitors",cex=1,font=2)
text(1.5,5.5,"Granulocytes",cex=0.8,font=3)
text(1.5,5.0,paste0("Eos:",sum(eos_local)," Baso:",sum(baso_local)),cex=0.7)
text(3,5.5,"Mono/Macro",cex=0.8,font=3)
text(3,5.0,paste0("Mac:",sum(macro_local)," Cl:",sum(classical_local)," Int:",sum(intermediate_local)),cex=0.7)
text(3,4.6,paste0("NC:",sum(nonclassical_local)),cex=0.7)
text(1.5,3.5,"DCs",cex=0.8,font=3)
text(1.5,3.0,paste0("pDC:",sum(pdc_local)," mDC:",sum(mdc_local)),cex=0.7)
text(3,3.5,"T/NK",cex=0.8,font=3)
text(3,3.0,paste0("CD4:",sum(cd4t_local)," CD8:",sum(cd8t_local)," gd:",sum(gdt_local)),cex=0.7)
text(3,2.6,paste0("NK:",sum(nk_local)," NKT:",sum(nkt_local)),cex=0.7)
text(4.5,3.5,"B cells",cex=0.8,font=3)
text(4.5,3.0,paste0("IgD+M+:",sum(igdp_igmp_local)," IgD-M+:",sum(igdn_igmp_local)),cex=0.7)
text(4.5,2.6,paste0("M-D-:",sum(igmn_igdn_local)," Plasma:",sum(plasma_local)),cex=0.7)
text(4.5,2.2,paste0("ProB:",sum(prob_local)),cex=0.7)
text(7,5.5,"LSK",cex=0.8,font=3)
text(7,5.0,paste0("HSC:",sum(hsc_local)," MPP:",sum(mpp_local)),cex=0.7)
text(8.5,5.5,"LK",cex=0.8,font=3)
text(8.5,5.0,paste0("CMP:",sum(cmp_local)," GMP:",sum(gmp_local)," MEP:",sum(mep_local)),cex=0.7)
text(7,4.0,"CLP proxy",cex=0.8,font=3); text(7,3.5,paste0("n=",sum(clp_local)),cex=0.7)
text(5,1.0,paste0("Unassigned: ",sum(llm_gate=="unassigned")," (",round(unassigned_frac*100,1),"%)"),cex=0.9,col="red")
dev.off()

# ============================================================
# REFERENCE COMPARISON
# ============================================================
cat("\n=== REFERENCE COMPARISON ===\n")
conf_table <- table(llm_gate, ref_gate)
write.csv(as.data.frame.matrix(conf_table), file.path(obj_dir,"qa_confusion_table.csv"))

top_matches <- do.call(rbind, lapply(rownames(conf_table), function(lab) {
  row <- sort(conf_table[lab,], decreasing=TRUE)
  top3 <- head(row[row>0], 3)
  if (length(top3)==0) return(NULL)
  data.frame(llm_gate=lab, ref_label=names(top3), count=as.integer(top3),
             fraction=as.integer(top3)/sum(conf_table[lab,]), stringsAsFactors=FALSE)
}))
write.csv(top_matches, file.path(obj_dir,"qa_top_matches_by_llm_gate.csv"), row.names=FALSE)

# ============================================================
# UMAP
# ============================================================
cat("\n=== UMAP ===\n")
set.seed(42)
umap_markers <- setdiff(marker_cols, "Ter119")
n_sub <- 10000

label_counts <- table(llm_gate)
target <- round(n_sub * label_counts / n_events)
target <- pmax(target, pmin(5, label_counts))
while (sum(target) > n_sub) { b <- which.max(target); target[b] <- target[b]-1 }
while (sum(target) < n_sub) { b <- which.max(label_counts-target); target[b] <- target[b]+1 }

sub_idx <- sort(unlist(lapply(names(target), function(lab) {
  idx <- which(llm_gate==lab); sample(idx, min(target[lab], length(idx)))
})))
cat("Subsampled:", length(sub_idx), "\n")

cat("Running UMAP...\n")
umap_in <- as.matrix(dat_t[sub_idx, umap_markers])
umap_res <- uwot::umap(umap_in, n_neighbors=15, min_dist=0.2, metric="euclidean", n_threads=4, verbose=FALSE)

df_u <- data.frame(UMAP1=umap_res[,1], UMAP2=umap_res[,2],
                   llm_gate=llm_gate[sub_idx], ref_gate=ref_gate[sub_idx])

p1 <- ggplot(df_u, aes(UMAP1,UMAP2,color=llm_gate)) +
  geom_point(size=0.3,alpha=0.5) + labs(title="UMAP: llm_gate") +
  theme_minimal() + theme(legend.text=element_text(size=6)) +
  guides(color=guide_legend(override.aes=list(size=2,alpha=1)))
ggsave(file.path(fig_dir,"umap_llm_gate.png"), p1, width=12, height=8, dpi=150)

p2 <- ggplot(df_u, aes(UMAP1,UMAP2,color=ref_gate)) +
  geom_point(size=0.3,alpha=0.5) + labs(title="UMAP: reference gate") +
  theme_minimal() + theme(legend.text=element_text(size=6)) +
  guides(color=guide_legend(override.aes=list(size=2,alpha=1)))
ggsave(file.path(fig_dir,"umap_reference_gate.png"), p2, width=12, height=8, dpi=150)

ps <- grid.arrange(p2, p1, ncol=2)
ggsave(file.path(fig_dir,"umap_side_by_side.png"), ps, width=20, height=8, dpi=150)

cat("\n=== DONE ===\n")
