suppressPackageStartupMessages({
  library(data.table) # for fread and data.table []
  library(ggplot2)
  library(grid) # for unit()
  library(openxlsx) # for reading the all-protein BH FDR results
})

source("decay_analysis_helpers.R")


# Use the protein-group matrix after sample-specific peptide-count filtering.
pg_matrix_path <- "report.pg_matrix.filtered_min2peptides.tsv"
if (!file.exists(pg_matrix_path)) {
  stop("Filtered pg_matrix not found. Run: Rscript filter_pg_by_peptide_count.R")
}

# normalization

df <- fread(pg_matrix_path)
df <- df[, c(1:4, 7:24, 25:42), with = FALSE]

#+SE

# ---- Cell-style theme: clean grid with clear axes and ticks ----
theme_cell <- function(base_size = 9, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major.x = element_line(linewidth = 0.5, colour = "grey80"),
      panel.grid.minor.x = element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.major.y = element_line(linewidth = 0.5, colour = "grey80"),
      panel.grid.minor.y = element_line(linewidth = 0.25, colour = "grey90"),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(2.5, "pt"),
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(t = 2, b = 6))
    )
}

# Time points and replicate indices (3 replicates per time point)
time_points <- c(0, 1, 2, 4, 6, 8)
rep_index <- list(1:3, 4:6, 7:9, 10:12, 13:15, 16:18)

# ---- Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Help message
if ("--help" %in% args || "-h" %in% args) {
  cat("Usage: Rscript 20250701CHX_decay_plot.R [gene_name] [--global]\n")
  cat("Arguments:\n")
  cat("  gene_name    The gene symbol to plot (default: IP6K2)\n")
  cat("  --global     Optionally generate the global median decay plot (default: FALSE)\n")
  cat("Output:\n")
  cat("  2026Figure_CHX_<gene_name>.svg (88 x 70 mm)\n")
  cat("Prerequisite:\n")
  cat("  Run 20250701CHX_p_value.R first to create HalfLife_AllProteins.xlsx.\n")
  q(save = "no")
}

# Defaults
gene_name <- "IP6K2"
run_global <- FALSE

# Check for --global or -g
if ("--global" %in% args || "-g" %in% args) {
  run_global <- TRUE
  args <- args[!args %in% c("--global", "-g")]
}

# The remaining first argument is the gene name
if (length(args) > 0) {
  gene_name <- args[1]
}

# Column 3 contains gene names
row_index <- which(df[[3]] == gene_name)
if (length(row_index) == 0) {
  stop(sprintf("Gene '%s' not found in dataset", gene_name))
}

otu_row <- df[row_index, ]
protein_group <- otu_row$Protein.Group[1]

# WT: columns 5–22; KO: columns 23–40
WT_vals <- as.numeric(otu_row[, 5:22])
KO_vals <- as.numeric(otu_row[, 23:40])

# A condition/time point must have at least two valid replicates. A lone
# replicate is excluded from the plotted mean, half-life fit, and slope test.
WT_summary <- summarize_timepoints(WT_vals, rep_index, min_replicates = 2L)
KO_summary <- summarize_timepoints(KO_vals, rep_index, min_replicates = 2L)
WT_means <- WT_summary$mean
KO_means <- KO_summary$mean
WT_sds <- WT_summary$sd
KO_sds <- KO_summary$sd

# Optional: calculate standard error instead of standard deviation
# WT_ses <- WT_sds / sqrt(3)
# KO_ses <- KO_sds / sqrt(3)

## ---- Normalization for plotting only (each group to its own 0h mean) ----
# Find which time point corresponds to 0 h
time0_pos <- which(time_points == 0)
if (length(time0_pos) != 1) stop("time_points must contain exactly one '0'")

# Indices of 0h replicates within each group (using your rep_index)
time0_idxs <- rep_index[[time0_pos]]

# Reference = each group's own 0h mean across replicates
ref_wt_0h_mean <- mean(WT_vals[time0_idxs], na.rm = TRUE)
ref_ko_0h_mean <- mean(KO_vals[time0_idxs], na.rm = TRUE)
if (!is.finite(ref_wt_0h_mean) || ref_wt_0h_mean <= 0) stop("Invalid WT 0h mean for normalization")
if (!is.finite(ref_ko_0h_mean) || ref_ko_0h_mean <= 0) stop("Invalid KO 0h mean for normalization")

# Scale each group's means and SDs by its own 0h reference
WT_means_norm <- WT_means / ref_wt_0h_mean
KO_means_norm <- KO_means / ref_ko_0h_mean
WT_sds_norm <- WT_sds / ref_wt_0h_mean
KO_sds_norm <- KO_sds / ref_ko_0h_mean

# Build plotting data with normalized intensities

df_plot_norm <- data.frame(
  Time = rep(time_points, 2),
  Intensity = c(WT_means_norm, KO_means_norm),
  SD = c(WT_sds_norm, KO_sds_norm),
  Group = factor(rep(c("WT", "KO"), each = length(time_points)),
    levels = c("WT", "KO")
  )
)

# ---- Half-life fitting based on group means; add a small pseudocount only for log2 transformation when values ≤ 0 ----
eps <- if (any(WT_means <= 0, na.rm = TRUE) || any(KO_means <= 0, na.rm = TRUE)) 1e-6 else 0
log2_WT <- log2(WT_means + eps)
log2_KO <- log2(KO_means + eps)

ok_wt <- is.finite(log2_WT) & is.finite(time_points)
ok_ko <- is.finite(log2_KO) & is.finite(time_points)
if (sum(ok_wt) < 2 || sum(ok_ko) < 2) stop("Not enough finite points for linear fit")

fit_wt <- lm(log2_WT[ok_wt] ~ time_points[ok_wt])
fit_ko <- lm(log2_KO[ok_ko] ~ time_points[ok_ko])

slope_wt <- coef(fit_wt)[2]
slope_ko <- coef(fit_ko)[2]
halflife_wt <- ifelse(slope_wt < 0, log2(2) / -slope_wt, Inf)
halflife_ko <- ifelse(slope_ko < 0, log2(2) / -slope_ko, Inf)

# ---- Generate fitted exponential curves, transform them back to the linear scale, and normalize to each group's own 0h ----
pred_t <- seq(min(time_points), max(time_points), length.out = 100)
pred_wt <- 2^(coef(fit_wt)[1] + coef(fit_wt)[2] * pred_t) / ref_wt_0h_mean
pred_ko <- 2^(coef(fit_ko)[1] + coef(fit_ko)[2] * pred_t) / ref_ko_0h_mean
df_fit <- data.frame(
  Time = c(pred_t, pred_t),
  Intensity = c(pred_wt, pred_ko),
  Group = factor(rep(c("WT", "KO"), each = length(pred_t)),
    levels = c("WT", "KO")
  )
)


# ----- Grid line breaks (major = 2h / 0.5, minor = 1h / 0.25) -----
x_major <- seq(min(time_points), max(time_points), by = 2)
x_minor <- setdiff(seq(min(time_points), max(time_points), by = 1), x_major)

# Generate Y-axis major/minor breaks
# ==== Adaptive Y-axis: symmetric around 1, step sizes scale with data range ====

# First determine the Y-axis limits from the error bars; if all values are NA, fall back to the mean values
y_low_cand <- df_plot_norm$Intensity - df_plot_norm$SD
y_high_cand <- df_plot_norm$Intensity + df_plot_norm$SD
if (all(!is.finite(y_low_cand)) || all(!is.finite(y_high_cand))) {
  y_low <- min(df_plot_norm$Intensity, na.rm = TRUE)
  y_high <- max(df_plot_norm$Intensity, na.rm = TRUE)
} else {
  y_low <- min(y_low_cand, na.rm = TRUE)
  y_high <- max(y_high_cand, na.rm = TRUE)
}

# Determine plotting range from actual data extent with 10% padding
pad <- (y_high - y_low) * 0.10
pad <- max(pad, 0.02) # minimum padding to avoid degenerate axis

ylim <- c(max(0, y_low - pad), y_high + pad)

# Choose step sizes adaptively based on the range
range_size <- ylim[2] - ylim[1]
y_step_major <- if (range_size > 1.5) 0.5 else if (range_size > 0.6) 0.25 else if (range_size > 0.3) 0.1 else 0.05
y_step_minor <- y_step_major / 2

# Major and minor tick marks (within the plotting range only)
y_major <- seq(floor(ylim[1] / y_step_major) * y_step_major,
  ceiling(ylim[2] / y_step_major) * y_step_major,
  by = y_step_major
)
y_minor <- setdiff(
  seq(floor(ylim[1] / y_step_minor) * y_step_minor,
    ceiling(ylim[2] / y_step_minor) * y_step_minor,
    by = y_step_minor
  ),
  y_major
)

# ----- Optional: use SE instead of SD for error bars by setting use_se <- TRUE -----
use_se <- TRUE
if (use_se) {
  n_wt <- WT_summary$n
  n_ko <- KO_summary$n
  WT_se <- WT_sds / sqrt(pmax(n_wt, 1))
  KO_se <- KO_sds / sqrt(pmax(n_ko, 1))
  df_plot_norm$SD[df_plot_norm$Group == "WT"] <- WT_se / ref_wt_0h_mean
  df_plot_norm$SD[df_plot_norm$Group == "KO"] <- KO_se / ref_ko_0h_mean
}

# ----- Raw replicate measurements (normalized using the same reference)-----
# build_raw_decay_data omits the whole condition/time point when fewer than
# two valid replicates remain.
df_wt_raw <- build_raw_decay_data(
  WT_vals, "WT", time_points, rep_index, min_replicates = 2L
)
df_ko_raw <- build_raw_decay_data(
  KO_vals, "KO", time_points, rep_index, min_replicates = 2L
)
df_wt_raw$Intensity <- df_wt_raw$Intensity / ref_wt_0h_mean
df_ko_raw$Intensity <- df_ko_raw$Intensity / ref_ko_0h_mean
df_raw <- rbind(df_wt_raw, df_ko_raw)
df_raw$Group <- factor(df_raw$Group, levels = c("WT", "KO"))

# ---- Calculate the 95% confidence intervals for the half-life (by transforming the slope CI using t½ = -1/slope) ----
# Confidence level (adjustable)
alpha <- 0.05 # 95% CI

get_slope_ci <- function(fit, level = 0.95) {
  if (fit$df.residual <= 0) {
    return(c(NA_real_, NA_real_))
  }
  ci <- suppressWarnings(try(confint(fit, level = level)[2, ], silent = TRUE))
  if (inherits(ci, "try-error") || any(!is.finite(ci))) {
    return(c(NA_real_, NA_real_))
  }
  ci
}

ci_wt <- get_slope_ci(fit_wt, level = 1 - alpha)
ci_ko <- get_slope_ci(fit_ko, level = 1 - alpha)

half_life_from_slope <- function(b) ifelse(b < 0, -1 / b, Inf)

hl_ci_from_slope_ci <- function(ci) {
  # ci: c(lower slope CI, upper slope CI)
  if (any(!is.finite(ci))) {
    c(NA_real_, NA_real_)
  } else if (all(ci < 0)) { # Negative slope (decay)
    sort(-1 / ci) # Finite two-sided confidence interval
  } else if (all(ci > 0)) { # Positive slope (accumulation)
    c(Inf, Inf) # Infinite half-life
  } else { # Confidence interval crosses zero
    c(-1 / ci[1], Inf) # Finite lower bound and infinite upper bound
  }
}

halflife_wt <- half_life_from_slope(coef(fit_wt)[2])
halflife_ko <- half_life_from_slope(coef(fit_ko)[2])

t12_wt_CI <- hl_ci_from_slope_ci(ci_wt)
t12_ko_CI <- hl_ci_from_slope_ci(ci_ko)

# ---- Test for differences in decay slopes between groups (interaction term in the log2 linear model) ----
eps <- if (any(c(WT_vals, KO_vals) <= 0, na.rm = TRUE)) 1e-6 else 0
df_raw_log2 <- transform(df_raw, log2I = log2(Intensity + eps))
fit_int <- lm(log2I ~ Time * Group, data = df_raw_log2)

coefs <- summary(fit_int)$coefficients
int_row <- grep("^Time:Group", rownames(coefs)) # Automatically identify the interaction term
p_int <- if (length(int_row) == 1) coefs[int_row, "Pr(>|t|)"] else NA_real_

# The displayed FDR must come from the all-protein Benjamini-Hochberg
# correction; it cannot be calculated from a single plotted protein.
half_life_results_path <- "HalfLife_AllProteins.xlsx"
if (!file.exists(half_life_results_path)) {
  stop(
    "HalfLife_AllProteins.xlsx not found. Run: ",
    "Rscript 20250701CHX_p_value.R"
  )
}

half_life_results <- read.xlsx(half_life_results_path)
required_result_cols <- c("Protein.Group", "Gene", "FDR")
missing_result_cols <- setdiff(required_result_cols, names(half_life_results))
if (length(missing_result_cols) > 0) {
  stop(
    "HalfLife_AllProteins.xlsx is missing required column(s): ",
    paste(missing_result_cols, collapse = ", ")
  )
}

fdr_match <- half_life_results[
  half_life_results$Protein.Group == protein_group &
    half_life_results$Gene == gene_name,
  , drop = FALSE
]
if (nrow(fdr_match) != 1) {
  stop(sprintf(
    "Expected one FDR result for gene '%s' / protein group '%s'; found %d. Rerun 20250701CHX_p_value.R.",
    gene_name, protein_group, nrow(fdr_match)
  ))
}
fdr_value <- as.numeric(fdr_match$FDR[1])

# Significance marker for the trajectory-wide slope-interaction FDR. The
# bracket spans the time course 
significance_label <- if (is.na(fdr_value)) {
  "NA"
} else if (fdr_value < 0.0001) {
  "****"
} else if (fdr_value < 0.001) {
  "***"
} else if (fdr_value < 0.01) {
  "**"
} else if (fdr_value < 0.05) {
  "*"
} else {
  "ns"
}

data_y_range <- diff(ylim)
bracket_y <- ylim[2] + 0.09 * data_y_range
bracket_tick_y <- bracket_y - 0.035 * data_y_range
bracket_label_y <- bracket_y + 0.025 * data_y_range
ylim[2] <- bracket_label_y + 0.06 * data_y_range
bracket_xmin <- min(time_points) + 0.35
bracket_xmax <- max(time_points) - 0.35
bracket_xmid <- mean(c(bracket_xmin, bracket_xmax))

# ---- Add WT/KO labels at the end of the fitted curves (no legend required) ----
end_x <- max(time_points)

end_wt <- df_fit$Intensity[df_fit$Group == "WT" & df_fit$Time == end_x]
end_ko <- df_fit$Intensity[df_fit$Group == "KO" & df_fit$Time == end_x]

end_wt <- end_wt[1]
end_ko <- end_ko[1]

# update

# ---- Error bar data (mean ± SE) ----
cap_w <- 0.35
lab_dx <- 0.25

df_lr <- transform(df_plot_norm,
  ymin = Intensity - SD,
  ymax = Intensity + SD
)

df_lr <- subset(df_lr, is.finite(ymin) & is.finite(ymax) & is.finite(Intensity))

# FULL caps (do NOT truncate with pmax/pmin)
df_cap <- transform(df_lr,
  x1 = Time - cap_w / 2,
  x2 = Time + cap_w / 2
)

# ---- End labels position (make sure finite) ----
end_x <- max(time_points)
end_wt <- tail(df_fit$Intensity[df_fit$Group == "WT"], 1)
end_ko <- tail(df_fit$Intensity[df_fit$Group == "KO"], 1)

# ==========================
# PLOT
# ==========================
# Keep only the two half-lives in the subtitle. Statistical significance is
# communicated by the trajectory-wide bracket and its marker in the panel.
subtitle_txt <- sprintf(
  "Parental t1/2 = %.2f h; KO t1/2 = %.2f h",
  halflife_wt, halflife_ko
)


p <- ggplot(df_plot_norm, aes(x = Time, y = Intensity, color = Group)) +
  geom_hline(
    yintercept = 1, linetype = "dotted",
    linewidth = 0.45, colour = "grey55"
  ) +
  # Bracket denotes the trajectory-wide Parental-versus-KO slope test.
  annotate("segment",
    x = bracket_xmin, xend = bracket_xmax,
    y = bracket_y, yend = bracket_y,
    linewidth = 0.45, colour = "black"
  ) +
  annotate("segment",
    x = bracket_xmin, xend = bracket_xmin,
    y = bracket_y, yend = bracket_tick_y,
    linewidth = 0.45, colour = "black"
  ) +
  annotate("segment",
    x = bracket_xmax, xend = bracket_xmax,
    y = bracket_y, yend = bracket_tick_y,
    linewidth = 0.45, colour = "black"
  ) +
  annotate("text",
    x = bracket_xmid, y = bracket_label_y,
    label = significance_label,
    size = 3.1, fontface = "bold", colour = "black"
  ) +
  geom_point(size = 1.0) +
  geom_linerange(
    data = df_lr,
    aes(x = Time, ymin = ymin, ymax = ymax, color = Group),
    linewidth = 0.38, inherit.aes = FALSE
  ) +
  geom_segment(
    data = df_cap,
    aes(x = x1, xend = x2, y = ymax, yend = ymax, color = Group),
    linewidth = 0.38, inherit.aes = FALSE
  ) +
  geom_segment(
    data = df_cap,
    aes(x = x1, xend = x2, y = ymin, yend = ymin, color = Group),
    linewidth = 0.38, inherit.aes = FALSE
  ) +
  geom_line(data = df_fit, linewidth = 0.6, linetype = "dashed", alpha = 0.5) +

  # Inline labels: text shows "Parental" but the data group is still WT
  annotate("text",
    x = end_x - lab_dx, y = end_wt,
    label = "Parental",
    hjust = 0, vjust = -0.6,
    size = 3.1, fontface = "bold",
    colour = "#377EB8"
  ) +
  annotate("text",
    x = end_x - lab_dx, y = end_ko,
    label = "KO",
    hjust = 0, vjust = -0.6,
    size = 3.1, fontface = "bold",
    colour = "#E41A1C"
  ) +
  labs(
    title = NULL,
    subtitle = subtitle_txt,
    x = "Time (h)",
    y = paste(gene_name, "relative intensity")
  ) +

  # Color mapping stays WT/KO
  scale_color_manual(values = c(
    "WT" = "#377EB8",
    "KO" = "#E41A1C"
  )) +
  scale_x_continuous(
    breaks = x_major,
    minor_breaks = x_minor,
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_y_continuous(
    breaks = y_major,
    minor_breaks = y_minor
  ) +
  coord_cartesian(
    xlim = c(0, max(time_points)),
    ylim = ylim,
    clip = "off"
  ) +
  theme_cell(base_size = 9) +
  theme(plot.margin = margin(5.5, 30, 5.5, 5.5))

ggsave(sprintf("2026Figure_CHX_%s.svg", gene_name),
  plot = p,
  device = "svg",
  width = 88, height = 70, units = "mm"
)

if (run_global) {
  message("Generating global change as control plot...")

  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
  })

# -----------------------------
# 1) Load matrix and select columns (same as your workflow)
# -----------------------------
df <- fread(pg_matrix_path)
df <- df[, c(1:4, 7:24, 25:42), with = FALSE] # keep your exact selection

# Time points and replicate indices (6 time points × 3 reps = 18 values per group)
time_points <- c(0, 1, 2, 4, 6, 8)
rep_index <- list(1:3, 4:6, 7:9, 10:12, 13:15, 16:18)

# WT columns and KO columns in your df (after subsetting above)
wt_cols <- 5:22
ko_cols <- 23:40

# -----------------------------
# 2) Helper: compute per-protein timepoint means from 18 values
# -----------------------------
tp_means_from_18 <- function(x18, rep_index) {
  summarize_timepoints(x18, rep_index, min_replicates = 2L)$mean
}

# -----------------------------
# 3) Build per-protein WT/KO timepoint means (N proteins × 6 time points)
# -----------------------------
WT_mat <- as.matrix(df[, ..wt_cols])
KO_mat <- as.matrix(df[, ..ko_cols])

WT_means_mat <- t(apply(WT_mat, 1, tp_means_from_18, rep_index = rep_index))
KO_means_mat <- t(apply(KO_mat, 1, tp_means_from_18, rep_index = rep_index))

colnames(WT_means_mat) <- paste0("t", time_points)
colnames(KO_means_mat) <- paste0("t", time_points)

# -----------------------------
# 4) Global normalization: each group normalized to its OWN 0h (recommended for "global kinetics")
#    This focuses on decay dynamics rather than baseline abundance differences.
# -----------------------------
t0_col <- which(time_points == 0)
if (length(t0_col) != 1) stop("time_points must contain exactly one 0")

WT_ref <- WT_means_mat[, t0_col]
KO_ref <- KO_means_mat[, t0_col]

# Keep proteins with valid positive 0h reference in each group
keep_wt_ref <- is.finite(WT_ref) & WT_ref > 0
keep_ko_ref <- is.finite(KO_ref) & KO_ref > 0
keep_ref <- keep_wt_ref & keep_ko_ref

WT_means_mat <- WT_means_mat[keep_ref, , drop = FALSE]
KO_means_mat <- KO_means_mat[keep_ref, , drop = FALSE]
WT_ref <- WT_ref[keep_ref]
KO_ref <- KO_ref[keep_ref]

WT_norm <- WT_means_mat / WT_ref
KO_norm <- KO_means_mat / KO_ref

# Optional: require enough finite time points per protein (prevents late-time missingness bias)
min_tp <- 6 # at least 6 time points present per group
keep_tp <- (rowSums(is.finite(WT_norm)) >= min_tp) & (rowSums(is.finite(KO_norm)) >= min_tp)

WT_norm <- WT_norm[keep_tp, , drop = FALSE]
KO_norm <- KO_norm[keep_tp, , drop = FALSE]

# -----------------------------
# 5) Global summary across proteins at each time point (median is robust)
# -----------------------------
wt_med <- apply(WT_norm, 2, median, na.rm = TRUE)
ko_med <- apply(KO_norm, 2, median, na.rm = TRUE)

df_global <- data.frame(
  Time      = rep(time_points, 2),
  Intensity = c(wt_med, ko_med),
  Group     = factor(rep(c("WT", "KO"), each = length(time_points)), levels = c("WT", "KO"))
)

# -----------------------------
# 6) (Optional) bootstrap 95% CI across proteins for the median curve
# -----------------------------
bootstrap_median_ci <- function(mat, B = 500, seed = 1) {
  set.seed(seed)
  n <- nrow(mat)
  out <- matrix(NA_real_, nrow = B, ncol = ncol(mat))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    out[b, ] <- apply(mat[idx, , drop = FALSE], 2, median, na.rm = TRUE)
  }
  lo <- apply(out, 2, quantile, probs = 0.025, na.rm = TRUE)
  hi <- apply(out, 2, quantile, probs = 0.975, na.rm = TRUE)
  data.frame(lo = lo, hi = hi)
}

ci_wt <- bootstrap_median_ci(WT_norm, B = 500, seed = 1)
ci_ko <- bootstrap_median_ci(KO_norm, B = 500, seed = 2)

df_ci <- rbind(
  data.frame(Time = time_points, lo = ci_wt$lo, hi = ci_wt$hi, Group = "WT"),
  data.frame(Time = time_points, lo = ci_ko$lo, hi = ci_ko$hi, Group = "KO")
)
df_ci$Group <- factor(df_ci$Group, levels = c("WT", "KO"))

# -----------------------------
# 7) Cell-like theme 
# -----------------------------
theme_cell <- function(base_size = 9, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey85"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.6),
      axis.ticks = element_line(colour = "black", linewidth = 0.6),
      axis.ticks.length = unit(2.5, "pt"),
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(t = 2, b = 6))
    )
}

# Grid breaks
x_major <- seq(min(time_points), max(time_points), by = 2)

# -----------------------------
# 8) Plot: global median decay curve (with optional CI ribbon)
# -----------------------------
p_global <- ggplot(df_global, aes(Time, Intensity, color = Group)) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.45, colour = "grey55") +

  # Optional CI ribbon (comment out if you want an even cleaner plot)
  geom_ribbon(
    data = df_ci,
    aes(x = Time, ymin = lo, ymax = hi, fill = Group),
    alpha = 0.18, color = NA, inherit.aes = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  annotate("text",
    x = max(time_points) - 0.2, y = df_global$Intensity[df_global$Group == "WT"][length(time_points)],
    label = "Parental", hjust = 0, vjust = -0.6, size = 3.1, fontface = "bold", colour = "#377EB8"
  ) +
  annotate("text",
    x = max(time_points) - 0.2, y = df_global$Intensity[df_global$Group == "KO"][length(time_points)],
    label = "KO", hjust = 0, vjust = -0.6, size = 3.1, fontface = "bold", colour = "#E41A1C"
  ) +
  labs(
    title = "Global CHX decay (median across proteins)",
    subtitle = sprintf("Normalization: each group normalized to its own 0h (proteins kept: %d)", nrow(WT_norm)),
    x = "Time (h)",
    y = "Median relative intensity"
  ) +
  scale_color_manual(values = c("WT" = "#377EB8", "KO" = "#E41A1C")) +
  scale_fill_manual(values = c("WT" = "#377EB8", "KO" = "#E41A1C")) +
  scale_x_continuous(breaks = x_major, expand = expansion(mult = c(0, 0.06))) +
  coord_cartesian(
    xlim = c(0, max(time_points)),
    ylim = (function() {
      # Dynamic y-axis: use CI ribbon bounds and median values
      y_all <- c(df_global$Intensity, df_ci$lo, df_ci$hi)
      y_all <- y_all[is.finite(y_all)]
      y_lo <- min(y_all)
      y_hi <- max(y_all)
      pad <- (y_hi - y_lo) * 0.10 # 10% padding
      # Round to nearest 0.05 for clean axis breaks
      c(
        floor((y_lo - pad) / 0.05) * 0.05,
        ceiling((y_hi + pad) / 0.05) * 0.05
      )
    })(),
    clip = "off"
  ) +
  theme_cell(base_size = 9) +
  theme(plot.margin = margin(5.5, 16, 5.5, 5.5))

  ggsave("Figure_CHX_GlobalMedian.svg",
    plot = p_global,
    device = "svg",
    width = 85, height = 65, units = "mm"
  )
}
