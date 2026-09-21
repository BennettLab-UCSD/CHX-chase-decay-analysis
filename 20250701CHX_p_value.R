# step1

library(data.table)
library(reshape2)
library(stats)
library(openxlsx)

source("decay_analysis_helpers.R")

# ============================
# 1. Load and filter the data
# ============================

pg_matrix_path <- "report.pg_matrix.filtered_min2peptides.tsv"
if (!file.exists(pg_matrix_path)) {
  stop("Filtered pg_matrix not found. Run: Rscript filter_pg_by_peptide_count.R")
}
df <- fread(pg_matrix_path)

# Keep only unique protein groups (remove entries with multiple proteins)
df <- df[!grepl(";", Protein.Group), ]

# Remove cRAP contaminants before model fitting and BH FDR correction.
contaminant_rows <- is_contaminant_protein_group(df$Protein.Group)
message(sprintf(
  "Excluding %d cRAP contaminant protein groups before FDR calculation.",
  sum(contaminant_rows)
))
df <- df[!contaminant_rows, ]

# Keep proteins with >= 2 proteotypic peptides
df <- df[N.Proteotypic.Sequences >= 2, ]

# Drop columns 5 and 6 (as in your original script)
df <- df[, -c(5, 6)]

# WT and KO columns (based on your file structure AFTER removing col 5–6)
WT_cols <- 5:22 # 18 columns = 6 time points × 3 replicates
KO_cols <- 23:40 # next 18 columns = KO samples

# Time points and replicate indices (fixed)
time_points <- c(0, 1, 2, 4, 6, 8)
rep_index <- list(1:3, 4:6, 7:9, 10:12, 13:15, 16:18)


# ======================================================
# 2. Function to compute half-life and slope statistics
# ======================================================

calc_one <- function(wt_vals, ko_vals, gene, pg_id) {
  # --- (1) Compute mean intensity at each time point ---
  # Exclude the whole condition/time point when fewer than two valid
  # replicates remain after sample-level peptide filtering.
  WT_summary <- summarize_timepoints(wt_vals, rep_index, min_replicates = 2L)
  KO_summary <- summarize_timepoints(ko_vals, rep_index, min_replicates = 2L)
  WT_means <- WT_summary$mean
  KO_means <- KO_summary$mean

  # Calculate 0h abundance fold change (KO / WT)
  wt_0h <- WT_means[1]
  ko_0h <- KO_means[1]
  fc_0h <- ko_0h / wt_0h

  # Add a small epsilon to avoid log2(0)
  eps <- if (any(c(WT_means, KO_means) <= 0, na.rm = TRUE)) 1e-6 else 0

  # Log2-transformed means
  log2_WT <- log2(WT_means + eps)
  log2_KO <- log2(KO_means + eps)

  ok_wt <- is.finite(log2_WT)
  ok_ko <- is.finite(log2_KO)

  # Need at least 2 valid time points for a regression line
  if (sum(ok_wt) < 2 | sum(ok_ko) < 2) {
    return(data.table(
      Protein.Group = pg_id, Gene = gene,
      slope_WT = NA_real_, slope_KO = NA_real_,
      t12_WT = NA_real_, t12_KO = NA_real_,
      t12_WT_low = NA_real_, t12_WT_high = NA_real_,
      t12_KO_low = NA_real_, t12_KO_high = NA_real_,
      p_slope = NA_real_,
      FC_KO_WT_0h = fc_0h
    ))
  }

  # --- (2) Linear regression in log2 space: log2(I) = a + b * time ---
  fit_wt <- lm(log2_WT[ok_wt] ~ time_points[ok_wt])
  fit_ko <- lm(log2_KO[ok_ko] ~ time_points[ok_ko])

  b_wt <- coef(fit_wt)[2] # WT slope
  b_ko <- coef(fit_ko)[2] # KO slope

  # Half-life: t1/2 = -1 / slope (valid only if slope < 0)
  t12_wt <- ifelse(b_wt < 0, -1 / b_wt, Inf)
  t12_ko <- ifelse(b_ko < 0, -1 / b_ko, Inf)

  # --- (3) Confidence intervals for slopes and half-lives (robust) ---

  # Safe wrapper for slope CI: returns c(NA, NA) if CI cannot be computed
  get_slope_ci <- function(fit) {
    # If there is no residual degree of freedom, CI cannot be estimated
    if (fit$df.residual <= 0) {
      return(c(NA_real_, NA_real_))
    }
    ci <- try(confint(fit)[2, ], silent = TRUE)
    if (inherits(ci, "try-error") || any(!is.finite(ci))) {
      return(c(NA_real_, NA_real_))
    }
    ci
  }

  ci_wt <- get_slope_ci(fit_wt)
  ci_ko <- get_slope_ci(fit_ko)

  # Convert slope CI → half-life CI
  hl_ci_from_slope <- function(ci) {
    # If CI is not available, return NA for both bounds
    if (any(!is.finite(ci))) {
      return(c(NA_real_, NA_real_))
    }
    # Proper decay: both bounds < 0
    if (all(ci < 0)) {
      return(sort(-1 / ci))
    }
    # Signal clearly increasing: both bounds > 0
    if (all(ci > 0)) {
      return(c(Inf, Inf))
    }
    # CI crosses zero: provide finite lower bound and infinite upper bound
    c(-1 / ci[1], Inf)
  }

  t12_wt_ci <- hl_ci_from_slope(ci_wt)
  t12_ko_ci <- hl_ci_from_slope(ci_ko)

  # --- (4) Interaction model using all replicates ---
  # This tests whether WT slope differs from KO slope

  df_wt_raw <- build_raw_decay_data(
    wt_vals, "WT", time_points, rep_index, min_replicates = 2L
  )
  df_ko_raw <- build_raw_decay_data(
    ko_vals, "KO", time_points, rep_index, min_replicates = 2L
  )
  df_raw <- rbind(df_wt_raw, df_ko_raw)
  df_raw$log2I <- log2(df_raw$Intensity + eps)
  df_raw$Group <- factor(df_raw$Group, levels = c("WT", "KO"))

  # Interaction model: log2(I) ~ Time * Group
  fit_int <- try(lm(log2I ~ Time * Group, data = df_raw), silent = TRUE)

  if (inherits(fit_int, "try-error")) {
    p_slope <- NA_real_
  } else {
    sm <- summary(fit_int)$coefficients
    row_int <- grep("Time:Group", rownames(sm))
    if (length(row_int) == 1) {
      p_slope <- sm[row_int, "Pr(>|t|)"]
    } else {
      p_slope <- NA_real_
    }
  }

  # Return results for this protein
  data.table(
    Protein.Group = pg_id, Gene = gene,
    slope_WT = b_wt, slope_KO = b_ko,
    t12_WT = t12_wt,
    t12_WT_low = t12_wt_ci[1], t12_WT_high = t12_wt_ci[2],
    t12_KO = t12_ko,
    t12_KO_low = t12_ko_ci[1], t12_KO_high = t12_ko_ci[2],
    p_slope = p_slope,
    FC_KO_WT_0h = fc_0h
  )
}


# ===============================================
# 3. Loop through all proteins and combine results
# ===============================================

results <- rbindlist(
  lapply(1:nrow(df), function(i) {
    wt_vals <- as.numeric(df[i, ..WT_cols])
    ko_vals <- as.numeric(df[i, ..KO_cols])

    calc_one(
      wt_vals, ko_vals,
      gene = df$Genes[i],
      pg_id = df$Protein.Group[i]
    )
  }),
  use.names = TRUE, fill = TRUE
)

# Multiple testing correction (Benjamini–Hochberg)
results[, FDR := p.adjust(p_slope, method = "BH")]

# Save output to Excel and remove old CSV if exists
write.xlsx(results, "HalfLife_AllProteins.xlsx", overwrite = TRUE)
if (file.exists("HalfLife_AllProteins.csv")) {
  file.remove("HalfLife_AllProteins.csv")
}

results
