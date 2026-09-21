suppressPackageStartupMessages(library(data.table))

filter_pg_by_peptide_count <- function(
    pr_path = "report.pr_matrix.tsv",
    pg_path = "report.pg_matrix.tsv",
    output_path = "report.pg_matrix.filtered_min2peptides.tsv",
    min_peptides = 2L) {
  min_peptides <- as.integer(min_peptides)
  if (length(min_peptides) != 1L || is.na(min_peptides) || min_peptides < 1L) {
    stop("min_peptides must be one positive integer")
  }

  pr <- fread(pr_path)
  pg <- fread(pg_path)

  required_pr <- c("Protein.Group", "Proteotypic", "Stripped.Sequence")
  required_pg <- c(
    "Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
    "N.Sequences", "N.Proteotypic.Sequences"
  )
  missing_pr <- setdiff(required_pr, names(pr))
  missing_pg <- setdiff(required_pg, names(pg))
  if (length(missing_pr) > 0L) {
    stop("Missing required pr_matrix columns: ", paste(missing_pr, collapse = ", "))
  }
  if (length(missing_pg) > 0L) {
    stop("Missing required pg_matrix columns: ", paste(missing_pg, collapse = ", "))
  }

  sample_cols <- setdiff(names(pg), required_pg)
  missing_samples <- setdiff(sample_cols, names(pr))
  if (length(missing_samples) > 0L) {
    stop(
      "These pg_matrix sample columns are absent from pr_matrix: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  proteotypic <- toupper(as.character(pr$Proteotypic)) %chin% c("1", "TRUE", "T")
  valid_sequence <- !is.na(pr$Stripped.Sequence) & nzchar(pr$Stripped.Sequence)
  filtered_by_sample <- integer(length(sample_cols))
  affected <- rep(FALSE, nrow(pg))

  for (j in seq_along(sample_cols)) {
    sample_col <- sample_cols[j]
    precursor_abundance <- suppressWarnings(as.numeric(pr[[sample_col]]))
    detected <- pr[
      proteotypic & valid_sequence & is.finite(precursor_abundance) & precursor_abundance > 0,
      .(N.Peptides = uniqueN(Stripped.Sequence)),
      by = Protein.Group
    ]

    peptide_count <- detected$N.Peptides[match(pg$Protein.Group, detected$Protein.Group)]
    peptide_count[is.na(peptide_count)] <- 0L

    protein_abundance <- suppressWarnings(as.numeric(pg[[sample_col]]))
    to_filter <- is.finite(protein_abundance) & peptide_count < min_peptides
    pg[[sample_col]][to_filter] <- NA_real_
    filtered_by_sample[j] <- sum(to_filter)
    affected <- affected | to_filter
  }

  fwrite(pg, output_path, sep = "\t", na = "", quote = FALSE)

  summary <- list(
    proteins = nrow(pg),
    samples = length(sample_cols),
    abundances_filtered = sum(filtered_by_sample),
    proteins_affected = sum(affected),
    min_peptides = min_peptides,
    filtered_by_sample = setNames(filtered_by_sample, sample_cols),
    output_path = normalizePath(output_path, winslash = "/", mustWork = FALSE)
  )
  invisible(summary)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(args %in% c("-h", "--help"))) {
    cat(
      "Usage: Rscript filter_pg_by_peptide_count.R [pr_matrix.tsv] [pg_matrix.tsv] [output.tsv] [min_peptides]\n",
      "Defaults: report.pr_matrix.tsv report.pg_matrix.tsv report.pg_matrix.filtered_min2peptides.tsv 2\n",
      sep = ""
    )
    return(invisible(NULL))
  }
  if (length(args) > 4L) {
    stop("Expected no more than four arguments; use --help for usage")
  }

  defaults <- c(
    "report.pr_matrix.tsv",
    "report.pg_matrix.tsv",
    "report.pg_matrix.filtered_min2peptides.tsv",
    "2"
  )
  defaults[seq_along(args)] <- args

  summary <- filter_pg_by_peptide_count(
    pr_path = defaults[1],
    pg_path = defaults[2],
    output_path = defaults[3],
    min_peptides = defaults[4]
  )
  cat(sprintf(
    "Wrote %s\nFiltered %d protein-sample abundances across %d proteins (%d samples; minimum %d unique proteotypic peptides).\n",
    summary$output_path,
    summary$abundances_filtered,
    summary$proteins_affected,
    summary$samples,
    summary$min_peptides
  ))
}

if (sys.nframe() == 0L) {
  main()
}
