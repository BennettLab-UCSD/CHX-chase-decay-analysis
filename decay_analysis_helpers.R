is_contaminant_protein_group <- function(protein_groups) {
  !is.na(protein_groups) & startsWith(as.character(protein_groups), "cRAP-")
}

summarize_timepoints <- function(values, rep_index, min_replicates = 2L) {
  min_replicates <- as.integer(min_replicates)
  if (length(min_replicates) != 1L || is.na(min_replicates) || min_replicates < 1L) {
    stop("min_replicates must be one positive integer")
  }
  if (length(values) != sum(lengths(rep_index))) {
    stop("values length must equal the total number of replicate indices")
  }

  n <- vapply(
    rep_index,
    function(idxs) sum(is.finite(values[idxs])),
    integer(1)
  )
  valid <- n >= min_replicates
  means <- rep(NA_real_, length(rep_index))
  sds <- rep(NA_real_, length(rep_index))

  for (i in which(valid)) {
    observed <- values[rep_index[[i]]]
    observed <- observed[is.finite(observed)]
    means[i] <- mean(observed)
    sds[i] <- sd(observed)
  }

  list(mean = means, sd = sds, n = n, valid = valid)
}

build_raw_decay_data <- function(
    values,
    group,
    time_points,
    rep_index,
    min_replicates = 2L) {
  if (length(time_points) != length(rep_index)) {
    stop("time_points and rep_index must have the same length")
  }

  summary <- summarize_timepoints(values, rep_index, min_replicates)
  rows <- lapply(which(summary$valid), function(i) {
    idxs <- rep_index[[i]]
    keep <- is.finite(values[idxs])
    data.frame(
      Time = rep(time_points[i], sum(keep)),
      Intensity = values[idxs][keep],
      Group = rep(group, sum(keep)),
      stringsAsFactors = FALSE
    )
  })

  if (length(rows) == 0L) {
    return(data.frame(
      Time = numeric(),
      Intensity = numeric(),
      Group = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}
