# Manual reference-comparison script for SDT implementations
# Run interactively when optional reference material is available locally.

devtools::load_all(".")

compare_to_bhsdtr2 <- function() {
  if (!requireNamespace("bhsdtr2", quietly = TRUE)) {
    message("bhsdtr2 not installed; skipping cross-package checks.")
    return(invisible(NULL))
  }

  dat <- rsdt_rating(
    n_per_cell = 200, n_subjects = 1,
    dprime = 1.4, criterion = 0.1,
    sdratio = 1.2, n_ratings = 6, spacing = 0.4
  )
  counts <- as.integer(dat[1, paste0("r", 1:6)])
  thresholds <- .sdt_make_thresholds(
    criterion = 0.1, n_ratings = 6,
    threshold_type = "parsimonious", spacing = 0.4
  )

  local_density <- dsdt_rating(
    counts = counts, stimulus = 1,
    dprime = 1.4, thresholds = thresholds,
    sdratio = 1.2, log = TRUE
  )

  message("Local rating log-density: ", signif(local_density, 6))
  message("Add package-specific bhsdtr2 calls here once the reference setup is available.")
  invisible(local_density)
}

compare_to_meyer_grant <- function() {
  ref_dir <- file.path("bmm-dev", "SDT_Gumbel")
  if (!dir.exists(ref_dir)) {
    message("Local Meyer-Grant reference directory not found at ", ref_dir)
    return(invisible(NULL))
  }

  message("Reference directory found: ", normalizePath(ref_dir))
  message("Use the scripts in that directory to compare:")
  message("  - gumbel binary densities against dsdt_binary(dist = 'gumbel_min')")
  message("  - aggregated rating fits against dsdt_rating(dist = 'gumbel_min')")
  message("  - ranking probabilities against dsdt_ranking(dist = 'gumbel_min')")
  invisible(ref_dir)
}

compare_to_bhsdtr2()
compare_to_meyer_grant()
