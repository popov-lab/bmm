# One warning for everything a single user call deprecates. A caller that
# translates many rows (priors, formula terms) collects the old and new forms
# first and calls this once, also when it found none; the deprecation policy
# in CONTRIBUTING.md ("Deprecating and removing features") gives the rules this
# implements.
# `what` and `with` are code, written as the user would type them; `with` is
# NULL when nothing replaces the old form.
warn_deprecated <- function(what, with, since, removed_in = next_minor(since),
                            details = NULL) {
  if (length(what) == 0) {
    return(invisible())
  }
  replacement <- if (length(with) == 0) "" else glue(" Use `{with}` instead.", .trim = FALSE)
  statements <- unique(glue(
    "`{what}` is deprecated since bmm {since} and will be removed in bmm ",
    "{removed_in} at the earliest.{replacement}"
  ))
  max_listed <- 20L
  not_listed <- length(statements) - max_listed
  if (not_listed > 0) {
    statements <- c(statements[seq_len(max_listed)], glue("Not listed: {not_listed} more."))
  }
  msg <- paste(c(statements, details), collapse = "\n")

  # R cuts a warning at 1000 characters by default, and `details` comes last
  old_options <- options(warning.length = 8170L)
  on.exit(options(old_options))
  warning2("{msg}")
}

next_minor <- function(version) {
  parts <- unlist(package_version(version))
  paste(parts[1], parts[2] + 1L, 0L, sep = ".")
}
