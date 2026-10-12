# One warning for everything a single user call deprecates. A caller that
# translates many rows (priors, formula terms) collects the old and new forms
# first and calls this once; the deprecation policy in CONTRIBUTING.md
# ("Deprecating and removing features") gives the rules this implements.
# `what` and `with` are code, written as the user would type them.
warn_deprecated <- function(what, with, since, removed_in = next_minor(since),
                            details = NULL) {
  statements <- unique(glue(
    "`{what}` is deprecated since bmm {since} and will be removed in bmm ",
    "{removed_in} at the earliest. Use `{with}` instead."
  ))
  msg <- paste(c(statements, details), collapse = "\n")
  warning2("{msg}")
}

next_minor <- function(version) {
  parts <- unlist(package_version(version))
  paste(parts[1], parts[2] + 1L, 0L, sep = ".")
}
