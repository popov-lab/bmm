# Local tests (ADR 0004): tests that sample with a real backend. They run on a
# developer machine and in the weekly local-tests workflow, never in R CMD
# check (the folder is build-ignored). From the package root:
#
#   Rscript tests/local/run.R
#
# (run.R turns the package's parallel workers off and loads the package.)
# Only files named test-*.R run; scripts without the prefix (make-fixtures.R,
# make-references.R) are run by hand.
library(testthat)

options(
  brms.backend = "cmdstanr",
  mc.cores = 2,
  bmm.sort_data = TRUE,
  bmm.silent = 2
)

# every value inside [lower, upper]; the brms local tests compare estimates
# to reference values this way
expect_range <- function(object, lower = -Inf, upper = Inf, ...) {
  testthat::expect_true(all(object >= lower & object <= upper), ...)
}

SW <- suppressWarnings
SM <- suppressMessages
