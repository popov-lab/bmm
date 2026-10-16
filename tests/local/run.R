# Runs the local tests (ADR 0004) from the package root:
#   Rscript tests/local/run.R
# Parallel workers are off: the package sets Config/testthat/parallel for
# tests/testthat, and the sampling tests here are better off one at a time.
Sys.setenv(NOT_CRAN = "true", TESTTHAT_PARALLEL = "false")
devtools::load_all(quiet = TRUE)
testthat::test_dir("tests/local", stop_on_failure = TRUE)
