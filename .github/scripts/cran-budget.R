# Runs the tests CRAN runs and enforces the CRAN-tier cap of ADR 0004 (R10):
# the built tarball's tests against the installed tarball, NOT_CRAN unset,
# one worker, no CmdStan and no fixtures, as on CRAN's Windows machine.
# Called by the tests-source workflow after the time budget. Writes
# cran-test-timings.csv, lists the slowest tests in the job summary, and
# exits with status 1 when a test fails or the cap is breached.
#
# CRAN returns a package whose check takes more than 10 min. The incoming
# pretest of 1.4.0 (2026-10-15) took 16 min on win-builder's Windows, 12 of
# them tests, which ran 4 to 5 times slower there than sequentially on a
# developer machine; the non-test steps take about 5 min there. A test
# over 0.4 s here on a developer machine (about 2 s on that machine) opens
# with skip_on_cran() (R10). The cap is provisional until the first run on
# the runner is measured (ADR 0004 § Changes). Seconds throughout.
cap <- 120

Sys.unsetenv("NOT_CRAN")
Sys.setenv(TESTTHAT_PARALLEL = "false", BMM_TEST_STAN = "false")

work <- tempfile("cran-budget")
dir.create(work)
tarball <- pkgbuild::build(".", dest_path = work, vignettes = FALSE, quiet = TRUE)
lib <- file.path(work, "lib")
dir.create(lib)
utils::install.packages(tarball, lib = lib, repos = NULL, type = "source", quiet = TRUE)
utils::untar(tarball, exdir = work)
withr::local_libpaths(lib, action = "prefix")

started <- Sys.time()
results <- testthat::test_dir(
  file.path(work, "bmm", "tests", "testthat"),
  package = "bmm", load_package = "installed",
  reporter = testthat::ListReporter$new(), stop_on_failure = FALSE
)
wall <- as.numeric(difftime(Sys.time(), started, units = "secs"))

timings <- as.data.frame(results)
timings <- timings[, c("file", "test", "nb", "failed", "skipped", "error", "warning", "real")]
write.csv(timings, "cran-test-timings.csv", row.names = FALSE)
slowest <- head(timings[order(-timings$real), c("file", "test", "real")], 20)

summary_path <- Sys.getenv("GITHUB_STEP_SUMMARY", "")
if (nzchar(summary_path)) {
  lines <- c(
    "## CRAN-tier timings (one worker, no fixtures, no CmdStan)",
    "",
    sprintf("Wall time %.0f s (cap %g s); %d tests, %d failed, %d errors, %d skipped.",
            wall, cap, nrow(timings), sum(timings$failed), sum(timings$error),
            sum(timings$skipped)),
    "",
    "| file | test | seconds |",
    "| --- | --- | --- |",
    sprintf("| %s | %s | %.1f |", slowest$file, substr(slowest$test, 1, 80), slowest$real)
  )
  cat(lines, file = summary_path, sep = "\n", append = TRUE)
}

cat(sprintf("\nCRAN-tier wall time %.0f s (cap %g s); %d tests, %d failed, %d errors, %d skipped\n",
            wall, cap, nrow(timings), sum(timings$failed), sum(timings$error),
            sum(timings$skipped)))
print(slowest, row.names = FALSE)

failed <- sum(timings$failed) + sum(timings$error)
if (failed > 0) {
  cat(sprintf("\n%d test(s) failed or errored in CRAN mode\n", failed))
}
if (wall > cap) {
  cat(sprintf("\nCRAN-tier wall time %.0f s exceeds the cap of %g s (ADR 0004, R10)\n", wall, cap))
}
if (failed > 0 || wall > cap) quit(status = 1)
