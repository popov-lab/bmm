# Runs the test suite from the source tree and enforces the time budget of
# ADR 0004. Called by the tests-source workflow on Ubuntu with
# BMM_TEST_STAN=true and four workers, so the fixtures in tests/testthat/assets
# and the Stan compile tests are included. Writes test-timings.csv, lists the
# slowest tests in the job summary, and exits with status 1 when a test fails
# or a cap is breached.
#
# The file and test caps, and the allowlist caps below, are 1.5x the maxima of
# the first measured run on the runner (ADR 0004), rounded: the slowest file
# took 221 s, the slowest test outside the allowlist 23 s, the two allowlisted
# tests 44 s and 28.5 s. That run had one worker, so its wall time cannot set
# the four-worker suite cap, which stays provisional until a four-worker run.
# Every later change goes through the ADR. Seconds throughout.
caps <- list(suite = 600, file = 330, test = 35)

# known offenders with their own caps, emptied by the #422 consolidation
allowed <- c(
  "test-update.R | update() takes a slice of the model frame as newdata (kfold(), reloo())" = 66,
  "test-update.R | only a slice of a fit's model frame counts as stored data" = 43
)

Sys.setenv(NOT_CRAN = "true")
# setup-r writes options(Ncpus = 1) to the site Rprofile and testthat reads
# Ncpus before TESTTHAT_CPUS, so without this the suite runs on one worker
options(Ncpus = NULL)
started <- Sys.time()
results <- devtools::test(
  reporter = testthat::ListReporter$new(),
  stop_on_failure = FALSE
)
wall <- as.numeric(difftime(Sys.time(), started, units = "secs"))

timings <- as.data.frame(results)
timings <- timings[, c("file", "test", "nb", "failed", "skipped", "error", "warning", "real")]
write.csv(timings, "test-timings.csv", row.names = FALSE)

by_file <- aggregate(real ~ file, timings, sum)
slowest <- head(timings[order(-timings$real), c("file", "test", "real")], 20)

summary_path <- Sys.getenv("GITHUB_STEP_SUMMARY", "")
if (nzchar(summary_path)) {
  lines <- c(
    "## Test timings",
    "",
    sprintf("Suite wall time %.0f s (cap %g s); %d tests, %d failed, %d errors, %d skipped.",
            wall, caps$suite, nrow(timings), sum(timings$failed), sum(timings$error),
            sum(timings$skipped)),
    "",
    "| file | test | seconds |",
    "| --- | --- | --- |",
    sprintf("| %s | %s | %.1f |", slowest$file, substr(slowest$test, 1, 80), slowest$real)
  )
  cat(lines, file = summary_path, sep = "\n", append = TRUE)
}

breaches <- character()
if (wall > caps$suite) {
  breaches <- c(breaches, sprintf("suite wall time %.0f s exceeds %g s", wall, caps$suite))
}
over_file <- by_file[by_file$real > caps$file, ]
breaches <- c(breaches, sprintf("%s takes %.0f s, cap %g s", over_file$file, over_file$real, caps$file))
key <- paste(timings$file, timings$test, sep = " | ")
test_cap <- ifelse(key %in% names(allowed), allowed[key], caps$test)
# which() drops the row of a file that errored outside test_that(), whose
# time is NA; the error itself fails the job below
slow <- which(timings$real > test_cap)
breaches <- c(breaches, sprintf("%s | %s takes %.1f s, cap %g s", timings$file[slow],
                                timings$test[slow], timings$real[slow], test_cap[slow]))

cat(sprintf("\nSuite wall time %.0f s; %d tests, %d failed, %d errors, %d skipped\n",
            wall, nrow(timings), sum(timings$failed), sum(timings$error), sum(timings$skipped)))
cat("\nSlowest tests:\n")
print(slowest, row.names = FALSE)
if (length(breaches)) {
  cat("\nBudget breaches:\n", paste0(" - ", breaches, collapse = "\n"), "\n", sep = "")
}
failed <- sum(timings$failed) + sum(timings$error) > 0
if (failed || length(breaches)) quit(status = 1)
