# Runs the test suite from the source tree and enforces the time budget of
# ADR 0004. Called by the tests-source workflow on Ubuntu with
# BMM_TEST_STAN=true and four workers, so the fixtures in tests/testthat/assets
# and the Stan compile tests are included. Writes test-timings.csv, lists the
# slowest tests in the job summary, and exits with status 1 when a test fails
# or a cap is breached.
#
# The caps are provisional: set from local runs scaled to the runner, to be
# replaced by 1.5x the maxima of the first measured run (ADR 0004, step 7).
# Every later change goes through the ADR. Seconds throughout.
caps <- list(suite = 600, file = 300, test = 60)

# known offenders with their own caps, emptied by the #422 consolidation
allowed <- c(
  "test-update.R | update() takes a slice of the model frame as newdata (kfold(), reloo())" = 90,
  "test-update.R | only a slice of a fit's model frame counts as stored data" = 90
)

Sys.setenv(NOT_CRAN = "true")
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
    sprintf("Suite wall time %.0f s (cap %d s); %d tests, %d failed, %d errors, %d skipped.",
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
  breaches <- c(breaches, sprintf("suite wall time %.0f s exceeds %d s", wall, caps$suite))
}
over_file <- by_file[by_file$real > caps$file, ]
breaches <- c(breaches, sprintf("%s takes %.0f s, cap %d s", over_file$file, over_file$real, caps$file))
key <- paste(timings$file, timings$test, sep = " | ")
test_cap <- ifelse(key %in% names(allowed), allowed[key], caps$test)
over_test <- timings[timings$real > test_cap, ]
breaches <- c(breaches, sprintf("%s | %s takes %.1f s, cap %d s", over_test$file, over_test$test,
                                over_test$real, test_cap[timings$real > test_cap]))

cat(sprintf("\nSuite wall time %.0f s; %d tests, %d failed, %d errors, %d skipped\n",
            wall, nrow(timings), sum(timings$failed), sum(timings$error), sum(timings$skipped)))
cat("\nSlowest tests:\n")
print(slowest, row.names = FALSE)
if (length(breaches)) {
  cat("\nBudget breaches:\n", paste0(" - ", breaches, collapse = "\n"), "\n", sep = "")
}
failed <- sum(timings$failed) + sum(timings$error) > 0
if (failed || length(breaches)) quit(status = 1)
