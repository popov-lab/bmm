# Runs the test suite from the source tree and enforces the time budget of
# ADR 0004. Called by the tests-source workflow on Ubuntu with
# BMM_TEST_STAN=true and four workers, so the fixtures in tests/testthat/assets
# and the Stan compile tests are included. Writes test-timings.csv, lists the
# slowest tests in the job summary, and exits with status 1 when a test fails
# or a cap is breached.
#
# The caps are 1.5x the maxima of the first four-worker run on the runner
# (ADR 0004), rounded: suite 532 s, slowest file 450 s, slowest test outside
# the allowlist 44 s. Four workers on
# the 4-vCPU runner double each test's time against one worker, so caps from
# a one-worker run do not carry over. Every later change goes through the ADR.
# Seconds throughout.
caps <- list(suite = 800, file = 680, test = 66)

# tests with a cap of their own, "file | test" = seconds; empty since #422,
# and an entry needs one in ADR 0004
allowed <- stats::setNames(numeric(), character())

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
# the list reporter keeps the messages out of the log, so print them here
bad <- which(timings$failed > 0 | timings$error)
if (length(bad)) {
  messages <- vapply(results[bad], function(test) {
    hits <- Filter(function(x) inherits(x, c("expectation_failure", "expectation_error")),
                   test$results)
    paste(vapply(hits, conditionMessage, character(1)), collapse = "\n")
  }, character(1))
  cat("\nFailing tests:\n",
      paste0("- ", timings$file[bad], " | ",
             ifelse(is.na(timings$test[bad]), "(code outside test_that())", timings$test[bad]),
             "\n", messages, collapse = "\n\n"),
      "\n", sep = "")
}
cat("\nSlowest tests:\n")
print(slowest, row.names = FALSE)
if (length(breaches)) {
  cat("\nBudget breaches:\n", paste0(" - ", breaches, collapse = "\n"), "\n", sep = "")
}
failed <- sum(timings$failed) + sum(timings$error) > 0
if (failed || length(breaches)) quit(status = 1)
