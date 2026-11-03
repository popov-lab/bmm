# Runs the test suite from the source tree and enforces the time budget of
# ADR 0004. Called by the tests-source workflow on Ubuntu with
# BMM_TEST_STAN=true and four workers, so the fixtures in tests/testthat/assets
# and the Stan compile tests are included. Writes test-timings.csv, lists the
# slowest tests in the job summary, and exits with status 1 when a test fails
# or a cap is breached.
#
# The runners of one image differ in speed, the slowest taking about twice
# as long as the fastest on every file alike (2026-11-03), so the file and
# test caps are shares of the summed time of all tests, which the runner
# leaves in place.
# Over 39 runs (2026-10-30 to 2026-11-03), with test-update.R counted as the
# three files it was split into, the largest file's share ranged 0.144 to
# 0.156 (test-conditional_effects.R) and the slowest test's 0.017 to 0.022;
# the caps are 1.5x those maxima, rounded. The suite cap stays in seconds of
# wall time, the one guard on the whole suite slowing down: 1.5x the maximum
# of twelve runs from 2026-10-19 and 2026-10-20 (668 s). The slowest runners
# reached 862 s by 2026-11-03, so read a breach against the runner's CPU in
# the job summary. A bad start order adds to that: the first run after the
# split took 906 s because of the order, not the runner (ADR 0004,
# 2026-11-03). Four workers on the 4-vCPU runner double each test's
# time against one worker, so caps from a one-worker run do not carry over.
# Every later change goes through the ADR.
caps <- list(suite = 1000, file = 0.23, test = 0.033)

# tests with a cap of their own, "file | test" = share; empty since #422, and
# an entry needs one in ADR 0004
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

# a file that errored outside test_that() has no time (see which() below)
total <- sum(timings$real, na.rm = TRUE)
# with no time every share is NaN and a breach would read "NA takes NA s"
if (!is.finite(total) || total <= 0) {
  stop("no test time recorded (summed time ", total, " s); see test-timings.csv")
}
timings$share <- timings$real / total
by_file <- aggregate(real ~ file, timings, sum)
by_file$share <- by_file$real / total
largest <- by_file[which.max(by_file$share), ]
slowest <- head(timings[order(-timings$real), c("file", "test", "real", "share")], 20)

summary_path <- Sys.getenv("GITHUB_STEP_SUMMARY", "")
if (nzchar(summary_path)) {
  lines <- c(
    "## Test timings",
    "",
    sprintf("Suite wall time %.0f s (cap %g s); %d tests, %d failed, %d errors, %d skipped.",
            wall, caps$suite, nrow(timings), sum(timings$failed), sum(timings$error),
            sum(timings$skipped)),
    sprintf("Summed test time %.0f s; largest file %s, share %.3f (cap %g); slowest test share %.4f (cap %g).",
            total, largest$file, largest$share, caps$file, max(timings$share, na.rm = TRUE), caps$test),
    "",
    "| file | test | seconds | share |",
    "| --- | --- | --- | --- |",
    sprintf("| %s | %s | %.1f | %.4f |", slowest$file, substr(slowest$test, 1, 80), slowest$real,
            slowest$share)
  )
  cat(lines, file = summary_path, sep = "\n", append = TRUE)
}

breaches <- character()
if (wall > caps$suite) {
  breaches <- c(breaches, sprintf("suite wall time %.0f s exceeds %g s", wall, caps$suite))
}
over_file <- by_file[by_file$share > caps$file, ]
breaches <- c(breaches, sprintf("%s takes %.0f s, share %.3f of %.0f s, cap %g", over_file$file,
                                over_file$real, over_file$share, total, caps$file))
key <- paste(timings$file, timings$test, sep = " | ")
test_cap <- ifelse(key %in% names(allowed), allowed[key], caps$test)
# which() drops the row of a file that errored outside test_that(), whose
# time is NA; the error itself fails the job below
slow <- which(timings$share > test_cap)
breaches <- c(breaches, sprintf("%s | %s takes %.1f s, share %.4f of %.0f s, cap %g",
                                timings$file[slow], timings$test[slow], timings$real[slow],
                                timings$share[slow], total, test_cap[slow]))

cat(sprintf("\nSuite wall time %.0f s; summed test time %.0f s; %d tests, %d failed, %d errors, %d skipped\n",
            wall, total, nrow(timings), sum(timings$failed), sum(timings$error), sum(timings$skipped)))
cat(sprintf("Largest file %s, share %.3f\n", largest$file, largest$share))
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
