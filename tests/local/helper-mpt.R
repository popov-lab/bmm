# the MPT fits use the tree and data builders of tests/testthat; testthat runs
# this file from tests/local, so the path is relative to it
source(file.path("..", "testthat", "helper-mpt.R"), local = environment())

# The known-value fits (test-mpt-known-values.R) simulate with rmultinom() on
# branch probabilities written out by hand, never through rmpt() or the tree
# objects, so an error in bmm's tree algebra cannot cancel out.
draw_counts <- function(probs, n) {
  counts <- stats::rmultinom(1, n, probs)[, 1]
  stats::setNames(as.integer(counts), names(probs))
}

# A population estimate is the posterior median of native_parameters() at
# re_formula = NA; it passes when it lies within `z` posterior SDs (MAD) of
# the value the data were simulated with, written as a ratio to that value.
# The SD itself must stay below `max_rel_sd` of the value, so that a loose
# posterior cannot pass by being loose. `truth` names one vector per parameter
# with one value per row of `newdata`, or a single value without newdata.
expect_known_values <- function(fit, truth, newdata = NULL, re_formula = NA,
                                z = 4, max_rel_sd = 0.05) {
  est <- native_parameters(
    fit, newdata = newdata, pars = names(truth), re_formula = re_formula,
    summary = TRUE, robust = TRUE
  )
  for (par in names(truth)) {
    rows <- est[est$parameter == par, ]
    expect_equal(nrow(rows), length(truth[[par]]), label = paste(par, "grid cells"))
    rel_sd <- rows$Est.Error / truth[[par]]
    expect_lt(max(rel_sd), max_rel_sd, label = paste(par, "posterior SD / true value"))
    expect_lt(
      max(abs(rows$Estimate / truth[[par]] - 1) / rel_sd), z,
      label = paste(par, "|median / true - 1| in posterior SDs")
    )
  }
  invisible(est)
}
