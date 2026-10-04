# A stand-in for the brmsprep object brms hands a custom family's
# posterior_epred and posterior_predict: posterior_epred.brmsprep() evaluates
# every dpar to a draws x observations matrix before it calls the family, and
# brms::get_dpar() slices such a matrix by observation
epred_prep <- function(dpars, data = list()) {
  dpars <- lapply(dpars, as.matrix)
  structure(
    list(dpars = dpars, data = data,
         ndraws = nrow(dpars[[1]]), nobs = ncol(dpars[[1]])),
    class = "brmsprep"
  )
}

# For each parameter set (a row of `sets`, one column per dpar), the family's
# posterior_epred on a single draw next to the mean of n_draws responses from
# its posterior_predict, which simulates once per draw of the observation
epred_vs_predict <- function(sets, epred, predict, data = list(), n_draws = 20000) {
  out <- t(vapply(seq_len(nrow(sets)), function(j) {
    data_j <- lapply(data, `[`, j)
    one <- lapply(sets[j, , drop = FALSE], function(v) matrix(v, 1, 1))
    many <- lapply(sets[j, , drop = FALSE], function(v) matrix(v, n_draws, 1))
    c(epred = epred(epred_prep(one, data_j))[1, 1],
      mc = mean(predict(1, epred_prep(many, data_j)), na.rm = TRUE))
  }, numeric(2)))
  cbind(out, rel_error = out[, "mc"] / out[, "epred"] - 1)
}

# posterior_epred on a draws x observations prep must match the same function
# evaluated one cell at a time; a data column recycled down the draws instead
# of across the observations fails this
expect_epred_by_cell <- function(epred, dpars, data = list()) {
  prep <- epred_prep(dpars, data)
  out <- epred(prep)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(prep$ndraws, prep$nobs))
  by_cell <- matrix(NA_real_, prep$ndraws, prep$nobs)
  for (s in seq_len(prep$ndraws)) {
    for (i in seq_len(prep$nobs)) {
      cell <- lapply(prep$dpars, function(x) x[s, i, drop = FALSE])
      by_cell[s, i] <- epred(epred_prep(cell, lapply(data, `[`, i)))[1, 1]
    }
  }
  expect_equal(unname(out), by_cell)
}

load_fixture_fit <- function(name) {
  path <- test_path("assets", name)
  skip_if_not(file.exists(path), "fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}
