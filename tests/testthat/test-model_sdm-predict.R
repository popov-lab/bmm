test_that("sdm posterior_predict draws each prediction from its own posterior draw", {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  withr::local_seed(445)
  yrep <- brms::posterior_predict(fit)
  c_draws <- brms::posterior_linpred(fit, dpar = "c", transform = TRUE)
  kappa_draws <- brms::posterior_linpred(fit, dpar = "kappa", transform = TRUE)
  set_size <- fit$data$set_size

  # mu is fixed at 0, and c and kappa are constant within a set size, so each
  # draw's mean cos(yrep) within a set size estimates that draw's E[cos(y)]
  cells <- expand.grid(draw = seq_len(nrow(yrep)), set_size = levels(set_size))
  obs_i <- match(cells$set_size, set_size)
  cells$ref <- vapply(seq_len(nrow(cells)), function(j) {
    d <- cells$draw[j]
    i <- obs_i[j]
    integrate(function(x) cos(x) * dsdm(x, 0, c_draws[d, i], kappa_draws[d, i]), -pi, pi)$value
  }, numeric(1))
  cells$obs <- vapply(seq_len(nrow(cells)), function(j) {
    mean(cos(yrep[cells$draw[j], set_size == cells$set_size[j]]))
  }, numeric(1))

  slope <- stats::coef(stats::lm(obs ~ ref + set_size, cells))[["ref"]]
  expect_gt(slope, 0.75)
  expect_lt(slope, 1.25)
})
