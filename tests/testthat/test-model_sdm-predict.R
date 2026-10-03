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

  # the slope's seed-to-seed SD is about 0.1 and pooled or mispaired samplers
  # stay below 0.4, so bounds 4 SD from 1 survive changes in RNG use
  slope <- stats::coef(stats::lm(obs ~ ref + set_size, cells))[["ref"]]
  expect_gt(slope, 0.6)
  expect_lt(slope, 1.4)
})

# posterior_epred (#475) ------------------------------------------------------

test_that("an sdm fit refuses posterior_epred() with a message naming the model", {
  skip_on_cran()
  dat <- data.frame(y = rsdm(50))
  fit <- bmm(bmf(c ~ 1, kappa ~ 1), dat, sdm(resp_error = "y"),
             backend = "mock", mock = 1, rename = FALSE)
  expect_true(is.function(fit$formula$family$posterior_epred))
  expect_error(fit$formula$family$posterior_epred(NULL),
               "The expected response is not defined for the sdm model")
})

test_that("posterior_epred() on an sdm fit saved without the function gives the message", {
  # the fixture's model predates the circular class, so the family function,
  # which restructure() adds, is what refuses here
  fit <- load_fixture_fit("bmmfit_example1.rds")
  expect_error(brms::posterior_epred(fit, ndraws = 5),
               "The expected response is not defined for the sdm model")
  # the model parameters stay available
  expect_equal(dim(brms::posterior_epred(fit, dpar = "c", ndraws = 5)),
               c(5L, nrow(fit$data)))
})
