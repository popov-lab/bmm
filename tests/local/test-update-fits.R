# update() on a real fit (moved from tests/testthat/test-update.R in #422)

skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")

test_that("brms::combine_models() merges a fit with its update() refit", {
  fit_args <- list(chains = 1, iter = 200, warmup = 100, refresh = 0, silent = 2)
  fit1 <- suppressWarnings(suppressMessages(brms::do_call(bmm, c(list(
    bmf(thetat ~ 1, kappa ~ 1),
    oberauer_lin_2017[oberauer_lin_2017$ID == 1, ],
    mixture2p("dev_rad"),
    backend = "cmdstanr", seed = 1
  ), fit_args))))
  fit2 <- suppressWarnings(suppressMessages(
    brms::do_call(update, c(list(fit1, seed = 2), fit_args))
  ))

  combined <- brms::combine_models(fit1, fit2)
  expect_s3_class(combined, "bmmfit")
  expect_equal(brms::nchains(combined), 2)
})
