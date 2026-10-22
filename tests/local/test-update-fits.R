# update() on real fits: brms::combine_models() on a refit (moved from
# tests/testthat/test-update.R in #422) and the no-recompile route of #520

test_that("brms::combine_models() merges a fit with its update() refit", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")
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

test_that("update() without recompiling stores the init it sampled with (#520)", {
  # tests/testthat reaches brms's no-recompile route only by mocking
  # compiled_model(); here rstan reuses its compiled model the way a user's fit
  # does, so a change in how brms picks that route shows up here first
  skip_if_not_installed("rstan")
  withr::local_seed(520)
  dat <- data.frame(y = rsdm(120, c = 4, kappa = 5), set_size = factor(rep(1:3, each = 40)))
  fit <- suppressMessages(suppressWarnings(bmm(
    bmf(c ~ 0 + set_size, kappa ~ 1), dat, sdm("y"),
    backend = "rstan", chains = 1, iter = 40, warmup = 20, refresh = 0, silent = 2, rename = FALSE
  )))
  new_data <- droplevels(dat[dat$set_size %in% 1:2, ])
  update_rstan <- function(...) suppressMessages(suppressWarnings(
    update(fit, newdata = new_data, ..., chains = 1, iter = 40, warmup = 20, refresh = 0,
           silent = 2, rename = FALSE)
  ))

  up <- update_rstan()
  # the same compiled model means brms took the route that leaves stan_args alone
  expect_identical(up$fit@stanmodel, fit$fit@stanmodel)
  expect_length(fit$stan_args$init()$b_c, 3)
  expect_length(up$stan_args$init()$b_c, brms::standata(up)$K_c)

  up <- update_rstan(init = 0)
  expect_identical(up$fit@stanmodel, fit$fit@stanmodel)
  expect_identical(up$stan_args$init, 0)
})
