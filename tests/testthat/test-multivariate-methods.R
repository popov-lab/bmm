mvm_dat_vwm <- data.frame(
  id = factor(rep(1:8, each = 10)),
  error = runif(80, -pi, pi)
)
mvm_dat_rt <- data.frame(
  id = factor(rep(1:8, each = 10)),
  rt = rlnorm(80, meanlog = -0.5, sdlog = 0.3),
  cond = rep(c("a", "b"), 40)
)

mvm_mock_fit <- function() {
  joint <- bmm_component(
    bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + (1 | p | id)),
    model = mixture2p(resp_error = "error"), data = mvm_dat_vwm
  ) +
    bmm_component(
      bmf(rt ~ 1 + (1 | p | id), sigma ~ cond),
      family = brms::lognormal(), data = mvm_dat_rt
    )
  bmm(joint, backend = "mock", mock_fit = 1, rename = FALSE)
}

test_that("unsupported post-processing methods give informative errors", {
  fit <- mvm_mock_fit()
  expect_error(update(fit, chains = 1), "not yet supported for multivariate")
  expect_error(conditional_effects(fit), "not yet supported for multivariate")
})

test_that("pp_check() requires a valid resp argument", {
  fit <- mvm_mock_fit()
  expect_error(pp_check(fit), "requires the 'resp' argument")
  expect_error(pp_check(fit), "'error', 'rt'")
  expect_error(pp_check(fit, resp = "bogus"), "Unknown response 'bogus'")
})

test_that("pp_check() rejects multinomial components", {
  m3dat <- oberauer_lewandowsky_2019_e1
  rt_dat <- data.frame(
    ID = unique(m3dat$ID),
    rt = rlnorm(length(unique(m3dat$ID)))
  )
  joint <- suppressMessages(
    bmm_component(
      bmf(c ~ 1 + (1 | p | ID), a ~ 1 + (1 | p | ID)),
      model = m3(
        resp_cats = c("corr", "other", "npl"),
        num_options = c("n_corr", "n_other", "n_npl"),
        choice_rule = "simple", version = "ss"
      ),
      data = m3dat
    ) +
      bmm_component(bmf(rt ~ 1 + (1 | p | ID)), family = brms::lognormal(), data = rt_dat)
  )
  fit <- suppressMessages(bmm(joint, backend = "mock", mock_fit = 1, rename = FALSE))
  expect_error(pp_check(fit, resp = "Y"), "not yet supported for multinomial")
})

test_that("parameters() reports all component parameters with their response", {
  fit <- mvm_mock_fit()
  pars <- parameters(fit)

  expect_s3_class(pars, "bmm_parameters")
  expect_true(all(c("parameter", "response") %in% names(pars)))
  expect_true(all(c("mu1", "kappa", "thetat") %in% pars$parameter[pars$response == "error"]))
  expect_equal(pars$parameter[pars$response == "rt"], "sigma")
  expect_match(attr(pars, "model_name"), "2 components")
  expect_output(print(pars), "sigma")
})

test_that("summary() of a multivariate fit shows components and correlations", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example_mv.rds")
  skip_if_not(file.exists(path), "multivariate fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  smry <- suppressWarnings(summary(fit))
  expect_s3_class(smry, "mvbmmsummary")
  expect_length(smry$components, 2)
  expect_equal(
    vapply(smry$components, function(x) x$resp_name, character(1)),
    c("error", "rt")
  )
  expect_true("cor(error_kappa_Intercept,rt_Intercept)" %in% rownames(smry$random$id))

  out <- capture.output(print(smry, color = FALSE))
  expect_true(any(grepl("Multivariate bmm model with 2 components", out)))
  expect_true(any(grepl("Component 1 \\[error\\]", out)))
  expect_true(any(grepl("Regression Coefficients \\[error\\]", out)))
  expect_true(any(grepl("Regression Coefficients \\[rt\\]", out)))
  expect_true(any(grepl("Constant Parameters \\[error\\]", out)))
  expect_true(any(grepl("error_kappa_Intercept", out)))

  brms_smry <- suppressWarnings(summary(fit, backend = "brms"))
  expect_s3_class(brms_smry, "brmssummary")
})

test_that("pp_check() works on a component of a multivariate fit", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example_mv.rds")
  skip_if_not(file.exists(path), "multivariate fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  plot <- pp_check(fit, resp = "rt", ndraws = 10)
  expect_s3_class(plot, "ggplot")
})

test_that("parameters() skips family components without predicted parameters", {
  joint <- bmm_component(
    bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + (1 | p | id)),
    model = mixture2p(resp_error = "error"), data = mvm_dat_vwm
  ) +
    bmm_component(bmf(rt ~ 1 + (1 | p | id)), family = brms::lognormal(), data = mvm_dat_rt)
  fit <- bmm(joint, backend = "mock", mock_fit = 1, rename = FALSE)
  pars <- parameters(fit)
  expect_equal(unique(pars$response), "error")
})
