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

test_that("helpers that read a single model refuse a multivariate fit", {
  fit <- mvm_mock_fit()
  expect_error(prior_info(fit), "not yet supported for multivariate")
  expect_error(native_parameters(fit), "not yet supported for multivariate")
  expect_error(model_citation(fit), "not yet supported for multivariate")
})

test_that("the expected response of a component is refused as for the model alone", {
  fit <- mvm_mock_fit()
  expect_error(posterior_epred(fit, resp = "error"), "not defined for the mixture2p model")
  expect_error(fitted(fit, resp = "error"), "not defined for the mixture2p model")
  expect_error(posterior_epred(fit), "not defined for the mixture2p model")
})

test_that("parameter_info() reports all component parameters with their response", {
  fit <- mvm_mock_fit()
  pars <- parameter_info(fit)

  expect_s3_class(pars, "bmm_parameters")
  expect_true(all(c("parameter", "response") %in% names(pars)))
  expect_true(all(c("mu1", "kappa", "thetat") %in% pars$parameter[pars$response == "error"]))
  expect_equal(pars$parameter[pars$response == "rt"], "sigma")
  expect_match(attr(pars, "model_name"), "2 components")
  expect_output(print(pars), "sigma")
})

test_that("parameter_info() skips family components without predicted parameters", {
  joint <- bmm_component(
    bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + (1 | p | id)),
    model = mixture2p(resp_error = "error"), data = mvm_dat_vwm
  ) +
    bmm_component(bmf(rt ~ 1 + (1 | p | id)), family = brms::lognormal(), data = mvm_dat_rt)
  fit <- bmm(joint, backend = "mock", mock_fit = 1, rename = FALSE)
  pars <- parameter_info(fit)
  expect_equal(unique(pars$response), "error")
})
