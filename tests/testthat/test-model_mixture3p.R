dat <- small_data()
model <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")

test_that("mixture3p works when set_size is not predicted and there is set_size 1", {
  skip_on_cran()
  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("mixture3p gives an error if set_size is a predictor but there is an intercept", {
  skip_on_cran()
  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ set_size)
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )

  formula <- bmf(kappa ~ 1 + session, thetat ~ 1, thetant ~ 0 + set_size)
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("constant priors are correct for mixture3p with set_size1 fixed effect predictor", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size)
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "thetant", "prior"], "constant(-100)")
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "thetat", "prior"], character(0))
})

test_that("constant priors are correct for mixture3p with set_size1 RANDOM-EFFECT predictor", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size + (0 + set_size | ID))
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "thetat", "prior"], character(0))
  expect_equal(
    prior[prior$coef == "set_size1" & prior$nlpar == "thetant", "prior"],
    c("constant(-100)", "constant(1e-8)")
  )
})
