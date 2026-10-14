dat <- small_data()

test_that("imm works when set_size is not predicted and there is set_size 1", {
  skip_on_cran()
  formula <- bmf(kappa ~ 1, a ~ 1, c ~ 1, s ~ 1)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size"
  )
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("imm_abc works when set_size is not predicted and there is set_size 1", {
  skip_on_cran()
  formula <- bmf(kappa ~ 1, a ~ 1, c ~ 1)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("imm_bsc works when set_size is not predicted and there is set_size 1", {
  skip_on_cran()
  formula <- bmf(kappa ~ 1, c ~ 1, s ~ 1)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "bsc"
  )
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("IMM models give an error if set_size is a predictor but there is an intercept", {
  skip_on_cran()
  formula <- bmf(kappa ~ 1, c ~ 1, a ~ set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )

  formula <- bmf(kappa ~ 1, c ~ 1, s ~ set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "bsc"
  )
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )

  formula <- bmf(kappa ~ 1, a ~ 1, c ~ 1, s ~ set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size"
  )
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )

  formula <- bmf(kappa ~ 1, c ~ 0 + set_size, a ~ set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )

  formula <- bmf(kappa ~ 0 + set_size, c ~ 1, s ~ set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "bsc"
  )
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )

  formula <- bmf(kappa ~ 1, a ~ 1, c ~ 0 + set_size, s ~ set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size"
  )
  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )
})

test_that("IMM models run when set_size is a predictor and intercept is supressed", {
  skip_on_cran()
  formula <- bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))

  formula <- bmf(kappa ~ 1, c ~ 1, s ~ 0 + set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "bsc"
  )
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))

  formula <- bmf(kappa ~ 1, a ~ 1, c ~ 1, s ~ 0 + set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size"
  )
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("constant priors are correct for IMM_abc with set_size1 fixed effect predictor", {
  formula <- bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "a", "prior"], "constant(0)")
})

test_that("constant priors are correct for IMM_bsc with set_size1 fixed effect predictor", {
  formula <- bmf(kappa ~ 1, c ~ 1, s ~ 0 + set_size)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "bsc"
  )
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "s", "prior"], "constant(0)")
})

test_that("constant priors are correct for IMM_full with set_size1 fixed effect predictor", {
  withr::local_options(bmm.silent = 2)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "full"
  )
  formula <- bmf(kappa ~ 1, c ~ 1, s ~ 0 + set_size)
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "s", "prior"], "constant(0)")
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "a", "prior"], character(0))

  formula <- bmf(kappa ~ 1, c ~ 1, s ~ 0 + set_size, a ~ 0 + set_size)
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "s", "prior"], "constant(0)")
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "a", "prior"], "constant(0)")
})

test_that("constant priors are correct for IMM_full with set_size1 RANDOM-EFFECT predictor", {
  withr::local_options(bmm.silent = 2)
  model <- imm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "full"
  )
  formula <- bmf(kappa ~ 1, c ~ 1, s ~ 0 + set_size + (0 + set_size | ID))
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "a", "prior"], character(0))
  expect_equal(
    prior[prior$coef == "set_size1" & prior$nlpar == "s", "prior"],
    c("constant(0)", "constant(1e-8)")
  )

  formula <- bmf(kappa ~ 1, c ~ 1, s ~ 0 + set_size, a ~ 0 + set_size + (0 + set_size | ID))
  fit <- bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  prior <- brms::prior_summary(fit)
  expect_equal(prior[prior$coef == "set_size1" & prior$nlpar == "s", "prior"], "constant(0)")
  expect_equal(
    prior[prior$coef == "set_size1" & prior$nlpar == "a", "prior"],
    c("constant(0)", "constant(1e-8)")
  )
})
