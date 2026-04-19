dat <- oberauer_lin_2017

test_that("sdm_simple works without non-target inputs", {
  withr::local_options(bmm.silent = 2)
  formula <- bmf(c ~ 1, kappa ~ 1)
  model <- sdm(resp_error = "dev_rad", version = "simple")
  expect_silent(bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("non-simple SDM versions run with mock backend", {
  withr::local_options(bmm.silent = 2)

  formula_abc <- bmf(c ~ 1, a ~ 1, kappa ~ 1)
  model_abc <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )
  expect_silent(bmm(formula_abc, dat, model_abc, backend = "mock", mock = 1, rename = FALSE))

  formula_bsc <- bmf(c ~ 1, s ~ 1, kappa ~ 1)
  model_bsc <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "bsc"
  )
  expect_silent(bmm(formula_bsc, dat, model_bsc, backend = "mock", mock = 1, rename = FALSE))

  formula_full <- bmf(c ~ 1, a ~ 1, s ~ 1, kappa ~ 1)
  model_full <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "full"
  )
  expect_silent(bmm(formula_full, dat, model_full, backend = "mock", mock = 1, rename = FALSE))
})

test_that("SDM non-simple versions reuse non-target preprocessing", {
  simple_checked <- check_data(
    sdm(resp_error = "dev_rad", version = "simple"),
    dat[1:8, ],
    bmf(c ~ 1, kappa ~ 1)
  )
  expect_false("ss_numeric" %in% names(simple_checked))
  expect_false(any(grepl("^LureIdx", names(simple_checked))))

  abc_checked <- check_data(
    sdm(
      resp_error = "dev_rad",
      nt_features = paste0("col_nt", 1:7),
      set_size = "set_size",
      version = "abc"
    ),
    dat[1:8, ],
    bmf(c ~ 1, a ~ 1, kappa ~ 1)
  )
  expect_true("ss_numeric" %in% names(abc_checked))
  expect_true(any(grepl("^LureIdx", names(abc_checked))))
})

test_that("SDM non-simple versions support regex expansion", {
  explicit <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "full"
  )
  regex_model <- sdm(
    resp_error = "dev_rad",
    nt_features = "col_nt",
    nt_distances = "dist_nt",
    set_size = "set_size",
    regex = TRUE,
    version = "full"
  )

  checked_explicit <- check_model(explicit, dat)
  checked_regex <- check_model(regex_model, dat)
  attributes(checked_explicit) <- NULL
  attributes(checked_regex) <- NULL
  expect_equal(checked_explicit, checked_regex)
})

test_that("SDM non-simple versions require intercept suppression for set_size predictors", {
  withr::local_options(bmm.silent = 2)
  formula <- bmf(c ~ 1, a ~ set_size, kappa ~ 1)
  model <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )

  expect_error(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE),
    "This model requires that the intercept is supressed when set_size is used as predictor."
  )
})

test_that("constant priors are correct for non-simple SDM versions with set_size1 predictors", {
  withr::local_options(bmm.silent = 2)

  model_abc <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    set_size = "set_size",
    version = "abc"
  )
  fit_abc <- bmm(
    bmf(c ~ 1, a ~ 0 + set_size, kappa ~ 1),
    dat,
    model_abc,
    backend = "mock",
    mock = 1,
    rename = FALSE
  )
  prior_abc <- brms::prior_summary(fit_abc)
  expect_equal(prior_abc[prior_abc$coef == "set_size1" & prior_abc$dpar == "a", "prior"], "constant(0)")

  model_bsc <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "bsc"
  )
  fit_bsc <- bmm(
    bmf(c ~ 1, s ~ 0 + set_size, kappa ~ 1),
    dat,
    model_bsc,
    backend = "mock",
    mock = 1,
    rename = FALSE
  )
  prior_bsc <- brms::prior_summary(fit_bsc)
  expect_equal(prior_bsc[prior_bsc$coef == "set_size1" & prior_bsc$dpar == "s", "prior"], "constant(0)")

  model_full <- sdm(
    resp_error = "dev_rad",
    nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size",
    version = "full"
  )
  fit_full <- bmm(
    bmf(c ~ 1, a ~ 0 + set_size, s ~ 0 + set_size, kappa ~ 1),
    dat,
    model_full,
    backend = "mock",
    mock = 1,
    rename = FALSE
  )
  prior_full <- brms::prior_summary(fit_full)
  expect_equal(prior_full[prior_full$coef == "set_size1" & prior_full$dpar == "a", "prior"], "constant(0)")
  expect_equal(prior_full[prior_full$coef == "set_size1" & prior_full$dpar == "s", "prior"], "constant(0)")
})
