test_that("Available mock models run without errors", {
  withr::local_options("bmm.silent" = 2)
  skip_on_cran()
  dat <- data.frame(
    resp_error = rimm(n = 5),
    Item2_rel = 2,
    Item3_rel = -1.5,
    spaD2 = 0.5,
    spaD3 = 2
  )

  # two-parameter model mock fit
  f <- bmf(kappa ~ 1, thetat ~ 1)
  mock_fit <- bmm(f, dat, mixture2p(resp_error = "resp_error"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")

  # three-parameter model mock fit
  f <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  model <- mixture3p(
    resp_error = "resp_error", set_size = 3,
    nt_features = paste0("Item", 2:3, "_rel")
  )
  mock_fit <- bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")

  # imm_abc model mock fit
  f <- bmf(kappa ~ 1, c ~ 1, a ~ 1)
  model <- imm(
    resp_error = "resp_error", set_size = 3,
    nt_features = paste0("Item", 2:3, "_rel"),
    version = "abc"
  )
  mock_fit <- bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")

  # imm_bsc model mock fit
  f <- bmf(kappa ~ 1, c ~ 1, s ~ 1)
  model <- imm(
    resp_error = "resp_error", set_size = 3,
    nt_features = paste0("Item", 2:3, "_rel"),
    nt_distances = paste0("spaD", 2:3),
    version = "bsc"
  )
  mock_fit <- bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")

  # imm_full model mock fit
  f <- bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1)
  model <- imm(
    resp_error = "resp_error", set_size = 3,
    nt_features = paste0("Item", 2:3, "_rel"),
    nt_distances = paste0("spaD", 2:3)
  )
  mock_fit <- bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")
})

test_that("Available models produce expected errors", {
  withr::local_options("bmm.silent" = 2)
  skip_on_cran()
  dat <- data.frame(
    resp_error = rimm(n = 5),
    Item2_rel = 2,
    Item3_rel = -1.5,
    spaD2 = -0.5,
    spaD3 = 2
  )

  okmodels <- c("mixture3p", "imm")
  for (model in okmodels) {
    model1 <- get_model(model)(resp_error = "resp_error",
      nt_features = "Item2_rel",
      set_size = 5,
      nt_distances = "spaD2")
    expect_error(
      bmm(bmf(kappa ~ 1), dat, model1,
        backend = "mock",
        mock_fit = 1, rename = FALSE
      ),
      "'nt_features' should equal max\\(set_size\\)-1"
    )

    model2 <- get_model(model)(resp_error = "resp_error",
      nt_features = "Item2_rel",
      set_size = TRUE,
      nt_distances = "spaD2")
    expect_error(
      bmm(bmf(kappa ~ 1), dat, model2,
        backend = "mock",
        mock_fit = 1, rename = FALSE
      ),
      "must be either a variable in your data or "
    )
  }

  for (version in c("bsc", "full")) {
    model1 <- imm(
      resp_error = "resp_error",
      nt_features = paste0("Item", 2:3, "_rel"),
      set_size = 3,
      nt_distances = paste0("spaD", 2:3),
      version = version
    )
    expect_error(
      bmm(bmf(kappa ~ 1), dat, model1,
        backend = "mock",
        mock_fit = 1, rename = FALSE
      ),
      "All non-target distances to the target need to be postive."
    )
  }
})

test_that("bmm() starts the step-size search at the package default unless the user sets it", {
  withr::local_options(bmm.step_size = 0.02)
  dat <- oberauer_lin_2017
  formula <- bmf(c ~ 1 + (1 | ID), kappa ~ 1 + (1 | ID))
  mock <- function(...) {
    bmm(formula, dat, sdm("dev_rad"), backend = "mock", mock_fit = 1, rename = FALSE, ...)
  }
  expect_equal(mock()$stan_args$control, list(step_size = 0.02))
  expect_equal(
    mock(control = list(adapt_delta = 0.95))$stan_args$control,
    list(adapt_delta = 0.95, step_size = 0.02)
  )
  expect_equal(mock(control = list(step_size = 0.5))$stan_args$control, list(step_size = 0.5))
  expect_equal(mock(threads = brms::threading(2))$stan_args$control, list(step_size = 0.02))
  # only the sampler has a step size; cmdstanr's other methods reject the argument
  expect_null(mock(algorithm = "meanfield")$stan_args$control)

  withr::local_options(bmm.step_size = FALSE)
  expect_null(mock()$stan_args$control)
})

test_that("bmm() builds the default prior and the inits from the data2 brm() gets", {
  dat <- oberauer_lin_2017
  ids <- levels(factor(dat$ID))
  A <- diag(length(ids))
  dimnames(A) <- list(ids, ids)
  formula <- bmf(kappa ~ 1 + (1 | gr(ID, cov = A)), thetat ~ 1)
  mock <- function() {
    bmm(formula, dat, mixture2p("dev_rad"),
      data2 = list(A = A), backend = "mock", mock_fit = 1, rename = FALSE
    )
  }
  expect_true(is.function(mock()$stan_args$init))
  expect_true("sd_1" %in% names(mock()$stan_args$init()))
  expect_match(stancode(formula, dat, mixture2p("dev_rad"), data2 = list(A = A)), "Lcov_1")
  expect_s3_class(default_prior(formula, dat, mixture2p("dev_rad"), data2 = list(A = A)), "brmsprior")

  withr::local_options(bmm.default_priors = FALSE)
  expect_true(is.function(mock()$stan_args$init))
})
