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

test_that("bmm() stores the name of the data the user passed", {
  withr::local_options("bmm.silent" = 2)
  my_data <- data.frame(y = rimm(n = 5))
  fit <- bmm(bmf(kappa ~ 1, thetat ~ 1), my_data, mixture2p("y"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_equal(attr(fit$data, "data_name"), "my_data")

  # m3's check_data() rebuilds the data frame after the name is first stored
  counts <- data.frame(corr = rpois(5, 10), other = rpois(5, 3), npl = rpois(5, 2))
  fit <- bmm(bmf(c ~ 1, a ~ 1), counts,
    m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5), version = "ss"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_equal(attr(fit$data, "data_name"), "counts")
})

test_that("check_data() methods see the name of the user's data", {
  seen <- NULL
  local_mocked_bindings(order_data_query = function(model, data, formula) {
    seen <<- attr(data, "data_name")
    data
  })
  my_sdm_data <- data.frame(y = rsdm(10))
  stancode(bmf(c ~ 1, kappa ~ 1), my_sdm_data, sdm("y"))
  expect_equal(seen, "my_sdm_data")

  my_matrix <- as.matrix(data.frame(y = rsdm(10)))
  seen <- NULL
  bmm(bmf(kappa ~ 1, c ~ 1), my_matrix, sdm("y"), backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(seen, "my_matrix")

  my_list <- list(y = rsdm(10))
  seen <- NULL
  bmm(bmf(kappa ~ 1, c ~ 1), my_list, sdm("y"), backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(seen, "my_list")
})

test_that("bmm() says when the formula is not a bmmformula", {
  dat <- data.frame(y = rimm(n = 5))
  for (f in list(y ~ 1, brms::bf(y ~ 1), list(kappa ~ 1), NULL, "kappa ~ 1")) {
    expect_error(
      bmm(f, dat, mixture2p("y"), backend = "mock", mock_fit = 1),
      "The provided formula is not a bmm formula"
    )
  }

  # problems with the model and the data are reported before the formula
  expect_error(
    bmm(y ~ 1, data.frame(z = 1:3), mixture2p("y"), backend = "mock", mock_fit = 1),
    "The response variable 'y' is not present in the data."
  )
  expect_error(
    bmm(y ~ 1, model = mixture2p("y"), backend = "mock", mock_fit = 1),
    "Data must be specified using the 'data' argument."
  )
})

test_that("bmm() says when data is not coercible to a data frame, whatever it is", {
  f <- bmf(kappa ~ 1, thetat ~ 1)
  my_env <- new.env()
  for (bad_data in list(sum, quote(abc), my_env)) {
    expect_error(
      bmm(f, bad_data, mixture2p("y"), backend = "mock", mock_fit = 1),
      "Argument 'data' must be coercible to a data.frame."
    )
  }
  expect_null(attr(my_env, "data_name"))
})

test_that("bmm() and the extractors say when data is missing", {
  f <- bmf(c ~ 1, kappa ~ 1)
  msg <- "Data must be specified using the 'data' argument."
  expect_error(bmm(f, model = sdm("y"), backend = "mock", mock_fit = 1), msg)
  expect_error(standata(f, model = sdm("y")), msg)
  expect_error(stancode(f, model = sdm("y")), msg)
  expect_error(default_prior(f, model = sdm("y")), msg)
})

test_that("bmm() stamps the bmm version into the Stan code of the fit", {
  fit <- bmm(bmf(c ~ 1, kappa ~ 1), data.frame(y = rsdm(10)), sdm("y"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_match(fit$model, paste0("brms [0-9.]+ and bmm ", utils::packageVersion("bmm")))
})

test_that("standata() does not configure the prior, which the Stan data does not use", {
  local_mocked_bindings(configure_prior = function(...) stop2("configure_prior() called"))
  expect_type(standata(bmf(c ~ 1, kappa ~ 1), data.frame(y = rsdm(10)), sdm("y")), "list")
})

test_that("stancode() and default_prior() do not build the inits, which they do not return", {
  local_mocked_bindings(create_initfun = function(...) stop2("create_initfun() called"))
  f <- bmf(c ~ 1, kappa ~ 1)
  dat <- data.frame(y = rsdm(10))
  expect_type(stancode(f, dat, sdm("y")), "character")
  expect_s3_class(default_prior(f, dat, sdm("y")), "brmsprior")
})
