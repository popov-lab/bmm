# =============================================================================
# Tests for cswald model (model-specific tests)
# Distribution function tests are in test-distributions.R
# =============================================================================

# -----------------------------------------------------------------------------
# Model construction tests
# -----------------------------------------------------------------------------

test_that("cswald() creates model with correct structure", {
  model <- cswald(rt = "rt", response = "response", version = "simple")

  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "cswald")
  expect_s3_class(model, "cswald_simple")
  expect_equal(model$resp_vars$rt, "rt")
  expect_equal(model$resp_vars$response, "response")
  expect_equal(model$version, "simple")
})

test_that("cswald() creates crisk version correctly", {
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  expect_s3_class(model, "cswald_crisk")
  expect_equal(model$version, "crisk")
  expect_true("zr" %in% names(model$parameters))
})

test_that("cswald simple version has correct parameters", {
  model <- cswald(rt = "rt", response = "response", version = "simple")

  expect_true(all(c("drift", "bound", "ndt", "s") %in% names(model$parameters)))
  expect_false("zr" %in% names(model$parameters))
})

test_that("cswald crisk version has correct parameters", {
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  expect_true(all(c("drift", "bound", "ndt", "zr", "s") %in% names(model$parameters)))
})

test_that("cswald accepts custom links", {
  model <- cswald(rt = "rt", response = "response",
                  links = list(drift = "identity"), version = "simple")
  expect_equal(model$links$drift, "identity")
})

# -----------------------------------------------------------------------------
# Data validation tests (check_data.cswald)
# -----------------------------------------------------------------------------

test_that("check_data.cswald errors when required variables missing", {
  model <- cswald(rt = "rt", response = "response")

  expect_error(
    check_data(model, data.frame(x = 1), bmf(drift ~ 1)),
    "RT variable 'rt' is not present"
  )

  expect_error(
    check_data(model, data.frame(rt = 1), bmf(drift ~ 1)),
    "response variable 'response' is not present"
  )
})

test_that("check_data.cswald errors when RT contains NA", {
  model <- cswald(rt = "rt", response = "response")
  dat <- data.frame(rt = c(0.5, NA, 0.8), response = c(1, 1, 0))

  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "RT variable 'rt' contains.*NA"
  )
})

test_that("check_data.cswald errors when response contains NA", {
  model <- cswald(rt = "rt", response = "response")
  dat <- data.frame(rt = c(0.5, 0.6, 0.8), response = c(1, NA, 0))

  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "response variable 'response' contains.*NA"
  )
})

test_that("check_data.cswald errors when RT contains negative values", {
  model <- cswald(rt = "rt", response = "response")
  dat <- data.frame(rt = c(-0.5, 0.6, 0.8), response = c(1, 1, 0))

  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "reaction times are lower than zero"
  )
})

test_that("check_data.cswald warns when RT > 10 seconds", {
  # Use crisk version and large dataset to avoid other warnings
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(
    rt = c(runif(99, 0.4, 1.5), 15),
    response = rep(c(1, 0), 50)
  )

  expect_warning(
    check_data(model, dat, bmf(drift ~ 1)),
    "larger than 10 seconds"
  )
})

test_that("check_data.cswald warns when RT < 0.1 seconds", {
  # Use crisk version and large dataset to avoid other warnings
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(
    rt = c(0.05, runif(99, 0.4, 1.5)),
    response = rep(c(1, 0), 50)
  )

  expect_warning(
    check_data(model, dat, bmf(drift ~ 1)),
    "smaller than 0.100 seconds"
  )
})

test_that("check_data.cswald handles different response formats", {
  # Use crisk version to avoid error rate warning, and use larger datasets
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  # integer 0/1 - should work silently
  dat_int <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c(1L, 0L), 50)
  )
  expect_silent(check_data(model, dat_int, bmf(drift ~ 1)))

  # logical - should warn and convert
  dat_logical <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c(TRUE, FALSE), 50)
  )
  expect_warning(
    check_data(model, dat_logical, bmf(drift ~ 1)),
    "boolean"
  )

  # character upper/lower - should warn and convert
  dat_char <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c("upper", "lower"), 50)
  )
  expect_warning(
    check_data(model, dat_char, bmf(drift ~ 1)),
    "character"
  )

  # factor - should warn and convert
  dat_factor <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = factor(rep(c("upper", "lower"), 50))
  )
  expect_warning(
    check_data(model, dat_factor, bmf(drift ~ 1)),
    "character"
  )
})

test_that("check_data.cswald errors on invalid response values", {
  # Use crisk version to avoid error rate warning
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  # invalid integer values
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep(1, 50), rep(2, 50))
  )
  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "values other than 0 and 1"
  )

  # invalid character values
  dat_char <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep("upper", 50), rep("invalid", 50))
  )
  expect_error(
    check_data(model, dat_char, bmf(drift ~ 1)),
    "invalid character values"
  )
})

test_that("check_data.cswald warns about high error rate for simple version", {
  model <- cswald(rt = "rt", response = "response", version = "simple")

  # 30% error rate
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep(0, 30), rep(1, 70))
  )

  expect_warning(
    check_data(model, dat, bmf(drift ~ 1)),
    "error rate"
  )
})

test_that("check_data.cswald does not warn about error rate for crisk version", {
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  # 30% error rate - no warning for crisk
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep(0, 30), rep(1, 70))
  )

  expect_silent(check_data(model, dat, bmf(drift ~ 1)))
})

test_that("check_data.cswald warns about small sample size", {
  # Use crisk version to avoid error rate warning
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(rt = runif(30, 0.4, 1.5), response = rep(c(0, 1), 15))

  expect_warning(
    check_data(model, dat, bmf(drift ~ 1)),
    "only.*observations"
  )
})

test_that("check_data.cswald returns a data.frame", {
  # Use crisk version to avoid error rate warning
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(rt = runif(100, 0.4, 1.5), response = rep(c(0, 1), 50))

  result <- check_data(model, dat, bmf(drift ~ 1))
  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Formula conversion tests (bmf2bf.cswald)
# -----------------------------------------------------------------------------

test_that("bmf2bf.cswald creates correct brms formula", {
  model <- cswald(rt = "rt", response = "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  bf <- bmf2bf(model, formula)

  expect_s3_class(bf, "brmsformula")
  # check that response includes dec() term
  expect_true(grepl("dec", deparse(bf$formula)))
})

# -----------------------------------------------------------------------------
# Model configuration tests
# -----------------------------------------------------------------------------

test_that("configure_model.cswald_simple returns correct components", {
  skip_on_cran()

  model <- cswald(rt = "rt", response = "response", version = "simple")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = sample(c(0, 1), 100, replace = TRUE)
  )
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  config <- configure_model(model, dat, formula)

  expect_true(all(c("formula", "data", "stanvars") %in% names(config)))
  expect_s3_class(config$formula, "brmsformula")
  expect_s3_class(config$formula$family, "customfamily")
  expect_equal(config$formula$family$name, "cswald")
})

test_that("configure_model.cswald_crisk returns correct components", {
  skip_on_cran()

  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = sample(c(0, 1), 100, replace = TRUE)
  )
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  config <- configure_model(model, dat, formula)

  expect_true(all(c("formula", "data", "stanvars") %in% names(config)))
  expect_s3_class(config$formula, "brmsformula")
  expect_s3_class(config$formula$family, "customfamily")
  expect_equal(config$formula$family$name, "cswald_crisk")
})

# -----------------------------------------------------------------------------
# Integration tests with mock backend
# -----------------------------------------------------------------------------

test_that("cswald simple version runs with mock backend", {
  skip_on_cran()

  dat <- rcswald(n = 100, drift = 2, bound = 1.5, ndt = 0.3)
  model <- cswald(rt = "rt", response = "response", version = "simple")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald crisk version runs with mock backend", {
  skip_on_cran()

  dat <- rcswald(n = 100, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5)
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald simple version with predictors runs with mock backend", {
  skip_on_cran()

  dat <- rcswald(n = 200, drift = 2, bound = 1.5, ndt = 0.3)
  dat$condition <- rep(c("A", "B"), each = 100)
  model <- cswald(rt = "rt", response = "response", version = "simple")
  formula <- bmf(drift ~ condition, bound ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald crisk version allows negative drift (more lower responses)", {
  skip_on_cran()

  # Generate data with negative drift (biased toward lower boundary)
  dat <- rcswald(n = 100, drift = -1.5, bound = 1.5, ndt = 0.3, zr = 0.5)

  # Most responses should be lower (0) with negative drift
  expect_true(mean(dat$response == 0) > 0.5)

  model <- cswald(rt = "rt", response = "response", version = "crisk")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  # Model should run without error
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald handles all-correct responses", {
  skip_on_cran()

  # Data with all correct responses
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(1, 100)
  )

  model <- cswald(rt = "rt", response = "response", version = "simple")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  # Should run without error (though may warn about 0% error rate)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald handles high error rate data with crisk version", {
  skip_on_cran()

  # Data with 50% error rate
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c(0, 1), 50)
  )

  model <- cswald(rt = "rt", response = "response", version = "crisk")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})