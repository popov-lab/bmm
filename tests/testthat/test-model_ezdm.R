# Test EZDM model specification and integration

test_that("ezdm model can be created with both versions", {
  expect_silent(ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par"))
  expect_silent(ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "4par"))
})

test_that("ezdm model has correct class structure", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "ezdm")
  expect_s3_class(model, "ezdm_3par")
})

test_that("ezdm model parameters are correctly defined for 3par version", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  expect_true("drift" %in% names(model$parameters))
  expect_true("bound" %in% names(model$parameters))
  expect_true("ndt" %in% names(model$parameters))
  expect_true("s" %in% names(model$parameters))
  # s is fixed to 0 (will be exponentiated to 1 in Stan)
  expect_equal(model$fixed_parameters$s, 0)
})

test_that("ezdm model parameters are correctly defined for 4par version", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "4par")
  expect_true("drift" %in% names(model$parameters))
  expect_true("bound" %in% names(model$parameters))
  expect_true("ndt" %in% names(model$parameters))
  expect_true("zr" %in% names(model$parameters))
  expect_true("s" %in% names(model$parameters))
  # s is fixed to 0 (will be exponentiated to 1 in Stan)
  expect_equal(model$fixed_parameters$s, 0)
})

test_that("ezdm model has correct link functions", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "4par")
  expect_equal(model$links$drift, "log")
  expect_equal(model$links$bound, "log")
  expect_equal(model$links$ndt, "log")
  expect_equal(model$links$zr, "logit")
  expect_equal(model$links$s, "log")
})

test_that("ezdm model accepts custom links", {
  custom_links <- list(drift = "identity")
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par", links = custom_links)
  expect_equal(model$links$drift, "identity")
  expect_equal(model$links$bound, "log") # unchanged
})

test_that("ezdm check_data validates mean_rt variable", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Valid data
  valid_data <- data.frame(
    mean_rt = c(0.5, 0.6, 0.7),
    var_rt = c(0.02, 0.03, 0.025),
    n_upper = c(80, 85, 75),
    n_trials = c(100, 100, 100)
  )
  expect_silent(check_data(model, valid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)))

  # Negative mean RTs should error
  invalid_data <- data.frame(
    mean_rt = c(-0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(80, 85),
    n_trials = c(100, 100)
  )
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "Mean RT values must be positive"
  )

  # Mean RT values > 10 should warn (likely milliseconds)
  ms_data <- data.frame(
    mean_rt = c(500, 600),
    var_rt = c(2000, 3000),
    n_upper = c(80, 85),
    n_trials = c(100, 100)
  )
  expect_warning(
    check_data(model, ms_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "milliseconds"
  )
})

test_that("ezdm check_data validates var_rt variable", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Negative or zero variance should error
  invalid_data <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(-0.02, 0.03),
    n_upper = c(80, 85),
    n_trials = c(100, 100)
  )
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "Variance of RT must be positive"
  )

  invalid_data2 <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0, 0.03),
    n_upper = c(80, 85),
    n_trials = c(100, 100)
  )
  expect_error(
    check_data(model, invalid_data2, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "Variance of RT must be positive"
  )
})

test_that("ezdm check_data validates n_trials variable", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Zero n_trials should error
  invalid_data <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(80, 85),
    n_trials = c(0, 100)
  )
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "must larger than two"
  )

  # n_trials = 1 should error (must be > 2)
  invalid_data_1 <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(1, 85),
    n_trials = c(1, 100)
  )
  expect_error(
    check_data(model, invalid_data_1, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "must larger than two"
  )

  # n_trials = 2 should error (must be larger than 2)
  invalid_data_2 <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(1, 85),
    n_trials = c(2, 100)
  )
  expect_error(
    check_data(model, invalid_data_2, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "must larger than two"
  )

  # Non-integer n_trials should warn
  invalid_data_nonint <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(80, 85),
    n_trials = c(100.5, 100)
  )
  expect_warning(
    check_data(model, invalid_data_nonint, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "whole numbers"
  )
})

test_that("ezdm check_data validates n_upper variable", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Negative n_upper should error
  invalid_data <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(-5, 85),
    n_trials = c(100, 100)
  )
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "needs to be positive"
  )

  # n_upper > n_trials should error
  invalid_data2 <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(120, 85),
    n_trials = c(100, 100)
  )
  expect_error(
    check_data(model, invalid_data2, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "cannot exceed total trials"
  )

  # Non-integer n_upper should warn
  invalid_data3 <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(80.5, 85),
    n_trials = c(100, 100)
  )
  expect_warning(
    check_data(model, invalid_data3, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "whole numbers"
  )
})

test_that("ezdm check_data handles missing values", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Missing values in required variables should not cause errors
  # (they may be handled by brms later in the fitting process)
  data_with_na <- data.frame(
    mean_rt = c(0.5, NA, 0.7),
    var_rt = c(0.02, 0.03, 0.025),
    n_upper = c(80, 85, 75),
    n_trials = c(100, 100, 100)
  )
  # check_data should complete without error
  expect_silent(
    check_data(model, data_with_na, bmf(drift ~ 1, bound ~ 1, ndt ~ 1))
  )
})

test_that("ezdm works with mock backend - 3par version", {
  skip_on_cran()

  # Simulate summary statistics
  sim_data <- rezdm(10, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3, version = "3par")
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ezdm works with mock backend - 4par version", {
  skip_on_cran()

  # For 4par, create data with separate upper/lower variables to match model expectations
  sim_data <- data.frame(
    mean_rt_upper = runif(10, 0.4, 0.6),
    mean_rt_lower = runif(10, 0.5, 0.7),
    var_rt_upper = runif(10, 0.01, 0.05),
    var_rt_lower = runif(10, 0.01, 0.05),
    n_upper = sample(30:70, 10, replace = TRUE),
    n_trials = rep(100, 10)
  )

  model <- ezdm(c("mean_rt_upper", "mean_rt_lower"), c("var_rt_upper", "var_rt_lower"), "n_upper", "n_trials", version = "4par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ezdm formula conversion works correctly for 3par", {
  skip_on_cran()

  sim_data <- rezdm(10, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3, version = "3par")
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)

  # Check that formula was converted properly
  expect_s3_class(fit$formula, "brmsformula")
  expect_s3_class(fit$formula$family, "customfamily")
  expect_equal(fit$formula$family$name, "ezdm_3par")
})

test_that("ezdm formula conversion works correctly for 4par", {
  skip_on_cran()

  # For 4par with separate upper/lower variables
  sim_data <- data.frame(
    mean_rt_upper = c(0.5, 0.55, 0.52, 0.48, 0.51),
    mean_rt_lower = c(0.55, 0.60, 0.57, 0.53, 0.56),
    var_rt_upper = c(0.02, 0.025, 0.022, 0.019, 0.021),
    var_rt_lower = c(0.024, 0.030, 0.026, 0.023, 0.025),
    n_upper = c(80, 85, 82, 78, 81),
    n_trials = c(100, 100, 100, 100, 100)
  )

  model <- ezdm(
    mean_rt = c("mean_rt_upper", "mean_rt_lower"),
    var_rt = c("var_rt_upper", "var_rt_lower"),
    n_upper = "n_upper",
    n_trials = "n_trials",
    version = "4par"
  )
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)

  # Check that formula was converted properly
  expect_s3_class(fit$formula, "brmsformula")
  expect_s3_class(fit$formula$family, "customfamily")
  expect_equal(fit$formula$family$name, "ezdm_4par")
})

test_that("ezdm with condition effects works", {
  skip_on_cran()

  # Simulate data with condition effects
  n_per_cond <- 10
  data_a <- rezdm(n_per_cond, n_trials = 100, drift = 2.5, bound = 1.5, ndt = 0.3, version = "3par")
  data_a$condition <- "A"
  data_b <- rezdm(n_per_cond, n_trials = 100, drift = 1.5, bound = 1.5, ndt = 0.3, version = "3par")
  data_b$condition <- "B"
  sim_data <- rbind(data_a, data_b)
  sim_data$condition <- factor(sim_data$condition)

  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  formula <- bmf(drift ~ 0 + condition, bound ~ 1, ndt ~ 1)

  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ezdm with hierarchical structure works", {
  skip_on_cran()

  # Simulate hierarchical data
  n_subjects <- 3
  n_per_subject <- 5

  data_list <- lapply(1:n_subjects, function(i) {
    d <- rezdm(n_per_subject,
      n_trials = 100,
      drift = rnorm(1, 2, 0.3),
      bound = 1.5,
      ndt = 0.3,
      version = "3par"
    )
    d$id <- paste0("S", i)
    d
  })
  sim_data <- do.call(rbind, data_list)
  sim_data$id <- factor(sim_data$id)

  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1, ndt ~ 1)

  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ezdm allows missing parameters with message", {
  skip_on_cran()

  sim_data <- rezdm(10, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3, version = "3par")
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Missing ndt parameter should work with message (not error)
  formula_incomplete <- bmf(drift ~ 1, bound ~ 1)
  expect_message(
    bmm(formula_incomplete, sim_data, model, backend = "mock", mock = 1, rename = FALSE),
    "No formula for parameter ndt"
  )
})

test_that("ezdm default priors are correctly set for 3par", {
  skip_on_cran()

  sim_data <- rezdm(10, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3, version = "3par")
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)
  prior_summary <- brms::prior_summary(fit)

  # Check that priors are set for main parameters (in dpar column)
  expect_true(any(grepl("drift", prior_summary$dpar)))
  expect_true(any(grepl("bound", prior_summary$dpar)))
  expect_true(any(grepl("ndt", prior_summary$dpar)))
  expect_true(any(grepl("^s$", prior_summary$dpar)))
})

test_that("ezdm default priors are correctly set for 4par", {
  skip_on_cran()

  # For 4par with separate upper/lower variables
  sim_data <- data.frame(
    mean_rt_upper = c(0.5, 0.55, 0.52, 0.48, 0.51),
    mean_rt_lower = c(0.55, 0.60, 0.57, 0.53, 0.56),
    var_rt_upper = c(0.02, 0.025, 0.022, 0.019, 0.021),
    var_rt_lower = c(0.024, 0.030, 0.026, 0.023, 0.025),
    n_upper = c(80, 85, 82, 78, 81),
    n_trials = c(100, 100, 100, 100, 100)
  )

  model <- ezdm(
    mean_rt = c("mean_rt_upper", "mean_rt_lower"),
    var_rt = c("var_rt_upper", "var_rt_lower"),
    n_upper = "n_upper",
    n_trials = "n_trials",
    version = "4par"
  )
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)
  prior_summary <- brms::prior_summary(fit)

  # Check that priors are set for all parameters (in dpar column)
  expect_true(any(grepl("drift", prior_summary$dpar)))
  expect_true(any(grepl("bound", prior_summary$dpar)))
  expect_true(any(grepl("ndt", prior_summary$dpar)))
  expect_true(any(grepl("zr", prior_summary$dpar)))
  expect_true(any(grepl("^s$", prior_summary$dpar)))
})

test_that("ezdm stanvars are correctly added", {
  skip_on_cran()

  sim_data <- rezdm(10, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3, version = "3par")
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)

  # Check that custom Stan functions were added
  expect_true(!is.null(fit$stanvars))
  expect_s3_class(fit$stanvars, "stanvars")
})

test_that("ezdm 3par requires single mean_rt and var_rt variables", {
  model_3par <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Check that resp_vars has correct structure
  expect_length(model_3par$resp_vars$mean_rt, 1)
  expect_length(model_3par$resp_vars$var_rt, 1)
  expect_equal(model_3par$resp_vars$mean_rt, "mean_rt")
  expect_equal(model_3par$resp_vars$var_rt, "var_rt")
})

test_that("ezdm 4par can accept vector or scalar for mean_rt and var_rt", {
  # Single variable (same for both boundaries)
  model_4par_single <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "4par")
  expect_length(model_4par_single$resp_vars$mean_rt, 1)
  expect_length(model_4par_single$resp_vars$var_rt, 1)

  # Separate variables for upper/lower boundaries
  model_4par_vector <- ezdm(
    mean_rt = c("mean_rt_upper", "mean_rt_lower"),
    var_rt = c("var_rt_upper", "var_rt_lower"),
    n_upper = "n_upper",
    n_trials = "n_trials",
    version = "4par"
  )
  expect_length(model_4par_vector$resp_vars$mean_rt, 2)
  expect_length(model_4par_vector$resp_vars$var_rt, 2)
  expect_equal(model_4par_vector$resp_vars$mean_rt, c("mean_rt_upper", "mean_rt_lower"))
})

test_that("ezdm model has void_mu flag set", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  expect_true(model$void_mu)
})

test_that("ezdm check_data validates all required variables exist", {
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  # Missing required variable
  incomplete_data <- data.frame(
    mean_rt = c(0.5, 0.6),
    var_rt = c(0.02, 0.03),
    n_upper = c(80, 85)
    # missing n_trials
  )

  expect_error(
    check_data(model, incomplete_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "missing from the data"
  )
})

test_that("ezdm posterior_predict function is defined for 3par", {
  # Test that the posterior_predict function exists and has correct structure
  # Note: Cannot test with mock backend as it doesn't create proper brmsfit objects
  # This would need actual sampling which is too slow for regular tests

  # Just verify the function is exported and callable
  expect_true(exists("posterior_predict_ezdm_3par", where = asNamespace("bmm")))
})

test_that("ezdm posterior_predict function is defined for 4par", {
  # Test that the posterior_predict function exists and has correct structure
  # Note: Cannot test with mock backend as it doesn't create proper brmsfit objects
  # This would need actual sampling which is too slow for regular tests

  # Just verify the function is exported and callable
  expect_true(exists("posterior_predict_ezdm_4par", where = asNamespace("bmm")))
})
