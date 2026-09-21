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
  expect_equal(model$links$drift, "identity")  # Changed to identity to allow negative drift
  expect_equal(model$links$bound, "log")
  expect_equal(model$links$ndt, "log")
  expect_equal(model$links$zr, "logit")
  expect_equal(model$links$s, "log")
})

test_that("ezdm model accepts custom links", {
  custom_links <- list(bound = "identity")  # Changed from drift since identity is now default
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par", links = custom_links)
  expect_equal(model$links$drift, "identity")  # default
  expect_equal(model$links$bound, "identity")  # custom
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
    "must be larger than two"
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
    "must be larger than two"
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
    "must be larger than two"
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


# R <-> Stan parity of the likelihood ------------------------------------------

# Drives the installed Stan chunks through a fixed_param generated-quantities
# run rather than cmdstanr::expose_functions(), which compiles through Rcpp and
# links against RcppParallel's libtbb -- that fails inside the test suite after
# the dependency chain is loaded, and silently skips, which is the worst failure
# mode for a parity test. sig_figs = 17 round-trips a double; the default 6
# would measure cmdstan's output precision instead of the code.
ezdm_stan_lpdf <- function(version, data) {
  sc_path <- system.file("stan_chunks", package = "bmm")
  declarations <- if (version == "3par") {
    "vector[N] mean_rt; vector[N] var_rt; vector[N] s;"
  } else {
    paste(
      "vector[N] mean_rt_upper; vector[N] mean_rt_lower;",
      "vector[N] var_rt_upper; vector[N] var_rt_lower;",
      "vector[N] zr; vector[N] s;"
    )
  }
  call <- if (version == "3par") {
    "ezdm_3par_lpdf(mean_rt[i] | 0.0, drift[i], bound[i], ndt[i], s[i],
                    var_rt[i], n_upper[i], n_trials[i])"
  } else {
    "ezdm_4par_lpdf(mean_rt_upper[i] | 0.0, drift[i], bound[i], ndt[i], zr[i],
                    s[i], mean_rt_lower[i], var_rt_upper[i], var_rt_lower[i],
                    n_upper[i], n_trials[i])"
  }

  program <- paste0(
    "functions {\n",
    read_lines2(file.path(sc_path, "ezdm_cumulants.stan")), "\n",
    read_lines2(file.path(sc_path, paste0("ezdm_", version, "_functions.stan"))),
    "\n}\ndata {\n  int<lower=1> N;\n  vector[N] drift; vector[N] bound;",
    " vector[N] ndt;\n  array[N] int n_upper; array[N] int n_trials;\n  ",
    declarations,
    "\n}\ngenerated quantities {\n  vector[N] stan_lpdf;\n",
    "  for (i in 1:N) {\n    stan_lpdf[i] = ", call, ";\n  }\n}\n"
  )

  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
    iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  # base-R read: cmdstanr's reader misparses long plain-decimal numbers
  csv <- utils::read.csv(
    fit$output_files()[1],
    comment.char = "#", check.names = FALSE
  )
  as.numeric(csv[1, paste0("stan_lpdf.", seq_len(data$N))])
}

# drift spans the series branch (t < 0.7), the closed forms, the t > 30
# saturation, and the values at which the old code returned NaN; negative
# drift takes the negative-drift branch of ezdm_logit_pc, and 2e-5 its series
ezdm_parity_drift <- c(
  -1500, -5, -0.2, 0, 1e-300, 1e-8, 2e-5, 0.001, 0.05, 0.2, 0.5, 1, 2, 5, 20, 1500
)

test_that("ezdm_3par_lpdf in Stan matches dezdm() in R", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  grid <- expand.grid(
    drift = ezdm_parity_drift, bound = c(0.4, 1.5, 3), s = c(1, 1.4),
    n_trials = c(3L, 500L)
  )
  grid$ndt <- 0.25
  moments <- .ezdm_moments_3par(grid$drift, grid$bound, grid$s)
  # keep the summaries near the implied moments so the lpdf stays moderate
  grid$mean_rt <- grid$ndt + moments$MDT * 1.02
  grid$var_rt <- moments$VRT * 0.9
  grid$n_upper <- pmin(grid$n_trials, round(grid$n_trials * 0.6))

  r_lpdf <- dezdm(
    mean_rt = grid$mean_rt, var_rt = grid$var_rt, n_upper = grid$n_upper,
    n_trials = grid$n_trials, drift = grid$drift, bound = grid$bound,
    ndt = grid$ndt, s = grid$s, version = "3par"
  )
  stan_lpdf <- ezdm_stan_lpdf("3par", c(
    list(N = nrow(grid)),
    as.list(grid[c(
      "mean_rt", "var_rt", "drift", "bound", "ndt", "s", "n_upper", "n_trials"
    )])
  ))

  # the binomial is evaluated on the logit scale, so the density stays finite
  # even at |drift| = 1500, where pC rounds to 1
  expect_true(all(is.finite(r_lpdf)))
  expect_true(all(is.finite(stan_lpdf)))
  expect_lt(max(abs(stan_lpdf - r_lpdf) / abs(r_lpdf)), 1e-9)
})

test_that("ezdm_4par_lpdf in Stan matches dezdm() in R", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  grid <- expand.grid(
    drift = ezdm_parity_drift, bound = c(0.4, 1.5, 3),
    zr = c(0.2, 0.5, 0.7, 0.95), s = c(0.5, 1, 2), n_upper = c(0L, 1L, 3L, 30L)
  )
  grid$n_trials <- 30L
  grid$ndt <- 0.25
  # n_upper 0 and 1 exercise the lower gate, 30 and 29 the upper one
  grid <- rbind(grid, transform(grid[grid$n_upper == 30L, ], n_upper = 29L))
  moments <- .ezdm_moments_4par(grid$drift, grid$bound, grid$zr, grid$s)
  grid$mean_rt_upper <- grid$ndt + moments$mdt_upper * 1.02
  grid$mean_rt_lower <- grid$ndt + moments$mdt_lower * 0.98
  grid$var_rt_upper <- moments$vrt_upper * 0.9
  grid$var_rt_lower <- moments$vrt_lower * 1.1

  r_lpdf <- dezdm(
    mean_rt = as.matrix(grid[c("mean_rt_upper", "mean_rt_lower")]),
    var_rt = as.matrix(grid[c("var_rt_upper", "var_rt_lower")]),
    n_upper = grid$n_upper, n_trials = grid$n_trials, drift = grid$drift,
    bound = grid$bound, ndt = grid$ndt, zr = grid$zr, s = grid$s,
    version = "4par"
  )
  stan_lpdf <- ezdm_stan_lpdf("4par", c(
    list(N = nrow(grid)),
    as.list(grid[c(
      "mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower",
      "drift", "bound", "ndt", "zr", "s", "n_upper", "n_trials"
    )])
  ))

  expect_false(any(is.nan(r_lpdf)))
  expect_false(any(is.nan(stan_lpdf)))
  expect_equal(is.finite(stan_lpdf), is.finite(r_lpdf))

  finite <- is.finite(r_lpdf)
  expect_lt(max(abs(stan_lpdf[finite] - r_lpdf[finite]) / abs(r_lpdf[finite])), 1e-9)
})

test_that("ezdm Stan code calls the model lpdf and parses", {
  skip_on_cran()
  data <- data.frame(
    mean_rt = c(0.5, 0.55), var_rt = c(0.02, 0.03), n_upper = c(60L, 40L),
    n_trials = c(100L, 100L), mean_rt_upper = c(0.5, 0.55),
    mean_rt_lower = c(0.6, 0.58), var_rt_upper = c(0.02, 0.03),
    var_rt_lower = c(0.03, 0.04)
  )

  code3 <- stancode(
    bmf(drift ~ 1, bound ~ 1, ndt ~ 1), data = data,
    model = ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  )
  expect_match(code3, "real ezdm_3par_lpdf(", fixed = TRUE)
  expect_match(code3, "real ezdm_symmetric_lpdf(", fixed = TRUE)
  expect_match(code3, "ezdm_3par_lpdf(Y[n] |", fixed = TRUE)

  code4 <- stancode(
    bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1), data = data,
    model = ezdm(
      c("mean_rt_upper", "mean_rt_lower"), c("var_rt_upper", "var_rt_lower"),
      "n_upper", "n_trials",
      version = "4par"
    )
  )
  expect_match(code4, "real ezdm_4par_lpdf(", fixed = TRUE)
  expect_match(code4, "real ezdm_boundary_lpdf(", fixed = TRUE)

  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  for (code in list(code3, code4)) {
    model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
    expect_true(model$check_syntax(quiet = TRUE))
  }
})

test_that("ezdm_4par_lpdf has the right drift gradient at and around zero drift", {
  # Returning the limit of pC at drift = 0 as a constant gave the right value
  # and a zero gradient, which only autodiff can see. CmdStan's standalone
  # log_prob method is used because expose_functions() returns no gradients.
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  chunk <- function(file) {
    paste(readLines(system.file("stan_chunks", file, package = "bmm")), collapse = "\n")
  }
  program <- paste0(
    "functions {\n", chunk("ezdm_cumulants.stan"), "\n", chunk("ezdm_4par_functions.stan"), "\n}\n",
    "data { int N; vector[2] mrt; vector[2] vrt; }\n",
    "parameters { vector[N] drift; }\n",
    "model { for (i in 1:N) target += ezdm_4par_lpdf(mrt[1] | 0.0, drift[i], 1.5, 0.25, 0.3, 1.0,\n",
    "  mrt[2], vrt[1], vrt[2], 12, 60); }\n"
  )
  # summaries of a cell with negative drift: most responses at the lower boundary
  moments <- .ezdm_moments_4par(-1.2, 1.5, 0.3, 1)
  mrt <- 0.25 + c(moments$mdt_upper, moments$mdt_lower)
  vrt <- c(moments$vrt_upper, moments$vrt_lower)
  drift <- c(-1.2, -1e-12, 0, 1e-12, 1.2)

  dir <- withr::local_tempdir()
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program, dir = dir), quiet = TRUE)
  cmdstanr::write_stan_json(list(N = length(drift), mrt = mrt, vrt = vrt), file.path(dir, "data.json"))
  cmdstanr::write_stan_json(list(drift = drift), file.path(dir, "pars.json"))
  status <- system2(model$exe_file(), c(
    "method=log_prob", paste0("constrained_params=", file.path(dir, "pars.json")), "jacobian=0",
    "data", paste0("file=", file.path(dir, "data.json")),
    "output", paste0("file=", file.path(dir, "out.csv")), "sig_figs=17"
  ), stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0)
  stan_gradient <- as.numeric(utils::read.csv(file.path(dir, "out.csv"), comment.char = "#")[1, -1])

  lpdf <- function(drift) {
    dezdm(mrt, vrt, n_upper = 12, n_trials = 60, drift = drift, bound = 1.5,
          ndt = 0.25, zr = 0.3, version = "4par")
  }
  r_gradient <- vapply(drift, \(d) (lpdf(d + 1e-5) - lpdf(d - 1e-5)) / 2e-5, numeric(1))
  expect_equal(stan_gradient / r_gradient, rep(1, length(drift)), tolerance = 1e-6)
})
