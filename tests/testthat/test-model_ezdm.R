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
  expect_warning(
    model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par", links = custom_links),
    "allow values that the model's default"
  )
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

# Sparse boundaries in 4par (#430) ----------------------------------------------

ezdm4_model <- function() {
  ezdm(
    mean_rt = c("mean_rt_upper", "mean_rt_lower"),
    var_rt = c("var_rt_upper", "var_rt_lower"),
    n_upper = "n_upper", n_trials = "n_trials", version = "4par"
  )
}

ezdm4_formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

# every warning rather than the first matching one: the contract is exactly
# one warning, and a second one (brms dropping rows) is the bug itself
collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

test_that("ezdm 4par keeps the cells of a boundary below two responses (#430)", {
  dat <- withr::with_seed(1, rezdm(
    30, n_trials = 40, drift = 3, bound = 1.5, ndt = 0.25, zr = 0.6,
    version = "4par"
  ))
  n_lower <- dat$n_trials - dat$n_upper
  expect_identical(sum(n_lower < 2), 29L)

  out <- collect_warnings(standata(ezdm4_formula, dat, ezdm4_model()))
  sdata <- out$value

  expect_identical(sdata$N, 30L)
  expect_identical(as.vector(sdata$vint3), as.integer(dat$n_upper >= 2))
  expect_identical(as.vector(sdata$vint4), as.integer(n_lower >= 2))
  expect_true(all(sdata$vreal1[n_lower < 2] == -1))
  expect_true(all(sdata$vreal3[n_lower < 2] == -1))
  expect_length(out$warnings, 1L)
  expect_match(out$warnings, "29 of 30 cells", fixed = TRUE)
})

# rows: both boundaries used; no lower responses; one lower response with a
# mean but no variance (what ezdm_summary_stats() can return after its
# correction); enough upper responses but no upper summaries (min_trials);
# n_upper missing; no upper responses; n_trials missing
test_that("check_data() marks unused 4par boundaries and fills them with -1", {
  d0 <- data.frame(
    mean_rt_upper = c(0.50, 0.52, 0.48, NA, 0.51, NA, 0.49),
    mean_rt_lower = c(0.60, NA, 0.62, 0.58, NA, 0.61, 0.57),
    var_rt_upper = c(0.020, 0.022, 0.018, NA, 0.021, NA, 0.019),
    var_rt_lower = c(0.030, NA, NA, 0.028, 0.031, 0.032, 0.027),
    n_upper = c(30L, 50L, 49L, 5L, NA, 0L, 30L),
    n_trials = c(50L, 50L, 50L, 50L, 50L, 50L, NA)
  )

  out <- collect_warnings(check_data(ezdm4_model(), d0, ezdm4_formula))
  d1 <- out$value

  expect_identical(names(d1), c(names(d0), "rt_used_upper", "rt_used_lower"))
  expect_identical(d1$rt_used_upper, c(1L, 1L, 1L, 0L, NA, 0L, NA))
  expect_identical(d1$rt_used_lower, c(1L, 0L, 0L, 1L, NA, 1L, NA))
  expect_identical(d1$mean_rt_upper, c(0.50, 0.52, 0.48, -1, 0.51, -1, 0.49))
  expect_identical(d1$var_rt_upper, c(0.020, 0.022, 0.018, -1, 0.021, -1, 0.019))
  expect_identical(d1$mean_rt_lower, c(0.60, -1, -1, 0.58, NA, 0.61, 0.57))
  expect_identical(d1$var_rt_lower, c(0.030, -1, -1, 0.028, 0.031, 0.032, 0.027))
  expect_identical(d1$n_upper, d0$n_upper)
  expect_identical(d1$n_trials, d0$n_trials)
  expect_length(out$warnings, 1L)
  expect_match(out$warnings, "4 of 7 cells", fixed = TRUE)
})

# an unused boundary's summaries are never read, so implausible values there
# (a lower mean of 0 with no lower responses, an upper mean of 50 from a
# single response, a mean without a variance) must not stop the fit
test_that("check_data() validates only the summaries of used 4par boundaries", {
  d0 <- data.frame(
    mean_rt_upper = c(0.50, 0.52, 50, -5, 0.49, 0.47),
    mean_rt_lower = c(0.60, 0, 0.60, 0.58, 50, -5),
    var_rt_upper = c(0.020, 0.022, 2, NA, 0.019, 0.018),
    var_rt_lower = c(0.030, 0, 0.030, 0.028, NA, -1),
    n_upper = c(30L, 50L, 1L, 5L, 30L, 49L),
    n_trials = 50L
  )

  out <- collect_warnings(check_data(ezdm4_model(), d0, ezdm4_formula))
  expect_false(any(grepl("milliseconds", out$warnings)))
})

test_that("check_data() still rejects implausible summaries of used 4par boundaries", {
  base <- data.frame(
    mean_rt_upper = 0.5, mean_rt_lower = 0.6, var_rt_upper = 0.02,
    var_rt_lower = 0.03, n_upper = 30L, n_trials = 50L
  )
  model <- ezdm4_model()

  expect_error(
    check_data(model, transform(base, mean_rt_lower = 0), ezdm4_formula),
    "Mean RT values must be positive"
  )
  expect_error(
    check_data(model, transform(base, var_rt_upper = -0.01), ezdm4_formula),
    "Variance of RT must be positive"
  )
  expect_warning(
    check_data(model, transform(base, mean_rt_upper = 50), ezdm4_formula),
    "milliseconds"
  )
})

# update() re-runs check_data() on fit$data, which already holds placeholders
# and indicators. Row 4 is the case that needs the indicator: its -1
# placeholders sit at a boundary with 20 responses, so a second pass that
# looked at counts and NAs alone would take them for data
test_that("check_data() on its own 4par output changes nothing", {
  d0 <- data.frame(
    mean_rt_upper = c(0.50, NA, 0.48, NA),
    mean_rt_lower = c(0.60, 0.61, NA, 0.58),
    var_rt_upper = c(0.020, NA, 0.018, NA),
    var_rt_lower = c(0.030, 0.029, NA, 0.028),
    n_upper = c(30L, 1L, 45L, 20L),
    n_trials = 50L
  )
  model <- ezdm4_model()

  d1 <- suppressWarnings(check_data(model, d0, ezdm4_formula))
  expect_true(all(c("rt_used_upper", "rt_used_lower") %in% names(d1)))
  d2 <- suppressWarnings(check_data(model, d1, ezdm4_formula))
  expect_identical(d2, d1)

  switched_off <- suppressWarnings(check_data(
    model, transform(d0[1, ], rt_used_upper = 0L), ezdm4_formula
  ))
  expect_identical(switched_off$rt_used_upper, 0L)
  expect_identical(switched_off$mean_rt_upper, -1)
  expect_identical(switched_off$var_rt_upper, -1)
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
  declarations <- if (version == "3par") {
    "vector[N] mean_rt; vector[N] var_rt; vector[N] s;"
  } else {
    paste(
      "vector[N] mean_rt_upper; vector[N] mean_rt_lower;",
      "vector[N] var_rt_upper; vector[N] var_rt_lower;",
      "vector[N] zr; vector[N] s;",
      "array[N] int rt_upper; array[N] int rt_lower;"
    )
  }
  call <- if (version == "3par") {
    "ezdm_3par_lpdf(mean_rt[i] | 0.0, drift[i], bound[i], ndt[i], s[i],
                    var_rt[i], n_upper[i], n_trials[i])"
  } else {
    "ezdm_4par_lpdf(mean_rt_upper[i] | 0.0, drift[i], bound[i], ndt[i], zr[i],
                    s[i], mean_rt_lower[i], var_rt_upper[i], var_rt_lower[i],
                    n_upper[i], n_trials[i], rt_upper[i], rt_lower[i])"
  }

  program <- paste0(
    "functions {\n", .ezdm_stan_functions(version),
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
  # the indicators check_data() sets where a boundary's summaries exist
  grid$rt_upper <- as.integer(grid$n_upper >= 2L)
  grid$rt_lower <- as.integer(grid$n_trials - grid$n_upper >= 2L)

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
      "drift", "bound", "ndt", "zr", "s", "n_upper", "n_trials",
      "rt_upper", "rt_lower"
    )])
  ))

  expect_false(any(is.nan(r_lpdf)))
  expect_false(any(is.nan(stan_lpdf)))
  expect_equal(is.finite(stan_lpdf), is.finite(r_lpdf))

  finite <- is.finite(r_lpdf)
  expect_lt(max(abs(stan_lpdf[finite] - r_lpdf[finite]) / abs(r_lpdf[finite])), 1e-9)
})

# Unused boundaries in the likelihood (#430) -----------------------------------

# The ways a boundary's summaries go unused: indicator 0 at a boundary with
# enough responses (what check_data() sets where the summaries are missing),
# fewer than two responses, and fewer than two responses with indicator 1,
# which the count must still gate
ezdm_unused_variants <- data.frame(
  n_upper = c(15L, 15L, 15L, 1L, 29L, 0L, 30L, 1L, 29L),
  rt_upper = c(0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 1L),
  rt_lower = c(1L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L)
)

ezdm_unused_grid <- function(placeholder) {
  grid <- merge(
    expand.grid(
      drift = ezdm_parity_drift, bound = c(0.4, 1.5, 3), zr = c(0.2, 0.7),
      s = c(0.5, 1, 2)
    ),
    ezdm_unused_variants
  )
  grid$n_trials <- 30L
  grid$ndt <- 0.25
  moments <- .ezdm_moments_4par(grid$drift, grid$bound, grid$zr, grid$s)
  grid$mean_rt_upper <- grid$ndt + moments$mdt_upper * 1.02
  grid$mean_rt_lower <- grid$ndt + moments$mdt_lower * 0.98
  grid$var_rt_upper <- moments$vrt_upper * 0.9
  grid$var_rt_lower <- moments$vrt_lower * 1.1
  unused_upper <- grid$rt_upper == 0L | grid$n_upper < 2L
  unused_lower <- grid$rt_lower == 0L | grid$n_trials - grid$n_upper < 2L
  grid[unused_upper, c("mean_rt_upper", "var_rt_upper")] <- placeholder
  grid[unused_lower, c("mean_rt_lower", "var_rt_lower")] <- placeholder
  grid
}

ezdm_stan_columns_4par <- c(
  "mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower",
  "drift", "bound", "ndt", "zr", "s", "n_upper", "n_trials",
  "rt_upper", "rt_lower"
)

dezdm_4par_grid <- function(grid) {
  dezdm(
    mean_rt = as.matrix(grid[c("mean_rt_upper", "mean_rt_lower")]),
    var_rt = as.matrix(grid[c("var_rt_upper", "var_rt_lower")]),
    n_upper = grid$n_upper, n_trials = grid$n_trials, drift = grid$drift,
    bound = grid$bound, ndt = grid$ndt, zr = grid$zr, s = grid$s,
    version = "4par"
  )
}

# -1 and 0 make gamma_lpdf raise if they are read, which in generated
# quantities turns every value of the run into NaN; 0.3 (a plausible RT) and
# 1e300 would only move the density, which the comparison with dezdm() sees.
# All placeholders share one run, stacked by rows.
test_that("ezdm_4par_lpdf never reads the summaries of an unused boundary", {
  skip_on_cran()
  placeholders <- c(-1, 0, 0.3, 1e300)
  grids <- lapply(placeholders, ezdm_unused_grid)

  # dezdm() has no indicator, so below two responses its count alone keeps
  # the placeholder unread; NA marks the unused boundaries with enough
  # responses, as for a boundary without summaries
  r_count_gated <- vapply(grids, function(grid) {
    off_upper <- grid$rt_upper == 0L & grid$n_upper >= 2L
    off_lower <- grid$rt_lower == 0L & grid$n_trials - grid$n_upper >= 2L
    grid[off_upper, c("mean_rt_upper", "var_rt_upper")] <- NA
    grid[off_lower, c("mean_rt_lower", "var_rt_lower")] <- NA
    dezdm_4par_grid(grid)
  }, numeric(nrow(grids[[1]])))
  for (k in seq_along(placeholders)[-1]) {
    expect_identical(r_count_gated[, k], r_count_gated[, 1])
  }

  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  stacked <- do.call(rbind, grids)
  stan_lpdf <- matrix(
    ezdm_stan_lpdf("4par", c(
      list(N = nrow(stacked)), as.list(stacked[ezdm_stan_columns_4par])
    )),
    ncol = length(placeholders)
  )
  for (k in seq_along(placeholders)[-1]) {
    expect_identical(stan_lpdf[, k], stan_lpdf[, 1])
  }

  grid <- grids[[1]]
  grid[grid$rt_upper == 0L, c("mean_rt_upper", "var_rt_upper")] <- NA
  grid[grid$rt_lower == 0L, c("mean_rt_lower", "var_rt_lower")] <- NA
  r_lpdf <- dezdm_4par_grid(grid)

  expect_false(anyNA(r_lpdf))
  expect_false(anyNA(stan_lpdf[, 1]))
  expect_equal(is.finite(stan_lpdf[, 1]), is.finite(r_lpdf))
  finite <- is.finite(r_lpdf)
  expect_lt(max(abs(stan_lpdf[finite, 1] - r_lpdf[finite]) / abs(r_lpdf[finite])), 1e-9)
})

# The reference is written out here rather than taken from bmm: the
# probability that a Wiener process started zr * bound above the lower
# boundary is absorbed at the upper one, checked against rtdists::ddiffusion()
# integrated over time. |drift * bound| <= 3 keeps this direct form accurate
# far below the tolerance.
test_that("an ezdm 4par cell without a used boundary is a Wiener binomial", {
  skip_on_cran()
  grid <- expand.grid(
    drift = c(-2, -0.9, -0.3, 0.3, 0.9, 2), bound = c(0.5, 1.5),
    zr = c(0.2, 0.5, 0.8), s = c(1, 1.3), n_upper = c(3L, 12L, 27L)
  )
  grid$n_trials <- 30L
  grid$ndt <- 0.25
  start <- grid$zr * grid$bound
  p_upper <- (1 - exp(-2 * grid$drift * start / grid$s^2)) /
    (1 - exp(-2 * grid$drift * grid$bound / grid$s^2))
  reference <- stats::dbinom(grid$n_upper, grid$n_trials, p_upper, log = TRUE)

  grid[c("mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower")] <- NA_real_
  expect_lt(max(abs(dezdm_4par_grid(grid) / reference - 1)), 1e-9)

  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  grid[c("mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower")] <- -1
  grid$rt_upper <- 0L
  grid$rt_lower <- 0L
  stan_lpdf <- ezdm_stan_lpdf("4par", c(
    list(N = nrow(grid)), as.list(grid[ezdm_stan_columns_4par])
  ))
  expect_lt(max(abs(stan_lpdf / reference - 1)), 1e-9)
})

# The placeholder is -1 so that a summary read by mistake cannot pass
# unnoticed: no variance is negative. cmdstan raises inside generated
# quantities and writes NaN for every value of that iteration, so each
# boundary gets a run of its own, and no such row may join a batched run.
test_that("a placeholder read as data gives no finite ezdm 4par density", {
  skip_on_cran()
  moments <- .ezdm_moments_4par(1, 1.5, 0.5, 1)
  mean_rt <- 0.25 + c(moments$mdt_upper, moments$mdt_lower)
  var_rt <- c(moments$vrt_upper, moments$vrt_lower)
  lpdf <- function(mean_rt, var_rt) {
    dezdm(mean_rt, var_rt, n_upper = 20, n_trials = 30, drift = 1,
          bound = 1.5, ndt = 0.25, zr = 0.5, version = "4par")
  }
  expect_identical(lpdf(c(-1, mean_rt[2]), c(-1, var_rt[2])), -Inf)
  expect_identical(lpdf(c(mean_rt[1], -1), c(var_rt[1], -1)), -Inf)

  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  for (boundary in c("upper", "lower")) {
    data <- list(
      N = 2L, drift = c(1, 1), bound = c(1.5, 1.5), ndt = c(0.25, 0.25),
      zr = c(0.5, 0.5), s = c(1, 1), n_upper = c(20L, 20L),
      n_trials = c(30L, 30L), mean_rt_upper = rep(mean_rt[1], 2),
      mean_rt_lower = rep(mean_rt[2], 2), var_rt_upper = rep(var_rt[1], 2),
      var_rt_lower = rep(var_rt[2], 2), rt_upper = c(1L, 1L),
      rt_lower = c(1L, 1L)
    )
    data[[paste0("mean_rt_", boundary)]][1] <- -1
    data[[paste0("var_rt_", boundary)]][1] <- -1
    stan_lpdf <- suppressMessages(ezdm_stan_lpdf("4par", data))
    expect_false(is.finite(stan_lpdf[1]), info = boundary)
  }
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

test_that("ezdm 4par Stan code hands the RT indicators to the likelihood", {
  skip_on_cran()
  data <- data.frame(
    mean_rt_upper = c(0.5, 0.55), mean_rt_lower = c(0.6, NA),
    var_rt_upper = c(0.02, 0.03), var_rt_lower = c(0.03, NA),
    n_upper = c(60L, 99L), n_trials = c(100L, 100L)
  )
  code <- suppressWarnings(
    stancode(ezdm4_formula, data = data, model = ezdm4_model())
  )
  expect_match(code, "vint1[n], vint2[n], vint3[n], vint4[n]);", fixed = TRUE)
})

test_that("ezdm_4par_lpdf has the right gradient at and around zero drift", {
  # Returning the limit of pC at drift = 0 as a constant gave the right value
  # and a zero gradient, which only autodiff can see. CmdStan's standalone
  # log_prob method is used because expose_functions() returns no gradients.
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  program <- paste0(
    "functions {\n", .ezdm_stan_functions("4par"), "\n}\n",
    "data { int N; vector[2] mrt; vector[2] vrt; }\n",
    "parameters { vector[N] drift; vector[N] bound; vector[N] zr; vector[N] s; }\n",
    "model { for (i in 1:N) target += ezdm_4par_lpdf(mrt[1] | 0.0, drift[i], bound[i], 0.25,\n",
    "  zr[i], s[i], mrt[2], vrt[1], vrt[2], 12, 60, 1, 1); }\n"
  )
  # summaries of a cell with negative drift: most responses at the lower boundary
  moments <- .ezdm_moments_4par(-1.2, 1.5, 0.3, 1)
  mrt <- 0.25 + c(moments$mdt_upper, moments$mdt_lower)
  vrt <- c(moments$vrt_upper, moments$vrt_lower)
  # sw is the drift at which ezdm_logit_pc() leaves its series. The points inside
  # that band catch a wrong series coefficient, which no R-side test can see,
  # and the ones at 99 sw a switch that has moved outwards
  sw <- 1e-4 / (2 * 1.5)
  drift <- c(-1.2, -1e-12, 0, 1e-12, 1.2, -0.5 * sw, 0.99 * sw, -99 * sw, 99 * sw)
  pars <- list(drift = drift, bound = 1.5, zr = 0.3, s = 1)

  dir <- withr::local_tempdir()
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program, dir = dir), quiet = TRUE)
  cmdstanr::write_stan_json(list(N = length(drift), mrt = mrt, vrt = vrt), file.path(dir, "data.json"))
  cmdstanr::write_stan_json(lapply(pars, rep_len, length(drift)), file.path(dir, "pars.json"))
  # On Windows the executable loads tbb.dll at startup, and CmdStan's makefile
  # copies that dll next to it only when TBB was absent from PATH at link time.
  # cmdstanr has it on PATH whenever it builds or runs a model, so launching
  # the executable ourselves has to put it there too
  run_path <- if (.Platform$OS.type == "windows") {
    c(
      getFromNamespace("toolchain_PATH_env_var", "cmdstanr")(),
      getFromNamespace("tbb_path", "cmdstanr")()
    )
  }
  status <- withr::with_path(run_path, system2(model$exe_file(), c(
    "method=log_prob", paste0("constrained_params=", file.path(dir, "pars.json")), "jacobian=0",
    "data", paste0("file=", file.path(dir, "data.json")),
    "output", paste0("file=", file.path(dir, "out.csv")), "sig_figs=17"
  ), stdout = FALSE, stderr = file.path(dir, "stderr.txt")))
  expect_equal(status, 0, info = paste(
    readLines(file.path(dir, "stderr.txt"), warn = FALSE), collapse = "\n"
  ))
  stan_gradient <- matrix(
    as.numeric(utils::read.csv(file.path(dir, "out.csv"), comment.char = "#")[1, -1]),
    ncol = length(pars), dimnames = list(NULL, names(pars))
  )

  lpdf <- function(par, step, drift) {
    point <- list(drift = drift, bound = 1.5, zr = 0.3, s = 1)
    point[[par]] <- point[[par]] + step
    dezdm(mrt, vrt, n_upper = 12, n_trials = 60, drift = point$drift, bound = point$bound,
          ndt = 0.25, zr = point$zr, s = point$s, version = "4par")
  }
  # Richardson-extrapolated central difference: a plain one at a step small
  # enough for its truncation error is too noisy for the tolerance below
  central <- function(par, h, drift) (lpdf(par, h, drift) - lpdf(par, -h, drift)) / (2 * h)
  r_gradient <- vapply(colnames(stan_gradient), function(par) {
    vapply(drift, \(d) (4 * central(par, 5e-5, d) - central(par, 1e-4, d)) / 3, numeric(1))
  }, numeric(length(drift)))

  # max, not expect_equal(): its tolerance applies to the mean difference, so
  # one wrong point among nine would pass
  expect_lt(max(abs(stan_gradient / r_gradient - 1)), 1e-8)
})

test_that("the Stan series literals are the coefficients of log(sinh x / x)", {
  # the only check of the Stan arithmetic that needs no Stan: every literal of
  # the generated chunk against the rationals the R code uses
  txt <- readLines(system.file("stan_chunks", "ezdm_series.stan", package = "bmm"))
  literals <- function(name) {
    from <- grep(paste0("real ", name, "\\(real x\\)"), txt)
    to <- from + grep("return acc;", txt[-seq_len(from)])[1]
    lines <- txt[from:to]
    as.numeric(regmatches(lines, regexpr("-?[0-9.]+(e-?[0-9]+)?(?=\\);|;$)", lines, perl = TRUE)))
  }
  a <- .EZDM_LOG_SINHC_COEF
  j <- seq_along(a)
  for (n in 1:4) {
    for (kind in c("log_sinhc", "log_cosh")) {
      coef <- a * (if (kind == "log_cosh") 4^j - 1 else 1) *
        ifelse(j >= n, factorial(j) / factorial(pmax(j - n, 0)), 0)
      expected <- rev(coef[j >= n])
      got <- literals(paste0("ezdm_", kind, "_d", n))
      expect_length(got, length(expected))
      expect_lt(max(abs(got / expected - 1)), 1e-14)
    }
  }
})

# log_lik() with newdata (#430) --------------------------------------------------

test_that("log_lik() of an ezdm 4par fit never reads a placeholder", {
  skip_on_cran()
  path <- test_path("assets", "bmmfit_ezdm4_ppcheck.rds")
  skip_if_not(file.exists(path), "fixture not available (excluded by .Rbuildignore)")
  old <- readRDS(path)
  placeholders <- c(-1, 0, 1e300)

  # a fit saved before #430 has no indicators, and its stored log_lik closure
  # calls the dezdm() of the loaded bmm, whose count must keep a boundary
  # below two responses unread
  nd <- old$data
  nd$n_upper[1] <- nd$n_trials[1] - 1
  old_ll <- lapply(placeholders, function(p) {
    nd[1, c("mean_rt_lower", "var_rt_lower")] <- p
    brms::log_lik(old, newdata = nd)
  })
  expect_identical(old_ll[[2]], old_ll[[1]])
  expect_identical(old_ll[[3]], old_ll[[1]])

  # the fixture's draws, passed through the current pipeline by the mock
  # backend on its data with one unused boundary of each kind; every
  # parameter is intercept-only, so the draws do not depend on the data
  sparse <- old$data
  sparse$n_upper[2] <- sparse$n_trials[2] - 1
  sparse[2, c("mean_rt_lower", "var_rt_lower")] <- NA
  sparse[5, c("mean_rt_upper", "var_rt_upper")] <- NA
  fit <- suppressWarnings(suppressMessages(bmm(
    old$bmm$user_formula, sparse, old$bmm$model,
    backend = "mock", mock_fit = old$fit, rename = FALSE
  )))
  nd <- fit$data
  expect_identical(nrow(nd), nrow(sparse))
  expect_identical(nd$rt_used_lower[2], 0L)
  expect_identical(nd$rt_used_upper[5], 0L)

  new_ll <- lapply(placeholders, function(p) {
    nd[2, c("mean_rt_lower", "var_rt_lower")] <- p
    nd[5, c("mean_rt_upper", "var_rt_upper")] <- p
    brms::log_lik(fit, newdata = nd)
  })
  expect_identical(new_ll[[2]], new_ll[[1]])
  expect_identical(new_ll[[3]], new_ll[[1]])

  prep <- brms::prepare_predictions(fit)
  for (i in c(2L, 5L)) {
    reference <- dezdm(
      mean_rt = c(sparse$mean_rt_upper[i], sparse$mean_rt_lower[i]),
      var_rt = c(sparse$var_rt_upper[i], sparse$var_rt_lower[i]),
      n_upper = sparse$n_upper[i], n_trials = sparse$n_trials[i],
      drift = brms::get_dpar(prep, "drift", i = i),
      bound = brms::get_dpar(prep, "bound", i = i),
      ndt = brms::get_dpar(prep, "ndt", i = i),
      zr = brms::get_dpar(prep, "zr", i = i),
      s = brms::get_dpar(prep, "s", i = i),
      version = "4par"
    )
    expect_equal(new_ll[[1]][, i], reference, info = paste("row", i))
  }
})
