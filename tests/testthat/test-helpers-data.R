test_that("check_data() produces expected errors and warnings", {
  expect_error(
    check_data(.model_mixture2p(resp_error = "y")),
    "Data must be specified using the 'data' argument."
  )
  expect_error(
    check_data(.model_mixture2p(resp_error = "y"), data.frame(), bmf(kappa ~ 1)),
    "Argument 'data' does not contain observations."
  )
  expect_error(
    check_data(.model_mixture2p(resp_error = "y"), data.frame(x = 1), bmf(kappa ~ 1)),
    "The response variable 'y' is not present in the data."
  )
  expect_error(
    check_data(.model_mixture2p(resp_error = "y"), y ~ 1),
    "Argument 'data' must be coercible to a data.frame."
  )

  mls <- lapply(c("mixture2p", "mixture3p", "imm"), get_model)
  for (ml in mls) {
    model <- ml(resp_error = "y", nt_features = "x", set_size = 2, nt_distances = "z")
    expect_warning(
      check_data(model, data.frame(y = 12, x = 1, z = 2), bmf(kappa ~ 1)),
      "It appears your response variable is in degrees.\n"
    )
    model <- ml(resp_error = "y", nt_features = "x", set_size = 2, nt_distances = "z")
    expect_silent(check_data(model, data.frame(y = 1, x = 1, z = 2), bmf(y ~ 1)))
  }

  mls <- lapply(c("mixture3p", "imm"), get_model)
  for (ml in mls) {
    model <- ml(resp_error = "y", nt_features = "x", set_size = 5, nt_distances = "z")
    expect_error(
      check_data(model, data.frame(y = 1, x = 1, z = 2), bmf(kappa ~ 1)),
      "'nt_features' should equal max\\(set_size\\)-1"
    )

    model <- ml(resp_error = "y", nt_features = "x", set_size = 2, nt_distances = "z")
    expect_warning(
      check_data(model, data.frame(y = 1, x = 2 * pi + 1, z = 2), bmf(kappa ~ 1)),
      "at least one of your non_target variables are in degrees"
    )
  }

  for (version in c("bsc", "full")) {
    ml <- imm(
      resp_error = "y", nt_features = paste0("x", 1:3), set_size = 4,
      nt_distances = "z", version = version
    )
    expect_error(
      check_data(ml, data.frame(y = 1, x1 = 1, x2 = 2, x3 = 3, z = 2), bmf(kappa ~ 1)),
      "'nt_distances' should equal max\\(set_size\\)-1"
    )
  }
})

test_that("check_var_set_size accepts valid input", {
  # Simple numeric vector is valid
  dat <- data.frame(y = rep(c(1, 2, 3), each = 3))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(names(check_var_set_size("y", dat)), c("max_set_size", "ss_numeric"))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 3)
  all(is.numeric(check_var_set_size("y", dat)$ss_numeric), na.rm = T)

  # Factor with numeric levels is valid
  dat <- data.frame(y = factor(rep(c(1, 2, 3), each = 3)))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 3)
  all(is.numeric(check_var_set_size("y", dat)$ss_numeric), na.rm = T)

  # Character vector representing numbers is valid
  dat <- data.frame(y = rep(c("1", "2", "3"), each = 3))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 3)
  all(is.numeric(check_var_set_size("y", dat)$ss_numeric), na.rm = T)

  # Numeric vector with NA values is valid (assuming NA is treated correctly)
  dat <- data.frame(y = rep(c(1, 5, NA), each = 3))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 5)
  all(is.numeric(check_var_set_size("y", dat)$ss_numeric), na.rm = T)

  # Factor with NA and numeric levels is valid
  dat <- data.frame(y = factor(rep(c(1, 5, NA), each = 3)))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 5)
  all(is.numeric(check_var_set_size("y", dat)$ss_numeric), na.rm = T)
})

test_that("check_var_set_size rejects invalid input", {
  # Factor with non-numeric levels is invalid
  dat <- data.frame(y = factor(rep(c("A", "B", "C"), each = 3)))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Character vector with non-numeric values is invalid
  dat <- data.frame(y = rep(c("A", "B", "C"), each = 3))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Character vector with NA and non-numeric values is invalid
  dat <- data.frame(y = rep(c("A", NA, "C"), each = 3))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Factor with NA and non-numeric levels is invalid
  dat <- data.frame(y = factor(rep(c("A", NA, "C"), each = 3)))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Character vector with numeric and non-numeric values is invalid
  dat <- data.frame(y = rep(c("A", 5, "C"), each = 3))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Factor with numeric and non-numeric levels is invalid
  dat <- data.frame(y = factor(rep(c("A", 5, "C"), each = 3)))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Numeric vector with invalid set sizes (less than 1) is invalid
  dat <- data.frame(y = rep(c(0, 1, 5), each = 3))
  expect_error(check_var_set_size("y", dat), "must be positive whole numbers")

  # Factor with levels less than 1 are invalid
  dat <- data.frame(y = factor(rep(c(0, 4, 5), each = 3)))
  expect_error(check_var_set_size("y", dat), "must be positive whole numbers")

  # Character vector representing set sizes with text is invalid
  dat <- data.frame(y = rep(paste0("set_size ", c(2, 3, 8)), each = 3))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Factor representing set sizes with text is invalid
  dat <- data.frame(y = factor(rep(paste0("set_size ", c(2, 3, 8)), each = 3)))
  expect_error(check_var_set_size("y", dat), "must be coercible to a numeric vector")

  # Numeric vector with decimals is invalid
  dat <- data.frame(y = c(1:8, 1.3))
  expect_error(check_var_set_size("y", dat), "must be positive whole numbers")

  # Setsize must be of length 1
  dat <- data.frame(y = c(1, 2, 3), z = c(1, 2, 3))
  expect_error(check_var_set_size(c("z", "y"), dat), "You provided a vector")
  expect_error(check_var_set_size(list("z", "y"), dat), "You provided a vector")
  expect_error(
    check_var_set_size(set_size = TRUE, dat),
    "must be either a variable in your data or a single numeric value"
  )
})

test_that("check_data() returns a data.frame()", {
  mls <- lapply(supported_models(print_call = FALSE), get_model)
  # test data includes variables for all model types:
  # - y, x, z, w, l, s for circular/mixture models
  # - mean_rt, var_rt, n_upper, n_trials for ezdm 3par
  # - mean_rt_upper/lower, var_rt_upper/lower for ezdm 4par
  test_data <- data.frame(
    y = 1, x = 1, z = 2, w = 1, s = 2, l = 1,
    mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
    mean_rt_upper = 0.45, mean_rt_lower = 0.55,
    var_rt_upper = 0.018, var_rt_lower = 0.025
  )
  for (ml in mls) {
    model <- ml(
      resp_error = "y", nt_features = "x", set_size = 2,
      nt_distances = "z", resp_cats = c("w", "l"), num_options = c(1, 1),
      mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
      n_trials = "n_trials"
    )
    expect_s3_class(
      check_data(model, test_data, bmf(kappa ~ 1)),
      "data.frame"
    )
  }
})

test_that("wrap(x) returns the same for values between -pi and pi", {
  x <- runif(100, -pi, pi)
  expect_equal(wrap(x), x)
  expect_equal(wrap(rad2deg(x), radians = F), rad2deg(wrap(x)))
})

test_that("wrap(x) returns the correct value for values between (pi, 2*pi)", {
  x <- pi + 1
  expect_equal(wrap(x), -(pi - 1))
  expect_equal(wrap(rad2deg(x), radians = F), rad2deg(wrap(x)))
})

test_that("wrap(x) returns the correct value for values between (-2*pi, -pi)", {
  x <- -pi - 1
  expect_equal(wrap(x), pi - 1)
  expect_equal(wrap(rad2deg(x), radians = F), rad2deg(wrap(x)))
})

test_that("wrap(x) returns the correct value for values over 2*pi", {
  x <- 2 * pi + 1
  expect_equal(wrap(x), 1)
  expect_equal(wrap(rad2deg(x), radians = F), rad2deg(wrap(x)))
})

test_that("wrap(x) returns the correct value for values between (3*pi,4*pi)", {
  x <- 3 * pi + 1
  expect_equal(wrap(x), -(pi - 1))
  expect_equal(wrap(rad2deg(x), radians = F), rad2deg(wrap(x)))
})

test_that("deg2rad returns the correct values for 0, 180, 360", {
  x <- c(0, 90, 180)
  expect_equal(round(deg2rad(x), 2), c(0.00, 1.57, 3.14))
  expect_equal(wrap(rad2deg(x), radians = F), rad2deg(wrap(x)))
})

test_that("rad2deg returns the correct values for 0, pi/2, 2*pi", {
  x <- c(0, pi / 2, 2 * pi)
  expect_equal(round(rad2deg(x), 2), c(0, 90, 360))
  expect_equal(wrap(rad2deg(x), radians = F), rad2deg(wrap(x)))
})

test_that("standata() works with brmsformula", {
  ff <- brms::bf(count ~ zAge + zBase * Trt + (1 | patient))
  sd <- standata(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(sd)[1], "standata")
})

test_that("standata() works with formula", {
  ff <- count ~ zAge + zBase * Trt + (1 | patient)
  sd <- standata(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(sd)[1], "standata")
})

test_that("standata() works with bmf", {
  ff <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  dat <- oberauer_lin_2017
  sd <- standata(ff, dat, mixture3p(
    resp_error = "dev_rad",
    nt_features = "col_nt",
    set_size = "set_size", regex = T
  ))
  expect_equal(class(sd)[1], "standata")
})

test_that("standata() returns a standata class", {
  ff <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  dat <- data.frame(y = rmixture3p(n = 3), nt1_loc = 2, nt2_loc = -1.5)

  standata <- standata(ff, dat, mixture3p(
    resp_error = "y",
    nt_features = paste0("nt", 1, "_loc"),
    set_size = 2
  ))
  expect_equal(class(standata)[1], "standata")
})

# first draft of tests was written by ChatGPT4
test_that("has_nonconsecutive_duplicates works", {
  expect_false(has_nonconsecutive_duplicates(c("a", "a", "b", "b", "c", "c")))
  expect_true(has_nonconsecutive_duplicates(c("a", "b", "a", "c", "c", "b")))
  expect_false(has_nonconsecutive_duplicates(rep("a", 5)))
  expect_false(has_nonconsecutive_duplicates(letters[1:5]))
  expect_true(has_nonconsecutive_duplicates(c("a", "a", "b", "a", "b", "b")))
  expect_true(has_nonconsecutive_duplicates(c(1, 2, 3, 1, 4, 2)))
  expect_false(has_nonconsecutive_duplicates(numeric(0)))
  expect_false(has_nonconsecutive_duplicates(c("a")))
  expect_true(has_nonconsecutive_duplicates(c("a", "b", "b", "a")))
  expect_false(has_nonconsecutive_duplicates(c(NA, NA, NA)))
  expect_false(has_nonconsecutive_duplicates(c(NA, 1, NA)))
  expect_true(has_nonconsecutive_duplicates(c(1, "a", 2, "b", 1, NA, "a")))
  expect_true(has_nonconsecutive_duplicates(c("1", 2, "2", 1)))
})

test_that("is_data_ordered works", {
  # Test with a data frame that is ordered
  data1 <- expand.grid(y = 1:3, B = 1:3, C = 1:3)
  formula1 <- bmf(y ~ B + C)
  expect_true(is_data_ordered(data1, formula1))

  # Test with a data frame that is not ordered
  data2 <- rbind(data1, data1[1, ])
  expect_false(is_data_ordered(data2, formula1))

  # Test when irrelevant variables are not ordered but predictors are
  data3 <- data1
  data3$A <- c(3, 2, 2, 1, 2, 1, 3, 1, 3, 3, 1, 2, 2, 1, 1, 1, 3, 3, 1, 3, 2, 3, 1, 2, 3, 2, 2)
  formula2 <- bmf(y ~ A + B + C)
  expect_true(is_data_ordered(data3, formula1))
  expect_false(is_data_ordered(data3, formula2))

  # test with a complex formula with shared covariance structure across parameters
  data <- oberauer_lin_2017
  formula <- bmf(
    c ~ 0 + set_size + (0 + set_size | p1 | ID),
    kappa ~ 0 + set_size + (0 + set_size | p1 | ID)
  )
  expect_false(is_data_ordered(data, formula))

  data <- dplyr::arrange(data, set_size, ID)
  expect_true(is_data_ordered(data, formula))
})

test_that("is_data_ordered works when there is only one predictor", {
  # Test with a data frame that is ordered
  data <- data.frame(
    y = rep(1:3, each = 2),
    B = rep(1:3, each = 2),
    C = factor(rep(1:3, each = 2)),
    D = rep(1:3, times = 2),
    E = factor(rep(1:3, times = 2))
  )
  expect_true(is_data_ordered(data, y ~ B))

  # Test with a data frame that is not ordered
  expect_false(is_data_ordered(data, y ~ D))

  # Test with a data frame that is ordered and predictor is a factor
  expect_true(is_data_ordered(data, y ~ C))

  # Test with a data frame that is not ordered and predictor is a factor
  expect_false(is_data_ordered(data, y ~ E))
})

test_that("is_data_ordered works when there are no predictors", {
  # Test with a data frame that is ordered
  data <- data.frame(y = 1:3)
  expect_true(is_data_ordered(data, y ~ 1))
})

test_that("is_data_ordered works when there are non-linear predictors", {
  data <- data.frame(
    y = rep(1:3, each = 2),
    B = rep(1:3, each = 2),
    C = rep(1:3, times = 2)
  )
  # Test with a data frame that is ordered
  formula1 <- bmf(y ~ nlD, nlD ~ B)
  expect_true(is_data_ordered(data, formula1))

  # Test with a data frame that is not ordered
  formula2 <- bmf(y ~ nlD, nlD ~ C)
  expect_false(is_data_ordered(data, formula2))
})

############################################################################# !
# ezdm_summary_stats TESTS                                                ####
############################################################################# !

test_that("ezdm_summary_stats() returns correct structure for 3par version", {
  set.seed(123)
  test_data <- data.frame(
    subject = rep(1:3, each = 50),
    rt = rgamma(150, shape = 5, rate = 10) + 0.3,
    correct = rbinom(150, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                .by = "subject", method = "simple")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true(all(c("subject", "mean_rt", "var_rt", "n_upper",
                    "n_trials", "contaminant_prop") %in% names(result)))

  # Check types
  expect_type(result$mean_rt, "double")
  expect_type(result$var_rt, "double")
  expect_type(result$n_upper, "integer")
  expect_type(result$n_trials, "integer")
})

test_that("ezdm_summary_stats() returns correct structure for 4par version", {
  set.seed(123)
  test_data <- data.frame(
    subject = rep(1:2, each = 100),
    rt = rgamma(200, shape = 5, rate = 10) + 0.3,
    correct = rbinom(200, 1, 0.7)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                .by = "subject", version = "4par", method = "simple")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("subject", "mean_rt_upper", "mean_rt_lower",
                    "var_rt_upper", "var_rt_lower", "n_upper", "n_trials",
                    "contaminant_prop_upper", "contaminant_prop_lower")
                  %in% names(result)))
})

test_that("ezdm_summary_stats() validates required arguments", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(rt = "rt", response = "correct"),
    "Argument 'data' is required"
  )
  expect_error(
    ezdm_summary_stats(test_data, response = "correct"),
    "Argument 'rt' is required"
  )
  expect_error(
    ezdm_summary_stats(test_data, rt = "rt"),
    "Argument 'response' is required"
  )
})

test_that("ezdm_summary_stats() validates column names", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "reaction_time", response = "correct"),
    "RT variable 'reaction_time' not found in data"
  )
  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "accuracy"),
    "Response variable 'accuracy' not found in data"
  )
})

test_that("ezdm_summary_stats() validates parameter options", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       version = "5par"),
    "version must be '3par' or '4par'"
  )
  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       distribution = "normal"),
    "distribution must be 'exgaussian', 'lognormal', or 'invgaussian'"
  )
  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       method = "invalid"),
    "method must be 'mixture', 'simple', or 'robust'"
  )
  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       method = "robust", robust_scale = "invalid"),
    "robust_scale must be 'iqr' or 'mad'"
  )
})

test_that("ezdm_summary_stats() warns for potential data issues", {
  # RT values > 10 (likely in milliseconds)
  test_data_ms <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) * 1000,  # milliseconds
    correct = rbinom(100, 1, 0.8)
  )

  expect_warning(
    ezdm_summary_stats(test_data_ms, rt = "rt", response = "correct",
                       method = "simple"),
    "Some RT values > 10. Ensure RTs are in seconds"
  )

  # Non-positive RT values
  test_data_neg <- data.frame(
    rt = c(-0.1, 0, rgamma(98, shape = 5, rate = 10) + 0.3),
    correct = rbinom(100, 1, 0.8)
  )

  expect_warning(
    ezdm_summary_stats(test_data_neg, rt = "rt", response = "correct",
                       method = "simple"),
    "non-positive RT values will be excluded"
  )
})

test_that("ezdm_summary_stats() handles too few trials", {
  # Create data with one subject having very few trials
  set.seed(123)
  test_data <- data.frame(
    subject = c(rep(1, 50), rep(2, 5)),  # Subject 2 has only 5 trials
    rt = rgamma(55, shape = 5, rate = 10) + 0.3,
    correct = rbinom(55, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                .by = "subject", method = "simple", min_trials = 10)

  # Subject 2 should have NA for mean_rt and var_rt
  expect_true(is.na(result$mean_rt[result$subject == 2]))
  expect_true(is.na(result$var_rt[result$subject == 2]))

  # But n_trials should still be reported
  expect_equal(result$n_trials[result$subject == 2], 5)
})

test_that("ezdm_summary_stats() simple method matches mean() and var()", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                method = "simple")

  expect_equal(result$mean_rt, mean(test_data$rt), tolerance = 1e-10)
  expect_equal(result$var_rt, var(test_data$rt), tolerance = 1e-10)
  expect_equal(result$n_trials, nrow(test_data))
  expect_equal(result$n_upper, sum(test_data$correct == 1))
})

test_that("ezdm_summary_stats() robust method uses median and IQR/MAD", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  # Test IQR-based robust method
  result_iqr <- ezdm_summary_stats(
    test_data, rt = "rt", response = "correct",
    method = "robust", robust_scale = "iqr"
  )

  expect_equal(result_iqr$mean_rt, median(test_data$rt), tolerance = 1e-10)
  expected_var_iqr <- (IQR(test_data$rt) / 1.349)^2
  expect_equal(result_iqr$var_rt, expected_var_iqr, tolerance = 1e-10)
  expect_equal(result_iqr$n_trials, nrow(test_data))
  expect_equal(result_iqr$n_upper, sum(test_data$correct == 1))
  expect_true(is.na(result_iqr$contaminant_prop))

  # Test MAD-based robust method
  result_mad <- ezdm_summary_stats(
    test_data, rt = "rt", response = "correct",
    method = "robust", robust_scale = "mad"
  )

  expect_equal(result_mad$mean_rt, median(test_data$rt), tolerance = 1e-10)
  expected_var_mad <- mad(test_data$rt)^2
  expect_equal(result_mad$var_rt, expected_var_mad, tolerance = 1e-10)
})

test_that("ezdm_summary_stats() robust method is resistant to outliers", {
  set.seed(123)
  # Create data with outliers
  clean_rt <- rnorm(90, mean = 0.5, sd = 0.1)
  outliers <- c(0.05, 0.08, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0)
  test_data <- data.frame(
    rt = c(clean_rt, outliers),
    correct = c(rep(1, 80), rep(0, 20))
  )

  result_simple <- ezdm_summary_stats(
    test_data, rt = "rt", response = "correct", method = "simple"
  )
  result_robust <- ezdm_summary_stats(
    test_data, rt = "rt", response = "correct", method = "robust"
  )

  # Robust method should give mean_rt closer to true median of clean data
  true_median_clean <- median(clean_rt)
  expect_true(
    abs(result_robust$mean_rt - true_median_clean) <
      abs(result_simple$mean_rt - true_median_clean)
  )
})

test_that("ezdm_summary_stats() handles multiple grouping variables", {
  set.seed(123)
  test_data <- data.frame(
    subject = rep(1:2, each = 100),
    condition = rep(rep(c("A", "B"), each = 50), 2),
    rt = rgamma(200, shape = 5, rate = 10) + 0.3,
    correct = rbinom(200, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                .by = c("subject", "condition"), method = "simple")

  expect_equal(nrow(result), 4)  # 2 subjects x 2 conditions
  expect_true(all(c("subject", "condition") %in% names(result)))
})

test_that("ezdm_summary_stats() handles no grouping variables", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                method = "simple")

  expect_equal(nrow(result), 1)
  expect_false("subject" %in% names(result))
})

test_that("ezdm_summary_stats() mixture method works with different distributions", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(200, shape = 5, rate = 10) + 0.3,
    correct = rbinom(200, 1, 0.8)
  )

  # Test each distribution
  for (dist in c("exgaussian", "lognormal", "invgaussian")) {
    result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                  distribution = dist, method = "mixture")
    expect_s3_class(result, "data.frame")
    expect_true(is.numeric(result$mean_rt))
    expect_true(is.numeric(result$var_rt))
    expect_true(is.numeric(result$contaminant_prop))
  }
})

test_that("ezdm_summary_stats() 4par handles all correct or all errors", {
  set.seed(123)
  # All correct responses
  test_data_all_correct <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rep(1, 100)
  )

  result <- ezdm_summary_stats(test_data_all_correct, rt = "rt",
                                response = "correct", version = "4par",
                                method = "simple", min_trials = 10)

  expect_false(is.na(result$mean_rt_upper))
  expect_true(is.na(result$mean_rt_lower))  # No error trials

  # All error responses
  test_data_all_errors <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rep(0, 100)
  )

  result2 <- ezdm_summary_stats(test_data_all_errors, rt = "rt",
                                 response = "correct", version = "4par",
                                 method = "simple", min_trials = 10)

  expect_true(is.na(result2$mean_rt_upper))  # No correct trials
  expect_false(is.na(result2$mean_rt_lower))
})

test_that("ezdm_summary_stats() attributes are correctly set", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                version = "4par", distribution = "lognormal",
                                method = "simple")

  expect_equal(attr(result, "version"), "4par")
  expect_equal(attr(result, "distribution"), "lognormal")
  expect_equal(attr(result, "method"), "simple")
})

test_that("ezdm_summary_stats() validates grouping variables", {
  set.seed(123)
  test_data <- data.frame(
    subject = rep(1:3, each = 50),
    rt = rgamma(150, shape = 5, rate = 10) + 0.3,
    correct = rbinom(150, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       .by = "nonexistent_var"),
    "Grouping variable 'nonexistent_var' not found in data"
  )
})

test_that("ezdm_summary_stats() validates contaminant_bound", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       contaminant_bound = c(3.0, 0.1)),
    "contaminant_bound\\[1\\] must be less than"
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       contaminant_bound = c(0.1)),
    "contaminant_bound must be a vector of length 2"
  )

  # Invalid element type
  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       contaminant_bound = c("invalid", 3.0)),
    "contaminant_bound elements must be numeric or"
  )
})

test_that("ezdm_summary_stats() accepts 'min' and 'max' for contaminant_bound", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  # c("min", "max") should work

  expect_no_error(
    result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                  contaminant_bound = c("min", "max"))
  )
  expect_true(is.data.frame(result))
  expect_true("mean_rt" %in% names(result))

  # c(0.1, "max") should work (mixed)
  expect_no_error(
    result2 <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                   contaminant_bound = c(0.1, "max"))
  )
  expect_true(is.data.frame(result2))

  # c("min", 3.0) should work (mixed)
  expect_no_error(
    result3 <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                   contaminant_bound = c("min", 3.0))
  )
  expect_true(is.data.frame(result3))

  # Case insensitive
  expect_no_error(
    result4 <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                   contaminant_bound = c("MIN", "MAX"))
  )
  expect_true(is.data.frame(result4))
})

test_that("ezdm_summary_stats() data-driven bounds use actual data range", {
  set.seed(456)
  # Create data with known min/max
  test_data <- data.frame(
    rt = c(0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5),
    correct = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0)
  )
  # Min RT = 0.2, Max RT = 2.5

  # Using c("min", "max") should capture all data in bounds
  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                contaminant_bound = c("min", "max"),
                                method = "simple")

  # Should have all 10 trials
  expect_equal(result$n_trials, 10)
})

test_that("ezdm_summary_stats() data-driven bounds work per group", {
  set.seed(789)
  # Create data with different RT ranges for different subjects
  test_data <- data.frame(
    subject = rep(1:2, each = 50),
    rt = c(
      runif(50, 0.3, 1.0),  # Subject 1: range 0.3-1.0
      runif(50, 0.5, 2.0)   # Subject 2: range 0.5-2.0
    ),
    correct = rbinom(100, 1, 0.8)
  )

  # Using "min"/"max" should adapt to each group's range
  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                .by = "subject", contaminant_bound = c("min", "max"),
                                method = "simple")

  expect_equal(nrow(result), 2)
  expect_equal(result$n_trials, c(50, 50))
})

test_that("ezdm_summary_stats() validates init_contaminant", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       init_contaminant = 0),
    "init_contaminant must be between 0 and 1"
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       init_contaminant = 1),
    "init_contaminant must be between 0 and 1"
  )
})

test_that("ezdm_summary_stats() validates max_contaminant", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       max_contaminant = 0),
    "max_contaminant must be between 0"
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       max_contaminant = 1.5),
    "max_contaminant must be between 0"
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       max_contaminant = -0.1),
    "max_contaminant must be between 0"
  )

  # max_contaminant = 1 should be valid (no upper bound enforcement)
  expect_no_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       max_contaminant = 1, method = "simple")
  )
})

test_that("ezdm_summary_stats() validates init < max contaminant", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       init_contaminant = 0.3, max_contaminant = 0.2),
    "init_contaminant must be less than max_contaminant"
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       init_contaminant = 0.5, max_contaminant = 0.5),
    "init_contaminant must be less than max_contaminant"
  )
})

test_that("ezdm_summary_stats() clips contaminant proportion to max", {
  set.seed(42)
  # Generate data with many contaminants (uniform RTs)
  n_true <- 50
  n_contam <- 50
  true_rt <- rgamma(n_true, shape = 5, rate = 10) + 0.3
  contam_rt <- runif(n_contam, 0.1, 3.0)
  test_data <- data.frame(
    rt = c(true_rt, contam_rt),
    correct = rbinom(100, 1, 0.8)
  )

  # With default max_contaminant = 0.5, high contamination should be clipped
  result <- suppressWarnings(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       max_contaminant = 0.3)
  )
  expect_true(result$contaminant_prop <= 0.3)

  # With max_contaminant = 0.1, should be clipped even more
  result_strict <- suppressWarnings(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       max_contaminant = 0.1)
  )
  expect_true(result_strict$contaminant_prop <= 0.1)
})

test_that("ezdm_summary_stats() warns when contaminant proportion is clipped", {
  set.seed(42)
  # Generate data with many contaminants
  n_true <- 50
  n_contam <- 50
  true_rt <- rgamma(n_true, shape = 5, rate = 10) + 0.3
  contam_rt <- runif(n_contam, 0.1, 3.0)
  test_data <- data.frame(
    rt = c(true_rt, contam_rt),
    correct = rbinom(100, 1, 0.8)
  )

  # Should warn when clipping occurs
  expect_warning(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       max_contaminant = 0.1),
    "clipped to max_contaminant"
  )
})

test_that("ezdm_summary_stats() max_contaminant works with 4par version", {
  set.seed(42)
  # Generate data with contaminants for both boundary types
  n_per_group <- 60
  true_rt_upper <- rgamma(n_per_group, shape = 5, rate = 10) + 0.3
  true_rt_lower <- rgamma(n_per_group, shape = 5, rate = 10) + 0.3
  contam_upper <- runif(20, 0.1, 3.0)
  contam_lower <- runif(20, 0.1, 3.0)

  test_data <- data.frame(
    rt = c(true_rt_upper, contam_upper, true_rt_lower, contam_lower),
    correct = c(rep(1, n_per_group + 20), rep(0, n_per_group + 20))
  )

  result <- suppressWarnings(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       version = "4par", max_contaminant = 0.2)
  )

  # Both contaminant proportions should be clipped
  expect_true(is.na(result$contaminant_prop_upper) ||
                result$contaminant_prop_upper <= 0.2)
  expect_true(is.na(result$contaminant_prop_lower) ||
                result$contaminant_prop_lower <= 0.2)
})

test_that("ezdm_summary_stats() handles character 'upper'/'lower' responses", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    response = sample(c("upper", "lower"), 100, replace = TRUE, prob = c(0.8, 0.2))
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "response",
                                method = "simple")

  expect_equal(result$n_upper, sum(test_data$response == "upper"))
  expect_equal(result$n_trials, nrow(test_data))
  expect_true(is.numeric(result$mean_rt))
})

test_that("ezdm_summary_stats() handles factor 'upper'/'lower' responses", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    response = factor(sample(c("upper", "lower"), 100, replace = TRUE, prob = c(0.8, 0.2)))
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "response",
                                method = "simple")

  expect_equal(result$n_upper, sum(test_data$response == "upper"))
  expect_equal(result$n_trials, nrow(test_data))
})

test_that("ezdm_summary_stats() handles character 'correct'/'error' responses", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    response = sample(c("correct", "error"), 100, replace = TRUE, prob = c(0.8, 0.2))
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "response",
                                method = "simple")

  expect_equal(result$n_upper, sum(test_data$response == "correct"))
  expect_equal(result$n_trials, nrow(test_data))
})

test_that("ezdm_summary_stats() handles logical TRUE/FALSE responses", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    response = sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.8, 0.2))
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "response",
                                method = "simple")

  expect_equal(result$n_upper, sum(test_data$response))
  expect_equal(result$n_trials, nrow(test_data))
})

test_that("ezdm_summary_stats() handles case-insensitive responses", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    response = sample(c("UPPER", "Lower", "UpPeR", "LOWER"), 100, replace = TRUE)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "response",
                                method = "simple")

  expect_equal(result$n_upper, sum(tolower(test_data$response) == "upper"))
  expect_equal(result$n_trials, nrow(test_data))
})

test_that("ezdm_summary_stats() errors on unrecognized response values", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    response = sample(c("fast", "slow"), 100, replace = TRUE)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "response"),
    "Unrecognized response values"
  )
})

test_that("ezdm_summary_stats() 4par version works with character responses", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(200, shape = 5, rate = 10) + 0.3,
    response = sample(c("upper", "lower"), 200, replace = TRUE, prob = c(0.7, 0.3))
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "response",
                                version = "4par", method = "simple")

  expect_true(all(c("mean_rt_upper", "mean_rt_lower", "var_rt_upper",
                    "var_rt_lower") %in% names(result)))
  expect_equal(result$n_upper, sum(test_data$response == "upper"))
  expect_false(is.na(result$mean_rt_upper))
  expect_false(is.na(result$mean_rt_lower))
})

############################################################################# !
# ezdm_summary_stats ACCURACY ADJUSTMENT TESTS                            ####
############################################################################# !

test_that("ezdm_summary_stats() adjust_accuracy adds columns when TRUE", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(200, shape = 5, rate = 10) + 0.3,
    correct = rbinom(200, 1, 0.8)
  )

  # Without adjustment
  result_no_adj <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                       adjust_accuracy = FALSE)
  expect_false("n_upper_adj" %in% names(result_no_adj))
  expect_false("n_trials_adj" %in% names(result_no_adj))

  # With adjustment
  result_adj <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                    adjust_accuracy = TRUE)
  expect_true("n_upper_adj" %in% names(result_adj))
  expect_true("n_trials_adj" %in% names(result_adj))
})

test_that("ezdm_summary_stats() adjust_accuracy produces integer values", {
  set.seed(42)
  # Create data with known properties
  n <- 200
  test_data <- data.frame(
    rt = rgamma(n, shape = 5, rate = 10) + 0.3,
    correct = rbinom(n, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                adjust_accuracy = TRUE, seed = 123)

  # Adjusted values should be integers
  expect_type(result$n_upper_adj, "integer")
  expect_type(result$n_trials_adj, "integer")

  # Adjusted values should be bounded
  expect_true(result$n_upper_adj >= 0)
  expect_true(result$n_upper_adj <= result$n_trials_adj)
  expect_true(result$n_trials_adj <= result$n_trials)
})

test_that("ezdm_summary_stats() seed ensures reproducibility", {
  set.seed(42)
  test_data <- data.frame(
    rt = rgamma(200, shape = 5, rate = 10) + 0.3,
    correct = rbinom(200, 1, 0.8)
  )

  # Same seed should give same results
  result1 <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                 adjust_accuracy = TRUE, seed = 999)
  result2 <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                 adjust_accuracy = TRUE, seed = 999)

  expect_equal(result1$n_upper_adj, result2$n_upper_adj)
  expect_equal(result1$n_trials_adj, result2$n_trials_adj)
})

test_that("ezdm_summary_stats() adjust_accuracy warns with simple method", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_warning(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       method = "simple", adjust_accuracy = TRUE),
    "adjust_accuracy has no effect"
  )
})

test_that("ezdm_summary_stats() adjust_accuracy validates guess_rate", {
  test_data <- data.frame(
    rt = rgamma(100, shape = 5, rate = 10) + 0.3,
    correct = rbinom(100, 1, 0.8)
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       guess_rate = -0.1),
    "guess_rate must be between 0 and 1"
  )

  expect_error(
    ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                       guess_rate = 1.5),
    "guess_rate must be between 0 and 1"
  )
})

test_that("ezdm_summary_stats() 4par with adjust_accuracy works", {
  set.seed(123)
  test_data <- data.frame(
    rt = rgamma(300, shape = 5, rate = 10) + 0.3,
    correct = rbinom(300, 1, 0.75)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                version = "4par", adjust_accuracy = TRUE,
                                seed = 456)

  expect_true("n_upper_adj" %in% names(result))
  expect_true("n_trials_adj" %in% names(result))
  expect_type(result$n_upper_adj, "integer")
  expect_type(result$n_trials_adj, "integer")
  expect_true(result$n_upper_adj >= 0)
  expect_true(result$n_upper_adj <= result$n_trials_adj)
})

test_that("ezdm_summary_stats() adjust_accuracy works with grouping", {
  set.seed(123)
  test_data <- data.frame(
    subject = rep(1:3, each = 100),
    rt = rgamma(300, shape = 5, rate = 10) + 0.3,
    correct = rbinom(300, 1, 0.8)
  )

  result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
                                .by = "subject", adjust_accuracy = TRUE,
                                seed = 789)

  expect_equal(nrow(result), 3)
  expect_true(all(c("n_upper_adj", "n_trials_adj") %in% names(result)))

  # Check types are integers
  expect_type(result$n_upper_adj, "integer")
  expect_type(result$n_trials_adj, "integer")

  # Check each group has valid adjusted values
  for (i in 1:3) {
    expect_true(result$n_upper_adj[i] >= 0)
    expect_true(result$n_upper_adj[i] <= result$n_trials_adj[i])
  }
})
