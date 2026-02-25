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

test_that("ezdm_summary_stats() returns 1-row data.frame for 3par version", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  result <- ezdm_summary_stats(rt, response, method = "simple")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c(
    "mean_rt", "var_rt", "n_upper",
    "n_trials", "contaminant_prop"
  ) %in% names(result)))
  expect_type(result$mean_rt, "double")
  expect_type(result$var_rt, "double")
})

test_that("ezdm_summary_stats() returns 1-row data.frame for 4par version", {
  set.seed(123)
  rt <- rgamma(200, shape = 5, rate = 10) + 0.3
  response <- rbinom(200, 1, 0.7)

  result <- ezdm_summary_stats(rt, response, version = "4par", method = "simple")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c(
    "mean_rt_upper", "mean_rt_lower",
    "var_rt_upper", "var_rt_lower", "n_upper", "n_trials",
    "contaminant_prop_upper", "contaminant_prop_lower"
  ) %in% names(result)))
})

test_that("ezdm_summary_stats() validates required arguments", {
  expect_error(
    ezdm_summary_stats(response = rbinom(10, 1, 0.5)),
    "required arguments are missing"
  )
  expect_error(
    ezdm_summary_stats(rt = rgamma(10, 5, 10)),
    "required arguments are missing"
  )
})

test_that("ezdm_summary_stats() validates rt is numeric", {
  expect_error(
    ezdm_summary_stats(rt = "not_numeric", response = c(1, 0)),
    "must be a numeric vector"
  )
  expect_error(
    ezdm_summary_stats(rt = numeric(0), response = numeric(0)),
    "has length 0"
  )
})

test_that("ezdm_summary_stats() validates rt and response have same length", {
  expect_error(
    ezdm_summary_stats(rt = c(0.5, 0.6), response = c(1, 0, 1)),
    "must have the same length"
  )
})

test_that("ezdm_summary_stats() validates parameter options", {
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  expect_error(ezdm_summary_stats(rt, response, version = "5par"), "should be one of")
  expect_error(ezdm_summary_stats(rt, response, distribution = "normal"), "should be one of")
  expect_error(ezdm_summary_stats(rt, response, method = "invalid"), "should be one of")
  expect_error(
    ezdm_summary_stats(rt, response, method = "robust", robust_scale = "invalid"),
    "should be one of"
  )
})

test_that("ezdm_summary_stats() warns for potential data issues", {
  rt_ms <- rgamma(100, shape = 5, rate = 10) * 1000
  response <- rbinom(100, 1, 0.8)

  expect_warning(
    ezdm_summary_stats(rt_ms, response, method = "simple"),
    "Some RT values > 10. Ensure RTs are in seconds"
  )
})

test_that("ezdm_summary_stats() errors for non-positive RT values", {
  rt <- c(-0.1, 0, rgamma(98, shape = 5, rate = 10) + 0.3)
  response <- rbinom(100, 1, 0.8)

  expect_error(
    ezdm_summary_stats(rt, response, method = "simple"),
    "Non-positive RT values found"
  )
})

test_that("ezdm_summary_stats() handles too few trials", {
  set.seed(123)
  rt <- rgamma(5, shape = 5, rate = 10) + 0.3
  response <- rbinom(5, 1, 0.8)

  result <- ezdm_summary_stats(rt, response, method = "simple", min_trials = 10)

  expect_true(is.na(result$mean_rt))
  expect_true(is.na(result$var_rt))
  expect_equal(result$n_trials, 5)
})

test_that("ezdm_summary_stats() simple method matches mean() and var()", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  result <- ezdm_summary_stats(rt, response, method = "simple")

  expect_equal(result$mean_rt, mean(rt), tolerance = 1e-10)
  expect_equal(result$var_rt, var(rt), tolerance = 1e-10)
  expect_equal(result$n_trials, 100L)
  expect_equal(result$n_upper, sum(response == 1))
})

test_that("ezdm_summary_stats() robust method uses median and IQR/MAD", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  result_iqr <- ezdm_summary_stats(rt, response,
    method = "robust", robust_scale = "iqr"
  )
  expect_equal(result_iqr$mean_rt, median(rt), tolerance = 1e-10)
  expect_equal(result_iqr$var_rt, (IQR(rt) / 1.349)^2, tolerance = 1e-10)
  expect_true(is.na(result_iqr$contaminant_prop))

  result_mad <- ezdm_summary_stats(rt, response,
    method = "robust", robust_scale = "mad"
  )
  expect_equal(result_mad$mean_rt, median(rt), tolerance = 1e-10)
  expect_equal(result_mad$var_rt, mad(rt)^2, tolerance = 1e-10)
})

test_that("ezdm_summary_stats() robust method is resistant to outliers", {
  set.seed(123)
  clean_rt <- rnorm(90, mean = 0.5, sd = 0.1)
  outliers <- c(0.05, 0.08, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0)
  rt <- c(clean_rt, outliers)
  response <- c(rep(1, 80), rep(0, 20))

  result_simple <- ezdm_summary_stats(rt, response, method = "simple")
  result_robust <- ezdm_summary_stats(rt, response, method = "robust")

  true_median_clean <- median(clean_rt)
  expect_true(
    abs(result_robust$mean_rt - true_median_clean) <
      abs(result_simple$mean_rt - true_median_clean)
  )
})

test_that("ezdm_summary_stats() mixture method works with different distributions", {
  set.seed(123)
  rt <- rgamma(200, shape = 5, rate = 10) + 0.3
  response <- rbinom(200, 1, 0.8)

  for (dist in c("exgaussian", "lognormal", "invgaussian")) {
    result <- ezdm_summary_stats(rt, response,
      distribution = dist, method = "mixture"
    )
    expect_s3_class(result, "data.frame")
    expect_true(is.numeric(result$mean_rt))
    expect_true(is.numeric(result$var_rt))
    expect_true(is.numeric(result$contaminant_prop))
  }
})

test_that("ezdm_summary_stats() 4par handles all correct or all errors", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3

  result <- ezdm_summary_stats(rt, rep(1, 100),
    version = "4par", method = "simple", min_trials = 10
  )
  expect_false(is.na(result$mean_rt_upper))
  expect_true(is.na(result$mean_rt_lower))

  result2 <- ezdm_summary_stats(rt, rep(0, 100),
    version = "4par", method = "simple", min_trials = 10
  )
  expect_true(is.na(result2$mean_rt_upper))
  expect_false(is.na(result2$mean_rt_lower))
})

test_that("ezdm_summary_stats() validates contaminant_bound", {
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  expect_error(
    ezdm_summary_stats(rt, response, contaminant_bound = c(3.0, 0.1)),
    "contaminant_bound\\[1\\] must be less than"
  )
  expect_error(
    ezdm_summary_stats(rt, response, contaminant_bound = c(0.1)),
    "contaminant_bound must be a vector of length 2"
  )
  expect_error(
    ezdm_summary_stats(rt, response, contaminant_bound = c("invalid", 3.0)),
    "contaminant_bound elements must be numeric or"
  )
})

test_that("ezdm_summary_stats() accepts 'min' and 'max' for contaminant_bound", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  expect_no_error(result <- ezdm_summary_stats(rt, response,
    contaminant_bound = c("min", "max")
  ))
  expect_true(is.data.frame(result))

  expect_no_error(ezdm_summary_stats(rt, response,
    contaminant_bound = c(0.1, "max")
  ))
  expect_no_error(ezdm_summary_stats(rt, response,
    contaminant_bound = c("min", 3.0)
  ))
  expect_no_error(ezdm_summary_stats(rt, response,
    contaminant_bound = c("MIN", "MAX")
  ))
})

test_that("ezdm_summary_stats() validates init_contaminant", {
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  expect_error(
    ezdm_summary_stats(rt, response, init_contaminant = 0),
    "init_contaminant must be between 0 and 1"
  )
  expect_error(
    ezdm_summary_stats(rt, response, init_contaminant = 1),
    "init_contaminant must be between 0 and 1"
  )
})

test_that("ezdm_summary_stats() validates max_contaminant", {
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  expect_error(ezdm_summary_stats(rt, response, max_contaminant = 0),
    "max_contaminant must be between 0")
  expect_error(ezdm_summary_stats(rt, response, max_contaminant = 1.5),
    "max_contaminant must be between 0")
  expect_error(ezdm_summary_stats(rt, response, max_contaminant = -0.1),
    "max_contaminant must be between 0")
  expect_no_error(ezdm_summary_stats(rt, response,
    max_contaminant = 1, method = "simple"
  ))
})

test_that("ezdm_summary_stats() validates init < max contaminant", {
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  expect_error(
    ezdm_summary_stats(rt, response, init_contaminant = 0.3, max_contaminant = 0.2),
    "init_contaminant must be less than max_contaminant"
  )
  expect_error(
    ezdm_summary_stats(rt, response, init_contaminant = 0.5, max_contaminant = 0.5),
    "init_contaminant must be less than max_contaminant"
  )
})

test_that("ezdm_summary_stats() clips contaminant proportion to max", {
  set.seed(42)
  true_rt <- rgamma(50, shape = 5, rate = 10) + 0.3
  contam_rt <- runif(50, 0.1, 3.0)
  rt <- c(true_rt, contam_rt)
  response <- rbinom(100, 1, 0.8)

  result <- suppressWarnings(
    ezdm_summary_stats(rt, response, max_contaminant = 0.3)
  )
  expect_true(result$contaminant_prop <= 0.3)

  result_strict <- suppressWarnings(
    ezdm_summary_stats(rt, response, max_contaminant = 0.1)
  )
  expect_true(result_strict$contaminant_prop <= 0.1)
})

test_that("ezdm_summary_stats() warns when contaminant proportion is clipped", {
  set.seed(42)
  true_rt <- rgamma(50, shape = 5, rate = 10) + 0.3
  contam_rt <- runif(50, 0.1, 3.0)
  rt <- c(true_rt, contam_rt)
  response <- rbinom(100, 1, 0.8)

  expect_warning(
    ezdm_summary_stats(rt, response, max_contaminant = 0.1),
    "clipped to max_contaminant"
  )
})

test_that("ezdm_summary_stats() max_contaminant works with 4par version", {
  set.seed(42)
  n_per_group <- 60
  rt <- c(
    rgamma(n_per_group, shape = 5, rate = 10) + 0.3,
    runif(20, 0.1, 3.0),
    rgamma(n_per_group, shape = 5, rate = 10) + 0.3,
    runif(20, 0.1, 3.0)
  )
  response <- c(rep(1, n_per_group + 20), rep(0, n_per_group + 20))

  result <- suppressWarnings(
    ezdm_summary_stats(rt, response, version = "4par", max_contaminant = 0.2)
  )

  expect_true(is.na(result$contaminant_prop_upper) ||
    result$contaminant_prop_upper <= 0.2)
  expect_true(is.na(result$contaminant_prop_lower) ||
    result$contaminant_prop_lower <= 0.2)
})

test_that("ezdm_summary_stats() handles character responses", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- sample(c("upper", "lower"), 100, replace = TRUE, prob = c(0.8, 0.2))

  result <- ezdm_summary_stats(rt, response, method = "simple")

  expect_equal(result$n_upper, sum(response == "upper"))
  expect_equal(result$n_trials, 100L)
})

test_that("ezdm_summary_stats() handles factor responses", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- factor(sample(c("upper", "lower"), 100, replace = TRUE, prob = c(0.8, 0.2)))

  result <- ezdm_summary_stats(rt, response, method = "simple")
  expect_equal(result$n_upper, sum(response == "upper"))
})

test_that("ezdm_summary_stats() handles correct/error and logical responses", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3

  response_ce <- sample(c("correct", "error"), 100, replace = TRUE, prob = c(0.8, 0.2))
  result_ce <- ezdm_summary_stats(rt, response_ce, method = "simple")
  expect_equal(result_ce$n_upper, sum(response_ce == "correct"))

  response_lgl <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.8, 0.2))
  result_lgl <- ezdm_summary_stats(rt, response_lgl, method = "simple")
  expect_equal(result_lgl$n_upper, sum(response_lgl))
})

test_that("ezdm_summary_stats() handles case-insensitive responses", {
  set.seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- sample(c("UPPER", "Lower", "UpPeR", "LOWER"), 100, replace = TRUE)

  result <- ezdm_summary_stats(rt, response, method = "simple")
  expect_equal(result$n_upper, sum(tolower(response) == "upper"))
})

test_that("ezdm_summary_stats() errors on unrecognized response values", {
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- sample(c("fast", "slow"), 100, replace = TRUE)

  expect_error(
    ezdm_summary_stats(rt, response),
    "Unrecognized response values"
  )
})

test_that("ezdm_summary_stats() 4par version works with character responses", {
  set.seed(123)
  rt <- rgamma(200, shape = 5, rate = 10) + 0.3
  response <- sample(c("upper", "lower"), 200, replace = TRUE, prob = c(0.7, 0.3))

  result <- ezdm_summary_stats(rt, response, version = "4par", method = "simple")

  expect_true(all(c(
    "mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower"
  ) %in% names(result)))
  expect_equal(result$n_upper, sum(response == "upper"))
  expect_false(is.na(result$mean_rt_upper))
  expect_false(is.na(result$mean_rt_lower))
})

test_that("ezdm_summary_stats() handles NAs in rt", {
  set.seed(123)
  rt <- c(rgamma(90, shape = 5, rate = 10) + 0.3, rep(NA, 10))
  response <- rbinom(100, 1, 0.8)

  result <- ezdm_summary_stats(rt, response, method = "simple")

  expect_equal(result$n_trials, 90L)
  expect_equal(result$mean_rt, mean(rt, na.rm = TRUE), tolerance = 1e-10)
})

test_that("ezdm_summary_stats() works with dplyr::reframe()", {
  skip_if_not_installed("dplyr")
  set.seed(123)
  test_data <- data.frame(
    subject = rep(1:3, each = 50),
    rt = rgamma(150, shape = 5, rate = 10) + 0.3,
    correct = rbinom(150, 1, 0.8)
  )

  result <- dplyr::group_by(test_data, subject) |>
    dplyr::reframe(ezdm_summary_stats(rt, correct, method = "simple"))

  expect_equal(nrow(result), 3)
  expect_true(all(c("subject", "mean_rt", "var_rt", "n_upper", "n_trials") %in%
    names(result)))
})

############################################################################# !
# adjust_ezdm_accuracy TESTS                                              ####
############################################################################# !

test_that("adjust_ezdm_accuracy() returns 1-row data.frame", {
  set.seed(42)
  result <- adjust_ezdm_accuracy(80, 100, 0.1)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c("n_upper_adj", "n_trials_adj") %in% names(result)))
})

test_that("adjust_ezdm_accuracy() returns integer values", {
  set.seed(42)
  result <- adjust_ezdm_accuracy(80, 100, 0.1)

  expect_type(result$n_upper_adj, "integer")
  expect_type(result$n_trials_adj, "integer")
  expect_true(result$n_upper_adj >= 0)
  expect_true(result$n_upper_adj <= result$n_trials_adj)
  expect_true(result$n_trials_adj <= 100)
})

test_that("adjust_ezdm_accuracy() handles NA contaminant_prop", {
  result <- adjust_ezdm_accuracy(80, 100, NA)

  expect_equal(result$n_upper_adj, 80L)
  expect_equal(result$n_trials_adj, 100L)
})

test_that("adjust_ezdm_accuracy() handles zero contaminant_prop", {
  result <- adjust_ezdm_accuracy(80, 100, 0)

  expect_equal(result$n_upper_adj, 80L)
  expect_equal(result$n_trials_adj, 100L)
})

test_that("adjust_ezdm_accuracy() validates inputs", {
  expect_error(adjust_ezdm_accuracy("a", 100, 0.1), "n_upper must be numeric")
  expect_error(adjust_ezdm_accuracy(80, "b", 0.1), "n_trials must be numeric")
  expect_error(adjust_ezdm_accuracy(80, 100, 0.1, guess_rate = -0.1),
    "guess_rate must be between 0 and 1")
  expect_error(adjust_ezdm_accuracy(80, 100, 0.1, guess_rate = 1.5),
    "guess_rate must be between 0 and 1")
})

############################################################################# !
# FLAG_CONTAMINANT_RTS TESTS                                              ####
############################################################################# !

# Helper to create test RT vector with known contaminants
.create_test_rt_vec <- function(n = 100, add_contaminants = TRUE, prop_contam = 0.1) {
  set.seed(123)
  if (add_contaminants) {
    n_legit <- floor((1 - prop_contam) * n)
    n_contam <- n - n_legit
    c(rgamma(n_legit, shape = 3, rate = 5), runif(n_contam, 0.05, 0.15))
  } else {
    rgamma(n, shape = 3, rate = 5)
  }
}

# Section 1: Return Structure Tests ------------------------------------------

test_that("flag_contaminant_rts() returns numeric vector", {
  rt <- .create_test_rt_vec(n = 100)

  result <- flag_contaminant_rts(rt)

  expect_type(result, "double")
  expect_equal(length(result), 100)
})

test_that("flag_contaminant_rts() attaches diagnostics attribute", {
  rt <- .create_test_rt_vec(n = 100)

  result <- flag_contaminant_rts(rt)
  diag <- attr(result, "diagnostics")

  expect_s3_class(diag, "data.frame")
  expect_equal(nrow(diag), 1)
  expect_true(all(
    c("mixture_params", "contaminant_prop", "converged") %in% names(diag)
  ))
  expect_type(diag$contaminant_prop, "double")
  expect_type(diag$converged, "logical")
})

# Section 2: Argument Validation Tests ---------------------------------------

test_that("flag_contaminant_rts() validates rt is numeric", {
  expect_error(
    flag_contaminant_rts(rt = "not_numeric"),
    "must be a numeric vector"
  )
})

test_that("flag_contaminant_rts() validates rt is non-empty", {
  expect_error(
    flag_contaminant_rts(rt = numeric(0)),
    "has length 0"
  )
})

test_that("flag_contaminant_rts() warns for RTs > 10", {
  rt <- c(15, 20, 25, 0.5, 0.6, 0.7, 0.8, 0.9)

  expect_warning(
    flag_contaminant_rts(rt),
    "Some RT values > 10"
  )
})

test_that("flag_contaminant_rts() errors for non-positive RTs", {
  rt <- c(-0.1, -0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9)

  expect_error(
    flag_contaminant_rts(rt),
    "Non-positive RT values found"
  )
})

# Section 3: NA Handling Tests ------------------------------------------------

test_that("flag_contaminant_rts() preserves NAs in output", {
  rt <- c(0.3, 0.4, NA, 0.5, 0.6, NA, 0.7, 0.8, 0.9, 1.0,
          0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 1.2, 1.3)

  result <- flag_contaminant_rts(rt)

  expect_equal(length(result), length(rt))
  expect_true(is.na(result[3]))
  expect_true(is.na(result[6]))
  expect_true(all(!is.na(result[!is.na(rt)])))
})

# Section 4: Distribution Tests -----------------------------------------------

test_that("flag_contaminant_rts() works with different distributions", {
  rt <- .create_test_rt_vec()

  result_exg <- flag_contaminant_rts(rt, distribution = "exgaussian")
  result_lnorm <- flag_contaminant_rts(rt, distribution = "lognormal")
  result_invg <- flag_contaminant_rts(rt, distribution = "invgaussian")

  diag_exg <- attr(result_exg, "diagnostics")
  diag_lnorm <- attr(result_lnorm, "diagnostics")
  diag_invg <- attr(result_invg, "diagnostics")

  expect_equal(diag_exg$distribution, "exgaussian")
  expect_equal(diag_lnorm$distribution, "lognormal")
  expect_equal(diag_invg$distribution, "invgaussian")

  expect_true(diag_exg$converged)
  expect_true(diag_lnorm$converged)
  expect_true(diag_invg$converged)
})

# Section 5: Detection Quality Tests ------------------------------------------

test_that("flag_contaminant_rts() detects simulated contaminants", {
  set.seed(456)
  n <- 200
  n_contam <- 20

  rt_legit <- rgamma(n - n_contam, shape = 3, rate = 5)
  rt_contam <- runif(n_contam, 0.05, 0.1)
  rt <- c(rt_legit, rt_contam)
  is_contam <- c(rep(FALSE, n - n_contam), rep(TRUE, n_contam))

  result <- flag_contaminant_rts(rt)

  mean_prob_contam <- mean(result[is_contam])
  mean_prob_legit <- mean(result[!is_contam])

  expect_true(mean_prob_contam > mean_prob_legit)
})

test_that("flag_contaminant_rts() contamination probabilities are in valid range", {
  rt <- .create_test_rt_vec()

  result <- flag_contaminant_rts(rt)

  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

# Section 6: Diagnostics Content Tests ----------------------------------------

test_that("flag_contaminant_rts() diagnostics contain mixture parameters", {
  rt <- .create_test_rt_vec()

  result <- flag_contaminant_rts(rt)
  diag <- attr(result, "diagnostics")

  expect_true("mixture_params" %in% names(diag))
  expect_type(diag$mixture_params, "list")

  params <- diag$mixture_params[[1]]
  expect_true(all(c("mu", "sigma", "tau") %in% names(params)))
})

test_that("flag_contaminant_rts() diagnostics report convergence correctly", {
  rt <- .create_test_rt_vec()

  result <- flag_contaminant_rts(rt)
  diag <- attr(result, "diagnostics")

  expect_true("converged" %in% names(diag))
  expect_type(diag$converged, "logical")
  expect_true("iterations" %in% names(diag))
  expect_type(diag$iterations, "integer")

  if (diag$converged) {
    expect_true("loglik" %in% names(diag))
    expect_false(is.na(diag$loglik))
  }
})

test_that("flag_contaminant_rts() diagnostics track trial counts", {
  rt <- .create_test_rt_vec(n = 150)

  result <- flag_contaminant_rts(rt)
  diag <- attr(result, "diagnostics")

  expect_equal(diag$n_trials, 150)
})

# Section 7: Small Input Tests ------------------------------------------------

test_that("flag_contaminant_rts() handles small inputs gracefully", {
  rt <- c(0.5, 0.6, 0.7)

  result <- suppressWarnings(flag_contaminant_rts(rt))

  expect_type(result, "double")
  expect_equal(length(result), 3)
})
