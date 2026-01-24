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
  # Use 50 rows to avoid small sample size warnings from cswald
  test_data <- data.frame(
    y = rep(1, 50), x = rep(1, 50), z = rep(2, 50), w = rep(1, 50), 
    s = rep(2, 50), l = rep(1, 50),
    mean_rt = rep(0.5, 50), var_rt = rep(0.02, 50), 
    n_upper = rep(80, 50), n_trials = rep(100, 50),
    mean_rt_upper = rep(0.45, 50), mean_rt_lower = rep(0.55, 50),
    var_rt_upper = rep(0.018, 50), var_rt_lower = rep(0.025, 50),
    rt = rep(0.6, 50), response = rep(1, 50)
  )
  for (ml in mls) {
    model <- ml(
      resp_error = "y", nt_features = "x", set_size = 2,
      nt_distances = "z", resp_cats = c("w", "l"), num_options = c(1, 1),
      mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
      n_trials = "n_trials", rt = "rt", response = "response"
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
