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
  expect_true(is.numeric(check_var_set_size("y", dat)$ss_numeric))

  # Factor with numeric levels is valid
  dat <- data.frame(y = factor(rep(c(1, 2, 3), each = 3)))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 3)
  expect_true(is.numeric(check_var_set_size("y", dat)$ss_numeric))

  # Character vector representing numbers is valid
  dat <- data.frame(y = rep(c("1", "2", "3"), each = 3))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 3)
  expect_true(is.numeric(check_var_set_size("y", dat)$ss_numeric))

  # Numeric vector with NA values is valid (assuming NA is treated correctly)
  dat <- data.frame(y = rep(c(1, 5, NA), each = 3))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 5)
  expect_true(is.numeric(check_var_set_size("y", dat)$ss_numeric))

  # Factor with NA and numeric levels is valid
  dat <- data.frame(y = factor(rep(c(1, 5, NA), each = 3)))
  expect_silent(check_var_set_size("y", dat))
  expect_equal(check_var_set_size("y", dat)$max_set_size, 5)
  expect_true(is.numeric(check_var_set_size("y", dat)$ss_numeric))
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
  # a model whose response the shared columns below cannot express (the
  # multi-column SDT responses) registers its own case in its provider
  models <- model_names()
  mls <- lapply(models, get_model)
  # test data includes variables for all model types:
  # - y, x, z, w, l, s for circular/mixture models
  # - mean_rt, var_rt, n_upper, n_trials for ezdm 3par
  # - mean_rt_upper/lower, var_rt_upper/lower for ezdm 4par
  # - stimulus (0/1) for sdt_yn (response counts come from `response`)
  # Use 50 rows to avoid small sample size warnings from cswald
  test_data <- data.frame(
    y = rep(1, 50), x = rep(1, 50), z = rep(2, 50), w = rep(1, 50),
    s = rep(2, 50), l = rep(1, 50),
    mean_rt = rep(0.5, 50), var_rt = rep(0.02, 50),
    n_upper = rep(80, 50), n_trials = rep(100, 50),
    mean_rt_upper = rep(0.45, 50), mean_rt_lower = rep(0.55, 50),
    var_rt_upper = rep(0.018, 50), var_rt_lower = rep(0.025, 50),
    rt = rep(0.6, 50), response = rep(1, 50),
    stimulus = rep(c(0L, 1L), 25)
  )
  registered <- model_test_cases("check_data")
  for (m in seq_along(models)) {
    case <- registered[[models[m]]]
    if (!is.null(case)) {
      expect_s3_class(check_data(case$model, case$data, bmf(kappa ~ 1)), "data.frame")
      next
    }
    ml <- mls[[m]]
    model <- ml(
      resp_error = "y", nt_features = "x", set_size = 2,
      nt_distances = "z", resp_cats = c("w", "l"), num_options = c(1, 1),
      mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
      n_trials = "n_trials", rt = "rt", response = "response",
      stimulus = "stimulus", m = 2,
      trees = mpt_tree("t", list(w = "p", l = "1 - p"))
    )
    expect_s3_class(
      check_data(model, test_data, bmf(kappa ~ 1)),
      "data.frame"
    )
  }
})

# A model's check_data() adds columns that pp_check(group = ) and
# conditional_effects() must not offer as the user's predictors, and that
# formula checks must not call unknown. built_data_columns() is the one
# declaration, so a column added without it fails here (#550). The reverse
# fails too: a declared name check_data() does not build hides a user column
# of that name and lets a formula naming a missing column pass the data
# check. It is compared with all returned columns, since a declared column may
# overwrite one the user supplied. Cases of the check_data kind are read too,
# for models whose stored-frame case cannot show a column (the ezdm 4par
# indicators).
test_that("built_data_columns() declares exactly the columns check_data() adds", {
  cases <- c(model_test_cases("stored_frame"), model_test_cases("check_data"))
  for (name in names(cases)) {
    case <- cases[[name]]
    checked <- suppressWarnings(suppressMessages(
      check_data(case$model, case$data, case$formula %||% bmf(kappa ~ 1))
    ))
    # Y is a matrix column that is no predictor, left out of the declaration
    added <- setdiff(names(checked), c(names(case$data), "Y"))
    declared <- built_data_columns(case$model)
    expect_identical(
      setdiff(added, declared), character(0),
      info = paste0(name, ": check_data() adds a column built_data_columns() does not declare")
    )
    expect_identical(
      setdiff(declared, names(checked)), character(0),
      info = paste0(name, ": built_data_columns() declares a column check_data() does not build")
    )
  }
})

test_that("built_data_columns() names what each model's check_data() builds", {
  cases <- model_test_cases("stored_frame")
  declared <- function(name) built_data_columns(cases[[name]]$model)

  expect_setequal(declared("mixture2p"), character(0))
  expect_setequal(declared("sdm"), character(0))
  expect_setequal(
    declared("mixture3p"),
    c(paste0("LureIdx", 1:7), "inv_ss", "ss_numeric")
  )
  expect_setequal(declared("imm"), c(paste0("LureIdx", 1:7), "inv_ss", "ss_numeric"))
  expect_setequal(declared("sdt_yn"), "dist_type")
  expect_setequal(declared("sdt_mafc"), c("m_afc", "dist_type"))
  expect_setequal(declared("sdt_ranking"), c("nTrials", "max_rank"))
  expect_setequal(declared("sdt_rating"), "nTrials")
  expect_setequal(declared("sdt_cdp"), "nTrials")
  expect_setequal(
    declared("mpt_impossible"),
    c("nTrials", "Idx_withdist", "Idx_nodist", "Poss_dist")
  )
  expect_setequal(declared("mpt_single_tree"), "nTrials")
  expect_setequal(
    declared("m3_num_options_by_category"),
    c("nTrials", paste0("Idx_", c("corr", "other", "dist", "npl")),
      paste0("n_opt_", c("corr", "other", "dist", "npl")))
  )

  ezdm_4par <- model_test_cases("check_data")$ezdm_4par$model
  expect_setequal(built_data_columns(ezdm_4par), c("rt_used_upper", "rt_used_lower"))
})

test_that("wrap(x) returns the same for values between -pi and pi", {
  withr::local_seed(445)
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

test_that("wrap() refuses a radians argument that is not TRUE or FALSE", {
  expect_error(wrap(1, radians = "yes"), "must be TRUE or FALSE")
  expect_error(wrap(1, radians = NA), "must be TRUE or FALSE")
  expect_error(wrap(1, radians = c(TRUE, FALSE)), "must be TRUE or FALSE")
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

test_that("standata() works with bmf", {
  ff <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  dat <- small_data()
  sd <- standata(ff, dat, mixture3p(
    resp_error = "dev_rad",
    nt_features = "col_nt",
    set_size = "set_size", regex = T
  ))
  expect_equal(class(sd)[1], "standata")
})

test_that("standata() returns a standata class", {
  withr::local_seed(445)
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
  data <- small_data()
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

test_that("calc_error_relative_to_nontargets() gives one row per non-target with the wrapped error", {
  data <- data.frame(id = 1:2, y = c(0.5, -3), nt1 = c(0.2, 3), nt2 = c(-0.4, 0))
  out <- calc_error_relative_to_nontargets(data, "y", c("nt1", "nt2"))
  expect_equal(nrow(out), 4)
  expect_false(any(c("nt1", "nt2") %in% names(out)))
  expect_equal(out$id, c(1, 2, 1, 2))
  # -3 - 3 = -6 lies outside [-pi, pi] and wraps to -6 + 2 * pi
  expect_equal(out$y_nt, c(0.3, -6 + 2 * pi, 0.9, -3))
})

# The RT models (ddm, cswald) take the same trial data, so they have to accept
# and refuse it alike (#557). One row per input; a row gives one expectation
# for all cases, or one per case where the models differ on purpose (the 0.100 s
# advice, cswald's error-rate warning). The messages are written out. The
# trial data have a 10% error rate, which cswald "simple" does not warn about.
rt_check_rt <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.55, 0.65, 0.75, 0.85, 0.95)
rt_check_response <- c(1, 0, 1, 1, 1, 1, 1, 1, 1, 1)
rt_check_chars <- c("upper", "LOWER", "Upper", rep("upper", 7))

rt_check_data <- function(rt = rt_check_rt, response = rt_check_response) {
  data.frame(rt = rt, response = response)
}

rt_check_error <- function(message) list(error = message)

# the warnings are given in the order they are raised
rt_check_ok <- function(response = rt_check_response, warnings = character(0),
                        integer = NULL) {
  list(response = response, warnings = warnings, integer = integer)
}

rt_check_rows <- function() {
  na_rt <- "Some values in rt were NA. These were removed from the analysis."
  na_response <- "Some values in response were NA. These were removed from the analysis."
  fast <- "Your data contains reaction times smaller than 0.100 seconds."
  fast_ddm <- paste(fast, "It is likely that the model will not be able to sample with the current settings of the initial values.")
  fast_cswald <- paste(fast, "Such fast responses are often anticipations and pull the estimate of `ndt` down; consider filtering them.")
  slow <- "Your data contains reaction times larger than 10 seconds."
  boolean <- "The response variable is boolean and will be internally transformed to an integer variable with values 0 for FALSE and 1 for TRUE."
  character <- "The response variable is a character variable and will be internally transformed to an integer variable with 0 for 'lower' and 1 for 'upper'."
  not_positive <- "Some reaction times are zero or negative, please check your data. The likelihood is zero for a response time that is not strictly greater than the non-decision time."
  not_numeric <- "The RT variable 'rt' needs to be of type double or integer."
  not_binary <- "The response variable 'response' contains values other than 0 and 1."
  no_rows <- "Every row has a missing value in 'rt' or 'response', so no observations are left."
  resp <- rt_check_response
  rt <- rt_check_rt

  list(
    valid = list(data = rt_check_data(), all = rt_check_ok(integer = TRUE)),
    no_rt_column = list(
      data = rt_check_data()["response"],
      all = rt_check_error("The RT variable 'rt' is not present in the data.")
    ),
    no_response_column = list(
      data = rt_check_data()["rt"],
      all = rt_check_error("The response variable 'response' is not present in the data.")
    ),
    rt_na = list(
      data = rt_check_data(rt = replace(rt, 1, NA)),
      all = rt_check_ok(resp[-1], na_rt)
    ),
    response_na = list(
      data = rt_check_data(response = replace(resp, 3, NA)),
      all = rt_check_ok(resp[-3], na_response)
    ),
    both_na = list(
      data = rt_check_data(rt = replace(rt, 1, NA), response = replace(resp, 3, NA)),
      all = rt_check_ok(resp[-c(1, 3)], c(na_rt, na_response))
    ),
    all_rt_na = list(data = rt_check_data(rt = NA_real_), all = rt_check_error(no_rows)),
    all_rt_na_logical = list(data = rt_check_data(rt = NA), all = rt_check_error(no_rows)),
    all_response_na = list(
      data = rt_check_data(response = NA_character_),
      all = rt_check_error(no_rows)
    ),
    logical_response_na = list(
      data = rt_check_data(response = replace(resp == 1, 3, NA)),
      all = rt_check_ok(resp[-3], c(na_response, boolean), integer = TRUE)
    ),
    rt_character = list(
      data = rt_check_data(rt = as.character(rt)),
      all = rt_check_error(not_numeric)
    ),
    rt_factor = list(
      data = rt_check_data(rt = factor(rt)),
      all = rt_check_error(not_numeric)
    ),
    rt_difftime = list(
      data = rt_check_data(rt = as.difftime(rt, units = "secs")),
      all = rt_check_error(not_numeric)
    ),
    rt_negative = list(
      data = rt_check_data(rt = replace(rt, 1, -0.1)),
      all = rt_check_error(not_positive)
    ),
    rt_zero = list(
      data = rt_check_data(rt = replace(rt, 1, 0)),
      all = rt_check_error(not_positive)
    ),
    rt_fast = list(
      data = rt_check_data(rt = replace(rt, 1, 0.05)),
      ddm = rt_check_ok(resp, fast_ddm),
      cswald = rt_check_ok(resp, fast_cswald), cswald_crisk = rt_check_ok(resp, fast_cswald)
    ),
    rt_infinite = list(
      data = rt_check_data(rt = replace(rt, 1, Inf)),
      all = rt_check_ok(resp, slow)
    ),
    rt_milliseconds = list(
      data = rt_check_data(rt = rt * 1000),
      all = rt_check_ok(resp, slow)
    ),
    response_2_3 = list(
      data = rt_check_data(response = replace(resp, 1:2, c(2, 3))),
      all = rt_check_error(not_binary)
    ),
    response_integer = list(
      data = rt_check_data(response = as.integer(resp)),
      all = rt_check_ok(resp, integer = TRUE)
    ),
    response_logical = list(
      data = rt_check_data(response = resp == 1),
      all = rt_check_ok(resp, boolean, integer = TRUE)
    ),
    response_character = list(
      data = rt_check_data(response = rt_check_chars),
      all = rt_check_ok(resp, character, integer = TRUE)
    ),
    response_character_invalid = list(
      data = rt_check_data(response = replace(rt_check_chars, 2, "up")),
      all = rt_check_error("The response variable 'response' contains invalid character values.")
    ),
    response_character_01 = list(
      data = rt_check_data(response = as.character(resp)),
      all = rt_check_error("The response variable 'response' contains invalid character values.")
    ),
    response_factor_labels = list(
      data = rt_check_data(response = factor(tolower(rt_check_chars))),
      all = rt_check_ok(resp, character, integer = TRUE)
    ),
    response_factor_01 = list(
      data = rt_check_data(response = factor(resp)),
      all = rt_check_error("The response variable 'response' contains invalid character values.")
    ),
    response_complex = list(
      data = rt_check_data(response = as.complex(resp)),
      all = rt_check_error("The response variable 'response' is of type 'complex'.")
    ),
    many_errors = list(
      data = rt_check_data(response = c(0, 0, 0, 0, rep(1, 6))),
      ddm = rt_check_ok(c(0, 0, 0, 0, rep(1, 6))),
      cswald = rt_check_ok(c(0, 0, 0, 0, rep(1, 6)), "Your data has an error rate of 40%."),
      cswald_crisk = rt_check_ok(c(0, 0, 0, 0, rep(1, 6)))
    )
  )
}

squish_text <- function(x) gsub("\\s+", " ", trimws(x))

# what check_data() returned, raised and signalled on one input
run_rt_check <- function(case, data) {
  warnings <- character(0)
  messages <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      check_data(case$model, data, case$formula),
      warning = function(w) {
        warnings <<- c(warnings, squish_text(conditionMessage(w)))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- c(messages, squish_text(conditionMessage(m)))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  list(value = value, warnings = warnings, messages = messages)
}

test_that("the RT models accept and refuse RTs and responses alike", {
  cases <- model_test_cases("rt_data")
  rows <- rt_check_rows()
  expect_setequal(names(cases), c("ddm", "cswald", "cswald_crisk"))
  for (case_name in names(cases)) {
    for (row_name in names(rows)) {
      row <- rows[[row_name]]
      want <- row[[case_name]] %||% row[["all"]]
      label <- paste(case_name, row_name)
      expect_false(is.null(want), label = paste(label, "has an expectation"))
      got <- run_rt_check(cases[[case_name]], row$data)
      expect_length(got$messages, 0)
      if (!is.null(want$error)) {
        expect_true(inherits(got$value, "error"), label = paste(label, "stops"))
        if (inherits(got$value, "error")) {
          expect_match(squish_text(conditionMessage(got$value)), want$error, fixed = TRUE,
                       info = label)
        }
        next
      }
      expect_false(inherits(got$value, "error"), label = paste(label, "passes"),
                   info = if (inherits(got$value, "error")) conditionMessage(got$value))
      if (inherits(got$value, "error")) next
      expect_length(got$warnings, length(want$warnings))
      for (i in seq_along(want$warnings)) {
        expect_match(got$warnings[i], want$warnings[i], fixed = TRUE, info = label)
      }
      expect_equal(got$value$response, want$response, info = label)
      expect_equal(nrow(got$value), length(want$response), info = label)
      if (!is.null(want$integer)) {
        expect_identical(is.integer(got$value$response), want$integer, info = label)
      }
    }
  }
})

# The table squishes whitespace, so the paragraph break in ddm's two-part
# advice is pinned here.
test_that("ddm() keeps the blank line in its 0.100 s warning", {
  case <- model_test_cases("rt_data")$ddm
  expect_warning(
    check_data(case$model, rt_check_data(rt = replace(rt_check_rt, 1, 0.05)), case$formula),
    paste0(
      "Your data contains reaction times smaller than 0.100 seconds.\n\n",
      "It is likely that the model will not be able to sample with the current ",
      "settings of the initial values.\n\n",
      "Either pass your own initial value function or consider filtering ",
      "reaction times below 0.100 seconds."
    ),
    fixed = TRUE
  )
})
