test_that("init argument is overwritten if the user supplies it", {
  config_args <- list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1)
  dots <- list(init = 2)
  out <- combine_args(nlist(config_args, dots))
  expect_equal(out, list(formula = "a", family = "b", data = "d", stanvars = "e", init = 2))
})

test_that("user cannot overwrite the custom family", {
  config_args <- list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1)
  dots <- list(family = "c")
  expect_error(combine_args(nlist(config_args, dots)), "Unsupported argument")
})

test_that("empty dots don't crash the function", {
  config_args <- list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1)
  out <- combine_args(nlist(config_args))
  expect_equal(out, list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1))
})

test_that("missing arguments in models are handled correctly", {
  expect_error(mixture2p(), "arguments are missing in mixture2p\\(\\)\\: resp_error")
  expect_error(sdm(), "arguments are missing in sdm\\(\\)\\: resp_error")
  expect_error(mixture3p("y"), "arguments are missing in mixture3p\\(\\)\\: nt_features, set_size")
  expect_error(mixture3p(set_size = "y"), "arguments are missing in mixture3p\\(\\)\\: resp_error, nt_features")
})

test_that("get_variables works", {
  expect_equal(get_variables("a", c("a", "b", "c")), "a")
  expect_equal(get_variables("a", c("a", "b", "c"), regex = TRUE), "a")
  expect_equal(get_variables("a", c("a", "b", "c"), regex = FALSE), "a")
  expect_equal(get_variables("a|b", c("a", "b", "c"), regex = TRUE), c("a", "b"))
  expect_equal(
    get_variables("abc", c("abc1", "abc2", "abc3", "other"), regex = TRUE),
    c("abc1", "abc2", "abc3")
  )
  expect_equal(
    get_variables("^abc", c("abc1", "abc2", "abc3", "other_abc4"), regex = TRUE),
    c("abc1", "abc2", "abc3")
  )
  expect_equal(
    get_variables("abc$", c("nt1_abc", "nt2_abc", "nt3_abc", "other_abc4"), regex = TRUE),
    c("nt1_abc", "nt2_abc", "nt3_abc")
  )
  expect_equal(
    get_variables("nt.*_abc", c("nt1_abc", "nt2_abc", "nt3_abc", "other_abc4"), regex = TRUE),
    c("nt1_abc", "nt2_abc", "nt3_abc")
  )
  expect_equal(get_variables("a|b", c("a", "b", "c"), regex = FALSE), "a|b")
  expect_error(get_variables("d", c("a", "b", "c"), regex = TRUE))
})

test_that("bmm_options works", {
  withr::defer(suppressMessages(bmm_options()))
  expect_message(bmm_options(), "Current bmm options")
  expect_message(bmm_options(sort_data = TRUE), "sort_data = TRUE")
  expect_equal(getOption("bmm.sort_data"), TRUE)
  op <- suppressMessages(bmm_options(sort_data = FALSE))
  expect_equal(getOption("bmm.sort_data"), FALSE)
  options(op)
  expect_equal(getOption("bmm.sort_data"), TRUE)
})


test_that("check_rds_file works", {
  good_files <- list("a.rds", "abc/a.rds", "a", "abc/a", "a.M")
  bad_files <- list(1, mean, c("a", "b"), TRUE)

  for (f in good_files) {
    expect_silent(res <- check_rds_file(f))
    expect_equal(fs::path_ext(res), "rds")
  }

  for (f in bad_files) {
    expect_error(check_rds_file(f))
  }

  expect_null(check_rds_file(NULL))
})

test_that("try_read_bmmfit works", {
  withr::local_options(bmm.sort_data = FALSE)
  mock_fit <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = F
  )
  file <- tempfile()
  mock_fit$file <- paste0(file, ".rds")
  saveRDS(mock_fit, paste0(file, ".rds"))
  expect_equal(try_read_bmmfit(file, FALSE), mock_fit,
    ignore_function_env = TRUE,
    ignore_formula_env = TRUE
  )

  x <- 1
  saveRDS(x, paste0(file, ".rds"))
  expect_error(try_read_bmmfit(file, FALSE), "not of class 'bmmfit'")
})

test_that("try_save_bmmfit works", {
  withr::local_options(bmm.sort_data = FALSE)
  file <- tempfile()
  mock_fit <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = F,
    file = file
  )
  rds_file <- paste0(file, ".rds")
  expect_true(file.exists(rds_file))
  expect_equal(readRDS(rds_file), mock_fit, ignore_function_env = TRUE, ignore_formula_env = TRUE)

  mock_fit2 <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 2, rename = F,
    file = file
  )
  expect_equal(mock_fit, mock_fit2, ignore_attr = TRUE)

  # they should not be the same if file_refit = TRUE
  mock_fit3 <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 3, rename = F,
    file = file, file_refit = TRUE
  )
  expect_error(expect_equal(mock_fit, mock_fit3))
})

test_that("is_namedlist works", {
  expect_true(is_namedlist(list(a = 1)))
  expect_true(is_namedlist(list(a = 1, b = 2)))
  expect_true(is_namedlist(nlist()))
  expect_true(is_namedlist(nlist(y ~ 1)))
  arg <- "hello"
  expect_true(is_namedlist(nlist(arg)))

  expect_false(is_namedlist(list(a = 1, 2)))
  expect_false(is_namedlist(list(1, 2)))
  expect_false(is_namedlist(list()))
  expect_false(is_namedlist(list(y ~ 1)))

})

test_that("softmax produces valid probability distributions", {
  # Basic test - output should sum to 1
  eta <- c(1, 2, 3)
  result <- softmax(eta)
  expect_equal(sum(result), 1)
  expect_true(all(result > 0))
  expect_true(all(result < 1))

  # Test with different lambda values
  result_lambda2 <- softmax(eta, lambda = 2)
  expect_equal(sum(result_lambda2), 1)
  expect_true(all(result_lambda2 > 0))

  # Higher lambda should increase difference between probabilities
  expect_true(max(result_lambda2) > max(result))
  expect_true(min(result_lambda2) < min(result))

  # Test with negative values
  eta_neg <- c(-2, -1, 0, 1, 2)
  result_neg <- softmax(eta_neg)
  expect_equal(sum(result_neg), 1)
  expect_true(all(result_neg > 0))

  # Test with all equal values - should produce uniform distribution
  eta_equal <- rep(5, 4)
  result_equal <- softmax(eta_equal)
  expect_equal(result_equal, rep(0.25, 4))

  # Test with single value
  result_single <- softmax(10)
  expect_equal(result_single, 1)

  # Test ordering is preserved
  eta_ordered <- 1:5
  result_ordered <- softmax(eta_ordered)
  expect_true(all(diff(result_ordered) > 0)) # should be monotonically increasing
})

test_that("softmax with extreme values doesn't overflow", {
  # Very large values
  eta_large <- c(100, 200, 300)
  result_large <- softmax(eta_large)
  expect_equal(sum(result_large), 1)
  expect_false(any(is.na(result_large)))
  expect_false(any(is.infinite(result_large)))

  # Very small values
  eta_small <- c(-300, -200, -100)
  result_small <- softmax(eta_small)
  expect_equal(sum(result_small), 1)
  expect_false(any(is.na(result_small)))
  expect_false(any(is.infinite(result_small)))
})

test_that("softmaxinv is the inverse of softmax", {
  # Test basic inverse relationship
  eta <- 5:7
  p <- softmax(eta)
  eta_recovered <- softmaxinv(p, ref_position = 1, ref_value = 5)
  expect_equal(eta_recovered, eta, tolerance = 1e-10)

  # Test with different reference positions
  eta2 <- c(2, 4, 6, 8)
  p2 <- softmax(eta2)

  # Reference at position 1
  eta_rec1 <- softmaxinv(p2, ref_position = 1, ref_value = 2)
  expect_equal(eta_rec1, eta2, tolerance = 1e-10)

  # Reference at position 2
  eta_rec2 <- softmaxinv(p2, ref_position = 2, ref_value = 4)
  expect_equal(eta_rec2, eta2, tolerance = 1e-10)

  # Reference at position 4
  eta_rec4 <- softmaxinv(p2, ref_position = 4, ref_value = 8)
  expect_equal(eta_rec4, eta2, tolerance = 1e-10)

  # Test with different lambda values
  eta3 <- c(1, 3, 5)
  lambda_val <- 2
  p3 <- softmax(eta3, lambda = lambda_val)
  eta_rec3 <- softmaxinv(p3, lambda = lambda_val, ref_position = 1, ref_value = 1)
  expect_equal(eta_rec3, eta3, tolerance = 1e-10)
})

test_that("softmaxinv with default parameters", {
  # Default reference position is length(p) and ref_value is 0
  eta <- c(1, 2, 3)
  p <- softmax(eta)

  # With defaults, last position should be 0
  eta_recovered <- softmaxinv(p)
  expect_equal(eta_recovered[length(eta_recovered)], 0, tolerance = 1e-10)

  # The differences should be preserved
  eta_shifted <- eta - eta[length(eta)]
  expect_equal(eta_recovered, eta_shifted, tolerance = 1e-10)
})

test_that("softmaxinv handles edge cases", {
  # Length 1 probability vector
  result <- softmaxinv(1)
  expect_equal(result, numeric(0))

  # Length 2 probability vector
  p2 <- c(0.3, 0.7)
  result2 <- softmaxinv(p2)
  expect_length(result2, 2)
  expect_equal(sum(softmax(result2) - p2), 0, tolerance = 1e-10)
})

test_that("softmaxinv validates inputs correctly", {
  # ref_position must be a single value
  expect_error(
    softmaxinv(c(0.2, 0.3, 0.5), ref_position = c(1, 2)),
    "single reference value"
  )

  # ref_position must be within valid range
  expect_error(
    softmaxinv(c(0.2, 0.3, 0.5), ref_position = 4),
    "less or equal than the length"
  )
})

test_that("softmax and softmaxinv work with example from documentation", {
  # Example from the documentation
  result <- softmax(5:7)
  recovered <- softmaxinv(result, ref_position = 1, ref_value = 5)
  expect_equal(recovered, 5:7, tolerance = 1e-10)
})
