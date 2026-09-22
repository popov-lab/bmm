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
  expect_equal(try_read_bmmfit(paste0(file, ".rds")), mock_fit,
    ignore_function_env = TRUE,
    ignore_formula_env = TRUE
  )

  x <- 1
  saveRDS(x, paste0(file, ".rds"))
  expect_error(try_read_bmmfit(paste0(file, ".rds")), "not of class 'bmmfit'")
})

test_that("validate_file_refit normalizes the user-facing values", {
  expect_equal(validate_file_refit(FALSE), "never")
  expect_equal(validate_file_refit(TRUE), "always")
  expect_equal(validate_file_refit("never"), "never")
  expect_equal(validate_file_refit("Always"), "always")
  expect_equal(validate_file_refit("On_Change"), "on_change")
  expect_error(validate_file_refit("sometimes"), "invalid option")
  expect_error(validate_file_refit(1), "invalid option")

  # a logical that is not a single TRUE/FALSE must error rather than degrade to
  # "never", which would silently return a stale fit
  expect_error(validate_file_refit(NA), "invalid option")
  expect_error(validate_file_refit(c(TRUE, TRUE)), "invalid option")
  expect_error(validate_file_refit(logical(0)), "invalid option")
  expect_error(validate_file_refit(NA_character_), "invalid option")
  expect_error(validate_file_refit(character(0)), "invalid option")
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

  # a mixed-case "always" must refit too, not silently read the cache
  mock_fit4 <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 4, rename = F,
    file = file, file_refit = "Always"
  )
  expect_equal(mock_fit4$fit, 4)
})

# `mock_fit` doubles as the identity of the fit: the cached object carries the
# value it was fitted with, so a returned fit whose $fit differs from the value
# passed to the current call is the cached one
test_that('file_refit = "on_change" returns the cached fit while nothing changes', {
  withr::local_options(bmm.sort_data = FALSE)
  file <- tempfile()
  cached <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = F,
    file = file, file_refit = "on_change"
  )
  expect_equal(cached$fit, 1)

  same <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 2, rename = F,
    file = file, file_refit = "on_change"
  )
  expect_equal(same$fit, 1)
  expect_equal(same$file, paste0(file, ".rds"))
})

test_that('file_refit = "on_change" refits when the model changes', {
  withr::local_options(bmm.sort_data = FALSE)
  fit_cache <- function(formula, data, ..., file, mock_fit) {
    bmm(formula, data, sdm("dev_rad"),
      backend = "mock", mock_fit = mock_fit, rename = F,
      file = file, file_refit = "on_change", ...
    )
  }
  data <- oberauer_lin_2017
  ff <- bmf(c ~ 1, kappa ~ 1)

  file <- tempfile()
  fit_cache(ff, data, file = file, mock_fit = 1)
  expect_equal(
    fit_cache(bmf(c ~ 0 + set_size, kappa ~ 1), data, file = file, mock_fit = 2)$fit,
    2
  )

  # a row subset leaves the Stan code untouched and is seen through sdata alone
  file <- tempfile()
  fit_cache(ff, data, file = file, mock_fit = 1)
  expect_equal(fit_cache(ff, data[1:100, ], file = file, mock_fit = 3)$fit, 3)

  # a prior change is the mirror case: identical sdata, different Stan code
  file <- tempfile()
  fit_cache(ff, data, file = file, mock_fit = 1)
  expect_equal(
    fit_cache(ff, data,
      file = file, mock_fit = 4,
      prior = brms::set_prior("normal(0, 0.1)", class = "Intercept", dpar = "kappa")
    )$fit,
    4
  )

  # renaming factor levels changes neither the Stan code nor the Stan data, so
  # only the data handed to brms::brmsfit_needs_refit() can catch it
  data$set_size <- factor(data$set_size)
  renamed <- data
  levels(renamed$set_size) <- paste0("ss", levels(renamed$set_size))
  ff_ss <- bmf(c ~ 0 + set_size, kappa ~ 1)
  file <- tempfile()
  fit_cache(ff_ss, data, file = file, mock_fit = 1)
  expect_equal(fit_cache(ff_ss, renamed, file = file, mock_fit = 5)$fit, 5)

  file <- tempfile()
  fit_cache(ff, data, file = file, mock_fit = 1)
  expect_equal(
    fit_cache(ff, data, file = file, mock_fit = 6, algorithm = "meanfield")$fit,
    6
  )

  # threading is the bmm-specific hazard brms does not have: a loop = FALSE
  # likelihood chunk is configured for the threading spec, so a cache hit across
  # a change of `threads` would return a fit whose likelihood is sliced
  # differently from the one the current call would compile
  file <- tempfile()
  fit_cache(ff, data, file = file, mock_fit = 1)
  expect_equal(fit_cache(ff, data, file = file, mock_fit = 7, threads = 2)$fit, 7)

  file <- tempfile()
  fit_cache(ff, data, file = file, mock_fit = 1, threads = 2)
  expect_equal(fit_cache(ff, data, file = file, mock_fit = 8)$fit, 8)
})

test_that('file_refit = "on_change" falls back when the cached fit has no algorithm', {
  withr::local_options(bmm.sort_data = FALSE)
  fit_cache <- function(mock_fit) {
    bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
      backend = "mock", mock_fit = mock_fit, rename = F,
      file = file, file_refit = "on_change"
    )
  }
  file <- tempfile()
  cached <- fit_cache(1)

  # brms compares the algorithm under a bare stopifnot(!is.null(fit$algorithm)),
  # so without the guard this errors instead of using the other three channels
  cached$algorithm <- NULL
  saveRDS(cached, paste0(file, ".rds"))
  expect_equal(fit_cache(2)$fit, 1)
})

test_that('file_refit = "on_change" refits a cached fit that cannot be restructured', {
  withr::local_options(bmm.sort_data = FALSE)
  fit_cache <- function(mock_fit) {
    bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
      backend = "mock", mock_fit = mock_fit, rename = F,
      file = file, file_refit = "on_change"
    )
  }
  file <- tempfile()
  cached <- fit_cache(1)

  # a pre-0.3.0 fit whose family carries no environment is what restructure()
  # gives up on; "on_change" must do the refit it asks for rather than abort
  cached$version$bmm_restructure <- NULL
  cached$version$bmm <- as.package_version("0.2.1")
  cached$family$env <- NULL
  saveRDS(cached, paste0(file, ".rds"))
  expect_error(restructure(readRDS(paste0(file, ".rds"))), "Unable to restructure")
  expect_equal(fit_cache(2)$fit, 2)
})

test_that('file_refit = "on_change" rejects a cached file that is not a bmmfit', {
  withr::local_options(bmm.sort_data = FALSE)
  file <- tempfile()
  saveRDS(1, paste0(file, ".rds"))
  expect_error(
    bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
      backend = "mock", mock_fit = 1, rename = F,
      file = file, file_refit = "on_change"
    ),
    "not of class 'bmmfit'"
  )
})

test_that('file_refit = "on_change" reports why it refits at silent = 0', {
  withr::local_options(bmm.sort_data = FALSE)
  file <- tempfile()
  suppressMessages(bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = F,
    file = file, file_refit = "on_change"
  ))
  expect_message(
    bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
      backend = "mock", mock_fit = 2, rename = F, silent = 0,
      file = file, file_refit = "on_change",
      prior = brms::set_prior("normal(0, 0.1)", class = "Intercept", dpar = "kappa")
    ),
    "Stan code has changed"
  )
})

test_that('bmm_options() accepts and applies file_refit = "on_change"', {
  withr::local_options(bmm.sort_data = FALSE)
  old_op <- suppressMessages(bmm_options(file_refit = "on_change"))
  withr::defer(options(old_op))
  expect_equal(getOption("bmm.file_refit"), "on_change")
  expect_error(suppressMessages(bmm_options(file_refit = "sometimes")), "invalid option")

  file <- tempfile()
  bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = F, file = file
  )
  same <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 2, rename = F, file = file
  )
  expect_equal(same$fit, 1)
  changed <- bmm(bmf(c ~ 0 + set_size, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 3, rename = F, file = file
  )
  expect_equal(changed$fit, 3)

  # bmm_options() stores the value as given rather than the normalized string,
  # so the logical forms have to survive the round trip through bmm() as well
  logical_op <- suppressMessages(bmm_options(file_refit = TRUE))
  withr::defer(options(logical_op))
  expect_true(getOption("bmm.file_refit"))
  refit <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 4, rename = F, file = file
  )
  expect_equal(refit$fit, 4)
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

test_that("configure_control() adds the starting step size under the user's control list", {
  withr::local_options(bmm.step_size = 0.1)
  expect_equal(configure_control(NULL, "cmdstanr"), list(step_size = 0.1))
  expect_equal(
    configure_control(list(adapt_delta = 0.95), "cmdstanr"),
    list(adapt_delta = 0.95, step_size = 0.1)
  )
  # the user's own value wins under either spelling
  expect_equal(configure_control(list(step_size = 0.5), "cmdstanr"), list(step_size = 0.5))
  expect_equal(configure_control(list(stepsize = 0.5), "cmdstanr"), list(stepsize = 0.5))
  # rstan spells the argument without the underscore
  expect_equal(configure_control(NULL, "rstan"), list(stepsize = 0.1))
  expect_equal(configure_control(list(stepsize = 0.5), "rstan"), list(stepsize = 0.5))

  withr::local_options(bmm.step_size = 0.02)
  expect_equal(configure_control(NULL, "cmdstanr"), list(step_size = 0.02))

  withr::local_options(bmm.step_size = FALSE)
  expect_null(configure_control(NULL, "cmdstanr"))
  expect_equal(configure_control(list(adapt_delta = 0.95), "cmdstanr"), list(adapt_delta = 0.95))
})

test_that("bmm_options(step_size = ) validates and applies the option", {
  withr::defer(suppressMessages(bmm_options(reset_options = TRUE)))
  expect_error(bmm_options(step_size = -1), "step_size")
  expect_error(bmm_options(step_size = "a"), "step_size")
  expect_error(bmm_options(step_size = c(0.1, 0.2)), "step_size")
  expect_message(bmm_options(step_size = 0.3), "step_size = 0.3")
  expect_equal(getOption("bmm.step_size"), 0.3)
  suppressMessages(bmm_options(step_size = FALSE))
  expect_false(getOption("bmm.step_size"))
  suppressMessages(bmm_options(reset_options = TRUE))
  expect_equal(getOption("bmm.step_size"), 0.01)
})
