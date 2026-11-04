############################################################################# !
# ezdm_summary_stats TESTS                                                ####
############################################################################# !

test_that("ezdm_summary_stats() returns 1-row data.frame for 3par version", {
  withr::local_seed(123)
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
  withr::local_seed(123)
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
  withr::local_seed(445)
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
  withr::local_seed(445)
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
  withr::local_seed(445)
  rt_ms <- rgamma(100, shape = 5, rate = 10) * 1000
  response <- rbinom(100, 1, 0.8)

  expect_warning(
    ezdm_summary_stats(rt_ms, response, method = "simple"),
    "Some RT values > 10. Ensure RTs are in seconds"
  )
})

test_that("ezdm_summary_stats() errors for non-positive RT values", {
  withr::local_seed(445)
  rt <- c(-0.1, 0, rgamma(98, shape = 5, rate = 10) + 0.3)
  response <- rbinom(100, 1, 0.8)

  expect_error(
    ezdm_summary_stats(rt, response, method = "simple"),
    "Non-positive RT values found"
  )
})

test_that("ezdm_summary_stats() handles too few trials", {
  withr::local_seed(123)
  rt <- rgamma(5, shape = 5, rate = 10) + 0.3
  response <- rbinom(5, 1, 0.8)

  result <- ezdm_summary_stats(rt, response, method = "simple", min_trials = 10)

  expect_true(is.na(result$mean_rt))
  expect_true(is.na(result$var_rt))
  expect_equal(result$n_trials, 5)
})

test_that("ezdm_summary_stats() simple method matches mean() and var()", {
  withr::local_seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- rbinom(100, 1, 0.8)

  result <- ezdm_summary_stats(rt, response, method = "simple")

  expect_equal(result$mean_rt, mean(rt), tolerance = 1e-10)
  expect_equal(result$var_rt, var(rt), tolerance = 1e-10)
  expect_equal(result$n_trials, 100L)
  expect_equal(result$n_upper, sum(response == 1))
})

test_that("ezdm_summary_stats() robust method uses median and IQR/MAD", {
  withr::local_seed(123)
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
  # Fixed data: 90 values tightly around 0.5 + 10 extreme outliers
  clean_rt <- seq(0.4, 0.6, length.out = 90)
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
  withr::local_seed(123)
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
  withr::local_seed(123)
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

# check_data() turns these NAs into placeholders for bmm(); the summaries
# themselves keep coding a boundary without summaries as NA (#430)
test_that("ezdm_summary_stats() 4par codes a boundary below min_trials as NA", {
  rt <- withr::with_seed(1, rgamma(100, shape = 5, rate = 10) + 0.3)
  response <- rep(c(1, 0), c(95, 5))

  result <- ezdm_summary_stats(rt, response, version = "4par",
                               method = "simple", min_trials = 10)

  expect_equal(result$n_trials - result$n_upper, 5)
  expect_true(is.na(result$mean_rt_lower))
  expect_true(is.na(result$var_rt_lower))
  expect_false(is.na(result$mean_rt_upper))
  expect_false(is.na(result$var_rt_upper))
})

test_that("ezdm_summary_stats() validates contaminant_bound", {
  withr::local_seed(445)
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
  withr::local_seed(123)
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
  withr::local_seed(445)
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
  withr::local_seed(445)
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
  withr::local_seed(445)
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
  # Fixed data: 50 legitimate RTs tightly clustered + 50 uniform contaminants
  # ensures the EM finds a large contaminant proportion that must be clipped
  rt <- c(rep(0.5, 25), rep(0.6, 25), seq(0.1, 5.0, length.out = 50))
  response <- rep(c(1, 0), 50)

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
  # Fixed data: tight cluster + spread-out values forces clipping at 0.1
  rt <- c(rep(0.5, 25), rep(0.6, 25), seq(0.1, 5.0, length.out = 50))
  response <- rep(c(1, 0), 50)

  expect_warning(
    ezdm_summary_stats(rt, response, max_contaminant = 0.1),
    "clipped to max_contaminant"
  )
})

test_that("ezdm_summary_stats() max_contaminant works with 4par version", {
  # Fixed data: tight clusters + spread-out contaminants per response group
  rt <- c(
    rep(0.5, 30), rep(0.6, 30), seq(0.1, 5.0, length.out = 20),
    rep(0.5, 30), rep(0.6, 30), seq(0.1, 5.0, length.out = 20)
  )
  response <- c(rep(1, 80), rep(0, 80))

  result <- suppressWarnings(
    ezdm_summary_stats(rt, response, version = "4par", max_contaminant = 0.2)
  )

  expect_true(is.na(result$contaminant_prop_upper) ||
    result$contaminant_prop_upper <= 0.2)
  expect_true(is.na(result$contaminant_prop_lower) ||
    result$contaminant_prop_lower <= 0.2)
})

test_that("ezdm_summary_stats() handles character responses", {
  withr::local_seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- sample(c("upper", "lower"), 100, replace = TRUE, prob = c(0.8, 0.2))

  result <- ezdm_summary_stats(rt, response, method = "simple")

  expect_equal(result$n_upper, sum(response == "upper"))
  expect_equal(result$n_trials, 100L)
})

test_that("ezdm_summary_stats() handles factor responses", {
  withr::local_seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- factor(sample(c("upper", "lower"), 100, replace = TRUE, prob = c(0.8, 0.2)))

  result <- ezdm_summary_stats(rt, response, method = "simple")
  expect_equal(result$n_upper, sum(response == "upper"))
})

test_that("ezdm_summary_stats() handles correct/error and logical responses", {
  withr::local_seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3

  response_ce <- sample(c("correct", "error"), 100, replace = TRUE, prob = c(0.8, 0.2))
  result_ce <- ezdm_summary_stats(rt, response_ce, method = "simple")
  expect_equal(result_ce$n_upper, sum(response_ce == "correct"))

  response_lgl <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.8, 0.2))
  result_lgl <- ezdm_summary_stats(rt, response_lgl, method = "simple")
  expect_equal(result_lgl$n_upper, sum(response_lgl))
})

test_that("ezdm_summary_stats() handles case-insensitive responses", {
  withr::local_seed(123)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- sample(c("UPPER", "Lower", "UpPeR", "LOWER"), 100, replace = TRUE)

  result <- ezdm_summary_stats(rt, response, method = "simple")
  expect_equal(result$n_upper, sum(tolower(response) == "upper"))
})

test_that("ezdm_summary_stats() errors on unrecognized response values", {
  withr::local_seed(445)
  rt <- rgamma(100, shape = 5, rate = 10) + 0.3
  response <- sample(c("fast", "slow"), 100, replace = TRUE)

  expect_error(
    ezdm_summary_stats(rt, response),
    "Unrecognized response values"
  )
})

test_that("ezdm_summary_stats() 4par version works with character responses", {
  withr::local_seed(123)
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
  withr::local_seed(123)
  rt <- c(rgamma(90, shape = 5, rate = 10) + 0.3, rep(NA, 10))
  response <- rbinom(100, 1, 0.8)

  result <- ezdm_summary_stats(rt, response, method = "simple")

  expect_equal(result$n_trials, 90L)
  expect_equal(result$mean_rt, mean(rt, na.rm = TRUE), tolerance = 1e-10)
})

test_that("ezdm_summary_stats() works with dplyr::reframe()", {
  skip_if_not_installed("dplyr")
  withr::local_seed(123)
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
# ezdm_summary_stats CONTAMINANT-FREE COUNT TESTS                         ####
############################################################################# !

# tight cognitive bulk + responses spread over the whole contaminant range, so
# the EM finds a contaminant proportion without any RNG in the test
.contaminated_rt <- function(n_clean, n_contam) {
  c(
    seq(0.50, 0.60, length.out = n_clean),
    seq(0.10, 4.00, length.out = n_contam)
  )
}

test_that("ezdm_summary_stats() removes the estimated contaminants from the counts", {
  rt <- .contaminated_rt(170, 30)
  response <- rep(c(1L, 1L, 1L, 0L), length.out = 200)

  result <- ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4))

  expect_gt(result$contaminant_prop, 0)
  expect_lt(result$n_trials, length(rt))
  expect_lt(result$n_upper, sum(response))
  expect_equal(result$n_trials, round(200 * (1 - result$contaminant_prop)))
  expect_equal(
    result$n_upper,
    round(sum(response) - 200 * result$contaminant_prop * 0.5)
  )
  expect_lte(result$n_upper, result$n_trials)
})

test_that("ezdm_summary_stats() keeps the raw counts without a contaminant proportion", {
  rt <- .contaminated_rt(170, 30)
  response <- rep(c(1L, 1L, 1L, 0L), length.out = 200)

  simple <- ezdm_summary_stats(rt, response, method = "simple")
  robust <- ezdm_summary_stats(rt, response, method = "robust")

  expect_true(is.na(simple$contaminant_prop))
  expect_equal(simple$n_trials, length(rt))
  expect_equal(simple$n_upper, sum(response))
  expect_equal(robust$n_trials, length(rt))
  expect_equal(robust$n_upper, sum(response))
})

test_that("ezdm_summary_stats() guess_rate moves n_upper but not n_trials", {
  rt <- .contaminated_rt(170, 30)
  response <- rep(c(1L, 1L, 1L, 0L), length.out = 200)

  none <- ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4),
                             guess_rate = 0)
  half <- ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4),
                             guess_rate = 0.5)

  expect_equal(none$n_upper, sum(response))
  expect_lt(half$n_upper, none$n_upper)
  expect_equal(none$n_trials, half$n_trials)
})

test_that("ezdm_summary_stats() validates guess_rate", {
  rt <- .contaminated_rt(90, 10)
  response <- rep(c(1L, 0L), length.out = 100)

  expect_error(ezdm_summary_stats(rt, response, guess_rate = -0.1),
    "guess_rate must be a single number between 0 and 1")
  expect_error(ezdm_summary_stats(rt, response, guess_rate = 1.5),
    "guess_rate must be a single number between 0 and 1")
  expect_error(ezdm_summary_stats(rt, response, guess_rate = c(0.2, 0.5)),
    "guess_rate must be a single number between 0 and 1")
})

test_that("ezdm_summary_stats() warns and degrades to an accuracy-neutral correction above the ceiling", {
  # this fixture's EM estimate is a free optimum rather than clipped, so the
  # counts are pinned against the returned proportion instead of as literals
  rt <- .contaminated_rt(170, 30)
  response <- c(rep(1L, 198), rep(0L, 2))

  expect_warning(
    result <- ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4)),
    "Observed accuracy \\(0.99\\) lies outside.*removed proportionally from both"
  )
  expect_gt(result$contaminant_prop, 0)
  expect_lt(result$n_upper, 198L)
  expect_lt(result$n_trials, 200L)
  expect_equal(result$n_upper, round(198 * (1 - result$contaminant_prop)))
  expect_equal(result$n_trials, round(200 * (1 - result$contaminant_prop)))
  expect_lte(result$n_upper, result$n_trials)
  expect_type(result$n_upper, "integer")
})

test_that("ezdm_summary_stats() warns and degrades to an accuracy-neutral correction below the floor", {
  rt <- .contaminated_rt(170, 30)
  response <- c(rep(1L, 2), rep(0L, 198))

  expect_warning(
    result <- ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4)),
    "Observed accuracy \\(0.01\\) lies outside"
  )
  expect_gt(result$contaminant_prop, 0)
  expect_lt(result$n_trials, 200L)
  expect_equal(result$n_upper, round(2 * (1 - result$contaminant_prop)))
  expect_equal(result$n_trials, round(200 * (1 - result$contaminant_prop)))
  expect_lte(result$n_upper, result$n_trials)
})

test_that("ezdm_summary_stats() corrects a cell that sits exactly on the ceiling", {
  # accuracy 0.90 with a contaminant proportion clipped to exactly 0.2 and a
  # guess rate of 0.5 is the boundary case: every non-contaminant response was
  # correct, which is consistent, so the counts are corrected rather than
  # refused. Only accuracy strictly above the ceiling is impossible.
  rt <- .contaminated_rt(152, 48)
  response <- c(rep(1L, 180), rep(0L, 20))

  result <- suppressWarnings(ezdm_summary_stats(
    rt, response, contaminant_bound = c(0.1, 4), max_contaminant = 0.2
  ))

  expect_equal(result$contaminant_prop, 0.2)
  expect_equal(result$n_upper, 160L)
  expect_equal(result$n_trials, 160L)
  expect_no_warning(.contaminant_free_counts(180L, 200L, 0.2, 0.5))
  expect_no_warning(.contaminant_free_counts(20L, 200L, 0.2, 0.5))
})

test_that("ezdm_summary_stats() rounding keeps n_upper <= n_trials at the ceiling", {
  # 197 trials, 177 correct, EM clipped to max_contaminant: the clean counts are
  # 157.3 and 157.6, a third of a count apart, so only a monotone rounding of
  # both keeps the constraint
  rt <- .contaminated_rt(150, 47)
  response <- c(rep(1L, 177), rep(0L, 20))

  result <- suppressWarnings(ezdm_summary_stats(
    rt, response, contaminant_bound = c(0.1, 4), max_contaminant = 0.2
  ))

  expect_equal(result$contaminant_prop, 0.2)
  expect_lte(result$n_upper, result$n_trials)
  expect_equal(result$n_upper, 157L)
  expect_equal(result$n_trials, 158L)
})

test_that("ezdm_summary_stats() rounds n_upper rather than truncating it", {
  # 184 trials, 165 correct, EM clipped to 0.2: clean_upper is 146.6, so only
  # rounding gives 147 -- truncation would drop a response
  rt <- .contaminated_rt(184 - 43, 43)
  response <- c(rep(1L, 165), rep(0L, 19))

  result <- suppressWarnings(ezdm_summary_stats(
    rt, response, contaminant_bound = c(0.1, 4), max_contaminant = 0.2
  ))

  expect_equal(result$contaminant_prop, 0.2)
  expect_equal(result$n_upper, 147L)
  expect_equal(result$n_trials, 147L)
  expect_type(result$n_upper, "integer")
})

test_that("ezdm_summary_stats() 4par consumes both boundary proportions", {
  rt <- c(
    .contaminated_rt(100, 10),
    seq(0.70, 0.80, length.out = 60), seq(0.10, 4.00, length.out = 30)
  )
  response <- c(rep(1L, 110), rep(0L, 90))

  result <- ezdm_summary_stats(rt, response, version = "4par",
                               contaminant_bound = c(0.1, 4))

  expect_gt(result$contaminant_prop_upper, 0)
  expect_gt(result$contaminant_prop_lower, 0)
  expect_false(isTRUE(all.equal(result$contaminant_prop_upper,
                                result$contaminant_prop_lower)))
  expect_equal(result$n_upper, round(110 * (1 - result$contaminant_prop_upper)))
  expect_equal(
    result$n_trials - result$n_upper,
    round(90 * (1 - result$contaminant_prop_lower))
  )
  expect_lte(result$n_upper, result$n_trials)
})

test_that("ezdm_summary_stats() 4par keeps n_upper <= n_trials with a near-empty boundary", {
  rt <- c(
    .contaminated_rt(160, 37),
    c(0.55, 0.58, 0.62)
  )
  response <- c(rep(1L, 197), rep(0L, 3))

  result <- ezdm_summary_stats(rt, response, version = "4par",
                               contaminant_bound = c(0.1, 4))

  expect_true(is.na(result$contaminant_prop_lower))
  expect_lte(result$n_upper, result$n_trials)
  expect_equal(result$n_trials - result$n_upper, 3L)
})

test_that("ezdm_summary_stats() 4par correction can take a boundary below two responses", {
  rt <- c(
    .contaminated_rt(150, 44),
    c(0.52, 0.55, 3.60, 3.90, 0.11, 3.99)
  )
  response <- c(rep(1L, 194), rep(0L, 6))

  result <- suppressWarnings(ezdm_summary_stats(
    rt, response, version = "4par", contaminant_bound = c(0.1, 4),
    min_trials = 2, max_contaminant = 0.9
  ))

  # the Stan likelihood gates the lower RT terms on this count, not on whether
  # the lower moments are NA, so the two can disagree after the correction
  expect_equal(result$n_trials - result$n_upper, 1L)
  expect_lt(result$n_trials - result$n_upper, 2)
  expect_false(is.na(result$mean_rt_lower))
  expect_lte(result$n_upper, result$n_trials)
})

test_that("ezdm_summary_stats() is deterministic across repeated calls", {
  rt <- .contaminated_rt(170, 30)
  response <- rep(c(1L, 1L, 1L, 0L), length.out = 200)

  three_par <- replicate(
    3, ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4)),
    simplify = FALSE
  )
  four_par <- replicate(
    3,
    ezdm_summary_stats(rt, response, version = "4par",
                       contaminant_bound = c(0.1, 4)),
    simplify = FALSE
  )

  expect_equal(three_par[[2]], three_par[[1]])
  expect_equal(three_par[[3]], three_par[[1]])
  expect_equal(four_par[[2]], four_par[[1]])
  expect_equal(four_par[[3]], four_par[[1]])
})

test_that("ezdm_summary_stats() counts pass check_data() unchanged", {
  rt <- .contaminated_rt(170, 30)
  response <- rep(c(1L, 1L, 1L, 0L), length.out = 200)
  dat <- ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4))

  expect_silent(check_data(
    ezdm(mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
         n_trials = "n_trials"),
    dat,
    bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  ))
})

############################################################################# !
# adjust_ezdm_accuracy TESTS                                              ####
############################################################################# !

test_that("adjust_ezdm_accuracy() signals its deprecation", {
  expect_warning(adjust_ezdm_accuracy(80, 100, 0.1), "deprecated")
  expect_warning(adjust_ezdm_accuracy(80, 100, 0.1), "ezdm_summary_stats")
})

test_that("adjust_ezdm_accuracy() returns 1-row data.frame", {
  withr::local_seed(42)
  result <- suppressWarnings(adjust_ezdm_accuracy(80, 100, 0.1))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c("n_upper_adj", "n_trials_adj") %in% names(result)))
})

test_that("adjust_ezdm_accuracy() returns integer values", {
  withr::local_seed(42)
  result <- suppressWarnings(adjust_ezdm_accuracy(80, 100, 0.1))

  expect_type(result$n_upper_adj, "integer")
  expect_type(result$n_trials_adj, "integer")
  expect_true(result$n_upper_adj >= 0)
  expect_true(result$n_upper_adj <= result$n_trials_adj)
  expect_true(result$n_trials_adj <= 100)
})

test_that("adjust_ezdm_accuracy() handles NA contaminant_prop", {
  result <- suppressWarnings(adjust_ezdm_accuracy(80, 100, NA))

  expect_equal(result$n_upper_adj, 80L)
  expect_equal(result$n_trials_adj, 100L)
})

test_that("adjust_ezdm_accuracy() handles zero contaminant_prop", {
  result <- suppressWarnings(adjust_ezdm_accuracy(80, 100, 0))

  expect_equal(result$n_upper_adj, 80L)
  expect_equal(result$n_trials_adj, 100L)
})

test_that("adjust_ezdm_accuracy() validates inputs", {
  suppressWarnings({
    expect_error(adjust_ezdm_accuracy("a", 100, 0.1), "n_upper must be numeric")
    expect_error(adjust_ezdm_accuracy(80, "b", 0.1), "n_trials must be numeric")
    expect_error(adjust_ezdm_accuracy(80, 100, 0.1, guess_rate = -0.1),
      "guess_rate must be between 0 and 1")
    expect_error(adjust_ezdm_accuracy(80, 100, 0.1, guess_rate = 1.5),
      "guess_rate must be between 0 and 1")
  })
})

test_that("adjust_ezdm_accuracy() rejects vector input and points to ezdm_summary_stats()", {
  suppressWarnings({
    expect_error(
      adjust_ezdm_accuracy(c(80, 90), c(100, 100), c(0.1, 0.2)),
      "takes one cell at a time.*ezdm_summary_stats\\(\\)"
    )
    expect_error(
      adjust_ezdm_accuracy(80, 100, c(0.1, 0.2)),
      "takes one cell at a time"
    )
  })
})

############################################################################# !
# FLAG_CONTAMINANT_RTS TESTS                                              ####
############################################################################# !

# Helper to create test RT vector with known contaminants (deterministic)
.create_test_rt_vec <- function(n = 100, add_contaminants = TRUE, prop_contam = 0.1) {
  if (add_contaminants) {
    n_legit <- floor((1 - prop_contam) * n)
    n_contam <- n - n_legit
    c(seq(0.3, 1.5, length.out = n_legit), seq(0.05, 0.15, length.out = n_contam))
  } else {
    seq(0.2, 1.5, length.out = n)
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
    expect_warning(
      flag_contaminant_rts(rt),
      "Some RT values > 10"
    ),
    "clipped to max_contaminant"
  )
})

test_that("flag_contaminant_rts() warns when lower bound > min RT", {
  rt <- c(0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

  expect_warning(
    flag_contaminant_rts(rt, contaminant_bound = c(0.3, "max")),
    "Lower contaminant bound.*greater than the minimum observed RT"
  )
})

test_that("flag_contaminant_rts() warns when upper bound < max RT", {
  rt <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0)

  expect_warning(
    flag_contaminant_rts(rt, contaminant_bound = c("min", 1.0)),
    "Upper contaminant bound.*less than the maximum observed RT"
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
  # Fixed data: tight legitimate cluster around 0.5-1.0, contaminants at 5-10
  # The separation is large enough that the mixture model must assign higher
  # contamination probability to the extreme values
  rt_legit <- seq(0.4, 1.2, length.out = 80)
  rt_contam <- seq(5, 10, length.out = 20)
  rt <- c(rt_legit, rt_contam)
  is_contam <- c(rep(FALSE, 80), rep(TRUE, 20))

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

# Without contaminants the EM approaches a contaminant proportion of zero
# slowly; this sample needs about 250 iterations
test_that("the EM converges by default where it needs more than 100 iterations", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(8)
  rt <- 0.3 + rnorm(60, 0.15, 0.03) + rexp(60, 10)
  response <- rep(c(1, 0), c(54, 6))

  expect_warning(
    capped <- ezdm_summary_stats(rt, response, maxit = 100),
    "did not converge within maxit = 100 iterations"
  )
  expect_true(is.na(capped$contaminant_prop))
  expect_warning(
    flag_contaminant_rts(rt, maxit = 100),
    "did not converge within maxit = 100 iterations"
  )

  stats <- expect_no_warning(ezdm_summary_stats(rt, response))
  expect_false(is.na(stats$contaminant_prop))
  probs <- expect_no_warning(flag_contaminant_rts(rt))
  expect_false(anyNA(probs))
  expect_gt(attr(probs, "diagnostics")$iterations, 100)
})

test_that("the EM warning suggests maxit only when the EM ran out of iterations", {
  expect_warning(
    flag_contaminant_rts(c(0.5, 0.6, 0.7)),
    "EM did not converge\\. Returning NA\\."
  )
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

# validate_fast_guesses(): three flagged trials below 0.5 s, two of them upper
# responses, update a Beta(1, 1) prior to Beta(3, 2); its density at 0.5 is
# 0.5^2 * 0.5 / B(3, 2) = 1.5, the Bayes factor for guessing
fast_guess_data <- function() {
  list(rt = c(0.1, 0.2, 0.3, 0.9, 1.0), flag = c(TRUE, TRUE, TRUE, TRUE, FALSE),
       response = c(1, 0, 1, 1, 0))
}

test_that("validate_fast_guesses() tests the flagged fast trials against guessing", {
  d <- fast_guess_data()
  out <- validate_fast_guesses(d$flag, d$rt, d$response, threshold_type = "absolute",
                               rt_threshold = 0.5)
  expect_equal(out$n_tested, 3)
  expect_equal(out$prop_upper, 2 / 3)
  expect_equal(c(out$posterior_alpha, out$posterior_beta), c(3, 2))
  expect_equal(out$bf_01, 1.5)
  expect_equal(out$bf_evidence, "anecdotal_for_guessing")
  expect_equal(out$mean_rt_tested, 0.2)
  # the highest-density interval: 95% of the mass, equal density at both ends
  expect_equal(pbeta(out$hdi_upper, 3, 2) - pbeta(out$hdi_lower, 3, 2), 0.95, tolerance = 1e-6)
  expect_equal(dbeta(out$hdi_lower, 3, 2), dbeta(out$hdi_upper, 3, 2), tolerance = 0.01)
  expect_lt(out$hdi_upper - out$hdi_lower, qbeta(0.975, 3, 2) - qbeta(0.025, 3, 2))
  expect_true(out$guess_in_hdi)

  # an unequal prior, so that each prior parameter has to reach its own side
  out <- validate_fast_guesses(d$flag, d$rt, d$response, threshold_type = "absolute",
                               rt_threshold = 0.5, prior_alpha = 2, prior_beta = 1)
  expect_equal(c(out$posterior_alpha, out$posterior_beta), c(4, 2))
  expect_equal(out$bf_01, 1.25)
})

test_that("validate_fast_guesses() reads a quantile threshold from the RTs", {
  d <- fast_guess_data()
  out <- validate_fast_guesses(d$flag, d$rt, d$response, rt_threshold = 0.5)
  expect_equal(out$rt_threshold, 0.3)
  expect_equal(out$n_tested, 2)
  expect_equal(out$threshold_type, "quantile")
})

test_that("validate_fast_guesses() falls back to the equal-tailed interval for a U-shaped posterior", {
  d <- fast_guess_data()
  out <- validate_fast_guesses(rep(FALSE, 5), d$rt, d$response, threshold_type = "absolute",
                               rt_threshold = 0.5, prior_alpha = 0.5, prior_beta = 0.5)
  expect_equal(out$n_tested, 0)
  expect_equal(c(out$hdi_lower, out$hdi_upper), qbeta(c(0.025, 0.975), 0.5, 0.5))
})

test_that("validate_fast_guesses() names the evidence on Jeffreys' scale", {
  expect_identical(
    vapply(c(20, 5, 2, 0.5, 0.2, 0.05), .categorize_bf, character(1)),
    c("strong_for_guessing", "moderate_for_guessing", "anecdotal_for_guessing",
      "anecdotal_against_guessing", "moderate_against_guessing", "strong_against_guessing")
  )
})

test_that("validate_fast_guesses() refuses priors, guess probabilities and thresholds out of range", {
  d <- fast_guess_data()
  call <- function(...) validate_fast_guesses(d$flag, d$rt, d$response, ...)
  expect_error(call(prior_alpha = 0), "Beta prior parameters must be positive")
  expect_error(call(guess_prob = 1), "guess_prob must be between 0 and 1")
  expect_error(call(rt_threshold = 1.5), "must be a quantile between 0 and 1")
  expect_error(call(threshold_type = "absolute", rt_threshold = 0), "must be positive")
})
