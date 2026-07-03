# =============================================================================
# Tests for lba model (model-specific tests)
# Distribution function tests are in test-distributions.R
# =============================================================================

# -----------------------------------------------------------------------------
# Model construction tests
# -----------------------------------------------------------------------------

test_that("lba() creates simple model with correct structure", {
  model <- lba(rt = "rt", response = "response", n_choices = 2)

  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "lba")
  expect_s3_class(model, "lba_simple")
  expect_equal(model$resp_vars$rt, "rt")
  expect_equal(model$resp_vars$response, "response")
  expect_equal(model$other_vars$n_choices, 2L)
  expect_equal(model$version, "simple")
  expect_equal(model$distribution, "normal")
})

test_that("lba simple version has correct parameters", {
  model <- lba(rt = "rt", response = "response", n_choices = 4)

  expect_true(all(c("driftc", "drifte", "gap", "sp", "ndt", "s") %in%
                    names(model$parameters)))
  expect_equal(model$other_vars$n_choices, 4L)
})

test_that("lba gamma distribution uses consistent drift parameter names", {
  model <- lba(rt = "rt", response = "response", n_choices = 2,
               distribution = "gamma")

  expect_true(all(c("driftc", "drifte", "gap", "sp", "ndt", "s") %in%
                    names(model$parameters)))
  expect_equal(model$distribution, "gamma")
  expect_equal(model$links$driftc, "log")
})

test_that("lba has correct link functions for gap and sp", {
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  expect_equal(model$links$gap, "log")
  expect_equal(model$links$sp, "log")
})

test_that("lba accepts custom links", {
  model <- lba(rt = "rt", response = "response", n_choices = 2,
               links = list(driftc = "log"))
  expect_equal(model$links$driftc, "log")
  expect_equal(model$links$drifte, "identity")
})

test_that("lba errors on invalid n_choices", {
  expect_error(lba(rt = "rt", response = "response", n_choices = 1))
  expect_error(lba(rt = "rt", response = "response", n_choices = 2.5))
})

test_that("lba errors on missing required arguments", {
  expect_error(lba(response = "response", n_choices = 2))
  expect_error(lba(rt = "rt", n_choices = 2))
})

test_that("lba() creates custom model with correct structure", {
  model <- lba(rt = "rt", response = "resp", version = "custom")

  expect_s3_class(model, "lba_custom")
  expect_equal(model$version, "custom")
  expect_null(model$other_vars$n_choices)
})

test_that("lba custom accepts accumulators", {
  model <- lba(rt = "rt", response = "resp", version = "custom",
               accumulators = c(correct = 1, lure = 3))
  expect_equal(model$other_vars$accumulators, c(correct = 1, lure = 3))
})

test_that("lba deprecates n_alternatives and num_alternatives", {
  expect_warning(
    model <- lba(rt = "rt", response = "response", n_alternatives = 2),
    "n_alternatives.*deprecated.*n_choices"
  )
  expect_equal(model$other_vars$n_choices, 2L)

  expect_warning(
    model <- lba(rt = "rt", response = "resp", version = "custom",
                 num_alternatives = c(correct = 1, lure = 3)),
    "num_alternatives.*deprecated.*accumulators"
  )
  expect_equal(model$other_vars$accumulators, c(correct = 1, lure = 3))
})


# -----------------------------------------------------------------------------
# Data validation tests — simple version
# -----------------------------------------------------------------------------

test_that("check_data.lba errors on missing RT variable", {
  dat <- data.frame(response = 1:10)
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  expect_error(check_data(model, dat, bmf(driftc ~ 1)))
})

test_that("check_data.lba errors on NA RT values", {
  dat <- data.frame(rt = c(0.5, NA, 0.6), response = c(1L, 2L, 1L))
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  expect_error(check_data(model, dat, bmf(driftc ~ 1)))
})

test_that("check_data.lba errors on negative RT", {
  dat <- data.frame(rt = c(0.5, -0.1, 0.6), response = c(1L, 2L, 1L))
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  expect_error(check_data(model, dat, bmf(driftc ~ 1)))
})

test_that("check_data.lba warns on RT > 10s", {
  dat <- data.frame(rt = c(0.5, 15, 0.6), response = c(1L, 2L, 1L))
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  expect_warning(check_data(model, dat, bmf(driftc ~ 1)))
})

test_that("check_data.lba_simple errors on out-of-range response", {
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), response = c(1L, 3L, 1L))
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  expect_error(check_data(model, dat, bmf(driftc ~ 1)))
})

test_that("check_data.lba_simple converts factor responses", {
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), response = factor(c(1, 2, 1)))
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  dat2 <- check_data(model, dat, bmf(driftc ~ 1))
  expect_true(is.integer(dat2$response))
})

test_that("check_data.lba_simple creates category columns", {
  dat <- data.frame(rt = c(0.5, 0.6, 0.7, 0.8),
                    response = c(1L, 2L, 3L, 1L))
  model <- lba(rt = "rt", response = "response", n_choices = 3)
  dat2 <- check_data(model, dat, bmf(driftc ~ 1))

  expect_true(".lba_cat" %in% names(dat2))
  expect_true(".lba_n1" %in% names(dat2))
  expect_true(".lba_n2" %in% names(dat2))
  expect_equal(dat2$.lba_cat[1], 1L)
  expect_equal(dat2$.lba_cat[2], 2L)
  expect_equal(dat2$.lba_n1[1], 1L)
  expect_equal(dat2$.lba_n2[1], 2L)
})


# -----------------------------------------------------------------------------
# Data validation tests — custom version
# -----------------------------------------------------------------------------

test_that("check_data.lba_custom maps character responses", {
  dat <- data.frame(rt = c(0.5, 0.6, 0.7),
                    response = c("correct", "wrong", "correct"))
  model <- lba(rt = "rt", response = "response", version = "custom")
  model$other_vars$resp_cats <- c("correct", "wrong")
  dat2 <- check_data(model, dat, bmf(correct ~ 1))
  expect_equal(dat2$.lba_cat[1], 1L)
  expect_equal(dat2$.lba_cat[2], 2L)
})

test_that("check_data.lba_custom errors on missing formula category", {
  dat <- data.frame(rt = c(0.5, 0.6), response = c("correct", "unknown"))
  model <- lba(rt = "rt", response = "response", version = "custom")
  model$other_vars$resp_cats <- c("correct", "wrong")
  expect_error(check_data(model, dat, bmf(correct ~ 1)))
})


# -----------------------------------------------------------------------------
# check_model tests — custom version
# -----------------------------------------------------------------------------

test_that("check_model.lba_custom discovers category params from formula", {
  model <- lba(rt = "rt", response = "response", version = "custom")
  formula <- bmf(fast ~ 1, slow ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  model <- check_model(model, data = NULL, formula = formula)

  expect_equal(model$other_vars$resp_cats, c("fast", "slow"))
  expect_true("fast" %in% names(model$parameters))
  expect_true("slow" %in% names(model$parameters))
})

test_that("check_model.lba_custom errors on Stan reserved words", {
  model <- lba(rt = "rt", response = "response", version = "custom")
  formula <- bmf(void ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_error(check_model(model, data = NULL, formula = formula))
})

test_that("check_model.lba_custom errors on category names ending in numbers", {
  model <- lba(rt = "rt", response = "response", version = "custom")
  formula <- bmf(correct ~ 1, error1 ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)

  expect_error(
    check_model(model, data = NULL, formula = formula),
    "cannot end in a number"
  )
})

test_that("check_model.lba_custom errors on category names containing underscores", {
  model <- lba(rt = "rt", response = "response", version = "custom")
  formula <- bmf(correct ~ 1, error_a ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)

  expect_error(
    check_model(model, data = NULL, formula = formula),
    "cannot contain underscores"
  )
})


# -----------------------------------------------------------------------------
# Formula conversion tests
# -----------------------------------------------------------------------------

test_that("bmf2bf.lba_simple creates correct formula", {
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  bf <- bmf2bf(model, bmf(driftc ~ 1))
  expect_true(grepl("vint\\(.lba_cat, .lba_n1, .lba_n2\\)", deparse(bf$formula)))
})

test_that("bmf2bf.lba_custom creates correct formula", {
  model <- lba(rt = "rt", response = "response", version = "custom")
  model$other_vars$resp_cats <- c("a", "b", "c")
  bf <- bmf2bf(model, bmf(a ~ 1))
  expect_true(grepl("vint\\(.lba_cat, .lba_n1, .lba_n2, .lba_n3\\)",
                    deparse(bf$formula)))
})


# -----------------------------------------------------------------------------
# Stan code generation tests
# -----------------------------------------------------------------------------

test_that("Stan code for normal LBA contains expected elements", {
  code <- bmm:::.lba_stan_code("lba_normal_simple", c("driftc", "drifte"), "normal")
  expect_true(grepl("lba_normal_single_lpdf", code))
  expect_true(grepl("lba_normal_single_lccdf", code))
  expect_true(grepl("std_normal_lcdf", code))
  expect_true(grepl("lba_log_clip", code))
  # single-pass race replaces the lba_race_loglik helper
  expect_false(grepl("lba_race_loglik", code))
  expect_true(grepl("reps", code))
  expect_false(grepl("Phi_approx", code))
})

test_that("Stan code for gamma LBA contains expected elements", {
  code <- bmm:::.lba_stan_code("lba_gamma_simple", c("driftc", "drifte"), "gamma")
  dat <- rlba(10, drift = c(2, 3), gap = 0.5, sp = 0.5, ndt = 0.2,
              distribution = "gamma")
  config <- configure_model(
    lba(rt = "rt", response = "response", n_choices = 2,
        distribution = "gamma"),
    dat,
    bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  )
  expect_true(grepl("lba_gamma_single_lpdf", code))
  expect_true(grepl("gamma_lcdf", code))
  expect_false(config$formula$family$loop)
})

test_that("Stan code for lognormal LBA contains expected elements", {
  code <- bmm:::.lba_stan_code("lba_lognormal_simple", c("driftc", "drifte"), "lognormal")
  expect_true(grepl("lba_lognormal_single_lpdf", code))
  expect_true(grepl("std_normal_lcdf", code))
  expect_true(grepl("log_diff_exp", code))
  expect_false(grepl("Phi_approx", code))
})

test_that("Stan code for frechet LBA contains expected elements", {
  code <- bmm:::.lba_stan_code("lba_frechet_simple", c("driftc", "drifte"), "frechet")
  expect_true(grepl("lba_frechet_single_lpdf", code))
  expect_true(grepl("array\\[16\\] real nodes", code))
  expect_true(grepl("log_sum_exp", code))
})

test_that("LBA generated Stan code uses vectorized custom likelihoods", {
  dat <- rlba(n = 20, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  simple_model <- lba(rt = "rt", response = "response", n_choices = 2)
  simple_formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  simple_code <- stancode(simple_formula, data = dat, model = simple_model,
                          backend = "cmdstanr")

  expect_true(grepl(
    "target \\+= lba_normal_simple_lpdf\\(Y \\| mu, driftc, drifte, gap, sp, ndt, s, vint1, vint2, vint3\\);",
    simple_code
  ))
  expect_false(grepl(
    "for \\(n in 1:N\\) \\{\\s+target \\+= lba_normal_simple_lpdf",
    simple_code,
    perl = TRUE
  ))

  dat$response <- ifelse(dat$response == 1, "correct", "wrong")
  custom_model <- lba(rt = "rt", response = "response", version = "custom")
  custom_formula <- bmf(correct ~ 1, wrong ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  custom_code <- stancode(custom_formula, data = dat, model = custom_model,
                          backend = "cmdstanr")

  expect_true(grepl("target \\+= lba_normal_custom_lpdf\\(Y \\|", custom_code))
  expect_false(grepl(
    "for \\(n in 1:N\\) \\{\\s+target \\+= lba_normal_custom_lpdf",
    custom_code,
    perl = TRUE
  ))
})


# -----------------------------------------------------------------------------
# Mock integration tests
# -----------------------------------------------------------------------------

test_that("lba normal simple version runs with mock backend (2-choice)", {
  dat <- rlba(n = 100, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("lba normal simple version runs with mock backend (4-choice)", {
  dat <- rlba(n = 100, drift = c(3, 1.5, 1, 0.8), gap = 0.5, sp = 0.5, ndt = 0.2)
  model <- lba(rt = "rt", response = "response", n_choices = 4)
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("lba gamma simple version runs with mock backend", {
  dat <- rlba(n = 100, drift = c(2, 3), gap = 0.5, sp = 0.5, ndt = 0.2,
              distribution = "gamma")
  model <- lba(rt = "rt", response = "response", n_choices = 2,
               distribution = "gamma")
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("lba lognormal simple version runs with mock backend", {
  dat <- rlba(n = 100, drift = c(0.5, 0.3), gap = 0.5, sp = 0.5, ndt = 0.2,
              distribution = "lognormal")
  model <- lba(rt = "rt", response = "response", n_choices = 2,
               distribution = "lognormal")
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("lba frechet simple version runs with mock backend", {
  dat <- rlba(n = 100, drift = c(2, 3), gap = 0.5, sp = 0.5, ndt = 0.2,
              distribution = "frechet")
  model <- lba(rt = "rt", response = "response", n_choices = 2,
               distribution = "frechet")
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("lba custom version runs with mock backend", {
  dat <- rlba(n = 100, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  dat$response <- ifelse(dat$response == 1, "correct", "wrong")
  model <- lba(rt = "rt", response = "response", version = "custom")
  formula <- bmf(correct ~ 1, wrong ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("lba simple with predictor runs with mock backend", {
  dat <- rlba(n = 200, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  dat$cond <- rep(c("a", "b"), each = 100)
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(driftc ~ cond, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  expect_no_error(
    suppressWarnings(
      bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  )
})


# -----------------------------------------------------------------------------
# Distribution function tests (dlba/rlba/plba use gap+sp interface)
# -----------------------------------------------------------------------------

test_that("dlba returns positive densities for valid inputs", {
  d <- dlba(c(0.5, 0.6), c(1, 2), drift = c(3, 1.5),
            gap = 0.5, sp = 0.5, ndt = 0.2)
  expect_true(all(d > 0))
})

test_that("rlba returns valid data.frame", {
  dat <- rlba(100, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  expect_s3_class(dat, "data.frame")
  expect_true(all(c("rt", "response") %in% names(dat)))
  expect_true(all(dat$rt > 0.2))
  expect_true(all(dat$response %in% 1:2))
})

test_that("plba is monotonically increasing", {
  rts <- seq(0.3, 1, by = 0.1)
  p <- plba(rts, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  expect_true(all(diff(p) >= 0))
})

test_that("dlba matches rtdists for normal distribution", {
  skip_if_not_installed("rtdists")
  set.seed(42)
  # rtdists: A = 0.5, b = 1 -> sp = A = 0.5, gap = b - A = 0.5
  ref <- rtdists::rLBA(5, A = 0.5, b = 1, t0 = 0.3,
                       mean_v = c(3, 1.5), sd_v = c(1, 1),
                       distribution = "norm")
  bmm_d <- dlba(ref$rt, ref$response, drift = c(3, 1.5),
                gap = 0.5, sp = 0.5, ndt = 0.3, distribution = "normal")
  rtd_d <- rtdists::dLBA(ref$rt, ref$response, A = 0.5, b = 1, t0 = 0.3,
                         mean_v = c(3, 1.5), sd_v = c(1, 1),
                         distribution = "norm", silent = TRUE)
  expect_equal(bmm_d, rtd_d, tolerance = 1e-10)
})

test_that("dlba matches rtdists for gamma distribution", {
  skip_if_not_installed("rtdists")
  set.seed(42)
  ref <- rtdists::rLBA(5, A = 0.5, b = 1, t0 = 0.3,
                       shape_v = c(2, 3), rate_v = c(1, 1),
                       distribution = "gamma")
  bmm_d <- dlba(ref$rt, ref$response, drift = c(2, 3),
                gap = 0.5, sp = 0.5, ndt = 0.3, s = 1,
                distribution = "gamma")
  rtd_d <- rtdists::dLBA(ref$rt, ref$response, A = 0.5, b = 1, t0 = 0.3,
                         shape_v = c(2, 3), rate_v = c(1, 1),
                         distribution = "gamma", silent = TRUE)
  expect_equal(bmm_d, rtd_d, tolerance = 1e-10)
})

test_that("dlba matches rtdists for lognormal distribution", {
  skip_if_not_installed("rtdists")
  set.seed(42)
  ref <- rtdists::rLBA(5, A = 0.5, b = 1, t0 = 0.3,
                       meanlog_v = c(0.5, 0.3), sdlog_v = c(1, 1),
                       distribution = "lnorm")
  bmm_d <- dlba(ref$rt, ref$response, drift = c(0.5, 0.3),
                gap = 0.5, sp = 0.5, ndt = 0.3, s = 1,
                distribution = "lognormal")
  rtd_d <- rtdists::dLBA(ref$rt, ref$response, A = 0.5, b = 1, t0 = 0.3,
                         meanlog_v = c(0.5, 0.3), sdlog_v = c(1, 1),
                         distribution = "lnorm", silent = TRUE)
  expect_equal(bmm_d, rtd_d, tolerance = 1e-10)
})

test_that("validate_lba_parameters catches invalid inputs", {
  expect_error(dlba(0.5, 1, drift = c(3, 1.5), gap = -1, sp = 0.3,
                    ndt = 0.2))
  expect_error(dlba(0.5, 1, drift = c(3, 1.5), gap = 0.5, sp = -0.1,
                    ndt = 0.2))
  expect_error(dlba(0.5, 1, drift = c(3, 1.5), gap = 0.5, sp = 0.3,
                    ndt = -1))
  expect_error(dlba(0.5, 1, drift = c(3, 1.5), gap = 0.5, sp = 0.3,
                    ndt = 0.2, s = 0))
})
