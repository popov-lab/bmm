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

test_that("lba refuses a link for a parameter it does not have", {
  expect_error(
    lba(rt = "rt", response = "response", n_choices = 2,
        links = list(typo = "log")),
    "Unrecognized link target\\(s\\): 'typo'"
  )
  expect_warning(
    model <- lba(rt = "rt", response = "response", n_choices = 2,
                 links = list(drifc = "log")),
    "'drifc' read as 'driftc'"
  )
  expect_equal(model$links$driftc, "log")
  expect_false("drifc" %in% names(model$links))
})

test_that("the custom version validates link targets once the formula names them", {
  model <- lba(rt = "rt", response = "response", version = "custom",
               links = list(fast = "log"))
  formula <- bmf(fast ~ 1, slow ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  model <- check_model(model, data = NULL, formula = formula)
  expect_equal(model$links$fast, "log")
  expect_equal(model$links$slow, "identity")
  # a second pass over the checked model must accept the category links
  expect_equal(check_model(model, data = NULL, formula = formula)$links, model$links)

  model <- lba(rt = "rt", response = "response", version = "custom",
               links = list(typo = "log"))
  expect_error(
    check_model(model, data = NULL, formula = formula),
    "Unrecognized link target\\(s\\): 'typo'"
  )
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

test_that("the vectorized family slices its vars only where brms really threads", {
  dat <- rlba(n = 20, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  model <- check_model(model, dat, formula)
  dat <- check_data(model, dat, formula)
  family_vars <- function() configure_model(model, dat, formula)$formula$family$vars

  expect_equal(family_vars(), c("vint1", "vint2", "vint3"))

  withr::with_options(list(brms.threads = brms::threading(2)), {
    expect_equal(family_vars(), paste0("vint", 1:3, "[start:end]"))
  })
  # brms compiles threading(force = TRUE) with threads but emits the serial
  # likelihood, where start/end do not exist
  withr::with_options(list(brms.threads = brms::threading(2, force = TRUE)), {
    expect_equal(family_vars(), c("vint1", "vint2", "vint3"))
  })
  withr::with_options(list(brms.threads = brms::threading(NULL)), {
    expect_equal(family_vars(), c("vint1", "vint2", "vint3"))
  })
})

test_that("the custom version slices every vint column when threading", {
  dat <- rlba(n = 99, drift = c(3, 2, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  dat$response <- c("fast", "mid", "slow")[dat$response]
  model <- lba(rt = "rt", response = "response", version = "custom",
               accumulators = c(fast = 1, mid = 2, slow = 3))
  formula <- bmf(fast ~ 1, mid ~ 1, slow ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  model <- check_model(model, dat, formula)
  dat <- check_data(model, dat, formula)

  withr::local_options(brms.threads = brms::threading(2))
  config <- configure_model(model, dat, formula)
  expect_equal(config$formula$family$vars, paste0("vint", 1:4, "[start:end]"))
})

test_that("stancode emits thread-safe slicing when threads is passed", {
  skip_on_cran()

  dat <- rlba(n = 100, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  serial_call <- paste0(
    "lba_normal_simple_lpdf(Y | mu, driftc, drifte, gap, sp, ndt, ",
    "s, vint1, vint2, vint3);"
  )

  code <- stancode(formula, dat, model = model, backend = "cmdstanr")
  expect_true(grepl(serial_call, code, fixed = TRUE))

  threaded <- stancode(formula, dat, model = model, backend = "cmdstanr",
                       threads = brms::threading(2))
  expect_true(grepl(
    paste0("lba_normal_simple_lpdf(Y[start:end] | mu, driftc, drifte, gap, ",
           "sp, ndt, s, vint1[start:end], vint2[start:end], ",
           "vint3[start:end]);"),
    threaded, fixed = TRUE
  ))
  expect_true(grepl("reduce_sum", threaded))

  # force = TRUE: brms keeps the serial likelihood, so sliced vars would
  # reference an undefined start/end and the model would not compile
  forced <- stancode(formula, dat, model = model, backend = "cmdstanr",
                     threads = brms::threading(2, force = TRUE))
  expect_true(grepl(serial_call, forced, fixed = TRUE))
  expect_false(grepl("start:end", forced, fixed = TRUE))

  # an explicit threads = NULL overrides a global option, as in brms
  withr::local_options(brms.threads = brms::threading(2))
  unthreaded <- stancode(formula, dat, model = model, backend = "cmdstanr",
                         threads = NULL)
  expect_true(grepl(serial_call, unthreaded, fixed = TRUE))
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


test_that("report_priors() omits the technical mu of the LBA family", {
  dat <- rlba(n = 60, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  fit <- bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)

  out <- report_priors(fit)
  expect_false("mu" %in% out$parameter)
  expect_true(all(c("driftc", "drifte", "gap", "sp", "ndt") %in% out$parameter))
})


# -----------------------------------------------------------------------------
# Posterior methods (brms contract: every draw of one observation in one call)
# -----------------------------------------------------------------------------

# A brmsprep with ndraws draws for one 2-choice observation whose starting
# point sits on both sides of the A ~ 0 branch of every kernel
lba_fake_prep <- function(distribution, driftc, drifte, sp, rt = 0.6,
                          response = 2L, gap = 0.5, ndt = 0.2, s = 1) {
  ndraws <- length(sp)
  structure(
    list(
      ndraws = ndraws, nobs = 1L,
      dpars = list(
        driftc = matrix(driftc, ndraws, 1), drifte = matrix(drifte, ndraws, 1),
        gap = matrix(gap, ndraws, 1), sp = matrix(sp, ndraws, 1),
        ndt = matrix(ndt, ndraws, 1), s = matrix(s, ndraws, 1)
      ),
      data = list(Y = rt, vint1 = response, vint2 = 1L, vint3 = 1L),
      family = list(name = paste0("lba_", distribution, "_simple"),
                    dpars = c("mu", "driftc", "drifte", "gap", "sp", "ndt", "s"))
    ),
    class = "brmsprep"
  )
}

test_that("log_lik evaluates all draws of an observation in one call", {
  drifts <- list(normal = c(3, 1.5), gamma = c(2, 3), frechet = c(2, 3),
                 lognormal = c(0.5, 0.3))
  sp <- c(1e-12, 0.3, 0.5)
  for (dist in names(drifts)) {
    d <- drifts[[dist]]
    prep <- lba_fake_prep(dist, driftc = d[1] + c(0, 0.2, -0.2),
                          drifte = d[2] + c(0, 0.1, -0.1), sp = sp)
    ll <- log_lik_lba_simple(1L, prep)
    expect_length(ll, 3)
    one_at_a_time <- vapply(1:3, function(k) {
      dlba(0.6, 2L, drift = c(d[1] + c(0, 0.2, -0.2)[k], d[2] + c(0, 0.1, -0.1)[k]),
           gap = 0.5, sp = sp[k], ndt = 0.2, distribution = dist, log = TRUE)
    }, numeric(1))
    expect_equal(ll, one_at_a_time, tolerance = 1e-10, info = dist)
  }
})

test_that("log_lik returns -Inf for draws whose ndt exceeds the response time", {
  prep <- lba_fake_prep("normal", driftc = 3, drifte = 1.5, sp = c(0.3, 0.3),
                        ndt = c(0.2, 0.7))
  ll <- log_lik_lba_simple(1L, prep)
  expect_true(is.finite(ll[1]))
  expect_equal(ll[2], -Inf)
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
  withr::local_seed(42)
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
  withr::local_seed(42)
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
  withr::local_seed(42)
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


# -----------------------------------------------------------------------------
# Density / rlba / pp_simulate joint-distribution consistency (Stage 3)
# -----------------------------------------------------------------------------

test_that("dlba integrates to 1 over rt and responses for every distribution", {
  gap <- 0.5; sp <- 0.3; ndt <- 0.2
  cases <- list(
    list(distribution = "normal", drift = c(3, 1.5)),
    list(distribution = "normal", drift = c(-0.5, 2)),
    list(distribution = "gamma", drift = c(2, 3)),
    list(distribution = "lognormal", drift = c(0.5, 0.3)),
    list(distribution = "lognormal", drift = c(-0.5, 0.3)),
    list(distribution = "frechet", drift = c(2, 3))
  )
  for (case in cases) {
    K <- length(case$drift)
    total <- sum(vapply(seq_len(K), function(resp) {
      stats::integrate(function(rt) {
        dlba(rt, resp, drift = case$drift, gap = gap, sp = sp, ndt = ndt,
             distribution = case$distribution)
      }, lower = ndt, upper = Inf, rel.tol = 1e-9, subdivisions = 500L)$value
    }, numeric(1)))
    expect_equal(total, 1, tolerance = 1e-6, label = case$distribution)
  }
})

test_that("rlba reproduces its own choice probability", {
  withr::local_seed(123)
  n <- 20000
  drift <- c(3, 1.5); gap <- 0.5; sp <- 0.5; ndt <- 0.2
  dat <- rlba(n, drift = drift, gap = gap, sp = sp, ndt = ndt)
  p_hat <- mean(dat$response == 2)
  p_true <- stats::integrate(
    function(rt) dlba(rt, 2, drift = drift, gap = gap, sp = sp, ndt = ndt),
    lower = ndt, upper = Inf, rel.tol = 1e-10
  )$value
  se <- sqrt(p_true * (1 - p_true) / n)
  expect_lt(abs(p_hat - p_true), 4 * se)
})

test_that("pp_simulate() draws rt and response jointly for the simple version", {
  withr::local_seed(11)
  n_draws <- 20000
  driftc <- 3; drifte <- 1.5; gap <- 0.5; sp <- 0.5; ndt <- 0.2
  prep <- lba_fake_prep("normal", driftc = driftc, drifte = drifte,
                        sp = rep(sp, n_draws), gap = gap, ndt = ndt)
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  sims <- pp_simulate(model, prep)

  drift <- c(driftc, drifte)
  p_resp2 <- mean(sims$response == 2)
  p_fast <- mean(sims$rt < 0.5)
  p_joint <- mean(sims$response == 2 & sims$rt < 0.5)

  int_resp2 <- stats::integrate(
    function(rt) dlba(rt, 2, drift = drift, gap = gap, sp = sp, ndt = ndt),
    ndt, Inf, rel.tol = 1e-10
  )$value
  int_fast <- sum(vapply(1:2, function(r) {
    stats::integrate(
      function(rt) dlba(rt, r, drift = drift, gap = gap, sp = sp, ndt = ndt),
      ndt, 0.5, rel.tol = 1e-10
    )$value
  }, numeric(1)))
  int_joint <- stats::integrate(
    function(rt) dlba(rt, 2, drift = drift, gap = gap, sp = sp, ndt = ndt),
    ndt, 0.5, rel.tol = 1e-10
  )$value

  se <- function(p) sqrt(p * (1 - p) / n_draws)
  expect_lt(abs(p_resp2 - int_resp2), 4 * se(int_resp2))
  expect_lt(abs(p_fast - int_fast), 4 * se(int_fast))
  expect_lt(abs(p_joint - int_joint), 4 * se(int_joint))
})

test_that("pp_simulate() draws rt and response jointly for the custom version", {
  withr::local_seed(12)
  n_draws <- 20000
  va <- 2; vb <- 1.2; gap <- 0.5; sp <- 0.5; ndt <- 0.2
  model <- check_model(
    lba(rt = "rt", response = "response", version = "custom",
        accumulators = c(a = 1, b = 2)),
    NULL, bmf(a ~ 1, b ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  )
  prep <- structure(
    list(
      ndraws = n_draws, nobs = 1L,
      dpars = list(
        a = matrix(va, n_draws, 1), b = matrix(vb, n_draws, 1),
        gap = matrix(gap, n_draws, 1), sp = matrix(sp, n_draws, 1),
        ndt = matrix(ndt, n_draws, 1), s = matrix(1, n_draws, 1)
      ),
      data = list(vint2 = 1L, vint3 = 2L)
    ),
    class = "brmsprep"
  )
  sims <- pp_simulate(model, prep)
  p_b <- mean(sims$response == 2)

  drift <- c(va, vb, vb)
  int_b <- sum(vapply(2:3, function(r) {
    stats::integrate(
      function(rt) dlba(rt, r, drift = drift, gap = gap, sp = sp, ndt = ndt),
      ndt, Inf, rel.tol = 1e-10
    )$value
  }, numeric(1)))
  se <- sqrt(int_b * (1 - int_b) / n_draws)
  expect_lt(abs(p_b - int_b), 4 * se)
})


# -----------------------------------------------------------------------------
# race_mean_decision_time() (deterministic posterior_epred integral)
# -----------------------------------------------------------------------------

test_that("race_mean_decision_time recovers the exact mean of an exponential race", {
  log_survivor <- function(t) -t
  mean_dt <- race_mean_decision_time(log_survivor, n_draws = 1)
  expect_equal(mean_dt, 1, tolerance = 1e-4)
})

test_that("race_mean_decision_time is deterministic across identical calls", {
  log_survivor <- function(t) -t
  a <- race_mean_decision_time(log_survivor, n_draws = 1)
  b <- race_mean_decision_time(log_survivor, n_draws = 1)
  expect_identical(a, b)
})

test_that("race_mean_decision_time is insensitive to t_max for a heavy-tailed normal race", {
  drift <- c(3, 1.5); b <- 1; A <- 0.5; s <- 1
  log_survivor <- function(t) {
    .lba_lsurv_single(t, drift[1], b, A, s, "normal") +
      .lba_lsurv_single(t, drift[2], b, A, s, "normal")
  }
  m1 <- race_mean_decision_time(log_survivor, n_draws = 1, t_max = 1e4)
  m2 <- race_mean_decision_time(log_survivor, n_draws = 1, t_max = 1e7)
  expect_lt(abs(m1 - m2), 1e-6)
})


# -----------------------------------------------------------------------------
# log_lik category-level convention (log(n[win]) for an error response)
# -----------------------------------------------------------------------------

test_that("log_lik_lba_simple adds log(n_win) for a category with several accumulators", {
  prep <- lba_fake_prep("normal", driftc = 3, drifte = 1.5, sp = 0.5,
                        rt = 0.6, response = 2L)
  prep$data$vint3 <- 2L

  ll <- log_lik_lba_simple(1L, prep)
  ref <- dlba(0.6, 2, drift = c(3, 1.5, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2,
             log = TRUE) + log(2)
  expect_equal(ll, ref, tolerance = 1e-10)
})

test_that("log_lik_lba_simple carries no log(n) term for a correct response", {
  prep <- lba_fake_prep("normal", driftc = 3, drifte = 1.5, sp = 0.5,
                        rt = 0.6, response = 1L)
  prep$data$vint3 <- 2L

  ll <- log_lik_lba_simple(1L, prep)
  ref <- dlba(0.6, 1, drift = c(3, 1.5, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2,
             log = TRUE)
  expect_equal(ll, ref, tolerance = 1e-10)
})


# -----------------------------------------------------------------------------
# check_data error messages (Stage 3)
# -----------------------------------------------------------------------------

test_that("check_data.lba_simple names the integer-coded convention on factor labels", {
  dat <- data.frame(rt = c(0.5, 0.6, 0.7),
                    response = factor(c("correct", "error", "correct")))
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "integer-coded")
})

test_that("check_data.lba_custom errors when a response category has zero accumulators on a trial", {
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), response = c("b", "a", "b"),
                    n_a = c(1L, 1L, 1L), n_b = c(0L, 1L, 1L))
  model <- lba(rt = "rt", response = "response", version = "custom",
               accumulators = c(a = "n_a", b = "n_b"))
  formula <- bmf(a ~ 1, b ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
  model <- check_model(model, dat, formula)
  expect_error(check_data(model, dat, formula), "no accumulator")
})


# -----------------------------------------------------------------------------
# check_model.lba warns when the formula frees the unidentified scale s
# -----------------------------------------------------------------------------

test_that("check_model.lba warns iff the formula estimates s", {
  model <- lba(rt = "rt", response = "response", n_choices = 2)
  free_s <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1, s ~ 1)
  fixed_s <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1, s = 1)
  no_s <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)

  expect_warning(check_model(model, NULL, free_s), "cannot identify")
  expect_no_warning(check_model(model, NULL, fixed_s))
  expect_no_warning(check_model(model, NULL, no_s))
})


# -----------------------------------------------------------------------------
# dlba() input validation and per-trial vector recycling (Stage 3)
# -----------------------------------------------------------------------------

test_that("dlba errors when response has neither length 1 nor the length of rt", {
  expect_error(
    dlba(c(0.5, 0.6, 0.7), c(1, 2), drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2),
    "length 1 or the length of rt"
  )
})

test_that("dlba errors when response indexes an accumulator that does not exist", {
  expect_error(
    dlba(0.5, 3, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2),
    "index the accumulators"
  )
})

test_that("dlba returns NA for an NA response time and a finite value elsewhere", {
  d <- dlba(c(0.5, NA), 1, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
  expect_true(is.finite(d[1]))
  expect_true(is.na(d[2]))
})

test_that("dlba recycles a per-trial gap vector like two separate scalar calls", {
  combined <- dlba(c(0.5, 0.6), 1, drift = c(3, 1.5), gap = c(0.5, 0.6),
                   sp = 0.5, ndt = 0.2)
  separate <- c(
    dlba(0.5, 1, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2),
    dlba(0.6, 1, drift = c(3, 1.5), gap = 0.6, sp = 0.5, ndt = 0.2)
  )
  expect_equal(combined, separate)
})

test_that("dlba accepts an unrestricted lognormal meanlog and rejects a non-positive gamma shape", {
  expect_gt(
    dlba(0.5, 1, drift = c(-0.5, 0.3), gap = 0.5, sp = 0.5, ndt = 0.2,
         distribution = "lognormal"),
    0
  )
  expect_error(
    dlba(0.5, 1, drift = c(-1, 2), gap = 0.5, sp = 0.5, ndt = 0.2,
         distribution = "gamma"),
    "positive"
  )
})
