# Test Binary SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_binary model can be created with default arguments", {
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials"))
})

test_that("sdt_binary model can be created with all distribution options", {
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "normal"))
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "logistic"))
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_min"))
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_max"))
})

test_that("sdt_binary model has correct class structure", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_binary")
})

test_that("sdt_binary model parameters are correctly defined", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_binary model has correct default link functions", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_equal(model$links$dprime, "identity")
  expect_equal(model$links$criterion, "identity")
  expect_equal(model$links$sdratio, "log")
})

test_that("sdt_binary fixes sdratio on the link scale and inits on the natural scale", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  # constant() prior on the Intercept, so exp(0) = 1 = equal variance
  expect_equal(model$fixed_parameters$sdratio, 0)
  # create_initfun() applies the forward link, so log() of these must be finite
  expect_true(all(model$init_ranges$sdratio > 0))
})

test_that("sdt_binary model accepts custom links", {
  custom_links <- list(dprime = "log")
  model <- sdt_binary("n_old", "stimulus", "n_trials", links = custom_links)
  expect_equal(model$links$dprime, "log")
  expect_equal(model$links$criterion, "identity")
})

test_that("sdt_binary model stores distribution info correctly", {
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_min")
  expect_equal(model$other_vars$dist, "gumbel_min")

  model2 <- sdt_binary("n_old", "stimulus", "n_trials", dist = "logistic")
  expect_equal(model2$other_vars$dist, "logistic")
})

test_that("sdt_binary model has default priors", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_true("dprime" %in% names(model$default_priors))
  expect_true("criterion" %in% names(model$default_priors))
  expect_true("main" %in% names(model$default_priors$dprime))
  expect_true("effects" %in% names(model$default_priors$dprime))
})

test_that("sdt_binary requires n_trials argument", {
  expect_error(sdt_binary("n_old", "stimulus"))
})

test_that("sdt_binary model has init_ranges for estimated SDT parameters", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_true(all(c("dprime", "criterion", "sdratio") %in% names(model$init_ranges)))
  for (range in model$init_ranges) {
    expect_length(range, 2)
    expect_true(range[1] < range[2])
  }
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdt_binary check_data validates required columns", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  valid_data <- data.frame(
    n_old = c(30, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  result <- check_data(model, valid_data, formula)
  expect_true("dist_type" %in% colnames(result))

  invalid_data <- data.frame(n_old = c(30, 40))
  expect_error(
    check_data(model, invalid_data, formula),
    "missing in the data"
  )
})

test_that("sdt_binary check_data validates stimulus coding", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  invalid_data <- data.frame(
    n_old = c(30, 40),
    stimulus = c(1, 2),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "must be coded as 0"
  )

  factor_data <- data.frame(
    n_old = c(30, 40),
    stimulus = factor(c(0, 1)),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, factor_data, formula),
    "must be coded as 0"
  )
})

test_that("sdt_binary check_data validates response counts", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  invalid_data <- data.frame(
    n_old = c(-5, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "non-negative"
  )

  invalid_data2 <- data.frame(
    n_old = c(60, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data2, formula),
    "must not exceed"
  )
})

test_that("sdt_binary check_data validates n_trials", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  invalid_data <- data.frame(
    n_old = c(0, 0),
    stimulus = c(0L, 1L),
    n_trials = c(-10, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "positive"
  )
})

test_that("sdt_binary check_data adds dist_type column", {
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  valid_data <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )

  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_min")
  result <- check_data(model, valid_data, formula)
  expect_equal(unique(result$dist_type), 2L)

  model2 <- sdt_binary("n_old", "stimulus", "n_trials", dist = "logistic")
  result2 <- check_data(model2, valid_data, formula)
  expect_equal(unique(result2$dist_type), 4L)
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_binary generates counts matching the design", {
  stimulus <- rep(c(0L, 1L), 5)
  n_old <- rsdt_binary(10, 50, stimulus, dprime = 1.5, criterion = 0)
  expect_length(n_old, 10)
  expect_true(all(n_old >= 0))
  expect_true(all(n_old <= 50))
})

test_that("rsdt_binary recycles vectorized parameters per observation", {
  n_old <- rsdt_binary(4, c(50, 100, 50, 100), c(0L, 1L, 0L, 1L),
                       dprime = c(0.5, 0.5, 2.5, 2.5), criterion = 0)
  expect_length(n_old, 4)
  expect_true(all(n_old <= c(50, 100, 50, 100)))
})

test_that("rsdt_binary generates different data for different distributions", {
  set.seed(42)
  n_norm <- rsdt_binary(2, 1000, c(0L, 1L), dprime = 2, criterion = 0)
  set.seed(42)
  n_logis <- rsdt_binary(2, 1000, c(0L, 1L), dprime = 2, criterion = 0,
                         dist = "logistic")
  expect_false(identical(n_norm, n_logis))
})

test_that("rsdt_binary with zero dprime produces equal hit and FA rates on average", {
  set.seed(123)
  n_old <- rsdt_binary(2, 10000, c(0L, 1L), dprime = 0, criterion = 0)
  expect_equal(n_old[1] / 10000, 0.5, tolerance = 0.05)
  expect_equal(n_old[2] / 10000, 0.5, tolerance = 0.05)
})

test_that("rsdt_binary validates input", {
  expect_error(rsdt_binary(c(2, 3), 100, c(0L, 1L), dprime = 1, criterion = 0),
               "single positive integer")
  expect_error(rsdt_binary(2, 100, c(0L, 2L), dprime = 1, criterion = 0),
               "0.*1")
})

test_that("dsdt_binary returns valid binomial density", {
  d <- dsdt_binary(n_old = 50, n_trials = 100, stimulus = 1,
                   dprime = 0, criterion = 0)
  expect_true(d > 0)
  expect_true(d <= 1)
})

test_that("dsdt_binary log density matches log of density", {
  d <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                   dprime = 1.5, criterion = 0.2)
  ld <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                    dprime = 1.5, criterion = 0.2, log = TRUE)
  expect_equal(log(d), ld, tolerance = 1e-10)
})

test_that("dsdt_binary is vectorized over observations", {
  d <- dsdt_binary(n_old = c(30, 80), n_trials = c(100, 100),
                   stimulus = c(0, 1), dprime = 1.5, criterion = 0.2)
  expect_length(d, 2)
  expect_true(all(d > 0))
})

test_that("dsdt_binary works for all distributions", {
  dists <- c("normal", "logistic", "gumbel_min", "gumbel_max")
  for (di in dists) {
    d <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 1,
                     dprime = 1.5, criterion = 0, dist = di)
    expect_true(d > 0, info = paste("dist:", di))
  }
})

test_that("dsdt_binary validates input", {
  expect_error(dsdt_binary(n_old = -1, n_trials = 100, stimulus = 1,
                           dprime = 1, criterion = 0), "non-negative")
  expect_error(dsdt_binary(n_old = 150, n_trials = 100, stimulus = 1,
                           dprime = 1, criterion = 0), "must not exceed")
  expect_error(dsdt_binary(n_old = 50, n_trials = 100, stimulus = 2,
                           dprime = 1, criterion = 0), "0.*1")
  expect_error(dsdt_binary(n_old = 50, n_trials = 100, stimulus = 1,
                           dprime = 1, criterion = 0, sdratio = 0),
               "sdratio must be positive")
  expect_error(dsdt_binary(n_old = 50, n_trials = 100, stimulus = 1,
                           dprime = 1, criterion = 0, dist = "cauchy"),
               "should be one of")
})

test_that("dsdt_binary stays finite where the probability scale underflows", {
  # naive dbinom(y, n, pnorm(eta)) returns -Inf here and poisons loo()
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    ld <- dsdt_binary(n_old = 1, n_trials = 100, stimulus = 1,
                      dprime = -80, criterion = 0, dist = d, log = TRUE)
    expect_true(is.finite(ld), info = d)
  }
})

test_that("dsdt_binary handles the y = 0 and y = n_trials edges", {
  # 0 * -Inf must contribute nothing rather than becoming NaN
  expect_equal(dsdt_binary(n_old = 0, n_trials = 10, stimulus = 1,
                           dprime = -80, criterion = 0, log = TRUE), 0)
  expect_equal(dsdt_binary(n_old = 10, n_trials = 10, stimulus = 1,
                           dprime = 80, criterion = 0, log = TRUE), 0)
})

test_that("dist names the evidence distribution, not the likelihood's CDF", {
  # bmm must reproduce the textbook generative story -- draw evidence from
  # `dist`, respond "old" when it exceeds the criterion -- for asymmetric
  # distributions too. Using F(eta) instead of S(-eta) silently fits the mirror
  # distribution and only shows up for the extreme-value families.
  dprime <- 1.5
  criterion_textbook <- 0.3
  criterion <- criterion_textbook - dprime / 2

  quantile_fun <- list(
    normal     = qnorm,
    gumbel_min = function(p) log(-log(1 - p)),
    gumbel_max = function(p) -log(-log(p)),
    logistic   = qlogis
  )

  for (d in names(quantile_fun)) {
    # inverse-transform sampling, independent of .sdt_dists
    u <- (seq_len(20000) - 0.5) / 20000
    p_noise <- mean(quantile_fun[[d]](u) > criterion_textbook)
    p_signal <- mean(dprime + quantile_fun[[d]](u) > criterion_textbook)

    expect_equal(
      dsdt_binary(1, 1, 0L, dprime, criterion, dist = d), p_noise,
      tolerance = 1e-3, info = paste(d, "noise")
    )
    expect_equal(
      dsdt_binary(1, 1, 1L, dprime, criterion, dist = d), p_signal,
      tolerance = 1e-3, info = paste(d, "signal")
    )
  }
})

test_that("sdt_dprime and sdt_criterion invert the model's decision rule", {
  dprime <- 1.2
  criterion <- -0.2
  for (d in c("normal", "gumbel_min", "gumbel_max", "logistic")) {
    fa <- dsdt_binary(1, 1, 0L, dprime, criterion, dist = d)
    hit <- dsdt_binary(1, 1, 1L, dprime, criterion, dist = d)
    expect_equal(sdt_dprime(hit, fa, dist = d), dprime,
                 tolerance = 1e-8, info = d)
    expect_equal(sdt_criterion(hit, fa, dist = d), criterion,
                 tolerance = 1e-8, info = d)
  }
})

test_that("dsdt_binary matches dbinom in the well-conditioned range", {
  eta <- (1.5 / 2 - 0.2) / 1.3
  expect_equal(
    dsdt_binary(n_old = 70, n_trials = 100, stimulus = 1, dprime = 1.5,
                criterion = 0.2, sdratio = 1.3, log = TRUE),
    dbinom(70, 100, pnorm(eta), log = TRUE)
  )
})

test_that("dsdt_binary returns one value per posterior draw for a scalar count", {
  # times_nonzero() must recycle: ifelse() alone would collapse to length 1
  d <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 1,
                   dprime = rep(1.5, 20), criterion = 0.2, log = TRUE)
  expect_length(d, 20)
})

test_that("dsdt_binary with sdratio produces different densities than EV-SDT", {
  d_ev <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                      dprime = 1.5, criterion = 0, sdratio = 1)
  d_uv <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                      dprime = 1.5, criterion = 0, sdratio = 1.3)
  expect_false(d_ev == d_uv)
})

test_that("dsdt_binary with sdratio = 1 matches equal-variance for noise trials", {
  d1 <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 0,
                    dprime = 1.5, criterion = 0.2, sdratio = 1)
  d2 <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 0,
                    dprime = 1.5, criterion = 0.2, sdratio = 1.5)
  expect_equal(d1, d2)
})

test_that("dsdt_binary is vectorized over sdratio", {
  d <- dsdt_binary(n_old = c(30, 80), n_trials = c(100, 100),
                   stimulus = c(0, 1), dprime = 1.5, criterion = 0.2,
                   sdratio = c(1.0, 1.3))
  expect_length(d, 2)
  expect_true(all(d > 0))
})


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

# brmsprep stand-in: get_dpar() returns the supplied per-draw vectors already on
# the natural scale, as brms does once it has applied the inverse link.
# dprime/criterion are held constant so the only across-draw variation comes
# from sdratio -- isolating the unequal-variance signal scaling that the
# .sdt_eta() collapse bug got wrong.
uv_prep <- function(ndraws = 50L) {
  draws <- list(
    dprime    = rep(1.2, ndraws),
    criterion = rep(0.0, ndraws),
    sdratio   = seq(1.1, 1.9, length.out = ndraws)
  )
  prep <- list(data = list(
    vint1  = c(0L, 1L),      # stimulus: noise, signal
    vint2  = c(1L, 1L),      # dist_type: normal
    trials = c(100L, 100L),
    Y      = c(30L, 70L)
  ))
  list(prep = prep, draws = draws)
}

test_that("log_lik_sdt_binary scales the signal per draw under UV-SDT", {
  mk <- uv_prep()
  local_mocked_bindings(
    get_dpar = function(prep, dpar, i = NULL, ...) mk$draws[[dpar]],
    .package = "brms"
  )
  # signal row (stimulus = 1): p = Phi((dprime/2 - criterion) / sdratio)
  p_expected <- pnorm((1.2 / 2 - 0) / mk$draws$sdratio)
  expect_equal(log_lik_sdt_binary(2L, mk$prep),
               dbinom(70L, 100L, p_expected, log = TRUE))
  # the collapse bug divided every draw by sdratio[1], giving a constant p
  expect_gt(length(unique(round(log_lik_sdt_binary(2L, mk$prep), 8))), 1L)
  # noise row (stimulus = 0) is scale-free and must ignore sdratio entirely
  expect_equal(log_lik_sdt_binary(1L, mk$prep),
               rep(dbinom(30L, 100L, pnorm(-1.2 / 2), log = TRUE),
                   length(mk$draws$dprime)))
})

test_that("posterior_predict_sdt_binary feeds per-draw signal probabilities", {
  mk <- uv_prep()
  captured <- NULL
  local_mocked_bindings(
    get_dpar = function(prep, dpar, i = NULL, ...) mk$draws[[dpar]],
    .package = "brms"
  )
  local_mocked_bindings(
    rbinom = function(n, size, prob) { captured <<- prob; rep(0L, n) },
    .package = "stats"
  )
  out <- posterior_predict_sdt_binary(2L, mk$prep)
  expect_length(out, length(mk$draws$dprime))
  expect_equal(captured, pnorm((1.2 / 2) / mk$draws$sdratio))
})


############################################################################# !
# FORMULA AND CONFIGURE_MODEL TESTS                                      ####
############################################################################# !

test_that("sdt_binary model produces valid stancode with brms", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    id = c(1, 1, 2, 2)
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_binary_lpmf", code))
  expect_true(grepl("sdt_log_cumprob", code))
})

test_that("sdt_binary model produces valid stancode for all distributions", {
  dat <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_binary("n_old", "stimulus", "n_trials", dist = d)
    formula <- bmf(dprime ~ 1, criterion ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0, info = paste("Failed for dist:", d))
  }
})

test_that("sdt_binary model handles predictors in formula", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    condition = c("A", "A", "B", "B")
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ condition, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_binary model handles random effects in formula", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    id = c(1, 1, 2, 2)
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1 + (1 | id), criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_binary default_prior returns valid prior object", {
  dat <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  prior <- default_prior(formula, data = dat, model = model)
  expect_s3_class(prior, "brmsprior")
  prior_strs <- prior$prior
  expect_true(any(grepl("normal", prior_strs)))
})
