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
  expect_true("d" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_binary model has correct default link functions", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_equal(model$links$d, "identity")
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
  custom_links <- list(d = "log")
  model <- sdt_binary("n_old", "stimulus", "n_trials", links = custom_links)
  expect_equal(model$links$d, "log")
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
  expect_true("d" %in% names(model$default_priors))
  expect_true("criterion" %in% names(model$default_priors))
  expect_true("main" %in% names(model$default_priors$d))
  expect_true("effects" %in% names(model$default_priors$d))
})

test_that("sdt_binary requires n_trials argument", {
  expect_error(sdt_binary("n_old", "stimulus"))
})

test_that("sdt_binary model has init_ranges for estimated SDT parameters", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_true(all(c("d", "criterion", "sdratio") %in% names(model$init_ranges)))
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
  formula <- bmf(d ~ 1, criterion ~ 1)

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
  formula <- bmf(d ~ 1, criterion ~ 1)

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
  formula <- bmf(d ~ 1, criterion ~ 1)

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
  formula <- bmf(d ~ 1, criterion ~ 1)

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
  formula <- bmf(d ~ 1, criterion ~ 1)
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
  n_old <- rsdt_binary(10, 50, stimulus, d = 1.5, criterion = 0)
  expect_length(n_old, 10)
  expect_true(all(n_old >= 0))
  expect_true(all(n_old <= 50))
})

test_that("rsdt_binary recycles vectorized parameters per observation", {
  n_old <- rsdt_binary(4, c(50, 100, 50, 100), c(0L, 1L, 0L, 1L),
                       d = c(0.5, 0.5, 2.5, 2.5), criterion = 0)
  expect_length(n_old, 4)
  expect_true(all(n_old <= c(50, 100, 50, 100)))
})

test_that("rsdt_binary generates different data for different distributions", {
  set.seed(42)
  n_norm <- rsdt_binary(2, 1000, c(0L, 1L), d = 2, criterion = 0)
  set.seed(42)
  n_logis <- rsdt_binary(2, 1000, c(0L, 1L), d = 2, criterion = 0,
                         dist = "logistic")
  expect_false(identical(n_norm, n_logis))
})

test_that("rsdt_binary with zero d produces equal hit and FA rates on average", {
  set.seed(123)
  n_old <- rsdt_binary(2, 10000, c(0L, 1L), d = 0, criterion = 0)
  expect_equal(n_old[1] / 10000, 0.5, tolerance = 0.05)
  expect_equal(n_old[2] / 10000, 0.5, tolerance = 0.05)
})

test_that("rsdt_binary validates input", {
  expect_error(rsdt_binary(c(2, 3), 100, c(0L, 1L), d = 1, criterion = 0),
               "single positive integer")
  expect_error(rsdt_binary(2, 100, c(0L, 2L), d = 1, criterion = 0),
               "0.*1")
})

test_that("dsdt_binary returns valid binomial density", {
  dens <- dsdt_binary(n_old = 50, n_trials = 100, stimulus = 1,
                   d = 0, criterion = 0)
  expect_true(dens > 0)
  expect_true(dens <= 1)
})

test_that("dsdt_binary log density matches log of density", {
  dens <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                   d = 1.5, criterion = 0.2)
  ld <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                    d = 1.5, criterion = 0.2, log = TRUE)
  expect_equal(log(dens), ld, tolerance = 1e-10)
})

test_that("dsdt_binary is vectorized over observations", {
  dens <- dsdt_binary(n_old = c(30, 80), n_trials = c(100, 100),
                   stimulus = c(0, 1), d = 1.5, criterion = 0.2)
  expect_length(dens, 2)
  expect_true(all(dens > 0))
})

test_that("dsdt_binary works for all distributions", {
  dists <- c("normal", "logistic", "gumbel_min", "gumbel_max")
  for (di in dists) {
    dens <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 1,
                     d = 1.5, criterion = 0, dist = di)
    expect_true(dens > 0, info = paste("dist:", di))
  }
})

test_that("dsdt_binary validates input", {
  expect_error(dsdt_binary(n_old = -1, n_trials = 100, stimulus = 1,
                           d = 1, criterion = 0), "non-negative")
  expect_error(dsdt_binary(n_old = 150, n_trials = 100, stimulus = 1,
                           d = 1, criterion = 0), "must not exceed")
  expect_error(dsdt_binary(n_old = 50, n_trials = 100, stimulus = 2,
                           d = 1, criterion = 0), "0.*1")
  expect_error(dsdt_binary(n_old = 50, n_trials = 100, stimulus = 1,
                           d = 1, criterion = 0, sdratio = 0),
               "sdratio must be positive")
  expect_error(dsdt_binary(n_old = 50, n_trials = 100, stimulus = 1,
                           d = 1, criterion = 0, dist = "cauchy"),
               "should be one of")
})

test_that("dsdt_binary stays finite where the probability scale underflows", {
  # naive dbinom(y, n, pnorm(eta)) returns -Inf here and poisons loo()
  for (dist_name in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    ld <- dsdt_binary(n_old = 1, n_trials = 100, stimulus = 1,
                      d = -80, criterion = 0, dist = dist_name,
                      log = TRUE)
    expect_true(is.finite(ld), info = dist_name)
  }
})

test_that("dsdt_binary handles the y = 0 and y = n_trials edges", {
  # 0 * -Inf must contribute nothing rather than becoming NaN
  expect_equal(dsdt_binary(n_old = 0, n_trials = 10, stimulus = 1,
                           d = -80, criterion = 0, log = TRUE), 0)
  expect_equal(dsdt_binary(n_old = 10, n_trials = 10, stimulus = 1,
                           d = 80, criterion = 0, log = TRUE), 0)
})

test_that("dist names the evidence distribution, not the likelihood's CDF", {
  # bmm must reproduce the textbook generative story -- draw evidence from
  # `dist`, respond "old" when it exceeds the criterion -- for asymmetric
  # distributions too. Using F(eta) instead of S(-eta) silently fits the mirror
  # distribution and only shows up for the extreme-value families.
  d <- 1.5
  criterion_textbook <- 0.3
  criterion <- criterion_textbook - d / 2

  quantile_fun <- list(
    normal     = qnorm,
    gumbel_min = function(p) log(-log(1 - p)),
    gumbel_max = function(p) -log(-log(p)),
    logistic   = qlogis
  )

  for (dist_name in names(quantile_fun)) {
    # inverse-transform sampling, independent of .sdt_dists
    u <- (seq_len(20000) - 0.5) / 20000
    p_noise <- mean(quantile_fun[[dist_name]](u) > criterion_textbook)
    p_signal <- mean(d + quantile_fun[[dist_name]](u) > criterion_textbook)

    expect_equal(
      dsdt_binary(1, 1, 0L, d, criterion, dist = dist_name), p_noise,
      tolerance = 1e-3, info = paste(dist_name, "noise")
    )
    expect_equal(
      dsdt_binary(1, 1, 1L, d, criterion, dist = dist_name), p_signal,
      tolerance = 1e-3, info = paste(dist_name, "signal")
    )
  }
})

test_that("sdt_d and sdt_criterion invert the model's decision rule", {
  d <- 1.2
  criterion <- -0.2
  for (dist_name in c("normal", "gumbel_min", "gumbel_max", "logistic")) {
    fa <- dsdt_binary(1, 1, 0L, d, criterion, dist = dist_name)
    hit <- dsdt_binary(1, 1, 1L, d, criterion, dist = dist_name)
    expect_equal(sdt_d(hit, fa, dist = dist_name), d,
                 tolerance = 1e-8, info = dist_name)
    expect_equal(sdt_criterion(hit, fa, dist = dist_name), criterion,
                 tolerance = 1e-8, info = dist_name)
  }
})

test_that("dsdt_binary matches dbinom in the well-conditioned range", {
  eta <- (1.5 * sqrt((1 + 1.3^2) / 2) / 2 - 0.2) / 1.3
  expect_equal(
    dsdt_binary(n_old = 70, n_trials = 100, stimulus = 1, d = 1.5,
                criterion = 0.2, sdratio = 1.3, log = TRUE),
    dbinom(70, 100, pnorm(eta), log = TRUE)
  )
})

test_that("dsdt_binary returns one value per posterior draw for a scalar count", {
  # times_nonzero() must recycle: ifelse() alone would collapse to length 1
  dens <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 1,
                      d = rep(1.5, 20), criterion = 0.2, log = TRUE)
  expect_length(dens, 20)
})

test_that("dsdt_binary with sdratio produces different densities than EV-SDT", {
  d_ev <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                      d = 1.5, criterion = 0, sdratio = 1)
  d_uv <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                      d = 1.5, criterion = 0, sdratio = 1.3)
  expect_false(d_ev == d_uv)
})

test_that("sdratio never divides the noise trials, it only widens the separation", {
  # Only the signal distribution carries the scale division. On a noise trial
  # sdratio therefore acts solely through d_a -> d_N = d * s, so the same
  # density must be reachable by an equal-variance call with the widened d.
  s <- sqrt((1 + 1.5^2) / 2)
  uv <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 0,
                    d = 1.5, criterion = 0.2, sdratio = 1.5)
  ev <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 0,
                    d = 1.5 * s, criterion = 0.2, sdratio = 1)
  expect_equal(uv, ev)
})

test_that("dsdt_binary is vectorized over sdratio", {
  dens <- dsdt_binary(n_old = c(30, 80), n_trials = c(100, 100),
                      stimulus = c(0, 1), d = 1.5, criterion = 0.2,
                      sdratio = c(1.0, 1.3))
  expect_length(dens, 2)
  expect_true(all(dens > 0))
})


############################################################################# !
# SENSITIVITY PARAMETERIZATION (d = d_a)                                 ####
############################################################################# !

test_that("d is a no-op relative to d' for every equal-variance call", {
  # sdratio = 1 gives s = 1, so nothing in the reparameterization may move
  grid <- expand.grid(d = c(0.4, 1.5, 3), criterion = c(-0.7, 0, 0.5),
                      stimulus = c(0L, 1L))
  for (dist_name in c("normal", "gumbel_min", "gumbel_max", "logistic")) {
    eta <- .sdt_eta(grid$d, grid$criterion, grid$stimulus, 1)
    textbook <- (grid$d / 2 * (2 * grid$stimulus - 1) - grid$criterion)
    expect_equal(eta, textbook, info = dist_name)
  }
})

test_that("d equals the noise-standardized separation divided by the RMS scale", {
  grid <- expand.grid(d = c(0.4, 1.5, 3), criterion = c(-0.7, 0, 0.5),
                      sdratio = c(0.7, 1, 1.25, 1.6, 2), stimulus = c(0L, 1L))
  s <- sqrt((1 + grid$sdratio^2) / 2)
  for (dist_name in c("normal", "gumbel_min", "gumbel_max", "logistic")) {
    from_da <- dsdt_binary(1, 1, grid$stimulus, grid$d, grid$criterion,
                           sdratio = grid$sdratio, dist = dist_name)
    # the same model written with the separation in noise units
    eta_dn <- ((grid$d * s) / 2 * (2 * grid$stimulus - 1) - grid$criterion) /
      grid$sdratio^(grid$stimulus == 1)
    expect_equal(from_da, exp(.sdt_log_p_old(eta_dn, dist_name)),
                 info = dist_name)
  }
})

test_that("d is constant across iso-discriminable unequal-variance models", {
  # the property the parameterization exists for: two models with the same
  # area under the ROC must report the same sensitivity, whatever sdratio is.
  # For the Gaussian, AUC = Phi(d_N / sqrt(1 + r^2)) = Phi(d / sqrt(2)).
  auc <- function(d, r) {
    stats::pnorm(d * sqrt((1 + r^2) / 2) / sqrt(1 + r^2))
  }
  expect_equal(auc(1.5, 1.0), auc(1.5, 1.6))
  expect_equal(auc(1.5, 1.0), stats::pnorm(1.5 / sqrt(2)))
  # and it is not vacuous: the noise-standardized separation does move
  expect_false(isTRUE(all.equal(1.5 * sqrt((1 + 1.6^2) / 2), 1.5)))
})

test_that("the sdt_yn sensitivity name is guarded against a data-column clash", {
  # `d` is short enough to collide with a real column; bmm must warn rather
  # than silently emitting it as a non-linear term (#378)
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  dat <- data.frame(
    n_old = c(10, 40), stimulus = c(0L, 1L), n_trials = c(50, 50),
    d = c(1, 2)
  )
  expect_warning(
    check_formula(model, dat, bmf(d ~ d, criterion ~ 1)),
    "used both as a predicted parameter and as a column"
  )
})


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

# brmsprep stand-in: get_dpar() returns the supplied per-draw vectors already on
# the natural scale, as brms does once it has applied the inverse link.
# d/criterion are held constant so the only across-draw variation comes
# from sdratio -- isolating the unequal-variance signal scaling that the
# .sdt_eta() collapse bug got wrong.
uv_prep <- function(ndraws = 50L) {
  draws <- list(
    d    = rep(1.2, ndraws),
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
  # `d` is d_a, so the separation in noise units is d * s, and the signal row
  # (stimulus = 1) is p = Phi((d * s / 2 - criterion) / sdratio)
  s <- sqrt((1 + mk$draws$sdratio^2) / 2)
  p_expected <- pnorm((1.2 * s / 2 - 0) / mk$draws$sdratio)
  expect_equal(log_lik_sdt_binary(2L, mk$prep),
               dbinom(70L, 100L, p_expected, log = TRUE))
  # the collapse bug divided every draw by sdratio[1], giving a constant p
  expect_gt(length(unique(round(log_lik_sdt_binary(2L, mk$prep), 8))), 1L)
  # the noise row carries no scale division, but it is not sdratio-free under
  # this parameterization: holding d_a fixed while sdratio grows widens the
  # separation d * s, which moves the noise mean away from the midpoint
  expect_equal(log_lik_sdt_binary(1L, mk$prep),
               dbinom(30L, 100L, pnorm(-1.2 * s / 2), log = TRUE))
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
  expect_length(out, length(mk$draws$d))
  s <- sqrt((1 + mk$draws$sdratio^2) / 2)
  expect_equal(captured, pnorm((1.2 * s / 2) / mk$draws$sdratio))
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
  formula <- bmf(d ~ 1, criterion ~ 1)
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
  for (dist_name in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_binary("n_old", "stimulus", "n_trials", dist = dist_name)
    formula <- bmf(d ~ 1, criterion ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0, info = paste("Failed for dist:", dist_name))
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
  formula <- bmf(d ~ condition, criterion ~ 1)
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
  formula <- bmf(d ~ 1 + (1 | id), criterion ~ 1)
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
  formula <- bmf(d ~ 1, criterion ~ 1)
  prior <- default_prior(formula, data = dat, model = model)
  expect_s3_class(prior, "brmsprior")
  prior_strs <- prior$prior
  expect_true(any(grepl("normal", prior_strs)))
})
