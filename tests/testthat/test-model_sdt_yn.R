# Test Binary SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_yn fixes sdratio on the link scale and inits on the natural scale", {
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  # constant() prior on the Intercept, so exp(0) = 1 = equal variance
  expect_equal(model$fixed_parameters$sdratio, 0)
  # create_initfun() applies the forward link, so log() of these must be finite
  expect_true(all(model$init_ranges$sdratio > 0))
})

test_that("the default sdratio reaches the Stan likelihood as 1, not 0", {
  # a 0 on the link scale is an SD ratio of 1; a 0 that reached the lpmf would
  # divide the signal trials by zero, so the generated code must carry both the
  # pin and the forward link
  dat <- data.frame(n_old = c(10, 40), stimulus = c(0L, 1L),
                    n_trials = c(50, 50))
  code <- stancode(bmf(d ~ 1, criterion ~ 1), data = dat,
                   model = sdt_yn("n_old", "stimulus", "n_trials"))
  expect_match(code, "Intercept_sdratio = 0;", fixed = TRUE)
  expect_match(code, "sdratio = exp(sdratio);", fixed = TRUE)
})

test_that("an sdratio of 0 makes the signal rows' evidence infinite", {
  # why sdt_yn refuses any link on sdratio but log (test-sdt-shared.R): the
  # fixed 0 read on the identity link is a ratio of zero
  expect_equal(.sdt_eta(1.5, 0, 1L, sdratio = 0), Inf)
})

test_that("sdt_yn model has default priors", {
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  expect_true("d" %in% names(model$default_priors))
  expect_true("criterion" %in% names(model$default_priors))
  expect_true("main" %in% names(model$default_priors$d))
  expect_true("effects" %in% names(model$default_priors$d))
})

test_that("sdt_yn requires n_trials argument", {
  expect_error(sdt_yn("n_old", "stimulus"))
})


############################################################################# !
# CHECK_FORMULA TESTS                                                    ####
############################################################################# !

# check_model() must run first: it is update_model_fixed_parameters() that
# strips sdratio from model$fixed_parameters when the user gives it a formula,
# and without that is_constant(formula)[["sdratio"]] is TRUE either way, so a
# direct check_formula() call would pass while testing nothing
sdt_yn_warns <- function(formula, data = NULL) {
  if (is.null(data)) {
    data <- data.frame(n_old = c(10, 40, 15, 35), stimulus = c(0L, 1L, 0L, 1L),
                       n_trials = c(50, 50, 50, 50),
                       condition = c("A", "A", "B", "B"), id = c(1, 1, 2, 2))
  }
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  model <- check_model(model, data, formula)
  check_formula(model, data, formula)
}

test_that("a free sdratio on an all-intercept design warns", {
  # provably three parameters for two sufficient statistics: the profile
  # log-likelihood over sdratio is flat to 8e-13 on this design
  expect_warning(
    sdt_yn_warns(bmf(d ~ 1, bias ~ 1, sdratio ~ 1)),
    "every parameter formula is intercept-only"
  )
  expect_warning(
    sdt_yn_warns(bmf(d ~ 1, criterion ~ 1, sdratio ~ 1)),
    "'criterion' is the same in every row"
  )
})

test_that("a free sdratio with a constant criterion warns whatever d does", {
  # measured from the noise mean, a constant criterion is read off the
  # false-alarm rate alone, and each hit rate then fixes only (d * s - criterion)
  # / sdratio, so sdratio trades against d in every cell
  unidentified <- list(
    bmf(d ~ condition, criterion ~ 1, sdratio ~ 1),
    bmf(d ~ condition, sdratio ~ 1),
    bmf(d ~ 1, criterion ~ 1, sdratio ~ condition),
    bmf(d ~ condition, criterion = 0.5, sdratio ~ 1)
  )
  for (formula in unidentified) {
    expect_warning(sdt_yn_warns(formula), "'criterion' is the same in every row")
  }
})

test_that("a constant criterion warns even when d's form identifies sdratio", {
  # d ~ 0 + xnum holds d proportional to xnum, so with the criterion away from
  # 0 (false-alarm rate .2) the cells cannot absorb the common shift a change
  # in sdratio needs, and a profile likelihood over sdratio is sharp. The
  # warning counts terms and cannot tell this from d ~ 1 + xnum, which is
  # unidentified, so it fires on purpose
  dat <- data.frame(n_old = c(10, 20, 10, 30, 10, 40), stimulus = rep(0:1, 3),
                    n_trials = rep(50, 6), xnum = rep(c(1, 2, 3), each = 2))
  expect_warning(
    sdt_yn_warns(bmf(d ~ 0 + xnum, criterion ~ 1, sdratio ~ 1), dat),
    "'criterion' is the same in every row"
  )
})

test_that("the sdratio warning sees through the 0 + Intercept spelling", {
  # rhs_vars(d ~ 0 + Intercept) is "Intercept", so a naive length() == 0 test
  # would miss the same unidentified design
  expect_warning(
    sdt_yn_warns(bmf(d ~ 0 + Intercept, bias ~ 0 + Intercept,
                     sdratio ~ 1)),
    "every parameter formula is intercept-only"
  )
  expect_warning(
    sdt_yn_warns(bmf(d ~ condition, criterion ~ 0 + Intercept, sdratio ~ 1)),
    "'criterion' is the same in every row"
  )
})

test_that("the sdratio warning stays silent on every identified design", {
  # a criterion that moves identifies sdratio, through a predictor or a random
  # effect; with bias held constant instead, so does anything that moves d or
  # sdratio, because the anchor then moves the criterion with them
  identified <- list(
    bmf(d ~ 1, criterion ~ condition, sdratio ~ 1),
    bmf(d ~ 1, criterion ~ 1 + (1 | id), sdratio ~ 1),
    bmf(d ~ condition, bias ~ 1, sdratio ~ 1),
    bmf(d ~ 1, bias ~ 1, sdratio ~ condition),
    bmf(d ~ 1, bias ~ 1 + (1 | id), sdratio ~ 1)
  )
  for (formula in identified) {
    expect_no_warning(sdt_yn_warns(formula))
  }
  # and the default, where sdratio is fixed, never warns
  expect_no_warning(sdt_yn_warns(bmf(d ~ 1, criterion ~ 1)))
})

test_that("the sdratio warning reaches the user through bmm()", {
  dat <- data.frame(n_old = c(10, 40), stimulus = c(0L, 1L),
                    n_trials = c(50, 50))
  expect_warning(
    bmm(bmf(d ~ 1, criterion ~ 1, sdratio ~ 1), dat,
        sdt_yn("n_old", "stimulus", "n_trials"),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "'criterion' is the same in every row"
  )
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdt_yn check_data coerces an unambiguous stimulus coding only", {
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  formula <- bmf(d ~ 1, criterion ~ 1)

  # factor, logical and character are the other shapes an unambiguous 0/1 column
  # arrives in, and all three are coerced
  for (coding in list(factor(c(0, 1)), c(FALSE, TRUE), c("0", "1"))) {
    coercible_data <- data.frame(
      n_old = c(30, 40),
      stimulus = coding,
      n_trials = c(50, 50)
    )
    result <- check_data(model, coercible_data, formula)
    expect_identical(result$stimulus, c(0L, 1L), info = class(coding)[1])
  }

  # the one coding bmm must not guess at: which level is the signal?
  labelled_data <- data.frame(
    n_old = c(30, 40),
    stimulus = factor(c("noise", "signal")),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, labelled_data, formula),
    "must not guess which level is the signal"
  )
})

test_that("the coerced stimulus column is what reaches Stan", {
  # coercing in check_data is pointless if standata() still sees the factor
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  dat <- data.frame(n_old = c(30, 40), stimulus = factor(c("0", "1")),
                    n_trials = c(50, 50))
  sdata <- standata(bmf(d ~ 1, criterion ~ 1), data = dat, model = model)
  expect_equal(as.integer(sdata$vint1), c(0L, 1L))
})

test_that("sdt_yn check_data says which stimulus problem it found", {
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  formula <- bmf(d ~ 1, criterion ~ 1)
  base <- data.frame(n_old = c(30, 40), stimulus = c(0L, 1L),
                     n_trials = c(50, 50))

  # one message for every problem told the user nothing about theirs
  value_data <- base
  value_data$stimulus <- c(2, 3)
  expect_error(check_data(model, value_data, formula), "found '2', '3'")

  list_data <- base
  list_data$stimulus <- I(list(1, 2))
  expect_error(check_data(model, list_data, formula),
               "found AsIs with values")
})

test_that("sdt_yn check_data adds dist_type column", {
  formula <- bmf(d ~ 1, criterion ~ 1)
  valid_data <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )

  model <- sdt_yn("n_old", "stimulus", "n_trials", dist = "gumbel_min")
  result <- check_data(model, valid_data, formula)
  expect_equal(unique(result$dist_type), 2L)

  model2 <- sdt_yn("n_old", "stimulus", "n_trials", dist = "logistic")
  result2 <- check_data(model2, valid_data, formula)
  expect_equal(unique(result2$dist_type), 4L)
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_yn recycles vectorized parameters per observation", {
  n_old <- rsdt_yn(4, c(50, 100, 50, 100), c(0L, 1L, 0L, 1L),
                       d = c(0.5, 0.5, 2.5, 2.5), criterion = 0)
  expect_length(n_old, 4)
  expect_true(all(n_old <= c(50, 100, 50, 100)))
})

test_that("rsdt_yn generates different data for different distributions", {
  withr::local_seed(42)
  n_norm <- rsdt_yn(2, 1000, c(0L, 1L), d = 2, criterion = 0)
  withr::local_seed(42)
  n_logis <- rsdt_yn(2, 1000, c(0L, 1L), d = 2, criterion = 0,
                         dist = "logistic")
  expect_false(identical(n_norm, n_logis))
})

test_that("rsdt_yn with zero d produces equal hit and FA rates on average", {
  withr::local_seed(123)
  n_old <- rsdt_yn(2, 10000, c(0L, 1L), d = 0, criterion = 0)
  expect_equal(n_old[1] / 10000, 0.5, tolerance = 0.05)
  expect_equal(n_old[2] / 10000, 0.5, tolerance = 0.05)
})

test_that("d/rsdt_yn refuse a stimulus type they cannot read", {
  # %in% compares a factor as character, so the value check passed and .sdt_eta()
  # returned a silent NA that propagates into a hand-written pp check
  for (bad in list(factor("1"), "1")) {
    expect_error(
      dsdt_yn(68, 100, bad, d = 1.5, criterion = 0),
      "must be numeric or logical",
      info = class(bad)[1]
    )
    expect_error(
      rsdt_yn(2, 100L, bad, d = 1.5, criterion = 0),
      "must be numeric or logical",
      info = class(bad)[1]
    )
  }

  # logical stays accepted, and gives the value the integer call gives
  expect_equal(dsdt_yn(68, 100, TRUE, d = 1.5, criterion = 0),
               dsdt_yn(68, 100, 1L, d = 1.5, criterion = 0))
  expect_equal(dsdt_yn(30, 100, FALSE, d = 1.5, criterion = 0),
               dsdt_yn(30, 100, 0L, d = 1.5, criterion = 0))
  expect_length(rsdt_yn(2, 100L, c(FALSE, TRUE), d = 1.5, criterion = 0), 2)
})

test_that("dsdt_yn stays finite where the probability scale underflows", {
  # naive dbinom(y, n, pnorm(eta)) returns -Inf here and poisons loo()
  for (dist_name in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    ld <- dsdt_yn(n_old = 1, n_trials = 100, stimulus = 1,
                      d = -80, criterion = 0, dist = dist_name,
                      log = TRUE)
    expect_true(is.finite(ld), info = dist_name)
  }
})

test_that("dsdt_yn handles the y = 0 and y = n_trials edges", {
  # 0 * -Inf must contribute nothing rather than becoming NaN
  expect_equal(dsdt_yn(n_old = 0, n_trials = 10, stimulus = 1,
                           d = -80, criterion = 0, log = TRUE), 0)
  expect_equal(dsdt_yn(n_old = 10, n_trials = 10, stimulus = 1,
                           d = 80, criterion = 0, log = TRUE), 0)
})

test_that("dist names the evidence distribution, not the likelihood's CDF", {
  # bmm must reproduce the textbook generative story -- draw evidence from
  # `dist`, respond "old" when it exceeds the criterion -- for asymmetric
  # distributions too. Using F(eta) instead of S(-eta) silently fits the mirror
  # distribution and only shows up for the extreme-value families.
  d <- 1.5
  criterion <- 0.3
  bias <- criterion - d / 2

  quantile_fun <- list(
    normal     = qnorm,
    gumbel_min = function(p) log(-log(1 - p)),
    gumbel_max = function(p) -log(-log(p)),
    logistic   = qlogis
  )

  for (dist_name in names(quantile_fun)) {
    # inverse-transform sampling, independent of .sdt_dists
    u <- (seq_len(20000) - 0.5) / 20000
    p_noise <- mean(quantile_fun[[dist_name]](u) > criterion)
    p_signal <- mean(d + quantile_fun[[dist_name]](u) > criterion)

    expect_equal(
      dsdt_yn(1, 1, 0L, d, criterion, dist = dist_name), p_noise,
      tolerance = 1e-3, info = paste(dist_name, "noise")
    )
    expect_equal(
      dsdt_yn(1, 1, 1L, d, criterion, dist = dist_name), p_signal,
      tolerance = 1e-3, info = paste(dist_name, "signal")
    )
    expect_equal(
      dsdt_yn(1, 1, 0L, d, bias = bias, dist = dist_name), p_noise,
      tolerance = 1e-3, info = paste(dist_name, "noise, bias")
    )
    expect_equal(
      dsdt_yn(1, 1, 1L, d, bias = bias, dist = dist_name), p_signal,
      tolerance = 1e-3, info = paste(dist_name, "signal, bias")
    )
  }
})

test_that("sdt_d, sdt_criterion and sdt_bias invert the model's decision rule", {
  d <- 1.2
  criterion <- 0.4
  for (dist_name in c("normal", "gumbel_min", "gumbel_max", "logistic")) {
    fa <- dsdt_yn(1, 1, 0L, d, criterion, dist = dist_name)
    hit <- dsdt_yn(1, 1, 1L, d, criterion, dist = dist_name)
    expect_equal(sdt_d(hit, fa, dist = dist_name), d,
                 tolerance = 1e-8, info = dist_name)
    expect_equal(sdt_criterion(hit, fa, dist = dist_name), 0.4,
                 tolerance = 1e-8, info = dist_name)
    expect_equal(sdt_bias(hit, fa, dist = dist_name), -0.2,
                 tolerance = 1e-8, info = dist_name)
  }
})

test_that("dsdt_yn matches dbinom in the well-conditioned range", {
  # criterion from the noise mean, bias from the midpoint of the separation
  sep <- 1.5 * sqrt((1 + 1.3^2) / 2)
  expect_equal(
    dsdt_yn(n_old = 70, n_trials = 100, stimulus = 1, d = 1.5,
                criterion = 0.2, sdratio = 1.3, log = TRUE),
    dbinom(70, 100, pnorm((sep - 0.2) / 1.3), log = TRUE)
  )
  expect_equal(
    dsdt_yn(n_old = 70, n_trials = 100, stimulus = 1, d = 1.5,
                bias = 0.2, sdratio = 1.3, log = TRUE),
    dbinom(70, 100, pnorm((sep / 2 - 0.2) / 1.3), log = TRUE)
  )
  expect_equal(
    dsdt_yn(n_old = 30, n_trials = 100, stimulus = 0, d = 1.5,
                criterion = 0.2, sdratio = 1.3, log = TRUE),
    dbinom(30, 100, pnorm(-0.2), log = TRUE)
  )
})

test_that("dsdt_yn and rsdt_yn take exactly one of criterion and bias", {
  expect_error(dsdt_yn(70, 100, 1, d = 1.5), "exactly one of 'criterion'")
  expect_error(dsdt_yn(70, 100, 1, d = 1.5, criterion = 0.5, bias = 0),
               "exactly one of 'criterion'")
  expect_error(rsdt_yn(2, 100, 1, d = 1.5), "exactly one of 'criterion'")
})

test_that("dsdt_yn returns one value per posterior draw for a scalar count", {
  # times_nonzero() must recycle: ifelse() alone would collapse to length 1
  dens <- dsdt_yn(n_old = 70, n_trials = 100, stimulus = 1,
                      d = rep(1.5, 20), criterion = 0.2, log = TRUE)
  expect_length(dens, 20)
})

test_that("sdratio never divides the noise trials, it only widens the separation", {
  # Only the signal distribution carries the scale division. On a noise trial
  # sdratio therefore acts solely through d_a -> d_N = d * s, so the same
  # density must be reachable by an equal-variance call with the widened d.
  s <- sqrt((1 + 1.5^2) / 2)
  uv <- dsdt_yn(n_old = 70, n_trials = 100, stimulus = 0,
                    d = 1.5, criterion = 0.2, sdratio = 1.5)
  ev <- dsdt_yn(n_old = 70, n_trials = 100, stimulus = 0,
                    d = 1.5 * s, criterion = 0.2, sdratio = 1)
  expect_equal(uv, ev)
})

test_that("n_trials may vary within a design", {
  # every other test holds n_trials constant, so a per-row count that silently
  # recycled the first cell would go unnoticed; broeder_schuetz_2009_e3 has
  # 6, 15, 30, 45 and 54 trials per cell
  n_trials <- c(6, 15, 30)
  n_old <- c(1, 10, 22)
  eta <- 1.5 - 0.2
  expect_equal(
    dsdt_yn(n_old, n_trials, 1L, d = 1.5, criterion = 0.2, log = TRUE),
    dbinom(n_old, n_trials, pnorm(eta), log = TRUE)
  )

  sdata <- standata(bmf(d ~ 1, criterion ~ 0 + condition),
                    data = broeder_schuetz_2009_e3,
                    model = sdt_yn("n_old", "stimulus", "n_trials"))
  expect_equal(sort(unique(as.integer(sdata$trials))), c(6L, 15L, 30L, 45L, 54L))
})


############################################################################# !
# SENSITIVITY PARAMETERIZATION (d = d_a)                                 ####
############################################################################# !

test_that("d is a no-op relative to d' for every equal-variance call", {
  # sdratio = 1 gives s = 1, so nothing in the reparameterization may move
  grid <- expand.grid(d = c(0.4, 1.5, 3), criterion = c(-0.7, 0, 0.5),
                      stimulus = c(0L, 1L))
  eta <- .sdt_eta(grid$d, grid$criterion, grid$stimulus, 1)
  expect_equal(eta, grid$d * grid$stimulus - grid$criterion)
})

test_that("d equals the noise-standardized separation divided by the RMS scale", {
  grid <- expand.grid(d = c(0.4, 1.5, 3), criterion = c(-0.7, 0, 0.5),
                      sdratio = c(0.7, 1, 1.25, 1.6, 2), stimulus = c(0L, 1L))
  s <- sqrt((1 + grid$sdratio^2) / 2)
  for (dist_name in c("normal", "gumbel_min", "gumbel_max", "logistic")) {
    from_da <- dsdt_yn(1, 1, grid$stimulus, grid$d, grid$criterion,
                           sdratio = grid$sdratio, dist = dist_name)
    # the same model written with the separation in noise units
    eta_dn <- ((grid$d * s) * grid$stimulus - grid$criterion) /
      grid$sdratio^(grid$stimulus == 1)
    expect_equal(from_da, exp(.sdt_log_p_old(eta_dn, dist_name)),
                 info = dist_name)
  }
})

test_that("d is constant across iso-discriminable unequal-variance models", {
  # the property the parameterization exists for: two models with the same
  # area under the ROC must report the same sensitivity, whatever sdratio is.
  # AUC is read off the model's own eta, so this goes red if .sdt_eta stops
  # dividing the separation by the RMS scale.
  d <- 1.5
  sdratio <- c(0.7, 1, 1.25, 1.6, 2)
  # a scalar stimulus against a vector of sdratio draws, the shape ROC points
  # arrive in; an ifelse() scale would apply sdratio[1] to all of them
  eta_noise <- .sdt_eta(d, 0, 0L, sdratio = sdratio)
  eta_signal <- .sdt_eta(d, 0, 1L, sdratio = sdratio)

  # undo the signal-scale division to get back the separation in noise units,
  # then AUC for two normals with SDs 1 and sdratio
  separation <- eta_signal * sdratio - eta_noise
  expect_equal(stats::pnorm(separation / sqrt(1 + sdratio^2)),
               rep(stats::pnorm(d / sqrt(2)), length(sdratio)))
  # not vacuous: the separation the model implies genuinely moves with sdratio
  expect_equal(separation, d * sqrt((1 + sdratio^2) / 2))
})

test_that("the sdt_yn sensitivity name is guarded against a data-column clash", {
  # `d` is short enough to collide with a real column; bmm must warn rather
  # than silently emitting it as a non-linear term (#378)
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  dat <- data.frame(
    n_old = c(10, 40), stimulus = c(0L, 1L), n_trials = c(50, 50),
    d = c(1, 2)
  )
  # the #378 warning is raised in check_formula.bmmodel and does not depend on fixed_parameters
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

test_that("log_lik_sdt_yn scales the signal per draw under UV-SDT", {
  mk <- uv_prep()
  local_mocked_bindings(
    get_dpar = function(prep, dpar, i = NULL, ...) mk$draws[[dpar]],
    .package = "brms"
  )
  # `d` is d_a, so the separation in noise units is d * s, and the signal row
  # (stimulus = 1) is p = Phi((d * s - criterion) / sdratio)
  s <- sqrt((1 + mk$draws$sdratio^2) / 2)
  p_expected <- pnorm((1.2 * s - 0) / mk$draws$sdratio)
  expect_equal(log_lik_sdt_yn(2L, mk$prep),
               dbinom(70L, 100L, p_expected, log = TRUE))
  # the collapse bug divided every draw by sdratio[1], giving a constant p
  expect_gt(length(unique(round(log_lik_sdt_yn(2L, mk$prep), 8))), 1L)
  # criterion is measured from the noise mean, so the noise row depends on it
  # alone, whatever d and sdratio do
  expect_equal(log_lik_sdt_yn(1L, mk$prep),
               rep(dbinom(30L, 100L, 0.5, log = TRUE), length(s)))
})

test_that("log_lik_sdt_yn reads a bias from the midpoint of the separation", {
  mk <- uv_prep()
  names(mk$draws)[names(mk$draws) == "criterion"] <- "bias"
  mk$prep$dpars <- mk$draws
  local_mocked_bindings(
    get_dpar = function(prep, dpar, i = NULL, ...) mk$draws[[dpar]],
    .package = "brms"
  )
  # with bias = 0 the criterion sits halfway along the separation d * s, which
  # grows with sdratio, so the noise row moves with it
  s <- sqrt((1 + mk$draws$sdratio^2) / 2)
  expect_equal(log_lik_sdt_yn(1L, mk$prep),
               dbinom(30L, 100L, pnorm(-1.2 * s / 2), log = TRUE))
  expect_equal(log_lik_sdt_yn(2L, mk$prep),
               dbinom(70L, 100L, pnorm((1.2 * s / 2) / mk$draws$sdratio),
                      log = TRUE))
})

test_that("posterior_predict_sdt_yn feeds per-draw signal probabilities", {
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
  out <- posterior_predict_sdt_yn(2L, mk$prep)
  expect_length(out, length(mk$draws$d))
  s <- sqrt((1 + mk$draws$sdratio^2) / 2)
  expect_equal(captured, pnorm((1.2 * s) / mk$draws$sdratio))
})


############################################################################# !
# FORMULA AND CONFIGURE_MODEL TESTS                                      ####
############################################################################# !

test_that("sdt_yn model produces valid stancode with brms", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    id = c(1, 1, 2, 2)
  )
  model <- sdt_yn("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(d ~ 1, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_match(code, paste("target += sdt_yn_lpmf(Y[n] | mu[n], d[n], criterion[n],",
                           "sdratio[n], vint1[n], vint2[n], trials[n]);"), fixed = TRUE)
  expect_match(code, "real sdt_log_cumprob(", fixed = TRUE)
})

############################################################################# !
# R <-> STAN DIST_TYPE CONTRACT                                          ####
############################################################################# !

test_that("the Stan dist_type branches match the registry order", {
  # The R side dispatches by position in .sdt_dists; the Stan side hardcodes
  # the integers. Nothing else in the suite compares the two, so swapping two
  # branches in the chunk would leave every test green while turning every
  # dist = "normal" fit into a Gumbel model. Text-level, as in
  # test-model_cswald.R, so it costs no compilation.
  sc <- read_lines2(file.path(system.file("stan_chunks", package = "bmm"),
                              "sdt_dist_funs.stan"))
  dispatcher <- function(fun) {
    sub("(?s)\n\\}.*", "",
        sub(paste0("(?s).*real ", fun, "\\("), "", sc, perl = TRUE), perl = TRUE)
  }
  branches <- list(
    sdt_cumprob = c("Phi\\(eta\\)", "1 - exp\\(-exp\\(eta\\)\\)",
                    "exp\\(-exp\\(-eta\\)\\)", "inv_logit\\(eta\\)"),
    sdt_log_cumprob = c("std_normal_lcdf\\(eta\\)",
                        "log1m_exp\\(-exp\\(eta\\)\\)",
                        "-exp\\(-eta\\)", "-log1p_exp\\(-eta\\)"),
    sdt_log_one_minus_cumprob = c("std_normal_lcdf\\(-eta\\)", "-exp\\(eta\\)",
                                  "log1m_exp\\(-exp\\(-eta\\)\\)",
                                  "-eta - log1p_exp\\(-eta\\)"),
    sdt_quantile = c("inv_Phi\\(u\\)", "log\\(-log1m\\(u\\)\\)",
                     "-log\\(-log\\(u\\)\\)", "logit\\(u\\)")
  )

  expect_equal(names(.sdt_dists),
               c("normal", "gumbel_min", "gumbel_max", "logistic"))
  for (fun in names(branches)) {
    for (i in seq_along(.sdt_dists)) {
      expect_match(dispatcher(fun),
                   paste0("dist_type == ", i, "\\) return ",
                          branches[[fun]][i]),
                   info = paste(fun, names(.sdt_dists)[i]))
    }
  }
})

# posterior_epred (#475) ------------------------------------------------------

test_that("posterior_epred_sdt_yn() is the mean count posterior_predict_sdt_yn() simulates", {
  withr::local_seed(475)
  sets <- data.frame(
    d = c(1.5, 0.5, 2, 1),
    criterion = c(0.2, -0.3, 0.8, 0),
    sdratio = c(1, 1.3, 0.8, 1.5)
  )
  for (dist in .sdt_dist_names) {
    data <- list(vint1 = c(1, 0, 1, 0), vint2 = rep(.sdt_dist_id(dist), 4),
                 trials = c(50, 80, 20, 100))
    res <- epred_vs_predict(sets, posterior_epred_sdt_yn, posterior_predict_sdt_yn,
                            data = data)
    expect_lt(max(abs(res[, "rel_error"])), 0.02)
  }
})

test_that("posterior_epred_sdt_yn() returns one column per observation", {
  dpars <- list(d = matrix(c(1.5, 0.5, 2, 1, 0.7, 1.2), 2),
                criterion = matrix(c(0.2, -0.3, 0.8, 0, 0.1, -0.5), 2),
                sdratio = matrix(c(1, 1.3, 0.8, 1.5, 1.1, 0.9), 2))
  # a stimulus that differs across observations but not across draws is what
  # a data column recycled down the draws gets wrong
  data <- list(vint1 = c(1, 0, 1), vint2 = rep(.sdt_dist_id("gumbel_min"), 3),
               trials = c(50, 80, 20))
  expect_epred_by_cell(posterior_epred_sdt_yn, dpars, data)
})

test_that("an sdt_yn fit stores posterior_epred_sdt_yn() in its family", {
  dat <- data.frame(n_old = c(10, 40), stimulus = c(0L, 1L), n_trials = c(50, 50))
  fit <- bmm(bmf(d ~ 1, criterion ~ 1), dat, sdt_yn("n_old", "stimulus", "n_trials"),
             backend = "mock", mock_fit = 1, rename = FALSE)
  expect_identical(fit$formula$family$posterior_epred, posterior_epred_sdt_yn)
})

test_that("both Stan yes/no likelihoods match dsdt_yn", {
  skip_unless_stan_tests()
  skip_on_cran()

  sc <- system.file("stan_chunks", package = "bmm")
  program <- paste0(
    "functions {\n", read_lines2(file.path(sc, "sdt_dist_funs.stan")), "\n",
    read_lines2(file.path(sc, "sdt_yn_funs.stan")), "\n}\n",
    "data {\n  int N; array[N] int y; array[N] int trials; array[N] int stimulus;\n",
    "  array[N] int dist; vector[N] d; vector[N] crit; vector[N] sdratio;\n}\n",
    "generated quantities {\n  vector[N] lp_crit; vector[N] lp_bias;\n",
    "  for (i in 1:N) {\n",
    "    lp_crit[i] = sdt_yn_lpmf(y[i] | 0, d[i], crit[i], sdratio[i], stimulus[i], dist[i], trials[i]);\n",
    "    lp_bias[i] = sdt_yn_bias_lpmf(y[i] | 0, d[i], crit[i], sdratio[i], stimulus[i], dist[i], trials[i]);\n",
    "  }\n}\n"
  )
  grid <- expand.grid(stimulus = 0:1, dist = seq_along(bmm:::.sdt_dists),
                      d = c(0.4, 2), crit = c(-0.6, 0.3, 1.4),
                      sdratio = c(0.7, 1, 1.6))
  grid$y <- rep_len(c(3L, 17L, 29L), nrow(grid))
  grid$trials <- 30L
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(list(N = nrow(grid)), as.list(grid)), fixed_param = TRUE,
    chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
    show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#",
                         check.names = FALSE)
  idx <- seq_len(nrow(grid))
  dist_names <- names(bmm:::.sdt_dists)[grid$dist]
  r_crit <- vapply(idx, function(i) with(grid[i, ], dsdt_yn(
    y, trials, stimulus, d, criterion = crit, sdratio = sdratio,
    dist = dist_names[i], log = TRUE)), numeric(1))
  r_bias <- vapply(idx, function(i) with(grid[i, ], dsdt_yn(
    y, trials, stimulus, d, bias = crit, sdratio = sdratio,
    dist = dist_names[i], log = TRUE)), numeric(1))
  expect_equal(as.numeric(csv[1, paste0("lp_crit.", idx)]), r_crit, tolerance = 1e-10)
  expect_equal(as.numeric(csv[1, paste0("lp_bias.", idx)]), r_bias, tolerance = 1e-10)
})

test_that("a formula with bias selects the bias family and dpar", {
  dat <- data.frame(n_old = c(10, 40), stimulus = c(0, 1), n_trials = 50)
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  code_bias <- bmm::stancode(bmf(d ~ 1, bias ~ 1), data = dat, model = model)
  expect_match(code_bias, "sdt_yn_bias_lpmf(Y[n] | mu[n], d[n], bias[n]", fixed = TRUE)
  code_crit <- bmm::stancode(bmf(d ~ 1), data = dat, model = model)
  expect_match(code_crit, "sdt_yn_lpmf(Y[n] | mu[n], d[n], criterion[n]", fixed = TRUE)
  expect_error(bmm::stancode(bmf(d ~ 1, criterion ~ 1, bias ~ 1),
                             data = dat, model = model),
               "gives both 'criterion', 'bias'")
})

test_that("an unbiased observer is bias fixed at 0", {
  dat <- data.frame(n_old = c(10, 40), stimulus = c(0, 1), n_trials = 50)
  model <- sdt_yn("n_old", "stimulus", "n_trials")
  prior <- default_prior(bmf(d ~ 1, bias = 0), data = dat, model = model)
  expect_true(any(prior$dpar == "bias" & grepl("constant\\(0\\)", prior$prior)))
  expect_false(any(prior$dpar == "criterion"))
})
