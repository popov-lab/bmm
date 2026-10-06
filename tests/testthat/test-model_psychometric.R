# Test the psychometric function model

psy_sigmoids <- c("normal", "logistic", "gumbel_min", "gumbel_max",
                  "weibull", "lognormal", "loglogistic")

psy_trial_data <- function(n = 40) {
  data.frame(
    x = rep(c(-1, 0, 1, 2), length.out = n),
    y = rep(c(0L, 1L, 1L, 1L), length.out = n),
    id = factor(rep(1:4, length.out = n))
  )
}

psy_count_data <- function() {
  data.frame(
    contrast = rep(c(0.01, 0.02, 0.04, 0.08, 0.16), 4),
    n_trials = 40L,
    n_correct = rep(c(20L, 23L, 30L, 37L, 39L), 4),
    chance = rep(c(0.5, 0.25), each = 10)
  )
}

############################################################################# !
# MODEL CONSTRUCTOR                                                      ####
############################################################################# !

test_that("psychometric() estimates exactly the rates it is not given", {
  pars <- function(...) names(psychometric("y", "x", ...)$parameters)
  expect_equal(pars(), c("midpoint", "width", "guess", "lapse"))
  expect_equal(pars(guess = 0.5), c("midpoint", "width", "lapse"))
  expect_equal(pars(lapse = 0), c("midpoint", "width", "guess"))
  expect_equal(pars(guess = "chance", lapse = 0), c("midpoint", "width"))
})

test_that("default priors, links and inits cover the estimated parameters only", {
  model <- psychometric("y", "x", guess = 0.5)
  expect_named(model$links, names(model$parameters))
  expect_named(model$default_priors, names(model$parameters))
  expect_named(model$init_ranges, names(model$parameters))
})

test_that("log-scale sigmoids give the midpoint a log link", {
  for (sigmoid in psy_sigmoids) {
    model <- psychometric("y", "x", sigmoid = sigmoid)
    log_x <- sigmoid %in% c("weibull", "lognormal", "loglogistic")
    expect_equal(model$links$midpoint, if (log_x) "log" else "identity",
                 info = sigmoid)
    expect_equal(model$links$width, "log")
  }
})

test_that("psychometric() refuses rates outside [0, 1) and malformed values", {
  for (bad in list(1, -0.1, c(0.5, 0.25), NA_real_, TRUE, c("a", "b"))) {
    expect_error(psychometric("y", "x", guess = bad), "guess must be")
    expect_error(psychometric("y", "x", lapse = bad), "lapse must be")
  }
  expect_silent(psychometric("y", "x", guess = 0, lapse = 0))
})

test_that("psychometric() requires response and intensity", {
  expect_error(psychometric("y"))
  expect_error(psychometric(intensity = "x"))
})

test_that("psychometric() refuses a link changed after construction", {
  model <- psychometric("y", "x")
  model$links$width <- "identity"
  expect_error(check_links(model), "cannot be changed")
})


############################################################################# !
# DATA-SCALED DEFAULT PRIORS                                             ####
############################################################################# !

prior_of <- function(pr, dpar, class = "Intercept") {
  pr[pr$dpar == dpar & pr$class == class & pr$coef == "" & pr$group == "", "prior"]
}

test_that("midpoint and width priors are scaled to the intensity range", {
  # range -1 to 2: centre 0.5, R = 3
  pr <- default_prior(bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1),
                      psy_trial_data(), psychometric("y", "x"))
  expect_equal(prior_of(pr, "midpoint"), "normal(0.5, 1.5)")
  expect_equal(prior_of(pr, "width"), glue("normal({signif(log(1.5), 3)}, 1)"))
  expect_equal(prior_of(pr, "guess"), "normal(-3, 1.3)")
  expect_equal(prior_of(pr, "lapse"), "normal(-3, 1.3)")

  pr <- default_prior(bmf(midpoint ~ 0 + cond, width ~ 1, guess ~ 1, lapse ~ 1),
                      transform(psy_trial_data(), cond = factor(rep(1:2, 20))),
                      psychometric("y", "x"))
  expect_equal(prior_of(pr, "midpoint", "b"), "normal(0.5, 1.5)")
})

test_that("log-scale sigmoids scale their priors to the log range", {
  dat <- psy_count_data()
  model <- psychometric("n_correct", "contrast", n_trials = "n_trials",
                        sigmoid = "weibull", guess = 0.5)
  pr <- default_prior(bmf(midpoint ~ 1, width ~ 1, lapse ~ 1),
                      dat, model)
  log_range <- log(c(0.01, 0.16))
  expect_equal(prior_of(pr, "midpoint"),
               glue("normal({signif(mean(log_range), 3)}, {signif(diff(log_range) / 2, 3)})"))
  expect_equal(prior_of(pr, "width"),
               glue("normal({signif(log(diff(log_range) / 2), 3)}, 1)"))
})

test_that("a user prior on the midpoint replaces the data-scaled default", {
  user <- brms::set_prior("normal(0, 5)", class = "Intercept", dpar = "midpoint")
  fit <- bmm(bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1),
             psy_trial_data(), psychometric("y", "x"), prior = user,
             backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(prior_of(fit$prior, "midpoint"), "normal(0, 5)")
})

test_that("a single intensity level centres the priors on it", {
  pr <- suppressWarnings(default_prior(
    bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1),
    transform(psy_trial_data(), x = 40), psychometric("y", "x")
  ))
  expect_equal(prior_of(pr, "midpoint"), "normal(40, 20)")
})

test_that("count data are treated as aggregated by bmm_data_check()", {
  expect_true(uses_aggregate_data(psychometric("k", "x", n_trials = "n")))
  expect_false(uses_aggregate_data(psychometric("y", "x")))
})

test_that("inits for midpoint and width lie inside the data range", {
  model <- check_model(psychometric("y", "x"), psy_trial_data())
  expect_true(all(model$init_ranges$midpoint >= -1 & model$init_ranges$midpoint <= 2))
  expect_true(all(model$init_ranges$width > 0 & model$init_ranges$width <= 3))

  model <- check_model(
    psychometric("n_correct", "contrast", n_trials = "n_trials", sigmoid = "weibull"),
    psy_count_data()
  )
  expect_true(all(model$init_ranges$midpoint >= 0.01 & model$init_ranges$midpoint <= 0.16))
})


############################################################################# !
# CHECK_DATA AND CHECK_FORMULA                                           ####
############################################################################# !

psy_check <- function(model, data, formula = bmf(midpoint ~ 1, width ~ 1)) {
  model <- check_model(model, data, formula)
  check_data(model, data, formula)
}

test_that("trial-level responses must be 0/1, and logical ones are accepted", {
  dat <- psy_trial_data()
  dat$y[1] <- 2L
  expect_error(psy_check(psychometric("y", "x"), dat), "must be 0/1")

  dat <- psy_trial_data()
  dat$y <- dat$y == 1L
  out <- psy_check(psychometric("y", "x"), dat)
  expect_identical(out$y, as.integer(dat$y))
  expect_true(all(out$psy_trials == 1L))
})

test_that("counts are checked against n_trials", {
  dat <- psy_count_data()
  dat$n_correct[1] <- 41L
  expect_error(
    psy_check(psychometric("n_correct", "contrast", n_trials = "n_trials"), dat),
    "must not exceed"
  )
})

test_that("the intensity column is validated", {
  model <- psychometric("y", "x")
  expect_error(psy_check(model, psy_trial_data()[, c("y", "id")]), "missing in the data")
  expect_error(psy_check(model, transform(psy_trial_data(), x = as.character(x))),
               "must be numeric")
  dat <- psy_trial_data()
  dat$x[3] <- NA
  expect_error(psy_check(model, dat), "must be finite")
  dat$x[3] <- Inf
  expect_error(psy_check(model, dat), "must be finite")
  expect_warning(psy_check(model, transform(psy_trial_data(), x = 1)),
                 "single value")
})

test_that("log-scale sigmoids refuse non-positive intensities and name the linear twin", {
  expect_error(
    psy_check(psychometric("y", "x", sigmoid = "weibull"), psy_trial_data()),
    "must be positive.*'gumbel_min'"
  )
  expect_error(
    psy_check(psychometric("y", "x", sigmoid = "lognormal"), psy_trial_data()),
    "'normal'"
  )
})

test_that("a rate column must hold probabilities in [0, 1)", {
  model <- psychometric("n_correct", "contrast", n_trials = "n_trials",
                        guess = "chance")
  expect_silent(psy_check(model, psy_count_data()))
  expect_error(psy_check(model, transform(psy_count_data(), chance = 1)),
               "guess rate column 'chance'")
  expect_error(psy_check(model, psy_count_data()[, -4]),
               "column 'chance' is missing")
})

test_that("reserved helper columns are announced before they are overwritten", {
  expect_warning(
    psy_check(psychometric("y", "x"), transform(psy_trial_data(), psy_dist = 9)),
    "psy_dist.*reserved"
  )
})

test_that("fixed rates reach the data as values, estimated ones as placeholders", {
  out <- psy_check(psychometric("y", "x", guess = 0.5), psy_trial_data())
  expect_true(all(out$psy_guess == 0.5))
  expect_true(all(out$psy_lapse == 0))
  expect_true(all(out$psy_logx == 0L))

  out <- psy_check(
    psychometric("n_correct", "contrast", n_trials = "n_trials",
                 sigmoid = "weibull", guess = "chance"),
    psy_count_data()
  )
  expect_false("psy_guess" %in% names(out))
  expect_true(all(out$psy_logx == 1L))
  expect_true(all(out$psy_dist == .sdt_dist_id("gumbel_min")))
})

test_that("intensity is refused in the midpoint and width formulas only", {
  dat <- psy_trial_data()
  model <- psychometric("y", "x")
  expect_error(
    bmm(bmf(midpoint ~ x, width ~ 1, guess ~ 1, lapse ~ 1), dat, model,
        backend = "mock", mock_fit = 1, rename = FALSE),
    "'midpoint' uses the intensity variable 'x'"
  )
  expect_error(
    bmm(bmf(midpoint ~ 1, width ~ 1 + (1 | x), guess ~ 1, lapse ~ 1), dat, model,
        backend = "mock", mock_fit = 1, rename = FALSE),
    "'width' uses the intensity variable"
  )
  expect_s3_class(
    bmm(bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ x), dat, model,
        backend = "mock", mock_fit = 1, rename = FALSE),
    "bmmfit"
  )
})


############################################################################# !
# CONFIGURE_MODEL                                                        ####
############################################################################# !

test_that("each combination of fixed rates compiles to its own family", {
  dat <- psy_count_data()
  cases <- list(
    psychometric = psychometric("n_correct", "contrast", "n_trials"),
    psychometric_fixguess = psychometric("n_correct", "contrast", "n_trials", guess = 0.5),
    psychometric_fixlapse = psychometric("n_correct", "contrast", "n_trials", lapse = 0),
    psychometric_fixboth = psychometric("n_correct", "contrast", "n_trials",
                                        guess = "chance", lapse = 0)
  )
  for (family in names(cases)) {
    model <- cases[[family]]
    sc <- suppressMessages(stancode(bmf(midpoint ~ 1, width ~ 1), dat, model))
    rates <- intersect(c("guess", "lapse"), names(model$parameters))
    dpars <- paste0(c("midpoint", "width", rates), "[n]", collapse = ", ")
    expect_match(
      sc,
      glue("target += {family}_lpmf(Y[n] | mu[n], {dpars}, trials[n], vreal1[n], \\
           vreal2[n], vreal3[n], vint1[n], vint2[n], psy_zmid, psy_zspan);"),
      fixed = TRUE, info = family
    )
    skip_if_not_installed("rstan")
    expect_true(rstan::stanc(model_code = sc)$status, info = family)
  }
})

test_that("a column guess rate is read from the user's column", {
  model <- psychometric("n_correct", "contrast", "n_trials", guess = "chance")
  fit <- bmm(bmf(midpoint ~ 1, width ~ 1, lapse ~ 1), psy_count_data(), model,
             backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(as.numeric(brms::standata(fit)$vreal2), psy_count_data()$chance)
  expect_equal(brms::standata(fit)$psy_zspan,
               .psychometric_z_constants("normal")[["zspan"]])
})


############################################################################# !
# DISTRIBUTION FUNCTIONS                                                 ####
############################################################################# !

# Independent forms of each sigmoid, written from their textbook definitions
# rather than from .sdt_dists: location-scale CDFs on x (or log x), with the
# location and scale implied by midpoint and width
psy_oracle_cdf <- function(x, midpoint, width, sigmoid) {
  cdf <- switch(sigmoid,
    normal = , lognormal = stats::pnorm,
    logistic = , loglogistic = stats::plogis,
    gumbel_min = , weibull = function(q) 1 - exp(-exp(q)),
    gumbel_max = function(q) exp(-exp(-q))
  )
  qf <- switch(sigmoid,
    normal = , lognormal = stats::qnorm,
    logistic = , loglogistic = stats::qlogis,
    gumbel_min = , weibull = function(p) log(-log(1 - p)),
    gumbel_max = function(p) -log(-log(p))
  )
  scale <- width / (qf(0.95) - qf(0.05))
  location <- if (sigmoid %in% c("weibull", "lognormal", "loglogistic")) log(midpoint) else midpoint
  axis <- if (sigmoid %in% c("weibull", "lognormal", "loglogistic")) log(x) else x
  cdf(qf(0.5) + (axis - location) / scale)
}

test_that("dpsychometric() is the binomial of guess + (1 - guess)(1 - lapse) F", {
  x <- c(0.2, 0.5, 1, 2, 4)
  for (sigmoid in psy_sigmoids) {
    for (rates in list(c(0, 0), c(0.5, 0.04), c(0.1, 0))) {
      psi <- rates[1] + (1 - rates[1]) * (1 - rates[2]) *
        psy_oracle_cdf(x, midpoint = 1, width = 1.5, sigmoid)
      expect_equal(
        dpsychometric(c(3, 5, 6, 9, 10), x, midpoint = 1, width = 1.5,
                      guess = rates[1], lapse = rates[2], n_trials = 10,
                      sigmoid = sigmoid),
        stats::dbinom(c(3, 5, 6, 9, 10), 10, psi),
        tolerance = 1e-12, info = paste(sigmoid, rates[1], rates[2])
      )
    }
  }
})

test_that("the Weibull sigmoid is the classical 1 - exp(-(x / alpha)^beta)", {
  midpoint <- 0.04
  width <- 2
  k <- .psychometric_z_constants("weibull")
  beta <- k[["zspan"]] / width
  alpha <- exp(log(midpoint) - k[["zmid"]] / beta)
  x <- c(0.005, 0.01, 0.04, 0.1, 0.5)
  expect_equal(
    dpsychometric(rep(1, 5), x, midpoint, width, sigmoid = "weibull"),
    1 - exp(-(x / alpha)^beta),
    tolerance = 1e-12
  )
})

test_that("midpoint and width put F at 0.5, 0.05 and 0.95 where they say", {
  for (sigmoid in psy_sigmoids) {
    at <- .psychometric_x_at(c(0.05, 0.5, 0.95), midpoint = 2, width = 1.2, sigmoid)
    expect_equal(at[2], 2, info = sigmoid)
    spread <- if (sigmoid %in% c("weibull", "lognormal", "loglogistic")) log(at[3] / at[1]) else at[3] - at[1]
    expect_equal(spread, 1.2, info = sigmoid)
    expect_equal(
      dpsychometric(rep(1, 3), at, midpoint = 2, width = 1.2, sigmoid = sigmoid),
      c(0.05, 0.5, 0.95), info = sigmoid
    )
  }
})

test_that("dpsychometric() stays finite far in the tails and at zero rates", {
  # a width of zspan makes z the intensity itself. psi rounds to 0 at z = -30,
  # but the log density does not
  zspan <- .psychometric_z_constants("normal")[["zspan"]]
  lp <- dpsychometric(1, -30, midpoint = 0, width = zspan, log = TRUE)
  expect_true(is.finite(lp))
  expect_equal(lp, stats::pnorm(-30, log.p = TRUE), tolerance = 1e-12)
  expect_equal(dpsychometric(0, -30, 0, zspan, log = TRUE), 0)
  # a lapse of zero makes a miss impossible only in the limit, not NaN
  expect_true(is.finite(dpsychometric(0, 30, 0, zspan, lapse = 0, log = TRUE)))
  expect_false(anyNA(dpsychometric(0:1, c(-5, 5), 0, 1, guess = 0, lapse = 0, log = TRUE)))
})

test_that("rpsychometric() draws counts with mean n_trials * psi", {
  withr::local_seed(1234)
  y <- rpsychometric(20000, 1, midpoint = 1, width = 2, guess = 0.5,
                     lapse = 0.1, n_trials = 10)
  # F = 0.5 at the midpoint: psi = 0.5 + 0.5 * 0.9 * 0.5
  expect_equal(mean(y) / (10 * 0.725), 1, tolerance = 0.01)
  expect_true(all(y >= 0 & y <= 10))
})

test_that("the distribution functions validate their arguments", {
  expect_error(dpsychometric(2, 1, 1, 1, n_trials = 1), "must not exceed")
  expect_error(dpsychometric(1, 1, 1, 0), "width must be positive")
  expect_error(dpsychometric(1, 1, 1, 1, guess = 1), "guess must be")
  expect_error(dpsychometric(1, 0, 1, 1, sigmoid = "weibull"), "intensity must be positive")
  expect_error(rpsychometric(1, 1, -1, 1, sigmoid = "lognormal"), "midpoint must be positive")
})


############################################################################# !
# R AND STAN PARITY                                                      ####
############################################################################# !

test_that("psychometric_log_lik in Stan matches dpsychometric()", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr not available")
  skip_if_not(!is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
              "CmdStan not installed")

  sc_path <- system.file("stan_chunks", package = "bmm")
  code <- glue(
    "functions {{\n{read_lines2(file.path(sc_path, 'sdt_dist_funs.stan'))}\n",
    "{read_lines2(file.path(sc_path, 'psychometric_funs.stan'))}\n}}\n",
    "data {{\n  int N;\n  array[N] int y;\n  array[N] int trials;\n",
    "  vector[N] x;\n  vector[N] midpoint;\n  vector[N] width;\n",
    "  vector[N] guess;\n  vector[N] lapse;\n  int dist_type;\n  int log_x;\n",
    "  real zmid;\n  real zspan;\n}}\n",
    "generated quantities {{\n  vector[N] lp;\n  for (n in 1:N) lp[n] = ",
    "psychometric_log_lik(y[n], trials[n], x[n], midpoint[n], width[n], ",
    "guess[n], lapse[n], dist_type, log_x, zmid, zspan);\n}}\n"
  )
  file <- file.path(tempdir(), "bmm_psychometric_parity.stan")
  writeLines(code, file)
  stan_model <- cmdstanr::cmdstan_model(file, quiet = TRUE)

  # a grid crossing tails, the midpoint, and rates fixed at exactly zero
  grid <- expand.grid(
    x = c(0.01, 0.3, 1, 3, 40), y = c(0L, 4L, 10L), midpoint = c(0.5, 2),
    width = c(0.3, 3), guess = c(0, 0.5), lapse = c(0, 0.03)
  )
  for (sigmoid in psy_sigmoids) {
    spec <- .psychometric_sigmoids[[sigmoid]]
    k <- .psychometric_z_constants(sigmoid)
    sdata <- c(
      as.list(grid), N = nrow(grid), list(trials = rep(10L, nrow(grid))),
      dist_type = .sdt_dist_id(spec$dist), log_x = as.integer(spec$log_x),
      zmid = k[["zmid"]], zspan = k[["zspan"]]
    )
    fit <- stan_model$sample(
      data = sdata, chains = 1, iter_sampling = 1, fixed_param = TRUE,
      refresh = 0, show_messages = FALSE, sig_figs = 18, seed = 1
    )
    lp_stan <- as.numeric(fit$draws("lp", format = "draws_matrix")[1, ])
    lp_r <- with(grid, dpsychometric(y, x, midpoint, width, guess, lapse,
                                     n_trials = 10, sigmoid = sigmoid, log = TRUE))
    expect_identical(is.finite(lp_stan), is.finite(lp_r), info = sigmoid)
    # relative, because the extreme-value sigmoids reach log densities near
    # -1e224 on this grid, where an absolute gap of 1e211 is rounding
    finite <- is.finite(lp_r)
    rel_error <- abs(lp_stan - lp_r) / pmax(1, abs(lp_r))
    expect_lt(max(rel_error[finite]), 1e-12)
  }
})


############################################################################# !
# FULL PIPELINE                                                          ####
############################################################################# !

test_that("trial-level and count data run through bmm() for every sigmoid", {
  for (sigmoid in psy_sigmoids) {
    expect_s3_class(
      bmm(bmf(midpoint ~ 1, width ~ 1),
          psy_count_data(),
          psychometric("n_correct", "contrast", "n_trials", sigmoid = sigmoid,
                       guess = 0.5, lapse = 0),
          backend = "mock", mock_fit = 1, rename = FALSE),
      "bmmfit"
    )
  }
  expect_s3_class(
    bmm(bmf(midpoint ~ 1 + (1 | id), width ~ 1, guess ~ 1, lapse ~ 1),
        psy_trial_data(), psychometric("y", "x"),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "bmmfit"
  )
})


############################################################################# !
# PSYCHOMETRIC_THRESHOLD                                                 ####
############################################################################# !

# a bmmfit stand-in whose posterior_linpred() returns constant draws, given on
# the link scale (see mock_linpred_factory in helper-sdt-analysis.R)
fake_psychometric_fit <- function(model, data, formula) {
  structure(
    list(data = data, bmm = list(model = model, user_formula = formula)),
    class = c("bmmfit", "brmsfit")
  )
}

summary_value <- function(thr, measure, p = NA, column = "mean") {
  smry <- attr(thr, "summary")
  rows <- smry$measure == measure & (is.na(p) & is.na(smry$p) | smry$p %in% p)
  smry[rows, column]
}

test_that("thresholds and slopes match their closed form", {
  zspan <- .psychometric_z_constants("normal")[["zspan"]]
  fit <- fake_psychometric_fit(
    psychometric("y", "x", guess = 0.5), psy_trial_data(),
    bmf(midpoint ~ 1, width ~ 1, lapse ~ 1)
  )
  # width = zspan makes F the standard normal of x - midpoint
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(midpoint = 1, width = log(zspan),
                                                  lapse = qlogis(0.1))),
    .package = "brms"
  )
  thr <- psychometric_threshold(fit, p = c(0.6, 0.9))
  gain <- 0.5 * 0.9
  for (p in c(0.6, 0.9)) {
    f <- (p - 0.5) / gain
    expect_equal(summary_value(thr, "threshold", p), 1 + qnorm(f))
    expect_equal(summary_value(thr, "slope", p), gain * dnorm(qnorm(f)))
  }
  expect_equal(summary_value(thr, "guess"), 0.5)
  expect_equal(summary_value(thr, "lapse"), 0.1)
  expect_equal(summary_value(thr, "lapse_wh"), 0.05)
  expect_equal(nrow(thr), n_draws_mock * 7)
})

test_that("Weibull thresholds invert the classical Weibull and slopes match a finite difference", {
  fit <- fake_psychometric_fit(
    psychometric("n_correct", "contrast", "n_trials", sigmoid = "weibull",
                 guess = 0.5, lapse = 0),
    psy_count_data(), bmf(midpoint ~ 1, width ~ 1)
  )
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(midpoint = log(0.04), width = log(2))),
    .package = "brms"
  )
  thr <- psychometric_threshold(fit, p = 0.8)
  k <- .psychometric_z_constants("weibull")
  beta <- k[["zspan"]] / 2
  alpha <- exp(log(0.04) - k[["zmid"]] / beta)
  x80 <- alpha * (-log(1 - 0.6))^(1 / beta)
  expect_equal(summary_value(thr, "threshold", 0.8), x80)

  psi <- function(x) dpsychometric(1, x, 0.04, 2, guess = 0.5, sigmoid = "weibull")
  h <- 1e-7
  expect_equal(summary_value(thr, "slope", 0.8), (psi(x80 + h) - psi(x80 - h)) / (2 * h),
               tolerance = 1e-6)
})

test_that("a target outside the asymptotes gives undefined draws and a warning", {
  fit <- fake_psychometric_fit(
    psychometric("y", "x", guess = 0.5), psy_trial_data(),
    bmf(midpoint ~ 1, width ~ 1, lapse ~ 1)
  )
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(midpoint = 0, width = 0, lapse = -3)),
    .package = "brms"
  )
  expect_warning(thr <- psychometric_threshold(fit, p = c(0.4, 0.75)),
                 "50% of the threshold draws are undefined")
  expect_equal(summary_value(thr, "threshold", 0.4, "undefined"), 1)
  expect_equal(summary_value(thr, "threshold", 0.75, "undefined"), 0)
  expect_true(is.nan(summary_value(thr, "threshold", 0.4)))
})

test_that("a per-row guess column is a condition of its own", {
  fit <- fake_psychometric_fit(
    psychometric("n_correct", "contrast", "n_trials", guess = "chance", lapse = 0),
    psy_count_data(), bmf(midpoint ~ 1, width ~ 1)
  )
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(midpoint = 0.05, width = log(0.05))),
    .package = "brms"
  )
  thr <- psychometric_threshold(fit, p = 0.75)
  smry <- attr(thr, "summary")
  expect_setequal(smry$chance, c(0.5, 0.25))
  # 75% correct is F = 0.5 with two alternatives, and F = 2/3 with four
  thresholds <- smry[smry$measure == "threshold", ]
  expect_equal(thresholds$mean[thresholds$chance == 0.5], 0.05)
  expect_gt(thresholds$mean[thresholds$chance == 0.25], 0.05)
})

test_that("user conditions must name every column the model reads", {
  fit <- fake_psychometric_fit(
    psychometric("n_correct", "contrast", "n_trials", guess = "chance", lapse = 0),
    psy_count_data(), bmf(midpoint ~ 1, width ~ 1)
  )
  local_mocked_bindings(
    posterior_linpred = mock_linpred_factory(list(midpoint = 0.05, width = log(0.05))),
    .package = "brms"
  )
  expect_error(psychometric_threshold(fit, conditions = data.frame(n_trials = 40L)),
               "missing: 'chance'")
  thr <- psychometric_threshold(fit, conditions = data.frame(chance = 0.5))
  expect_equal(summary_value(thr, "guess"), 0.5)
})

test_that("psychometric_threshold() validates its arguments", {
  fit <- fake_psychometric_fit(psychometric("y", "x"), psy_trial_data(),
                               bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1))
  expect_error(psychometric_threshold(fit, p = 1), "strictly between 0 and 1")
  expect_error(psychometric_threshold(fit, p = NA), "strictly between 0 and 1")
  expect_error(psychometric_threshold(fit, ndraws = 10), "draw_ids")
  expect_error(psychometric_threshold(fake_mafc_fit()), "only available for psychometric")
})

test_that("summary() of a psychometric fit points to the Wichmann-Hill lapse rate", {
  note <- summary_notes(psychometric("y", "x"), NULL)
  expect_match(note, "lapse * (1 - guess)", fixed = TRUE)
  expect_match(note, "psychometric_threshold()", fixed = TRUE)
})


############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT AND POSTERIOR_EPRED                         ####
############################################################################# !

# a brmsprep stand-in: dpars as draws x observations matrices on the natural
# scale, data as brms passes the addition terms
fake_psychometric_prep <- function(dpars, data) {
  n_draws <- nrow(dpars[[1]])
  structure(list(dpars = dpars, data = data, ndraws = n_draws, nobs = length(data$Y)),
            class = "brmsprep")
}

test_that("log_lik, posterior_predict and posterior_epred agree with the distribution functions", {
  withr::local_seed(99)
  n_draws <- 6
  x <- c(0.01, 0.04, 0.1)
  draws <- function(lo, hi) matrix(stats::runif(n_draws * 3, lo, hi), n_draws, 3)
  dpars <- list(midpoint = draws(0.02, 0.06), width = draws(1, 3), lapse = draws(0.01, 0.1))
  data <- list(Y = c(10L, 25L, 38L), trials = rep(40L, 3), vreal1 = x,
               vreal2 = c(0.5, 0.25, 0.5), vreal3 = rep(0, 3),
               vint1 = rep(.sdt_dist_id("gumbel_min"), 3), vint2 = rep(1L, 3))
  prep <- fake_psychometric_prep(dpars, data)

  for (i in 1:3) {
    expected <- dpsychometric(data$Y[i], x[i], dpars$midpoint[, i], dpars$width[, i],
                              guess = data$vreal2[i], lapse = dpars$lapse[, i],
                              n_trials = 40, sigmoid = "weibull", log = TRUE)
    expect_equal(log_lik_psychometric(i, prep), expected)
    pp <- posterior_predict_psychometric(i, prep)
    expect_length(pp, n_draws)
    expect_true(all(pp >= 0 & pp <= 40))
  }

  epred <- posterior_epred_psychometric(prep)
  expect_equal(dim(epred), c(n_draws, 3))
  psi <- vapply(1:3, function(i) {
    dpsychometric(1, x[i], dpars$midpoint[, i], dpars$width[, i],
                  guess = data$vreal2[i], lapse = dpars$lapse[, i], sigmoid = "weibull")
  }, numeric(n_draws))
  expect_equal(epred, 40 * psi)
})

test_that("each sigmoid is recovered from its Stan ids", {
  for (sigmoid in psy_sigmoids) {
    spec <- .psychometric_sigmoids[[sigmoid]]
    expect_equal(.psychometric_sigmoid_name(.sdt_dist_id(spec$dist), as.integer(spec$log_x)),
                 sigmoid)
  }
})
