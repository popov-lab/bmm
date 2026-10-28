# Test EZ-CDM model specification and integration

ezcdm_test_model <- function(...) {
  ezcdm("mean_angle", "var_angle", "mean_rt", "var_rt", "n_trials", ...)
}

ezcdm_test_data <- function(n = 10, ...) {
  rezcdm(n, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3, ...)
}

ezcdm_test_formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

# columns passed as NULL are dropped, as in utils::modifyList()
ezcdm_valid_data <- function(...) {
  as.data.frame(utils::modifyList(
    list(
      mean_angle = c(0.1, 0.2),
      var_angle = c(0.2, 0.3),
      mean_rt = c(0.5, 0.6),
      var_rt = c(0.02, 0.03),
      n_trials = c(100, 100)
    ),
    list(...)
  ))
}

# brms::prepare_predictions() needs posterior draws, which mock fits lack
ezcdm_prep <- function(data, dpars) {
  structure(
    list(
      ndraws = nrow(dpars$drift),
      nobs = nrow(data),
      dpars = dpars,
      data = list(
        Y = data$mean_angle, vreal1 = data$var_angle, vreal2 = data$mean_rt,
        vreal3 = data$var_rt, trials = data$n_trials
      )
    ),
    class = "brmsprep"
  )
}

test_that("ezcdm model can be created with both versions", {
  expect_silent(ezcdm_test_model(version = "3par"))
  expect_silent(ezcdm_test_model(version = "4par"))
})

test_that("ezcdm requires all response variables", {
  expect_error(ezcdm("mean_angle", "var_angle", "mean_rt", "var_rt"), "n_trials")
})

test_that("ezcdm model has correct class structure", {
  model_3par <- ezcdm_test_model()
  expect_s3_class(model_3par, "bmmodel")
  expect_s3_class(model_3par, "ezcdm")
  expect_s3_class(model_3par, "ezcdm_3par")
  expect_false(inherits(model_3par, "circular"))

  expect_s3_class(ezcdm_test_model(version = "4par"), "ezcdm_4par")
})

test_that("ezcdm 3par fixes the drift angle and 4par estimates it", {
  model_3par <- ezcdm_test_model(version = "3par")
  model_4par <- ezcdm_test_model(version = "4par")
  expect_named(model_3par$parameters, c("drift", "driftangle", "bound", "ndt"))
  expect_identical(model_3par$parameters, model_4par$parameters)
  expect_identical(model_3par$fixed_parameters, list(mu = 0, driftangle = 0))
  expect_identical(model_4par$fixed_parameters, list(mu = 0))
})

test_that("check_model() fits a 3par ezcdm with a driftangle formula as 4par", {
  formula <- bmf(drift ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  expect_message(
    model <- check_model(ezcdm_test_model(), ezcdm_valid_data(), formula),
    "predicts driftangle, which ezcdm version '3par' fixes to 0, so the model is fitted as version '4par'"
  )
  expect_identical(model$version, "4par")
  expect_identical(class(model), c("bmmodel", "ezcdm", "ezcdm_4par"))
  expect_identical(attr(model, "call")$version, "4par")
  expect_identical(model$fixed_parameters, list(mu = 0))

  # a constant drift angle keeps the version, as does the default formula
  model <- expect_silent(check_model(ezcdm_test_model(), ezcdm_valid_data(),
                                     bmf(drift ~ 1, bound ~ 1, ndt ~ 1, driftangle = 0.2)))
  expect_s3_class(model, "ezcdm_3par")
  expect_identical(model$fixed_parameters, list(mu = 0, driftangle = 0.2))
  expect_s3_class(expect_silent(check_model(ezcdm_test_model(), ezcdm_valid_data(), ezcdm_test_formula)),
                  "ezcdm_3par")
})

test_that("ezcdm versions share the same default links", {
  expected <- list(drift = "log", driftangle = "identity", bound = "log", ndt = "log")
  expect_identical(ezcdm_test_model(version = "3par")$links, expected)
  expect_identical(ezcdm_test_model(version = "4par")$links, expected)
})

test_that("ezcdm model accepts custom links", {
  model <- ezcdm_test_model(links = list(drift = "softplus"))
  expect_equal(model$links$drift, "softplus")
  expect_equal(model$links$bound, "log")
})

test_that("ezcdm validates its links like the other models", {
  expect_warning(
    model <- ezcdm_test_model(links = list(bond = "softplus")),
    "'bond' read as 'bound'"
  )
  expect_equal(model$links$bound, "softplus")
  expect_error(ezcdm_test_model(links = list(ndt = "logg")), "Unknown link function")
})

test_that("ezcdm check_data accepts valid data", {
  valid_data <- ezcdm_valid_data(
    mean_angle = c(0.1, -0.2, 3),
    var_angle = c(0, 0.3, 1),
    mean_rt = c(0.5, 0.6, 0.7),
    var_rt = c(0.02, 0.03, 0.025),
    n_trials = c(100, 100, 3)
  )
  expect_silent(check_data(ezcdm_test_model(), valid_data, ezcdm_test_formula))
})

test_that("ezcdm check_data warns when mean_angle looks like degrees", {
  expect_warning(
    check_data(ezcdm_test_model(), ezcdm_valid_data(mean_angle = c(7, 0.1)), ezcdm_test_formula),
    "mean_angle variable \\('mean_angle'\\) is in degrees"
  )
  expect_silent(check_data(ezcdm_test_model(), ezcdm_valid_data(mean_angle = c(6.2, 0.1)), ezcdm_test_formula))
})

test_that("ezcdm check_data wraps circular means coded on [0, 2 * pi)", {
  data <- check_data(ezcdm_test_model(), ezcdm_valid_data(mean_angle = c(4.5, 0.2)), ezcdm_test_formula)
  expect_equal(data$mean_angle[1], 4.5 - 2 * pi)
  # means already in range are kept bit for bit, so update() sees the same data
  expect_identical(data$mean_angle[2], 0.2)
  expect_identical(
    check_data(ezcdm_test_model(), data, ezcdm_test_formula)$mean_angle,
    data$mean_angle
  )
})

test_that("ezcdm check_data names the columns that are not numeric", {
  data <- ezcdm_valid_data(mean_rt = c("0.5", "0.6"), n_trials = factor(c(100, 100)))
  expect_error(
    check_data(ezcdm_test_model(), data, ezcdm_test_formula),
    "These variables must be numeric: 'mean_rt', 'n_trials'"
  )
})

test_that("ezcdm check_data rejects invalid summary statistics", {
  cases <- data.frame(
    column = c("var_angle", "var_angle", "mean_rt", "var_rt", "var_rt", "n_trials", "n_trials"),
    value = c(1.2, -0.1, -0.5, -0.02, 0, 2, 100.5),
    error = c(
      "\\('var_angle'\\) must be between 0 and 1", "between 0 and 1", "Mean RT values must be positive",
      "Variance of RT must be positive", "Variance of RT must be positive",
      "must be larger than two", "'n_trials'.*whole numbers"
    )
  )
  for (k in seq_len(nrow(cases))) {
    data <- ezcdm_valid_data()
    data[[cases$column[k]]][1] <- cases$value[k]
    expect_error(
      check_data(ezcdm_test_model(), data, ezcdm_test_formula),
      cases$error[k],
      info = glue("{cases$column[k]} = {cases$value[k]}")
    )
  }
})

test_that("ezcdm check_data warns when mean_rt looks like milliseconds", {
  expect_warning(
    check_data(ezcdm_test_model(), ezcdm_valid_data(mean_rt = c(500, 600)), ezcdm_test_formula),
    "milliseconds"
  )
})

test_that("ezcdm check_data rejects several variable names for one argument", {
  model <- ezcdm(c("ma1", "ma2"), "var_angle", "mean_rt", "var_rt", c("n1", "n2"))
  data <- ezcdm_valid_data(ma1 = c(0.1, 0.2), ma2 = c(0.1, 0.2), n1 = c(100, 100), n2 = c(100, 100))
  expect_error(
    check_data(model, data, ezcdm_test_formula),
    "single variable name: 'mean_angle', 'n_trials'"
  )
})

test_that("ezcdm check_data validates all required variables exist", {
  expect_error(
    check_data(ezcdm_test_model(), ezcdm_valid_data(var_rt = NULL), ezcdm_test_formula),
    "missing from the data: 'var_rt'"
  )
})

test_that("ezcdm check_data handles missing values", {
  data_with_na <- ezcdm_valid_data(
    mean_angle = c(0.1, NA, 0.3),
    var_angle = c(0.2, 0.3, NA),
    mean_rt = c(0.5, NA, 0.7),
    var_rt = c(0.02, 0.03, 0.025),
    n_trials = c(100, 100, 100)
  )
  expect_silent(check_data(ezcdm_test_model(), data_with_na, ezcdm_test_formula))
})

test_that("ezcdm uses one family for both versions", {
  skip_on_cran()
  withr::local_seed(3)
  sim_data <- ezcdm_test_data()

  fit_3par <- expect_silent(bmm(ezcdm_test_formula, sim_data, ezcdm_test_model(),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  fit_4par <- expect_silent(bmm(bmf(drift ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1), sim_data,
    ezcdm_test_model(version = "4par"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))

  for (fit in list(fit_3par, fit_4par)) {
    family <- fit$formula$family
    expect_s3_class(fit$formula, "brmsformula")
    expect_s3_class(family, "customfamily")
    expect_equal(family$name, "ezcdm")
    expect_equal(family$dpars, c("mu", "drift", "driftangle", "bound", "ndt"))
    expect_identical(family$log_lik, log_lik_ezcdm)
    expect_identical(family$posterior_predict, posterior_predict_ezcdm)
    expect_error(family$posterior_epred(NULL), "not defined for the ezcdm model")
  }
})

test_that("ezcdm formula conversion passes all summaries to the family", {
  bf <- bmf2bf(ezcdm_test_model(), ezcdm_test_formula)
  expect_equal(
    deparse1(bf$formula),
    "mean_angle | vreal(var_angle, mean_rt, var_rt) + trials(n_trials) ~ 1"
  )
})

test_that("ezcdm 3par fixes driftangle to exactly zero", {
  skip_on_cran()
  withr::local_seed(4)
  fit <- bmm(ezcdm_test_formula, ezcdm_test_data(), ezcdm_test_model(),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  priors <- brms::prior_summary(fit)
  expect_true(any(priors$dpar == "driftangle" & priors$prior == "constant(0)"))
})

test_that("ezcdm 3par with a driftangle formula is fitted and stored as 4par", {
  skip_on_cran()
  withr::local_seed(14)
  formula <- bmf(drift ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  expect_message(
    fit <- bmm(formula, ezcdm_test_data(driftangle = 0.3), ezcdm_test_model(),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "fitted as version '4par'"
  )
  expect_s3_class(fit$bmm$model, "ezcdm_4par")
  expect_identical(fit$bmm$model$version, "4par")
  expect_match(construct_model_call(fit$bmm$model), 'version = "4par"', fixed = TRUE)
  expect_false("driftangle" %in% names(fit$bmm$model$fixed_parameters))
  priors <- brms::prior_summary(fit)
  expect_false(any(priors$dpar == "driftangle" & priors$prior == "constant(0)"))
  expect_true(any(priors$dpar == "driftangle" & priors$class == "Intercept" & priors$prior == "normal(0, 1)"))
})

test_that("ezcdm 4par asks for a driftangle formula when it is missing", {
  skip_on_cran()
  withr::local_seed(5)
  expect_message(
    bmm(ezcdm_test_formula, ezcdm_test_data(), ezcdm_test_model(version = "4par"),
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "No formula for parameter driftangle"
  )
})

test_that("ezcdm with condition effects works", {
  skip_on_cran()
  withr::local_seed(6)
  data_a <- rezcdm(10, n_trials = 100, drift = 2.5, bound = 1.5, ndt = 0.3)
  data_a$condition <- "A"
  data_b <- rezcdm(10, n_trials = 100, drift = 1.5, bound = 1.5, ndt = 0.3)
  data_b$condition <- "B"
  sim_data <- rbind(data_a, data_b)
  sim_data$condition <- factor(sim_data$condition)

  formula <- bmf(drift ~ 0 + condition, bound ~ 1, ndt ~ 1)
  expect_silent(bmm(formula, sim_data, ezcdm_test_model(), backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("ezcdm with hierarchical structure works", {
  skip_on_cran()
  withr::local_seed(7)
  sim_data <- do.call(rbind, lapply(1:3, function(i) {
    d <- rezcdm(5, n_trials = 100, drift = rnorm(1, 2, 0.3), bound = 1.5, ndt = 0.3)
    d$id <- paste0("S", i)
    d
  }))
  sim_data$id <- factor(sim_data$id)

  formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1 + (1 | id), ndt ~ 1)
  expect_silent(bmm(formula, sim_data, ezcdm_test_model(), backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("ezcdm default priors are correctly set", {
  skip_on_cran()
  withr::local_seed(8)
  sim_data <- ezcdm_test_data()
  sim_data$cond <- factor(rep(c("a", "b"), 5))
  formula <- bmf(drift ~ cond, driftangle ~ cond, bound ~ cond, ndt ~ cond)
  fit <- bmm(formula, sim_data, ezcdm_test_model(version = "4par"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  priors <- brms::prior_summary(fit)
  # bmm sets the priors on the intercept and on the whole b class of each dpar
  dpar_prior <- function(dpar, class) {
    priors$prior[priors$dpar == dpar & priors$class == class & priors$coef == ""]
  }

  expect_equal(dpar_prior("drift", "Intercept"), "normal(0.5, 1)")
  expect_equal(dpar_prior("bound", "Intercept"), "normal(0.5, 1)")
  expect_equal(dpar_prior("ndt", "Intercept"), "normal(-1.5, 0.5)")
  expect_equal(dpar_prior("driftangle", "Intercept"), "normal(0, 1)")
  expect_equal(dpar_prior("drift", "b"), "normal(0, 0.5)")
  expect_equal(dpar_prior("bound", "b"), "normal(0, 0.5)")
  expect_equal(dpar_prior("ndt", "b"), "normal(0, 0.3)")
  expect_equal(dpar_prior("driftangle", "b"), "normal(0, 0.5)")
})

# the Stan code of a 4par model with a condition effect, and of one with the
# tan_half drift angle link
ezcdm_stancode_both <- function() {
  sim_data <- withr::with_seed(10, ezcdm_test_data(driftangle = 3))
  sim_data$cond <- factor(rep(c("a", "b"), 5))
  formula <- bmf(drift ~ cond, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  list(
    identity = stancode(formula, data = sim_data, model = ezcdm_test_model(version = "4par")),
    tan_half = stancode(
      formula, data = sim_data,
      model = ezcdm_test_model(version = "4par", links = list(driftangle = "tan_half"))
    )
  )
}

test_that("ezcdm Stan code calls ezcdm_lpdf", {
  skip_on_cran()
  code <- ezcdm_stancode_both()
  expect_match(code$identity, "real ezcdm_lpdf(", fixed = TRUE)
  expect_match(code$identity, "ezcdm_lpdf(Y[n] |", fixed = TRUE)
  expect_match(code$tan_half, "driftangle = inv_tan_half(driftangle);", fixed = TRUE)
  expect_false(grepl("inv_tan_half", code$identity, fixed = TRUE))
})

test_that("ezcdm Stan code parses", {
  skip_unless_stan_tests()
  skip_on_cran()
  for (code in ezcdm_stancode_both()) {
    model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
    expect_true(model$check_syntax(quiet = TRUE))
  }
})

test_that("ezcdm 4par with a tan_half drift angle link runs through the pipeline", {
  skip_on_cran()
  withr::local_seed(11)
  model <- ezcdm_test_model(version = "4par", links = list(driftangle = "tan_half"))
  formula <- bmf(drift ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  expect_silent(bmm(formula, ezcdm_test_data(driftangle = 3), model,
                    backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("log_lik_ezcdm() maps the standata slots onto dezcdm()", {
  withr::local_seed(12)
  sim_data <- rezcdm(3, n_trials = c(20, 50, 200), drift = 2, driftangle = 0.2, bound = 1.5, ndt = 0.3)
  ndraws <- 4
  draws <- list(
    drift = matrix(runif(12, 1, 3), ndraws),
    driftangle = matrix(runif(12, -0.3, 0.3), ndraws),
    bound = matrix(runif(12, 1, 2), ndraws),
    ndt = matrix(runif(12, 0.2, 0.4), ndraws)
  )
  prep <- ezcdm_prep(sim_data, draws)

  for (i in seq_len(nrow(sim_data))) {
    expect_equal(
      log_lik_ezcdm(i, prep),
      dezcdm(
        sim_data$mean_angle[i], sim_data$var_angle[i], sim_data$mean_rt[i],
        sim_data$var_rt[i], sim_data$n_trials[i],
        drift = draws$drift[, i], driftangle = draws$driftangle[, i],
        bound = draws$bound[, i], ndt = draws$ndt[, i], log = TRUE
      )
    )
  }
})

# posterior_epred() and fitted() refuse on this, before brms reads the family;
# the refusals on a fit are in test-helpers-postprocess.R
test_that("ezcdm has no expected response", {
  expect_false(expected_response_defined(ezcdm_test_model()))
  expect_false(expected_response_defined(ezcdm_test_model(version = "4par")))
})

test_that("posterior_predict_ezcdm() returns one circular mean per draw", {
  withr::local_seed(13)
  sim_data <- ezcdm_valid_data()
  ndraws <- 50
  prep <- ezcdm_prep(sim_data, list(
    drift = matrix(20, ndraws, 2),
    driftangle = matrix(1, ndraws, 2),
    bound = matrix(2, ndraws, 2),
    ndt = matrix(0.3, ndraws, 2)
  ))

  pred <- posterior_predict_ezcdm(1, prep)
  expect_length(pred, ndraws)
  expect_true(all(pred > -pi & pred <= pi))
  # kappa = 40 over 100 trials concentrates the circular mean at the drift angle
  expect_lt(max(abs(pred - 1)), 0.1)
})

# -----------------------------------------------------------------------------
# R <-> Stan parity of the likelihood
# -----------------------------------------------------------------------------

# RT summaries near the implied moments keep the lpdf moderate at every kappa;
# kappa 1e3 to 1e6 cross the asymptotic Bessel branch of .ezcdm_moments(), and
# drift = 1e-300 the branch that avoids the overflow of bound / drift
ezcdm_parity_grid <- function() {
  grid <- expand.grid(
    kappa = c(1e-4, 5e-3, 0.2, 0.3, 1, 10, 100, 150, 800, 1e3, 1e4, 1e6),
    driftangle = c(0, -2.9),
    var_angle = c(0, 0.95),
    n_trials = c(3L, 500L),
    bound = c(0.7, 2.5)
  )
  grid$drift <- grid$kappa / grid$bound
  grid <- rbind(grid, data.frame(
    kappa = 1e-300, driftangle = 0, var_angle = 0.3, n_trials = 30L, bound = 1,
    drift = 1e-300
  ))
  grid$ndt <- 0.3
  grid$mean_angle <- 0.5
  moments <- .ezcdm_moments(grid$drift, grid$bound, grid$ndt)
  grid$mean_rt <- moments$MRT * 1.02
  grid$var_rt <- moments$VRT * 0.9
  grid
}

ezcdm_parity_stan_lpdf <- function(grid) {
  program <- paste0(
    "functions {\n",
    read_lines2(file.path(system.file("stan_chunks", package = "bmm"), "ezcdm_functions.stan")),
    "\n}\n", "
data {
  int<lower=1> N;
  vector[N] mean_angle;
  vector[N] var_angle;
  vector[N] mean_rt;
  vector[N] var_rt;
  array[N] int trials;
  vector[N] drift;
  vector[N] driftangle;
  vector[N] bound;
  vector[N] ndt;
}
generated quantities {
  vector[N] stan_lpdf;
  for (i in 1:N) {
    stan_lpdf[i] = ezcdm_lpdf(mean_angle[i] | 0.0, drift[i], driftangle[i],
                              bound[i], ndt[i], var_angle[i], mean_rt[i],
                              var_rt[i], trials[i]);
  }
}
"
  )
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(
      list(N = nrow(grid), trials = grid$n_trials),
      as.list(grid[c("mean_angle", "var_angle", "mean_rt", "var_rt", "drift", "driftangle", "bound", "ndt")])
    ),
    # 17 significant figures round-trip a double; 18 switches to plain
    # decimals that misparse
    fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0,
    sig_figs = 17, refresh = 0, show_messages = FALSE
  )
  # base-R CSV read: cmdstanr's reader misparses long plain-decimal numbers
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  as.numeric(csv[1, paste0("stan_lpdf.", seq_len(nrow(grid)))])
}

test_that("ezcdm_lpdf in Stan matches dezcdm() in R", {
  skip_unless_stan_tests()
  skip_on_cran()

  grid <- ezcdm_parity_grid()
  stan_lpdf <- ezcdm_parity_stan_lpdf(grid)
  r_lpdf <- with(grid, dezcdm(
    mean_angle, var_angle, mean_rt, var_rt, n_trials,
    drift = drift, driftangle = driftangle, bound = bound, ndt = ndt, log = TRUE
  ))

  expect_true(all(is.finite(r_lpdf)))
  expect_true(all(is.finite(stan_lpdf)))
  relative_diff <- abs(stan_lpdf - r_lpdf) / abs(r_lpdf)
  expect_lt(max(relative_diff), 1e-9)
})

test_that("ezcdm_lpdf has the right gradient on both sides of every regime switch", {
  # The chunk switches between series, closed forms and asymptotic expansions
  # at kappa = 0.01, 0.25, 100 and 1000. A series that is a term short still
  # matches dezcdm() to 1e-9 in value, but its gradient is off by up to 3e-5,
  # which only autodiff shows. The references are 60-digit mpmath derivatives
  # of the closed forms: near kappa = 0.01 the closed forms cancel too much
  # for R finite differences. CmdStan's standalone log_prob method gives the
  # gradients, because expose_functions() returns none.
  skip_unless_stan_tests()
  skip_on_cran()

  program <- paste0(
    "functions {\n",
    read_lines2(file.path(system.file("stan_chunks", package = "bmm"), "ezcdm_functions.stan")),
    "\n}\n",
    "data { int N; vector[N] var_angle; vector[N] mean_rt; vector[N] var_rt; }\n",
    "parameters { vector[N] drift; vector[N] bound; vector[N] ndt; vector[N] driftangle; }\n",
    "model { for (i in 1:N) target += ezcdm_lpdf(0.1 | 0.0, drift[i], driftangle[i], bound[i],\n",
    "  ndt[i], var_angle[i], mean_rt[i], var_rt[i], 60); }\n"
  )
  # bound = 1.2, so kappa = 0.00096, then 0.0099 and 0.01008, 0.249 and
  # 0.25104, 0.96, 99 and 101.04, 999 and 1001.04, and 4800; the summaries lie
  # near the implied moments
  points <- data.frame(
    drift = c(0.0008, 0.00825, 0.0084, 0.2075, 0.2092, 0.8, 82.5, 84.2, 832.5, 834.2, 4000),
    var_angle = c(
      0.999534, 0.995199, 0.995111, 0.880161, 0.879195, 0.580924,
      0.0349115, 0.0348121, 0.0304856, 0.0304846, 0.030101
    ),
    mean_rt = c(
      1.0506, 1.05059, 1.05059, 1.04491, 1.04482, 0.976497,
      0.323906, 0.323607, 0.310484, 0.310481, 0.309309
    ),
    var_rt = c(
      0.23328, 0.233272, 0.233272, 0.228535, 0.228458, 0.175633,
      1.90386e-6, 1.79123e-6, 1.86998e-9, 1.85857e-9, 1.68715e-11
    )
  )
  pars <- list(drift = points$drift, bound = 1.2, ndt = 0.3, driftangle = 0.4)
  reference <- cbind(
    drift = c(
      -3.985085047122548e-3, -4.1409283043692387e-2, -4.2113428251397891e-2,
      -1.0340796073967198, -1.0425289382784974, -3.6785899399477387,
      -9.9291545065453891e+1, -1.0205277071599188e+2, -4.7235423990908405e+3,
      -4.7418097610350061e+3, -9.9636260853584198e+4
    ),
    bound = c(
      1.2184261873900556e+1, 1.2183783457217963e+1, 1.2183884517724231e+1,
      1.2109472969641551e+1, 1.2109372936029e+1, 1.127488424300486e+1,
      4.1436250621114998e+3, 4.3378730801896028e+3, 1.3047310559010595e+6,
      1.3120687356795456e+6, 1.1565739876262507e+8
    ),
    ndt = c(
      1.8999729835289725e+1, 1.8999776913024024e+1, 1.8999877176769284e+1,
      1.924742663569219e+1, 1.9252730916479183e+1, 2.2755499628905862e+1,
      2.7983158253353205e+5, 2.9689877459860675e+5, 2.6230610199912953e+8,
      2.6391147220565566e+8, 2.8861504929228611e+10
    ),
    driftangle = c(
      -7.9322351791210123e-6, -8.4276275223556823e-4, -8.738140060141365e-4,
      -5.290977999285588e-1, -5.3773247836434803e-1, -7.1334965449041145,
      -1.6941067286209044e+3, -1.7291936757770677e+3, -1.717347508519773e+4,
      -1.7208561792982526e+4, -8.2547928841140459e+4
    )
  )

  dir <- withr::local_tempdir()
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program, dir = dir), quiet = TRUE)
  cmdstanr::write_stan_json(c(list(N = nrow(points)), as.list(points[-1])), file.path(dir, "data.json"))
  cmdstanr::write_stan_json(lapply(pars, rep_len, nrow(points)), file.path(dir, "pars.json"))
  # On Windows the executable loads tbb.dll at startup, and CmdStan's makefile
  # copies that dll next to it only when TBB was absent from PATH at link time.
  # cmdstanr has it on PATH whenever it builds or runs a model, so launching
  # the executable ourselves has to put it there too
  run_path <- if (.Platform$OS.type == "windows") {
    c(
      getFromNamespace("toolchain_PATH_env_var", "cmdstanr")(),
      getFromNamespace("tbb_path", "cmdstanr")()
    )
  }
  status <- withr::with_path(run_path, system2(model$exe_file(), c(
    "method=log_prob", paste0("constrained_params=", file.path(dir, "pars.json")), "jacobian=0",
    "data", paste0("file=", file.path(dir, "data.json")),
    "output", paste0("file=", file.path(dir, "out.csv")), "sig_figs=17"
  ), stdout = FALSE, stderr = file.path(dir, "stderr.txt")))
  expect_equal(status, 0, info = paste(
    readLines(file.path(dir, "stderr.txt"), warn = FALSE), collapse = "\n"
  ))
  stan_gradient <- matrix(
    as.numeric(utils::read.csv(file.path(dir, "out.csv"), comment.char = "#")[1, -1]),
    ncol = length(pars), dimnames = list(NULL, names(pars))
  )

  relative_error <- abs(stan_gradient / reference - 1)
  # Just above kappa = 0.01, R^2 - 1 + 2R/kappa cancels to 1e-4 of its terms,
  # which magnifies rounding in R and its derivative into a relative error
  # of about 5e-7 in the drift gradient, harmless for HMC
  above_switch <- which(points$drift == 0.0084)
  expect_lt(relative_error[above_switch, "drift"], 2e-6)
  relative_error[above_switch, "drift"] <- 0
  # max, not expect_equal(): its tolerance applies to the mean difference, so
  # one wrong entry among 44 would pass
  expect_lt(max(relative_error), 1e-8)
})
