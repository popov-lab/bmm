# Test EZ-CDM model specification and integration

ezcdm_test_model <- function(...) {
  ezcdm("mean_angle", "var_angle", "mean_rt", "var_rt", "n_trials", ...)
}

ezcdm_test_data <- function(n = 10, ...) {
  rezcdm(n, n_trials = 100, driftrate = 2, bound = 1.5, ndt = 0.3, ...)
}

ezcdm_test_formula <- bmf(driftrate ~ 1, bound ~ 1, ndt ~ 1)

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
      ndraws = nrow(dpars$driftrate),
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
  expect_error(ezcdm_test_model(version = "5par"), "should be one of")
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
  expect_named(model_3par$parameters, c("driftrate", "driftangle", "bound", "ndt"))
  expect_identical(model_3par$parameters, model_4par$parameters)
  expect_identical(model_3par$fixed_parameters, list(mu = 0, driftangle = 0))
  expect_identical(model_4par$fixed_parameters, list(mu = 0))
})

test_that("ezcdm versions share the same default links", {
  expected <- list(driftrate = "log", driftangle = "identity", bound = "log", ndt = "log")
  expect_identical(ezcdm_test_model(version = "3par")$links, expected)
  expect_identical(ezcdm_test_model(version = "4par")$links, expected)
})

test_that("ezcdm model accepts custom links", {
  model <- ezcdm_test_model(links = list(driftrate = "softplus"))
  expect_equal(model$links$driftrate, "softplus")
  expect_equal(model$links$bound, "log")
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
    check_data(ezcdm_test_model(), ezcdm_valid_data(mean_angle = c(10, -170)), ezcdm_test_formula),
    "degrees"
  )
})

test_that("ezcdm check_data rejects invalid summary statistics", {
  cases <- data.frame(
    column = c("var_angle", "var_angle", "mean_rt", "var_rt", "var_rt", "n_trials", "n_trials"),
    value = c(1.2, -0.1, -0.5, -0.02, 0, 2, 100.5),
    error = c(
      "between 0 and 1", "between 0 and 1", "Mean RT values must be positive",
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

test_that("ezcdm works with mock backend - 3par version", {
  skip_on_cran()
  withr::local_seed(1)
  expect_silent(
    bmm(ezcdm_test_formula, ezcdm_test_data(), ezcdm_test_model(), backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("ezcdm works with mock backend - 4par version", {
  skip_on_cran()
  withr::local_seed(2)
  formula <- bmf(driftrate ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  expect_silent(
    bmm(formula, ezcdm_test_data(driftangle = 0.3), ezcdm_test_model(version = "4par"),
      backend = "mock", mock_fit = 1, rename = FALSE
    )
  )
})

test_that("ezcdm uses one family for both versions", {
  skip_on_cran()
  withr::local_seed(3)
  sim_data <- ezcdm_test_data()

  fit_3par <- bmm(ezcdm_test_formula, sim_data, ezcdm_test_model(),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  fit_4par <- bmm(bmf(driftrate ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1), sim_data,
    ezcdm_test_model(version = "4par"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )

  for (fit in list(fit_3par, fit_4par)) {
    expect_s3_class(fit$formula, "brmsformula")
    expect_s3_class(fit$formula$family, "customfamily")
    expect_equal(fit$formula$family$name, "ezcdm")
    expect_equal(fit$formula$family$dpars, c("mu", "driftrate", "driftangle", "bound", "ndt"))
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

test_that("ezcdm 3par with a driftangle formula frees the drift angle", {
  skip_on_cran()
  withr::local_seed(14)
  formula <- bmf(driftrate ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  fit <- bmm(formula, ezcdm_test_data(driftangle = 0.3), ezcdm_test_model(),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
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
  data_a <- rezcdm(10, n_trials = 100, driftrate = 2.5, bound = 1.5, ndt = 0.3)
  data_a$condition <- "A"
  data_b <- rezcdm(10, n_trials = 100, driftrate = 1.5, bound = 1.5, ndt = 0.3)
  data_b$condition <- "B"
  sim_data <- rbind(data_a, data_b)
  sim_data$condition <- factor(sim_data$condition)

  formula <- bmf(driftrate ~ 0 + condition, bound ~ 1, ndt ~ 1)
  expect_silent(bmm(formula, sim_data, ezcdm_test_model(), backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("ezcdm with hierarchical structure works", {
  skip_on_cran()
  withr::local_seed(7)
  sim_data <- do.call(rbind, lapply(1:3, function(i) {
    d <- rezcdm(5, n_trials = 100, driftrate = rnorm(1, 2, 0.3), bound = 1.5, ndt = 0.3)
    d$id <- paste0("S", i)
    d
  }))
  sim_data$id <- factor(sim_data$id)

  formula <- bmf(driftrate ~ 1 + (1 | id), bound ~ 1 + (1 | id), ndt ~ 1)
  expect_silent(bmm(formula, sim_data, ezcdm_test_model(), backend = "mock", mock_fit = 1, rename = FALSE))
})

test_that("ezcdm default priors are correctly set", {
  skip_on_cran()
  withr::local_seed(8)
  formula <- bmf(driftrate ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  fit <- bmm(formula, ezcdm_test_data(), ezcdm_test_model(version = "4par"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  priors <- brms::prior_summary(fit)
  intercept_prior <- function(dpar) priors$prior[priors$dpar == dpar & priors$class == "Intercept"]

  expect_equal(intercept_prior("driftrate"), "normal(0.5, 1)")
  expect_equal(intercept_prior("bound"), "normal(0.5, 1)")
  expect_equal(intercept_prior("ndt"), "normal(-1.5, 0.5)")
  expect_equal(intercept_prior("driftangle"), "normal(0, 1)")
})

# create_initfun() builds its Stan code without the constant priors, so the mu
# intercept is still a sampled parameter there and needs an init range
test_that("ezcdm initial values cover every parameter create_initfun() sees", {
  skip_on_cran()
  withr::local_seed(9)
  sim_data <- ezcdm_test_data()
  formulas <- list(
    "3par" = ezcdm_test_formula,
    "4par" = bmf(driftrate ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  )
  for (version in names(formulas)) {
    model <- check_model(ezcdm_test_model(version = version), sim_data, formulas[[version]])
    data <- check_data(model, sim_data, formulas[[version]])
    formula <- check_formula(model, data, formulas[[version]])
    brms_formula <- configure_model(model, data, formula)$formula
    inits <- create_initfun(model, data, brms_formula)()
    expect_setequal(
      names(inits),
      c("Intercept", "Intercept_driftrate", "Intercept_driftangle", "Intercept_bound", "Intercept_ndt")
    )
    expect_true(all(vapply(inits, is.finite, logical(1))))
  }
})

test_that("ezcdm Stan code calls ezcdm_lpdf and parses", {
  skip_on_cran()
  withr::local_seed(10)
  sim_data <- ezcdm_test_data()
  sim_data$cond <- factor(rep(c("a", "b"), 5))

  code <- stancode(
    bmf(driftrate ~ cond, driftangle ~ 1, bound ~ 1, ndt ~ 1),
    data = sim_data,
    model = ezcdm_test_model(version = "4par")
  )
  expect_match(code, "real ezcdm_lpdf(", fixed = TRUE)
  expect_match(code, "ezcdm_lpdf(Y[n] |", fixed = TRUE)

  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
  expect_true(model$check_syntax(quiet = TRUE))
})

test_that("ezcdm 4par with a tan_half drift angle link wraps the angle in Stan", {
  skip_on_cran()
  withr::local_seed(11)
  model <- ezcdm_test_model(version = "4par", links = list(driftangle = "tan_half"))
  formula <- bmf(driftrate ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1)
  sim_data <- ezcdm_test_data(driftangle = 3)

  code <- stancode(formula, sim_data, model)
  expect_match(code, "driftangle = inv_tan_half(driftangle);", fixed = TRUE)
  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock_fit = 1, rename = FALSE))

  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  stan_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
  expect_true(stan_model$check_syntax(quiet = TRUE))
})

test_that("log_lik_ezcdm() maps the standata slots onto dezcdm()", {
  withr::local_seed(12)
  sim_data <- rezcdm(3, n_trials = c(20, 50, 200), driftrate = 2, driftangle = 0.2, bound = 1.5, ndt = 0.3)
  ndraws <- 4
  draws <- list(
    driftrate = matrix(runif(12, 1, 3), ndraws),
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
        driftrate = draws$driftrate[, i], driftangle = draws$driftangle[, i],
        bound = draws$bound[, i], ndt = draws$ndt[, i]
      )
    )
  }
})

test_that("posterior_predict_ezcdm() returns one circular mean per draw", {
  withr::local_seed(13)
  sim_data <- ezcdm_valid_data()
  ndraws <- 50
  prep <- ezcdm_prep(sim_data, list(
    driftrate = matrix(20, ndraws, 2),
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
# driftrate = 1e-300 the branch that avoids the overflow of bound / driftrate
ezcdm_parity_grid <- function() {
  grid <- expand.grid(
    kappa = c(1e-4, 5e-3, 0.2, 0.3, 1, 10, 100, 150, 800, 1e3, 1e4, 1e6),
    driftangle = c(0, -2.9),
    var_angle = c(0, 0.95),
    n_trials = c(3L, 500L),
    bound = c(0.7, 2.5)
  )
  grid$driftrate <- grid$kappa / grid$bound
  grid <- rbind(grid, data.frame(
    kappa = 1e-300, driftangle = 0, var_angle = 0.3, n_trials = 30L, bound = 1,
    driftrate = 1e-300
  ))
  grid$ndt <- 0.3
  grid$mean_angle <- 0.5
  moments <- .ezcdm_moments(grid$driftrate, grid$bound, grid$ndt)
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
  vector[N] driftrate;
  vector[N] driftangle;
  vector[N] bound;
  vector[N] ndt;
}
generated quantities {
  vector[N] stan_lpdf;
  for (i in 1:N) {
    stan_lpdf[i] = ezcdm_lpdf(mean_angle[i] | 0.0, driftrate[i], driftangle[i],
                              bound[i], ndt[i], var_angle[i], mean_rt[i],
                              var_rt[i], trials[i]);
  }
}
"
  )
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(
      list(N = nrow(grid), trials = grid$n_trials),
      as.list(grid[c("mean_angle", "var_angle", "mean_rt", "var_rt", "driftrate", "driftangle", "bound", "ndt")])
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
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan not installed")

  grid <- ezcdm_parity_grid()
  stan_lpdf <- ezcdm_parity_stan_lpdf(grid)
  r_lpdf <- with(grid, dezcdm(
    mean_angle, var_angle, mean_rt, var_rt, n_trials,
    driftrate = driftrate, driftangle = driftangle, bound = bound, ndt = ndt
  ))

  expect_true(all(is.finite(r_lpdf)))
  expect_true(all(is.finite(stan_lpdf)))
  relative_diff <- abs(stan_lpdf - r_lpdf) / abs(r_lpdf)
  expect_lt(max(relative_diff), 1e-9)
})
