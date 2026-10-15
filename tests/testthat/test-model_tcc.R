# Target confusability competition model: constructor, check_data(), dtcc() and
# rtcc(), the generated Stan code, and the Stan table against R

# exponential stand-in for a measured colour similarity, f(0) = 1, f(180 deg) = 0
tcc_test_similarity <- function(n_options = 360) {
  delta <- seq(0, 180, length.out = n_options / 2 + 1)
  (exp(-delta / 20) - exp(-9)) / (1 - exp(-9))
}

tcc_test_model <- function(n_options = 360, ...) {
  tcc("err", similarity = tcc_test_similarity(n_options), n_options = n_options, ...)
}

# log P(option at distance k) by adaptive quadrature over all n options of the
# wheel, written without the package's quadrature or its distance bookkeeping
ref_tcc_logp <- function(d, k, similarity) {
  n <- 2 * (length(similarity) - 1)
  offsets <- 0:(n - 1)
  f_all <- similarity[pmin(offsets, n - offsets) + 1]
  own <- k + 1
  logint <- function(x) {
    dnorm(x - d * f_all[own], log = TRUE) +
      rowSums(pnorm(outer(x, d * f_all[-own], "-"), log.p = TRUE))
  }
  # every integrand lies inside this range; split at its peak, which at large d
  # is too narrow for integrate() to find on an infinite range
  grid <- seq(min(0, d) - 12, max(0, d) + 12, by = 0.01)
  values <- logint(grid)
  peak <- max(values)
  at <- grid[which.max(values)]
  part <- function(lo, hi) {
    integrate(function(x) exp(logint(x) - peak), lo, hi,
              rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
  }
  peak + log(part(min(grid), at) + part(at, max(grid)))
}

############################################################################# !
# CONSTRUCTOR                                                            ####
############################################################################# !

test_that("tcc() builds an independent model with d on the log link", {
  model <- tcc_test_model()
  expect_s3_class(model, "tcc")
  expect_false(inherits(model, "circular"))
  expect_named(model$parameters, "d")
  expect_equal(model$links$d, "log")
  expect_equal(model$other_vars$n_options, 360)
  expect_equal(model$other_vars$units, "radians")
  expect_equal(model$other_vars$d_max, 30)
  expect_equal(tcc_test_model(180, units = "bins")$other_vars$units, "bins")
})

test_that("tcc() refuses a similarity vector the model cannot read", {
  f <- tcc_test_similarity(360)
  expect_error(tcc("err", similarity = f[-1]), "181 values")
  expect_error(tcc("err", similarity = f, n_options = 180), "91 values")
  expect_error(tcc("err", similarity = c(0.9, f[-1])), "1 at distance 0")
  expect_error(tcc("err", similarity = c(1, f[-c(1, 2)], 1.2)), "between 0 and 1")
  expect_error(tcc("err", similarity = replace(f, 5, NA)), "between 0 and 1")
  expect_error(tcc("err", similarity = as.character(f)), "numeric")
  expect_error(tcc("err"), "missing.*similarity")
})

test_that("tcc() refuses odd or tiny option counts and a non-positive d_max", {
  f <- tcc_test_similarity(360)
  expect_error(tcc("err", similarity = f[1:3], n_options = 5), "even number")
  expect_error(tcc("err", similarity = 1, n_options = 0), "even number")
  expect_error(tcc("err", similarity = f, d_max = 0), "positive number")
  expect_error(tcc("err", similarity = f, d_max = c(6, 12)), "positive number")
  expect_error(tcc("err", similarity = f, units = "degrees"), "should be one of")
})

############################################################################# !
# CHECK_DATA                                                             ####
############################################################################# !

test_that("check_data() bins radians to the nearest option, with both edges at +-pi in one bin", {
  w <- 2 * pi / 360
  dat <- data.frame(err = c(0, pi, -pi, 10 * w, -10 * w, 0.4999 * w, 0.5001 * w,
                            -0.5001 * w, 3 * pi / 2, NA))
  out <- check_data(tcc_test_model(), dat, bmf(d ~ 1))
  expect_identical(out$tcc_bin, c(0L, 180L, 180L, 10L, -10L, 0L, 1L, -1L, -90L, NA))

  w180 <- 2 * pi / 180
  out <- check_data(tcc_test_model(180), data.frame(err = c(pi, -pi, 0.51 * w180)), bmf(d ~ 1))
  expect_identical(out$tcc_bin, c(90L, 90L, 1L))
})

test_that("check_data() reduces integer bins modulo n_options", {
  dat <- data.frame(err = c(0, 1, -1, 180, -180, 359, 360, 541))
  out <- check_data(tcc_test_model(units = "bins"), dat, bmf(d ~ 1))
  expect_identical(out$tcc_bin, c(0L, 1L, -1L, 180L, 180L, -1L, 0L, -179L))
})

test_that("check_data() refuses responses it cannot bin", {
  expect_error(check_data(tcc_test_model(), data.frame(x = 1), bmf(d ~ 1)),
               "'err' is not present")
  expect_error(check_data(tcc_test_model(), data.frame(err = c(0, 45)), bmf(d ~ 1)),
               "radians")
  expect_error(check_data(tcc_test_model(units = "bins"), data.frame(err = c(1, 2.5)), bmf(d ~ 1)),
               "whole numbers")
  expect_error(check_data(tcc_test_model(), data.frame(err = c("a", "b")), bmf(d ~ 1)),
               "must be numeric")
  expect_warning(
    check_data(tcc_test_model(), data.frame(err = 0, tcc_bin = 3), bmf(d ~ 1)),
    "reserved"
  )
})

############################################################################# !
# dtcc() AND rtcc()                                                      ####
############################################################################# !

test_that("dtcc() matches adaptive quadrature where the error peaks", {
  # the largest d of the default table and the least likely options, plus
  # negative d (repulsion) at the low end of the table
  cells <- rbind(
    expand.grid(n = 360, d = c(-1.5, 0, 0.37, 12, 29.97, 30), k = c(0, 1, 179, 180)),
    expand.grid(n = 180, d = c(0.37, 30), k = c(0, 1, 89, 90))
  )
  cells$sharp <- FALSE
  # a far option under a 1-degree similarity peaks halfway between its own
  # mean and the target's, below the window's former lower end
  cells <- rbind(cells, data.frame(n = 360, d = c(12, 30), k = rep(c(0, 1, 5, 90, 180), each = 2),
                                   sharp = TRUE))
  for (i in seq_len(nrow(cells))) {
    cl <- cells[i, ]
    f <- tcc_test_similarity(cl$n)
    if (cl$sharp) f <- exp(-(0:180)) - exp(-180)
    x <- 2 * pi * cl$k / cl$n
    # absolute: testthat's tolerance is relative; at d = 30 the farthest option
    # (log P about -357) is off by 3.4e-7, a relative error in P
    logp <- dtcc(x, d = cl$d, similarity = f, n_options = cl$n, log = TRUE)
    error <- logp - ref_tcc_logp(cl$d, cl$k, f)
    expect_lt(abs(error), if (logp > -100) 1e-7 else 1e-6,
              label = paste("n =", cl$n, "d =", cl$d, "k =", cl$k, "sharp =", cl$sharp))
  }
})

test_that("dtcc() is a probability mass function over the options", {
  f <- tcc_test_similarity(360)
  x <- 2 * pi * (-179:180) / 360
  for (d in c(-2, 0, 1.3, 12, 20)) {
    expect_equal(sum(dtcc(x, d = d, similarity = f)), 1, tolerance = 1e-7, info = paste("d =", d))
  }
  expect_equal(dtcc(x, d = 0, similarity = f), rep(1 / 360, 360))
  # symmetric, vectorised over x and d, and the same anywhere inside a bin
  expect_equal(dtcc(-x, d = 2, similarity = f), dtcc(x, d = 2, similarity = f))
  expect_equal(dtcc(c(0.1, 0.1), d = c(1, 3), similarity = f),
               c(dtcc(0.1, 1, f), dtcc(0.1, 3, f)))
  expect_equal(dtcc(0.4 * 2 * pi / 360, 2, f), dtcc(0, 2, f))
})

test_that("dtcc() and rtcc() keep the length of empty input and return NA for a missing d", {
  f <- tcc_test_similarity(360)
  expect_length(dtcc(numeric(0), 2, f), 0)
  expect_equal(is.na(dtcc(c(0, 0.1), c(2, NA), f)), c(FALSE, TRUE))
  withr::local_seed(566)
  expect_equal(is.na(rtcc(3, c(1, NA, 2), f)), c(FALSE, TRUE, FALSE))
})

test_that("rtcc() matches the argmax of independent Gaussian memory-match signals", {
  withr::local_seed(566)
  n <- 20000
  d <- 1.6
  f <- tcc_test_similarity(360)
  offsets <- -179:180
  means <- d * f[abs(offsets) + 1]
  signals <- matrix(stats::rnorm(n * 360), n) + rep(means, each = n)
  direct <- abs(offsets[max.col(signals, ties.method = "first")])

  sim <- rtcc(n, d = d, similarity = f)
  expect_true(all(sim > -pi & sim <= pi))
  sim_dist <- abs(round(sim * 360 / (2 * pi)))
  # coarse distance classes, compared by a two-sample chi-square test
  breaks <- c(-1, 0, 2, 5, 10, 20, 45, 90, 180)
  counts <- rbind(table(cut(direct, breaks)), table(cut(sim_dist, breaks)))
  expect_gt(stats::chisq.test(counts)$p.value, 0.001)
})

############################################################################# !
# CONFIGURE_MODEL AND STAN CODE                                          ####
############################################################################# !

tcc_mock_fit <- function(model = tcc_test_model(), formula = bmf(d ~ 1), ...) {
  dat <- data.frame(err = c(0, 0.1, -0.2, 1, 3, -2), cond = rep(c("a", "b"), 3))
  suppressMessages(bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE, ...))
}

# #566 design decision (2026-10-15): the similarity argument selects the
# similarity function, the formula never does
test_that("a formula for lambda or kappa is refused when the similarity is a vector", {
  expect_error(tcc_mock_fit(formula = bmf(d ~ 1, lambda ~ 1)),
               "Unrecognized model parameters: 'lambda'")
  expect_error(tcc_mock_fit(formula = bmf(d ~ cond, kappa ~ cond)),
               "Unrecognized model parameters: 'kappa'")
})

test_that("d reaches the likelihood through exp() of its linear predictor", {
  code <- brms::stancode(tcc_mock_fit(formula = bmf(d ~ cond)))
  expect_match(code, "d = exp(d);", fixed = TRUE)
  expect_match(code, "target += tcc_lpmf(Y | mu, d, tcc_tab, tcc_grid_lo, tcc_grid_step);",
               fixed = TRUE)
})

test_that("the Stan table is built in transformed data and read by the lpmf", {
  code <- brms::stancode(tcc_mock_fit(formula = bmf(d ~ cond)))
  expect_match(code, "matrix[tcc_grid_rows, rows(tcc_similarity)] tcc_tab;", fixed = TRUE)
  expect_match(code, "target += tcc_lpmf(Y | mu, d, tcc_tab, tcc_grid_lo, tcc_grid_step);",
               fixed = TRUE)
  sdata <- brms::standata(tcc_mock_fit())
  expect_equal(sdata$tcc_grid_lo, -2)
  expect_equal(sdata$tcc_grid_step, 0.02)
  expect_equal(sdata$tcc_grid_rows, 1603)
  expect_equal(sdata$tcc_similarity, tcc_test_similarity(360), ignore_attr = TRUE)
})

test_that("a table from R reaches the same lpmf as data and skips the build", {
  model <- tcc_test_model(d_max = 6)
  model$other_vars$table <- "r"
  fit <- tcc_mock_fit(model)
  code <- brms::stancode(fit)
  expect_match(code, "matrix[403, 181] tcc_tab;", fixed = TRUE)
  expect_false(grepl("tcc_tab[g] = tcc_logp_by_distance(", code, fixed = TRUE))
  expect_match(code, "target += tcc_lpmf(Y | mu, d, tcc_tab, tcc_grid_lo, tcc_grid_step);",
               fixed = TRUE)
  expect_equal(dim(brms::standata(fit)$tcc_tab), c(403, 181))
})

test_that("threading passes the table to the sliced likelihood as data", {
  code <- brms::stancode(tcc_mock_fit(threads = brms::threading(2)))
  expect_match(code, "data matrix tcc_tab", fixed = TRUE)
  expect_match(code, "tcc_lpmf(Y[start:end] | mu, d, tcc_tab, tcc_grid_lo, tcc_grid_step)",
               fixed = TRUE)
})

test_that("draws of d beyond the end of the table are reported after the fit", {
  model <- tcc_test_model(d_max = 6)
  # draws x observations; d of draws 3 and 4 passes d_max in one observation
  inside <- epred_prep(list(d = cbind(c(1, 5.9, 2, 3), c(2, 2, 2, 2))))
  beyond <- epred_prep(list(d = cbind(c(1, 5.9, 7.2, 2), c(2, 2, 2, 6.5))))
  expect_silent(.tcc_check_grid(model, inside))
  expect_warning(.tcc_check_grid(model, beyond),
                 "In 2 of 4 posterior draws, d passes.*-2 to 6.*Refit with d_max = 12")
})

test_that("a fit without draws is not checked against the table", {
  expect_silent(postprocess_brm(tcc_test_model(), tcc_mock_fit()))
})

############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT AND PP_CHECK                                ####
############################################################################# !

# the two response conventions the stored-frame cases of test-update.R cover:
# radians on 360 options with the table built in Stan, whole option steps on
# 180 options with the table passed from R
tcc_postprocess_cases <- function() {
  bins_model <- tcc("err", similarity = tcc_test_similarity(180), n_options = 180,
                    units = "bins", d_max = 6)
  bins_model$other_vars$table <- "r"
  list(
    radians = list(model = tcc_test_model(d_max = 6), n = 360),
    bins = list(model = bins_model, n = 180,
                data = data.frame(err = c(0, 3, -7, 90, 181, 12), cond = rep(c("a", "b"), 3)))
  )
}

tcc_case_prep <- function(case, d) {
  dat <- case$data %||% data.frame(err = c(0, 0.1, -0.2, 1, 3, -2), cond = rep(c("a", "b"), 3))
  fit <- suppressMessages(bmm(bmf(d ~ 1), dat, case$model, backend = "mock", mock_fit = 1,
                              rename = FALSE))
  sdata <- brms::standata(fit)
  epred_prep(list(d = matrix(d, length(d), sdata$N)), sdata)
}

test_that("log_lik_tcc() reads the table the fit used, for either response convention", {
  d <- c(-1.234, 0.011, 2.345, 5.987)
  for (case in tcc_postprocess_cases()) {
    prep <- tcc_case_prep(case, d)
    f <- tcc_test_similarity(case$n)
    for (i in seq_len(prep$nobs)) {
      exact <- dtcc(2 * pi * prep$data$Y[i] / case$n, d, f, n_options = case$n, log = TRUE)
      expect_lt(max(abs(log_lik_tcc(i, prep) - exact)), 1e-6)
    }
  }
})

test_that("posterior_predict_tcc() draws signed offsets from the model's pmf", {
  withr::local_seed(567)
  n <- 20000
  for (case in tcc_postprocess_cases()) {
    half <- case$n / 2
    sims <- posterior_predict_tcc(1, tcc_case_prep(case, rep(1.6, n)))
    expect_true(all(sims > -half & sims <= half & sims == round(sims)))
    p_dist <- dtcc(2 * pi * (0:half) / case$n, 1.6, tcc_test_similarity(case$n),
                   n_options = case$n) * c(1, rep(2, half - 1), 1)
    breaks <- c(-1, 0, 1, 2, 5, 10, 20, half / 2, half)
    expected <- tapply(p_dist, cut(0:half, breaks), sum)
    observed <- table(cut(abs(sims), breaks))
    expect_gt(stats::chisq.test(observed, p = expected / sum(expected))$p.value, 0.001)
    off_axis <- sims[sims != 0 & sims != half]
    expect_equal(mean(off_axis > 0), 0.5, tolerance = 0.05)
  }
})

test_that("pp_simulate.tcc() draws each observation from its own d", {
  withr::local_seed(568)
  case <- tcc_postprocess_cases()$radians
  prep <- tcc_case_prep(case, rep(0, 40))
  prep$dpars$d[, 1:3] <- 8
  prep$dpars$d[, 4:6] <- -1.5
  sims <- pp_simulate(case$model, prep)$error
  expect_identical(dim(sims), c(40L, 6L))
  expect_lt(max(colMeans(abs(sims[, 1:3]))), 3)
  expect_gt(min(colMeans(abs(sims[, 4:6]))), 60)
})

test_that("the pp_check views put the response error on the user's scale", {
  offsets <- c(0, 45, -90, 180)
  checks <- pp_observables(tcc_test_model())$checks
  expect_equal(checks$error$compute(list(error = offsets)), c(0, pi / 4, -pi / 2, pi))
  expect_equal(checks$abs_error$compute(list(error = offsets)), c(0, pi / 4, pi / 2, pi))
  expect_match(checks$error$label, "radians")
  checks <- pp_observables(tcc_postprocess_cases()$bins$model)$checks
  expect_equal(checks$error$compute(list(error = offsets)), offsets)
  expect_match(checks$abs_error$label, "options")
})

test_that("a tcc fit stores its post-processing functions and has no expected response", {
  family <- tcc_mock_fit()$formula$family
  expect_identical(family$log_lik, log_lik_tcc)
  expect_identical(family$posterior_predict, posterior_predict_tcc)
  expect_error(family$posterior_epred(NULL), "not defined for the tcc model")
})

# the kernels run in a fixed_param program over the stanvars bmm generates, as
# the other parity tests do (expose_functions() fails inside the suite); `decl`
# and `gq` add data and generated quantities to the generated blocks, and the
# draw comes back as a named numeric vector
tcc_stan_run <- function(stanvars, decl, gq, data) {
  block <- function(name) {
    paste(vapply(stanvars, function(s) if (s$block == name) s$scode else "", ""), collapse = "\n")
  }
  program <- paste0(
    "functions {\n", block("functions"), "\n}\n",
    "data {\n", block("data"), "\n", decl, "\n}\n",
    "transformed data {\n", block("tdata"), "\n}\n",
    "generated quantities {\n", gq, "\n}\n"
  )
  sdata <- lapply(stanvars, function(s) s$sdata)
  sdata <- sdata[!vapply(sdata, is.null, TRUE)]
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(sdata, data), fixed_param = TRUE, chains = 1, iter_sampling = 1,
    iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  list(draw = unlist(csv[1, ]), sdata = sdata)
}

tcc_stan_lp <- function(model, rows) {
  stanvars <- bmm:::.tcc_stanvars(model, table = "stan")
  r_table <- bmm:::.tcc_table(bmm:::.tcc_similarity(model), bmm:::.tcc_grid(model))
  run <- tcc_stan_run(
    stanvars,
    decl = paste("int<lower=1> N;", "array[N] int y;", "vector[N] d;",
                 "matrix[tcc_grid_rows, rows(tcc_similarity)] r_tab;", sep = "\n"),
    gq = paste0(
      "  real table_diff = max(abs(tcc_tab - r_tab));\n",
      "  real lp_runs = tcc_lpmf(y | rep_vector(0, N), d, tcc_tab, tcc_grid_lo, tcc_grid_step);\n",
      "  real lp_runs_r = tcc_lpmf(y | rep_vector(0, N), d, r_tab, tcc_grid_lo, tcc_grid_step);\n",
      "  vector[N] lp;\n",
      "  for (i in 1:N) lp[i] = tcc_lpmf({y[i]} | rep_vector(0, 1), rep_vector(d[i], 1), tcc_tab,\n",
      "                                  tcc_grid_lo, tcc_grid_step);\n"
    ),
    data = list(N = nrow(rows), y = rows$y, d = rows$d, r_tab = r_table)
  )
  # log_lik_tcc() on the same rows, the R twin of the Stan lookup
  prep <- epred_prep(list(d = matrix(rows$d, 1)), c(run$sdata, list(Y = rows$y)))
  list(table_diff = run$draw[["table_diff"]], lp_runs = run$draw[["lp_runs"]],
       lp_runs_r = run$draw[["lp_runs_r"]], lp = unname(run$draw[paste0("lp.", seq_len(nrow(rows)))]),
       log_lik = vapply(seq_len(nrow(rows)), log_lik_tcc, numeric(1), prep = prep))
}

test_that("the Stan likelihood matches dtcc() on and off the d grid, from either table", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_test_model(d_max = 6)
  f <- tcc_test_similarity(360)
  # runs of equal d (the per-cell path) and single rows (the per-trial path)
  rows <- data.frame(
    d = c(rep(-1.234, 3), rep(0.011, 2), 2, 2.345, rep(5.987, 4), 4.4, 6),
    y = c(0L, 7L, -180L, 3L, -3L, 0L, 45L, 0L, 1L, -2L, 180L, -90L, 179L)
  )
  out <- tcc_stan_lp(model, rows)
  expect_lt(out$table_diff, 1e-10)
  expect_equal(out$lp_runs, out$lp_runs_r, tolerance = 1e-10)
  expect_equal(out$lp_runs, sum(out$lp), tolerance = 1e-10)
  expect_lt(max(abs(out$lp - dtcc(2 * pi * rows$y / 360, rows$d, f, log = TRUE))), 1e-6)
  expect_lt(max(abs(out$log_lik - out$lp)), 1e-9)
})

test_that("the Stan table stays finite where Phi() underflows, past d = 63", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_test_model(d_max = 70)
  f <- tcc_test_similarity(360)
  rows <- data.frame(d = c(62.9, 63.5, 69.9), y = c(0L, 0L, 3L))
  out <- tcc_stan_lp(model, rows)
  expect_lt(out$table_diff, 1e-9)
  expect_lt(max(abs(out$lp - dtcc(2 * pi * rows$y / 360, rows$d, f, log = TRUE))), 1e-6)
})

test_that("beyond the d grid the Stan likelihood extends the table linearly", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_test_model(d_max = 6)
  f <- tcc_test_similarity(360)
  rows <- data.frame(d = c(-3, 7.5), y = c(0L, 0L))
  out <- tcc_stan_lp(model, rows)
  edge <- function(at, d) {
    h <- 0.02
    val <- function(x) dtcc(0, x, f, log = TRUE)
    val(at) + (d - at) * (val(at + h) - val(at - h)) / (2 * h)
  }
  expect_lt(max(abs(out$lp - c(edge(-2, -3), edge(6, 7.5)))), 1e-9)
  expect_lt(max(abs(out$log_lik - out$lp)), 1e-9)
})

############################################################################# !
# PARAMETRIC SIMILARITY                                                  ####
############################################################################# !

# The Laplace and von Mises similarities in their direct form at distances 0,
# 1, ..., n / 2 options, the distance in radians, scaled to 1 at the target and
# 0 opposite it. Oberauer (2021, eq. 2 and footnote 2) gives the two shapes,
# Schurgin et al. (2020) scale a similarity to [0, 1]
ref_similarity <- function(shape, par, n_options = 360) {
  delta <- 2 * pi * (0:(n_options / 2)) / n_options
  g <- if (shape == "laplace") exp(-delta / par) else exp(par * (cos(delta) - 1))
  (g - g[length(g)]) / (1 - g[length(g)])
}

tcc_shape_model <- function(shape = "laplace", ...) {
  tcc("err", similarity = shape, ...)
}

# the name of the similarity parameter is the argument of dtcc() and rtcc()
tcc_shape_args <- function(shape, par) {
  stats::setNames(list(par), if (shape == "laplace") "lambda" else "kappa")
}

test_that("a named similarity adds its parameter on the log link with a proper prior", {
  lap <- tcc_shape_model("laplace")
  expect_named(lap$parameters, c("d", "lambda"))
  expect_equal(lap$links, list(d = "log", lambda = "log"))
  expect_match(lap$default_priors$lambda$main, "^normal\\(")
  vm <- tcc_shape_model("vonmises")
  expect_named(vm$parameters, c("d", "kappa"))
  expect_equal(vm$links, list(d = "log", kappa = "log"))
  expect_match(vm$default_priors$kappa$main, "^normal\\(")
  expect_named(tcc_test_model()$parameters, "d")
  expect_error(tcc("err", similarity = "gaussian"), "'laplace', 'vonmises' or a numeric vector")
  expect_error(tcc("err", similarity = c("laplace", "vonmises")),
               "'laplace', 'vonmises' or a numeric vector")
})

# #566 design decision (2026-10-15): the similarity argument selects the
# similarity function, the formula never does
test_that("the similarity argument decides which similarity parameter a formula may name", {
  lap <- tcc_shape_model("laplace", d_max = 2)
  vm <- tcc_shape_model("vonmises", d_max = 2)
  fit <- tcc_mock_fit(lap, bmf(d ~ cond, lambda ~ cond))
  expect_match(brms::stancode(fit), "lambda = exp(lambda);", fixed = TRUE)
  expect_message(
    bmm(bmf(d ~ 1), data.frame(err = c(0, 1)), vm, backend = "mock", mock_fit = 1, rename = FALSE),
    "No formula for parameter kappa"
  )
  # a fixed value is on the log scale, as for every bmm parameter on a log link
  expect_match(brms::stancode(tcc_mock_fit(lap, bmf(d ~ 1, lambda = -1.25))),
               "Intercept_lambda = -1.25;", fixed = TRUE)
  expect_error(tcc_mock_fit(lap, bmf(d ~ 1, kappa ~ 1)), "Unrecognized model parameters: 'kappa'")
  expect_error(tcc_mock_fit(vm, bmf(d ~ 1, lambda ~ 1)), "Unrecognized model parameters: 'lambda'")
})

test_that("the parametric similarities are the Laplace and von Mises shapes scaled to [0, 1]", {
  delta <- 2 * pi * (0:180) / 360
  for (lambda in c(0.004, 0.05, 0.35, 3, 21)) {
    expect_lt(max(abs(.tcc_shape_similarity("laplace", lambda, delta) -
                        ref_similarity("laplace", lambda))), 1e-13, label = paste("lambda", lambda))
  }
  for (kappa in c(0.02, 0.5, 4, 100, 5e4)) {
    expect_lt(max(abs(.tcc_shape_similarity("vonmises", kappa, delta) -
                        ref_similarity("vonmises", kappa))), 1e-13, label = paste("kappa", kappa))
  }
})

test_that("dtcc() with a parametric similarity matches adaptive quadrature where the error peaks", {
  # both ends of each similarity grid and a value in between, at small d, at
  # the old and the new top of the d grid. Where P > 1e-4 the 120 nodes hold to
  # 2e-9 up to d = 12 and 6e-7 at d = 30; on responses the model calls
  # impossible, at the broad end, to 4e-6 and 6e-3
  cells <- rbind(
    expand.grid(shape = "laplace", par = c(0.004, 0.0189, 0.5, 21), d = c(0.37, 11.95, 29.95),
                k = c(0, 1, 90, 180), stringsAsFactors = FALSE),
    expand.grid(shape = "vonmises", par = c(0.02, 4.1, 5e4), d = c(11.95, 29.95),
                k = c(0, 1, 90, 180), stringsAsFactors = FALSE)
  )
  for (i in seq_len(nrow(cells))) {
    cl <- cells[i, ]
    exact <- do.call(dtcc, c(list(2 * pi * cl$k / 360, d = cl$d, similarity = cl$shape, log = TRUE),
                             tcc_shape_args(cl$shape, cl$par)))
    tolerance <- if (exact > log(1e-4)) 1e-6 else if (cl$d < 12) 1e-5 else 1e-2
    expect_lt(abs(exact - ref_tcc_logp(cl$d, cl$k, ref_similarity(cl$shape, cl$par))), tolerance,
              label = paste(cl$shape, cl$par, "d =", cl$d, "k =", cl$k))
  }
})

test_that("dtcc() takes the similarity parameter per observation and sums to one", {
  x <- 2 * pi * (-179:180) / 360
  for (lambda in c(0.004, 0.5, 21)) {
    expect_equal(sum(dtcc(x, 2, "laplace", lambda = lambda)), 1, tolerance = 1e-7,
                 info = paste("lambda", lambda))
  }
  expect_equal(sum(dtcc(x, 29, "vonmises", kappa = 5e4)), 1, tolerance = 1e-7)
  expect_equal(dtcc(c(0.1, 0.1), d = 2, similarity = "laplace", lambda = c(0.2, 0.9)),
               c(dtcc(0.1, 2, "laplace", lambda = 0.2), dtcc(0.1, 2, "laplace", lambda = 0.9)))
  expect_equal(dtcc(x, 2, "vonmises", kappa = 4, log = TRUE),
               dtcc(x, 2, ref_similarity("vonmises", 4), log = TRUE), tolerance = 1e-12)
  expect_equal(dtcc(x, 2, "laplace", lambda = 0.3, n_options = 180),
               dtcc(x, 2, ref_similarity("laplace", 0.3, 180), n_options = 180), tolerance = 1e-12)
})

test_that("dtcc() and rtcc() refuse a similarity parameter that does not fit the similarity", {
  expect_error(dtcc(0, 1, "laplace"), "`lambda`")
  expect_error(dtcc(0, 1, "laplace", kappa = 2), "`lambda`")
  expect_error(dtcc(0, 1, "vonmises", lambda = 2), "`kappa`")
  expect_error(dtcc(0, 1, "laplace", lambda = 0), "positive")
  expect_error(dtcc(0, 1, "laplace", lambda = c(1, NA)), "positive")
  expect_error(dtcc(0, 1, tcc_test_similarity(), lambda = 2), "numeric vector")
  expect_error(rtcc(5, 1, "vonmises"), "`kappa`")
  expect_error(dtcc(0, 1, "gaussian", lambda = 1), "'laplace', 'vonmises' or a numeric vector")
})

test_that("rtcc() with a parametric similarity matches the argmax of Gaussian signals", {
  withr::local_seed(569)
  n <- 20000
  offsets <- -179:180
  breaks <- c(-1, 0, 2, 5, 10, 20, 45, 90, 180)
  for (cl in list(list(shape = "laplace", par = 0.5), list(shape = "vonmises", par = 4))) {
    means <- 1.6 * ref_similarity(cl$shape, cl$par)[abs(offsets) + 1]
    signals <- matrix(stats::rnorm(n * 360), n) + rep(means, each = n)
    direct <- abs(offsets[max.col(signals, ties.method = "first")])
    sim <- do.call(rtcc, c(list(n, d = 1.6, similarity = cl$shape), tcc_shape_args(cl$shape, cl$par)))
    expect_true(all(sim > -pi & sim <= pi))
    counts <- rbind(table(cut(direct, breaks)), table(cut(abs(round(sim * 360 / (2 * pi))), breaks)))
    expect_gt(stats::chisq.test(counts)$p.value, 0.001, label = cl$shape)
  }
  # one parameter per observation: P(target) is high for a sharp similarity
  # and low for a broad one
  withr::local_seed(570)
  sims <- rtcc(4000, d = 8, similarity = "laplace", lambda = rep(c(0.01, 20), 2000))
  for (j in 1:2) {
    p <- dtcc(0, 8, "laplace", lambda = c(0.01, 20)[j])
    hits <- mean(sims[seq(j, 4000, 2)] == 0)
    expect_lt(abs(hits - p), 4 * sqrt(p * (1 - p) / 2000), label = paste("lambda", c(0.01, 20)[j]))
  }
})

test_that("the similarity parameter reaches a 2-D table in transformed data through exp()", {
  fit <- tcc_mock_fit(tcc_shape_model("laplace"), bmf(d ~ cond, lambda ~ cond))
  code <- brms::stancode(fit)
  expect_match(code, "lambda = exp(lambda);", fixed = TRUE)
  expect_match(code, paste("target += tcc_parametric_lpmf(Y | mu, d, lambda, tcc_tab2, tcc_grid_lo,",
                           "tcc_grid_step, tcc_u_lo, tcc_u_step);"), fixed = TRUE)
  expect_match(code, "array[tcc_u_rows] matrix[tcc_grid_rows, rows(tcc_delta)] tcc_tab2;",
               fixed = TRUE)
  sdata <- brms::standata(fit)
  expect_equal(c(sdata$tcc_grid_lo, sdata$tcc_grid_step, sdata$tcc_grid_rows), c(-2, 0.1, 324))
  expect_equal(c(sdata$tcc_u_lo, sdata$tcc_u_step, sdata$tcc_u_rows, sdata$tcc_shape),
               c(log(0.004), 0.1, 89, 1))
  expect_equal(sdata$tcc_delta, 2 * pi * (0:180) / 360, ignore_attr = TRUE)

  vm <- tcc_shape_model("vonmises", n_options = 180, units = "bins")
  sdata <- brms::standata(suppressMessages(bmm(
    bmf(d ~ 1, kappa ~ cond), data.frame(err = c(0, 3, -7, 90), cond = c("a", "b")), vm,
    backend = "mock", mock_fit = 1, rename = FALSE
  )))
  expect_equal(c(sdata$tcc_u_lo, sdata$tcc_u_rows, sdata$tcc_shape), c(log(0.02), 151, 2))
  expect_equal(sdata$tcc_delta, 2 * pi * (0:90) / 180, ignore_attr = TRUE)

  code <- brms::stancode(tcc_mock_fit(tcc_shape_model("laplace"), threads = brms::threading(2)))
  expect_match(code, "data array[] matrix tcc_tab2", fixed = TRUE)
})

test_that("draws of lambda or kappa beyond the similarity grid are reported after the fit", {
  model <- tcc_shape_model("laplace", d_max = 6)
  d <- matrix(2, 4, 2)
  inside <- epred_prep(list(d = d, lambda = cbind(c(0.005, 0.3, 21, 1), 1)))
  beyond <- epred_prep(list(d = d, lambda = cbind(c(0.003, 0.3, 1, 1), c(1, 1, 1, 22))))
  expect_silent(.tcc_check_grid(model, inside))
  expect_warning(.tcc_check_grid(model, beyond),
                 "In 2 of 4 posterior draws, lambda passes the range from 0.004 to 21.7")
  vm <- tcc_shape_model("vonmises", d_max = 6)
  expect_warning(.tcc_check_grid(vm, epred_prep(list(d = d, kappa = cbind(c(1, 2, 3, 9e4), 1)))),
                 "In 1 of 4 posterior draws, kappa passes the range from 0.02 to")
})

# a mock prep of the parametric family: one draw per (d, parameter) pair, one
# observation per distance k
tcc_shape_prep <- function(model, d, par, k) {
  fit <- suppressMessages(bmm(bmf(d ~ 1), data.frame(err = 2 * pi * k / 360), model,
                              backend = "mock", mock_fit = 1, rename = FALSE))
  sdata <- brms::standata(fit)
  dpars <- list(d = matrix(d, length(d), sdata$N))
  dpars[[setdiff(names(model$parameters), "d")]] <- matrix(par, length(par), sdata$N)
  epred_prep(dpars, sdata)
}

test_that("log_lik_tcc() reads the 2-D table within its accuracy where the error peaks", {
  # between grid nodes in both directions, where cubic interpolation errs
  # most, from sharp to broad similarity and up to the top of the d range: up
  # to 6e-5 where P > 1e-4, and on any response 5e-4 to d = 12 and 3e-3 at 30
  u_mid <- function(lo, j) exp(lo + (j + 0.5) * 0.1)
  cells <- list(
    laplace = expand.grid(d = c(0.35, 2.05, 11.95, 29.85), par = u_mid(log(0.004), c(5, 15, 48, 84))),
    vonmises = expand.grid(d = c(0.35, 2.05, 11.95, 29.85), par = u_mid(log(0.02), c(5, 60, 120, 145)))
  )
  k <- c(0, 1, 2, 5, 20, 90, 180)
  for (shape in names(cells)) {
    cl <- cells[[shape]]
    prep <- tcc_shape_prep(tcc_shape_model(shape), cl$d, cl$par, k)
    for (i in seq_along(k)) {
      exact <- vapply(seq_len(nrow(cl)), function(j) {
        do.call(dtcc, c(list(2 * pi * k[i] / 360, cl$d[j], shape, log = TRUE),
                        tcc_shape_args(shape, cl$par[j])))
      }, numeric(1))
      error <- abs(log_lik_tcc(i, prep) - exact)
      expect_lt(max(error[cl$d < 12]), 6e-4, label = paste(shape, "k =", k[i]))
      expect_lt(max(error), 4e-3, label = paste(shape, "k =", k[i]))
      expect_lt(max(error[exact > log(1e-4)], 0), 1e-4, label = paste(shape, "k =", k[i], "P > 1e-4"))
    }
  }
})

test_that("beyond the similarity grid log_lik_tcc() extends the table linearly", {
  model <- tcc_shape_model("laplace", d_max = 6)
  # d = 2 is a grid node, so only the lambda direction interpolates
  top <- log(0.004) + 86 * 0.1
  lo <- log(0.004)
  edge <- function(at, u) {
    val <- function(v) dtcc(0, 2, "laplace", lambda = exp(v), log = TRUE)
    val(at) + (u - at) * (val(at + 0.1) - val(at - 0.1)) / 0.2
  }
  prep <- tcc_shape_prep(model, c(2, 2), exp(c(lo - 0.7, top + 0.9)), 0)
  expect_lt(max(abs(log_lik_tcc(1, prep) - c(edge(lo, lo - 0.7), edge(top, top + 0.9)))), 1e-9)
})

# With a broad similarity and no far responses, d and kappa trade off along a
# ridge that leaves the d grid; there a rising linear extension made the
# posterior improper and the chains ran to d = 1e5 (S5 recovery, 2026-10-15)
test_that("above the d grid a response's log probability may fall but never rise", {
  model <- tcc_shape_model("vonmises", d_max = 6)
  # kappa on a node of its grid and the last d node at d_max plus one step, so
  # that only the extension in d is read
  kappa <- 0.02 * exp(2)
  top <- 6.1
  val <- function(d, k) dtcc(2 * pi * k / 360, d, "vonmises", kappa = kappa, log = TRUE)
  expect_gt(val(top + 0.1, 0), val(top - 0.1, 0))
  expect_lt(val(top + 0.1, 90), val(top - 0.1, 90))
  prep <- tcc_shape_prep(model, c(top, 8, 20), rep(kappa, 3), c(0, 90))
  expect_lt(max(abs(log_lik_tcc(1, prep) - val(top, 0))), 1e-9)
  slope <- (val(top + 0.1, 90) - val(top - 0.1, 90)) / 0.2
  expect_lt(max(abs(log_lik_tcc(2, prep) - (val(top, 90) + (c(top, 8, 20) - top) * slope))), 1e-9)
})

test_that("posterior_predict_tcc() draws from the 2-D table at each draw's d and lambda", {
  withr::local_seed(571)
  n <- 20000
  model <- tcc_shape_model("laplace", d_max = 6)
  prep <- tcc_shape_prep(model, rep(1.6, n), rep(0.5, n), 0)
  sims <- posterior_predict_tcc(1, prep)
  p_dist <- dtcc(2 * pi * (0:180) / 360, 1.6, "laplace", lambda = 0.5) * c(1, rep(2, 179), 1)
  breaks <- c(-1, 0, 1, 2, 5, 10, 20, 90, 180)
  expected <- tapply(p_dist, cut(0:180, breaks), sum)
  expect_gt(stats::chisq.test(table(cut(abs(sims), breaks)), p = expected / sum(expected))$p.value,
            0.001)
})

test_that("the Stan 2-D table and lookup match R on and off both grids", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_shape_model("laplace", d_max = 2)
  stanvars <- bmm:::.tcc_stanvars(model, table = "stan")
  # runs of equal (d, lambda) and single rows, inside and beyond both grids
  rows <- data.frame(
    d = c(rep(0.37, 3), rep(0.37, 2), 1.234, 1.234, -0.5, 2.6, 1.1, 1.1),
    lambda = c(rep(0.31, 3), rep(0.62, 2), 0.0189, 21.5, 2.2, 0.4, 0.002, 30),
    y = c(0L, 7L, -180L, 3L, -3L, 0L, 45L, 1L, -2L, 0L, 90L)
  )
  # the last d row pads the grid and holds the smaller of its own value and
  # that two rows down; at lambda 3.6 the target's log P rises with d there
  nodes <- data.frame(u = c(1, 2, 45, 88, 89, 70, 70), d = c(1, 43, 20, 2, 30, 44, 44),
                      col = c(1, 181, 2, 90, 1, 1, 90))
  run <- tcc_stan_run(
    stanvars,
    decl = paste("int<lower=1> N;", "array[N] int y;", "vector[N] d;", "vector[N] lambda;",
                 "int<lower=1> K;", "array[K] int node_u;", "array[K] int node_d;",
                 "array[K] int node_col;", "vector[2] sim_par;", sep = "\n"),
    gq = paste0(
      "  real lp_runs = tcc_parametric_lpmf(y | rep_vector(0, N), d, lambda, tcc_tab2, tcc_grid_lo,\n",
      "                                     tcc_grid_step, tcc_u_lo, tcc_u_step);\n",
      "  vector[N] lp;\n",
      "  for (i in 1:N) lp[i] = tcc_parametric_lpmf({y[i]} | rep_vector(0, 1), rep_vector(d[i], 1),\n",
      "                         rep_vector(lambda[i], 1), tcc_tab2, tcc_grid_lo, tcc_grid_step,\n",
      "                         tcc_u_lo, tcc_u_step);\n",
      "  vector[K] tab;\n",
      "  for (k in 1:K) tab[k] = tcc_tab2[node_u[k]][node_d[k], node_col[k]];\n",
      "  vector[rows(tcc_delta)] sim_laplace = tcc_similarity(tcc_delta, sim_par[1], 1);\n",
      "  vector[rows(tcc_delta)] sim_vonmises = tcc_similarity(tcc_delta, sim_par[2], 2);\n"
    ),
    data = list(N = nrow(rows), y = rows$y, d = rows$d, lambda = rows$lambda, K = nrow(nodes),
                node_u = nodes$u, node_d = nodes$d, node_col = nodes$col, sim_par = c(0.37, 6.3))
  )
  draw <- run$draw
  sdata <- run$sdata
  delta <- 2 * pi * (0:180) / 360
  expect_lt(max(abs(draw[paste0("sim_laplace.", 1:181)] - ref_similarity("laplace", 0.37))), 1e-13)
  expect_lt(max(abs(draw[paste0("sim_vonmises.", 1:181)] - ref_similarity("vonmises", 6.3))), 1e-13)
  r_nodes <- vapply(seq_len(nrow(nodes)), function(j) {
    lambda <- exp(sdata$tcc_u_lo + (nodes$u[j] - 2) * sdata$tcc_u_step)
    d <- sdata$tcc_grid_lo + (nodes$d[j] - 2) * sdata$tcc_grid_step
    at <- function(x) dtcc(2 * pi * (nodes$col[j] - 1) / 360, x, "laplace", lambda = lambda, log = TRUE)
    if (nodes$d[j] == sdata$tcc_grid_rows) min(at(d), at(d - 0.2)) else at(d)
  }, numeric(1))
  expect_lt(max(abs(draw[paste0("tab.", seq_len(nrow(nodes)))] - r_nodes)), 1e-10)
  lp <- unname(draw[paste0("lp.", seq_len(nrow(rows)))])
  expect_equal(draw[["lp_runs"]], sum(lp), tolerance = 1e-10)
  prep <- epred_prep(list(d = matrix(rows$d, 1), lambda = matrix(rows$lambda, 1)),
                     c(sdata, list(Y = rows$y)))
  expect_lt(max(abs(vapply(seq_len(nrow(rows)), log_lik_tcc, numeric(1), prep = prep) - lp)), 1e-9)
  inside <- rows$d <= 2 & rows$d >= -2 & rows$lambda >= 0.004 & rows$lambda <= 21.7
  exact <- vapply(which(inside), function(j) {
    dtcc(2 * pi * rows$y[j] / 360, rows$d[j], "laplace", lambda = rows$lambda[j], log = TRUE)
  }, numeric(1))
  expect_lt(max(abs(lp[inside] - exact)), 1e-3)
})
