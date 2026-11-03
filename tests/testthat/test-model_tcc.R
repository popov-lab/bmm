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
  grid <- seq(min(0, d) - 12, max(0, d) + 12, by = 0.05)
  peak <- max(logint(grid))
  res <- integrate(function(x) exp(logint(x) - peak), -Inf, Inf,
                   rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)
  peak + log(res$value)
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
  expect_error(tcc("err", similarity = "vonmises"), "tcc_similarity\\(\\).*sdm\\(\\)")
  expect_error(tcc("err"), "missing.*similarity")
})

test_that("tcc_similarity() scales the Laplace and von Mises shapes to 1 at the target and 0 opposite", {
  delta <- (0:180) * pi / 180
  laplace <- exp(-2.5 * delta)
  vonmises <- exp(4 * (cos(delta) - 1))
  expect_equal(tcc_similarity("laplace", kappa = 2.5),
               (laplace - laplace[181]) / (1 - laplace[181]), tolerance = 1e-14)
  expect_equal(tcc_similarity("vonmises", kappa = 4),
               (vonmises - vonmises[181]) / (1 - vonmises[181]), tolerance = 1e-14)
  expect_equal(tcc_similarity("laplace", kappa = 1 / deg2rad(20)), tcc_test_similarity(),
               tolerance = 1e-14)
  expect_length(tcc_similarity("vonmises", kappa = 4, n_options = 180), 91)
  # four options lie at 0, pi / 2 and pi from the target
  expect_equal(tcc_similarity("laplace", kappa = 1, n_options = 4),
               c(1, (exp(-pi / 2) - exp(-pi)) / (1 - exp(-pi)), 0), tolerance = 1e-14)
})

test_that("tcc_similarity() stays exact where the scaled shapes approach a line or a cosine", {
  delta <- (0:180) * pi / 180
  expect_equal(tcc_similarity("laplace", kappa = 1e-9), 1 - delta / pi, tolerance = 1e-8)
  expect_equal(tcc_similarity("vonmises", kappa = 1e-9), (1 + cos(delta)) / 2, tolerance = 1e-8)
  sharp <- tcc_similarity("vonmises", kappa = 1e4)
  expect_equal(sharp[c(1, 181)], c(1, 0))
  expect_true(all(is.finite(sharp)) && all(sharp >= 0))
})

test_that("a vector from tcc_similarity() is accepted by tcc()", {
  expect_s3_class(tcc("err", similarity = tcc_similarity("vonmises", kappa = 8)), "tcc")
})

test_that("tcc_similarity() refuses an unknown shape, a non-positive kappa and an odd option count", {
  expect_error(tcc_similarity("gaussian", kappa = 1), "should be one of")
  expect_error(tcc_similarity("laplace", kappa = 0), "single positive finite number")
  expect_error(tcc_similarity("laplace", kappa = c(1, 2)), "single positive finite number")
  expect_error(tcc_similarity("laplace", kappa = Inf), "single positive finite number")
  expect_error(tcc_similarity("laplace", kappa = 1, n_options = 361), "even number")
  expect_error(tcc_similarity("laplace", kappa = 1, n_options = Inf), "even number")
  expect_error(tcc_similarity("laplace"), "missing.*kappa")
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
  expect_equal(sdata$tcc_grid_lo, 0)
  expect_equal(sdata$tcc_grid_step, 0.02)
  expect_equal(sdata$tcc_grid_rows, 1504)
  # d can be negative only under a link that allows it
  identity_fit <- suppressWarnings(tcc_mock_fit(model = tcc_test_model(links = list(d = "identity"))))
  expect_equal(brms::standata(identity_fit)$tcc_grid_lo, -2)
  expect_equal(sdata$tcc_similarity, tcc_test_similarity(360), ignore_attr = TRUE)
})

test_that("a table from R reaches the same lpmf as data and skips the build", {
  model <- tcc_test_model(d_max = 6)
  model$other_vars$table <- "r"
  fit <- tcc_mock_fit(model)
  code <- brms::stancode(fit)
  expect_match(code, "matrix[304, 181] tcc_tab;", fixed = TRUE)
  expect_false(grepl("tcc_tab[g] = tcc_logp_by_distance(", code, fixed = TRUE))
  expect_match(code, "target += tcc_lpmf(Y | mu, d, tcc_tab, tcc_grid_lo, tcc_grid_step);",
               fixed = TRUE)
  expect_equal(dim(brms::standata(fit)$tcc_tab), c(304, 181))
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
                 "In 2 of 4 posterior draws, d passes.*0 to 6.*Refit with d_max = 12")
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
  d <- c(1.234, 0.011, 2.345, 5.987)
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
# the other parity tests do (expose_functions() fails inside the suite)
# `r_table` is the table passed from R (by default the quadrature table);
# lp_r reads it row by row, for the tables of correlated or motor noise
tcc_stan_lp <- function(model, rows, r_table = NULL) {
  stanvars <- bmm:::.tcc_stanvars(model, table = "stan")
  r_table <- r_table %||% bmm:::.tcc_table(bmm:::.tcc_similarity(model), bmm:::.tcc_grid(model))
  data_decl <- paste(vapply(stanvars, function(s) if (s$block == "data") s$scode else "", ""),
                     collapse = "\n")
  tdata <- paste(vapply(stanvars, function(s) if (s$block == "tdata") s$scode else "", ""),
                 collapse = "\n")
  funs <- paste(vapply(stanvars, function(s) if (s$block == "functions") s$scode else "", ""),
                collapse = "\n")
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n", data_decl, "\n  int<lower=1> N;\n  array[N] int y;\n  vector[N] d;\n",
    "  matrix[tcc_grid_rows, rows(tcc_similarity)] r_tab;\n}\n",
    "transformed data {\n", tdata, "\n}\n",
    "generated quantities {\n",
    "  real table_diff = max(abs(tcc_tab - r_tab));\n",
    "  real lp_runs = tcc_lpmf(y | rep_vector(0, N), d, tcc_tab, tcc_grid_lo, tcc_grid_step);\n",
    "  real lp_runs_r = tcc_lpmf(y | rep_vector(0, N), d, r_tab, tcc_grid_lo, tcc_grid_step);\n",
    "  vector[N] lp;\n",
    "  vector[N] lp_r;\n",
    "  for (i in 1:N) lp[i] = tcc_lpmf({y[i]} | rep_vector(0, 1), rep_vector(d[i], 1), tcc_tab,\n",
    "                                  tcc_grid_lo, tcc_grid_step);\n",
    "  for (i in 1:N) lp_r[i] = tcc_lpmf({y[i]} | rep_vector(0, 1), rep_vector(d[i], 1), r_tab,\n",
    "                                    tcc_grid_lo, tcc_grid_step);\n",
    "}\n"
  )
  sdata <- lapply(stanvars, function(s) s$sdata)
  sdata <- sdata[!vapply(sdata, is.null, TRUE)]
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(sdata, list(N = nrow(rows), y = rows$y, d = rows$d, r_tab = r_table)),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
    show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  # log_lik_tcc() on the same rows, the R twin of the Stan lookup
  prep <- epred_prep(list(d = matrix(rows$d, 1)), c(sdata, list(Y = rows$y)))
  list(table_diff = csv$table_diff, lp_runs = csv$lp_runs, lp_runs_r = csv$lp_runs_r,
       lp = as.numeric(csv[1, paste0("lp.", seq_len(nrow(rows)))]),
       lp_r = as.numeric(csv[1, paste0("lp_r.", seq_len(nrow(rows)))]),
       log_lik = vapply(seq_len(nrow(rows)), log_lik_tcc, numeric(1), prep = prep),
       log_lik_r = bmm:::.tcc_interp(rows$d, abs(rows$y) + 1,
                                     list(lo = sdata$tcc_grid_lo, step = sdata$tcc_grid_step,
                                          tab = r_table)))
}

test_that("the Stan likelihood matches dtcc() on and off the d grid, from either table", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_test_model(d_max = 6)
  f <- tcc_test_similarity(360)
  # runs of equal d (the per-cell path) and single rows (the per-trial path)
  rows <- data.frame(
    d = c(rep(1.234, 3), rep(0.011, 2), 2, 2.345, rep(5.987, 4), 4.4, 6),
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

test_that("beyond the d grid the Stan likelihood extends the table linearly, never rising above it", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_test_model(d_max = 6)
  f <- tcc_test_similarity(360)
  # the target's log P rises with d, the opposite option's falls
  rows <- data.frame(d = c(-0.5, 7.5, 7.5), y = c(0L, 0L, 180L))
  out <- tcc_stan_lp(model, rows)
  h <- 0.02
  val <- function(x, y) dtcc(2 * pi * y / 360, x, f, log = TRUE)
  edge <- function(at, d, y) val(at, y) + (d - at) * (val(at + h, y) - val(at - h, y)) / (2 * h)
  # the last node is d_max + h
  expected <- c(edge(0, -0.5, 0), val(6 + h, 0), edge(6 + h, 7.5, 180))
  expect_lt(max(abs(out$lp - expected)), 1e-9)
  expect_lt(max(abs(out$log_lik - out$lp)), 1e-9)
})

test_that("the Stan likelihood reads a correlated table from R as log_lik_tcc() does, beyond the grid too", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_test_model(d_max = 6, correlation = exp(-(0:180) / 10), correlation_draws = 6000,
                          motor_sd = deg2rad(2))
  r_table <- bmm:::.tcc_noise_table(model, bmm:::.tcc_grid(model))$tab
  rows <- data.frame(
    d = c(rep(1.234, 3), 0.011, 2.345, 5.987, 6.01, 7.5, 7.5, 9),
    y = c(0L, 7L, -180L, 3L, 45L, 1L, 0L, 0L, 180L, -90L)
  )
  out <- tcc_stan_lp(model, rows, r_table)
  expect_lt(max(abs(out$lp_r - out$log_lik_r)), 1e-9)
  expect_equal(out$lp_runs_r, sum(out$lp_r), tolerance = 1e-10)
  # the correlated table is not the independent-noise one the program built
  expect_gt(max(abs(out$lp_r - out$lp)), 0.1)
})

############################################################################# !
# GUMBEL NOISE                                                           ####
############################################################################# !

# log P of one option at distance k under Gumbel noise, by integrate() over the
# argmax of all n options of the wheel, written without the closed form:
# P(j) = int g(m - d f_j) prod_{o != j} G(m - d f_o) dm, with g and G the
# density and CDF of the standard Gumbel distribution
ref_tcc_gumbel_logp <- function(d, k, similarity) {
  n <- 2 * (length(similarity) - 1)
  offsets <- 0:(n - 1)
  means <- d * similarity[pmin(offsets, n - offsets) + 1]
  own <- k + 1
  logint <- function(m) {
    z <- m - means[own]
    -z - exp(-z) - rowSums(exp(-outer(m, means[-own], "-")))
  }
  grid <- seq(-10, max(means) + 30, by = 0.01)
  values <- logint(grid)
  peak <- max(values)
  at <- grid[which.max(values)]
  part <- function(lo, hi) {
    integrate(function(m) exp(logint(m) - peak), lo, hi,
              rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
  }
  peak + log(part(-Inf, at) + part(at, Inf))
}

tcc_gumbel_model <- function(similarity = tcc_test_similarity(360), ...) {
  tcc("err", similarity = similarity, noise = "gumbel", ...)
}

test_that("tcc(noise = 'gumbel') keeps the parameters, links and priors and refuses other noise", {
  gaussian <- tcc_test_model()
  gumbel <- tcc_gumbel_model()
  expect_s3_class(gumbel, "tcc")
  expect_equal(gumbel$other_vars$noise, "gumbel")
  expect_equal(gaussian$other_vars$noise, "gaussian")
  expect_named(gumbel$parameters, "d")
  expect_equal(gumbel$links, gaussian$links)
  expect_equal(gumbel$default_priors, gaussian$default_priors)
  expect_match(gumbel$parameters$d, "Gumbel")
  expect_error(tcc("err", similarity = tcc_test_similarity(360), noise = "logistic"),
               "should be one of")
})

test_that("each noise names its model and cites its authors", {
  gaussian <- tcc_test_model()
  gumbel <- tcc_gumbel_model()
  expect_equal(gaussian$name, "Target Confusability Competition (TCC) model")
  expect_equal(gumbel$name, "TCC with Gumbel noise")
  expect_length(gaussian$citation, 1)
  expect_match(gaussian$citation, "^Schurgin, M. W., Wixted, J. T., & Brady, T. F. \\(2020\\)")
  expect_length(gumbel$citation, 2)
  expect_equal(gumbel$citation[1], gaussian$citation)
  expect_match(gumbel$citation[2], "^Oberauer, K. \\(2023\\).*130\\(3\\), 841-852.*rev0000328$")
})

test_that("model_citation() keeps the Oberauer reference of a Gumbel model and fit", {
  expect_length(model_citation(tcc_test_model()), 1)
  expect_length(model_citation(tcc_gumbel_model()), 2)
  refs <- model_citation(tcc_mock_fit(tcc_gumbel_model()))
  expect_length(refs, 2)
  expect_match(refs[2], "^Oberauer, K. \\(2023\\)")
  # a model stored before noise existed is the Gaussian one
  old <- tcc_test_model()
  old$other_vars$noise <- NULL
  expect_length(model_citation(old), 1)
})

test_that("dtcc(noise = 'gumbel') matches the argmax integral of Gumbel signals", {
  for (f in list(tcc_test_similarity(360), tcc_similarity("laplace", kappa = 1.7),
                 tcc_similarity("vonmises", kappa = 4))) {
    for (d in c(0.3, 2, 7)) {
      k <- c(0, 5, 60, 180)
      got <- dtcc(2 * pi * k / 360, d, f, noise = "gumbel", log = TRUE)
      ref <- vapply(k, ref_tcc_gumbel_logp, numeric(1), d = d, similarity = f)
      expect_lt(max(abs(got - ref)), 1e-9, label = paste("d", d))
    }
  }
})

test_that("dtcc(noise = 'gumbel') is the softmax arm of m3() with one category per distance", {
  f <- tcc_test_similarity(36)
  mult <- c(1, rep(2, 17), 1)
  cats <- paste0("k", 0:18)
  model <- m3(resp_cats = cats, num_options = mult, choice_rule = "softmax")
  act_funs <- do.call(bmf, lapply(0:18, function(k) stats::as.formula(paste0("k", k, " ~ a", k))))
  for (d in c(0.5, 3, 9)) {
    pars <- stats::setNames(d * f, paste0("a", 0:18))
    p_cat <- vapply(seq_along(cats), function(j) {
      exp(dm3(as.numeric(seq_along(cats) == j), pars, model, act_funs, log = TRUE))
    }, numeric(1))
    got <- dtcc(2 * pi * (0:18) / 36, d, f, n_options = 36, noise = "gumbel")
    expect_equal(got * mult, p_cat, tolerance = 1e-12, info = paste("d =", d))
  }
})

test_that("dtcc(noise = 'gumbel') sums to one and rtcc() matches the argmax of Gumbel signals", {
  x <- 2 * pi * (-179:180) / 360
  for (d in c(0, 1.3, 12, 60)) {
    expect_equal(sum(dtcc(x, d, tcc_test_similarity(360), noise = "gumbel")), 1, tolerance = 1e-12)
  }
  expect_equal(dtcc(x, 0, tcc_similarity("vonmises", kappa = 3), noise = "gumbel"), rep(1 / 360, 360))
  withr::local_seed(573)
  n <- 20000
  offsets <- -179:180
  breaks <- c(-1, 0, 2, 5, 10, 20, 45, 90, 180)
  for (f in list(tcc_test_similarity(360), tcc_similarity("vonmises", kappa = 4))) {
    gumbel <- -log(-log(matrix(stats::runif(n * 360), n)))
    direct <- abs(offsets[max.col(gumbel + rep(3 * f[abs(offsets) + 1], each = n),
                                  ties.method = "first")])
    sim <- rtcc(n, 3, f, noise = "gumbel")
    expect_true(all(sim > -pi & sim <= pi))
    counts <- rbind(table(cut(direct, breaks)), table(cut(abs(round(sim * 360 / (2 * pi))), breaks)))
    expect_gt(stats::chisq.test(counts)$p.value, 0.001)
  }
})

# The bridge between the two constructors (#566, design-note row 18): softmax
# of d f with the scaled von Mises similarity
# f = (exp(kappa (cos x - 1)) - exp(-2 kappa)) / (1 - exp(-2 kappa)) is the sdm
# density exp(c sqrt(kappa / (2 pi)) exp(kappa (cos x - 1))) / Z at the option
# centres, with c = d sqrt(2 pi / kappa) / (1 - exp(-2 kappa)) in bmm's sqrtexp
# parametrization. They differ by the log option width, and by the sum over
# options that replaces Z: the trapezoid rule on a smooth periodic function,
# far below the error of dsdm()
test_that("tcc() with Gumbel noise and a fixed von Mises similarity is sdm() on the option centres", {
  x <- 2 * pi * (-179:180) / 360
  for (kappa in c(0.5, 4, 30)) {
    f <- tcc_similarity("vonmises", kappa = kappa)
    for (d in c(0.5, 3, 12)) {
      c_sdm <- d * sqrt(2 * pi / kappa) / -expm1(-2 * kappa)
      expect_lt(max(abs(dtcc(x, d, f, noise = "gumbel", log = TRUE) -
                          dsdm(x, 0, c_sdm, kappa, log = TRUE) - log(2 * pi / 360))),
                1e-8, label = paste("kappa", kappa, "d", d))
    }
  }
})

test_that("the Gumbel likelihood reads no table and passes the similarity to its own lpmf", {
  fit <- tcc_mock_fit(tcc_gumbel_model(), bmf(d ~ cond))
  code <- brms::stancode(fit)
  expect_match(code, "d = exp(d);", fixed = TRUE)
  expect_match(code, "target += tcc_gumbel_lpmf(Y | mu, d, tcc_similarity);", fixed = TRUE)
  expect_false(grepl("tcc_tab", code, fixed = TRUE))
  expect_false(grepl("tcc_logp_by_distance(", code, fixed = TRUE))
  expect_equal(brms::standata(fit)$tcc_similarity, tcc_test_similarity(360), ignore_attr = TRUE)
  expect_null(brms::standata(fit)$tcc_grid_rows)

  code <- brms::stancode(tcc_mock_fit(tcc_gumbel_model(), threads = brms::threading(2)))
  expect_match(code, "tcc_gumbel_lpmf(Y[start:end] | mu, d, tcc_similarity)", fixed = TRUE)
})

# md5 of the Gaussian Stan programs at 566-tcc 355aca27, before Gumbel noise
# was added, without the brms version line and hashed through a binary
# connection (a text-mode file gets CRLF on Windows); taken under brms 2.23.0,
# other brms versions generate other code
test_that("adding Gumbel noise leaves the Gaussian Stan programs bit-identical", {
  skip_if_not(packageVersion("brms") == "2.23.0", "hashes taken under brms 2.23.0")
  md5_code <- function(model, formula, ...) {
    lines <- strsplit(brms::stancode(tcc_mock_fit(model, formula, ...)), "\n")[[1]]
    path <- withr::local_tempfile(fileext = ".stan")
    con <- file(path, "wb")
    writeBin(charToRaw(paste(lines[!grepl("generated with brms", lines)], collapse = "\n")), con)
    close(con)
    unname(tools::md5sum(path))
  }
  f <- exp(-(0:180) / 20)
  expect_equal(md5_code(tcc("err", similarity = f), bmf(d ~ cond)),
               "5ee0c72d0d4ace4d2cff6816bd208014")
  expect_equal(md5_code(tcc("err", similarity = f), bmf(d ~ cond), threads = brms::threading(2)),
               "9a1c977823b5fe0979b495cb6514c182")
})

# md5 of the Gumbel and table-from-R programs at 566-tcc-options ef664176,
# before correlated and motor noise were added; taken as above
test_that("correlated and motor noise leave the Gumbel and table-from-R programs bit-identical", {
  skip_on_cran()
  skip_if_not(packageVersion("brms") == "2.23.0", "hashes taken under brms 2.23.0")
  md5_code <- function(model, formula, ...) {
    lines <- strsplit(brms::stancode(tcc_mock_fit(model, formula, ...)), "\n")[[1]]
    path <- withr::local_tempfile(fileext = ".stan")
    con <- file(path, "wb")
    writeBin(charToRaw(paste(lines[!grepl("generated with brms", lines)], collapse = "\n")), con)
    close(con)
    unname(tools::md5sum(path))
  }
  f <- exp(-(0:180) / 20)
  gumbel <- tcc("err", similarity = f, noise = "gumbel")
  expect_equal(md5_code(gumbel, bmf(d ~ cond)), "9c11f0e60651bbbf79eb154023d59313")
  expect_equal(md5_code(gumbel, bmf(d ~ cond), threads = brms::threading(2)),
               "426e4125e03be789b7361e08b1415d74")
  from_r <- tcc("err", similarity = f, d_max = 6)
  from_r$other_vars$table <- "r"
  expect_equal(md5_code(from_r, bmf(d ~ cond)), "4ce5b52728d80d7247b022de3ef8a2c6")
  correlated <- tcc("err", similarity = f, d_max = 6, correlation = exp(-(0:180) / 10),
                    correlation_draws = 6000, motor_sd = 0.03)
  expect_equal(md5_code(correlated, bmf(d ~ cond)), "4ce5b52728d80d7247b022de3ef8a2c6")
})

test_that("a Gumbel fit is not checked against a table, whatever its draws", {
  prep <- epred_prep(list(d = cbind(c(1, 45, 300), 2)))
  expect_silent(.tcc_check_grid(tcc_gumbel_model(), prep))
  expect_silent(postprocess_brm(tcc_gumbel_model(), tcc_mock_fit(tcc_gumbel_model())))
})

test_that("log_lik_tcc_gumbel() is dtcc() at each draw, and each cell reads its own draws", {
  d <- c(0.4, 2.2, 9, 41)
  f <- tcc_test_similarity(360)
  fit <- suppressMessages(bmm(bmf(d ~ 1), data.frame(err = 2 * pi * c(0, 7, -180, 33) / 360),
                              tcc_gumbel_model(), backend = "mock", mock_fit = 1, rename = FALSE))
  sdata <- brms::standata(fit)
  prep <- epred_prep(list(d = matrix(d, 4, sdata$N)), sdata)
  for (i in seq_len(prep$nobs)) {
    exact <- dtcc(2 * pi * prep$data$Y[i] / 360, d, f, noise = "gumbel", log = TRUE)
    expect_lt(max(abs(log_lik_tcc_gumbel(i, prep) - exact)), 1e-12)
  }

  # d ~ cond: the observations of each cell share their draws of d, and the next
  # cell has other draws (180 options in bins, one cell per sign of the offset)
  f <- tcc_test_similarity(180)
  model <- tcc("bin", similarity = f, n_options = 180, units = "bins", noise = "gumbel")
  fit <- suppressMessages(bmm(bmf(d ~ 1), data.frame(bin = c(0, 4, -4, 90)), model,
                              backend = "mock", mock_fit = 1, rename = FALSE))
  sdata <- brms::standata(fit)
  d_obs <- cbind(c(1, 2, 3), c(1, 2, 3), c(6, 0.5, 30), c(6, 0.5, 30))
  prep <- epred_prep(list(d = d_obs), sdata)
  for (i in 1:4) {
    expect_equal(log_lik_tcc_gumbel(i, prep),
                 dtcc(2 * pi * sdata$Y[i] / 180, d_obs[, i], f, n_options = 180, noise = "gumbel",
                      log = TRUE),
                 info = paste("observation", i))
  }
  withr::local_seed(575)
  prep <- epred_prep(list(d = cbind(rep(40, 2000), rep(0.01, 2000))), sdata)
  expect_lt(mean(abs(posterior_predict_tcc_gumbel(1, prep))), 3)
  expect_gt(mean(abs(posterior_predict_tcc_gumbel(2, prep))), 30)

  family <- tcc_mock_fit(tcc_gumbel_model())$formula$family
  expect_identical(family$log_lik, log_lik_tcc_gumbel)
  expect_identical(family$posterior_predict, posterior_predict_tcc_gumbel)
  expect_error(family$posterior_epred(NULL), "not defined for the tcc model")
})

test_that("posterior_predict_tcc_gumbel() and pp_simulate() draw from the Gumbel pmf", {
  withr::local_seed(574)
  n <- 20000
  breaks <- c(-1, 0, 1, 2, 5, 10, 20, 90, 180)
  mult <- c(1, rep(2, 179), 1)
  f <- tcc_similarity("vonmises", kappa = 4)
  model <- tcc_gumbel_model(f)
  fit <- suppressMessages(bmm(bmf(d ~ 1), data.frame(err = 0), model,
                              backend = "mock", mock_fit = 1, rename = FALSE))
  prep <- epred_prep(list(d = matrix(3, n, 1)), brms::standata(fit))
  expected <- tapply(dtcc(2 * pi * (0:180) / 360, 3, f, noise = "gumbel") * mult, cut(0:180, breaks), sum)
  for (sims in list(posterior_predict_tcc_gumbel(1, prep), pp_simulate(model, prep)$error[, 1])) {
    expect_true(all(sims > -180 & sims <= 180 & sims == round(sims)))
    expect_gt(stats::chisq.test(table(cut(abs(sims), breaks)), p = expected)$p.value, 0.001)
  }
  # each observation from its own d
  model <- tcc_gumbel_model()
  fit <- suppressMessages(bmm(bmf(d ~ 1), data.frame(err = rep(0, 4)), model,
                              backend = "mock", mock_fit = 1, rename = FALSE))
  sdata <- brms::standata(fit)
  sims <- pp_simulate(model, epred_prep(list(d = matrix(c(40, 40, 0.01, 0.01), 50, 4, byrow = TRUE)),
                                        sdata))$error
  expect_lt(max(colMeans(abs(sims[, 1:2]))), 3)
  expect_gt(min(colMeans(abs(sims[, 3:4]))), 60)
})

test_that("the Stan Gumbel likelihood matches dtcc() for runs and single rows", {
  skip_on_cran()
  skip_unless_stan_tests()
  rows <- data.frame(
    d = c(rep(1.234, 3), rep(0.011, 2), 2, 2.345, rep(45, 2), 4.4),
    y = c(0L, 7L, -180L, 3L, -3L, 0L, 45L, 1L, -2L, 90L)
  )
  f <- tcc_test_similarity(360)
  stanvars <- bmm:::.tcc_stanvars(tcc_gumbel_model(f), table = "stan")
  funs <- paste(vapply(stanvars, function(s) if (s$block == "functions") s$scode else "", ""),
                collapse = "\n")
  data_decl <- paste(vapply(stanvars, function(s) if (s$block == "data") s$scode else "", ""),
                     collapse = "\n")
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n", data_decl, "\n  int<lower=1> N;\n  array[N] int y;\n  vector[N] d;\n}\n",
    "generated quantities {\n",
    "  real lp_runs = tcc_gumbel_lpmf(y | rep_vector(0, N), d, tcc_similarity);\n",
    "  vector[N] lp;\n",
    "  for (i in 1:N) lp[i] = tcc_gumbel_lpmf({y[i]} | rep_vector(0, 1), rep_vector(d[i], 1),\n",
    "                                         tcc_similarity);\n",
    "}\n"
  )
  sdata <- lapply(stanvars, function(s) s$sdata)
  sdata <- sdata[!vapply(sdata, is.null, TRUE)]
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(sdata, list(N = nrow(rows), y = rows$y, d = rows$d)),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
    show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  lp <- as.numeric(csv[1, paste0("lp.", seq_len(nrow(rows)))])
  expect_lt(max(abs(lp - dtcc(2 * pi * rows$y / 360, rows$d, f, noise = "gumbel", log = TRUE))), 1e-10)
  expect_equal(csv$lp_runs, sum(lp), tolerance = 1e-12)
})
