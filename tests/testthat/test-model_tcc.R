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

test_that("tcc() builds an independent model with d on the identity link", {
  model <- tcc_test_model()
  expect_s3_class(model, "tcc")
  expect_false(inherits(model, "circular"))
  expect_named(model$parameters, "d")
  expect_equal(model$links$d, "identity")
  expect_equal(model$other_vars$n_options, 360)
  expect_equal(model$other_vars$units, "radians")
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
    expand.grid(n = 360, d = c(-1.5, 0, 0.37, 11.97, 12), k = c(0, 1, 179, 180)),
    expand.grid(n = 180, d = c(0.37, 12), k = c(0, 1, 89, 90))
  )
  for (i in seq_len(nrow(cells))) {
    cl <- cells[i, ]
    f <- tcc_test_similarity(cl$n)
    x <- 2 * pi * cl$k / cl$n
    # absolute: testthat's tolerance is relative, and log P reaches -100 here
    error <- dtcc(x, d = cl$d, similarity = f, n_options = cl$n, log = TRUE) -
      ref_tcc_logp(cl$d, cl$k, f)
    expect_lt(abs(error), 1e-7, label = paste("n =", cl$n, "d =", cl$d, "k =", cl$k))
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

test_that("the Stan table is built in transformed data and read by the lpmf", {
  code <- brms::stancode(tcc_mock_fit(formula = bmf(d ~ cond)))
  expect_match(code, "matrix[tcc_grid_rows, rows(tcc_similarity)] tcc_tab;", fixed = TRUE)
  expect_match(code, "target += tcc_lpmf(Y | mu, d, tcc_tab, tcc_grid_lo, tcc_grid_step);",
               fixed = TRUE)
  sdata <- brms::standata(tcc_mock_fit())
  expect_equal(sdata$tcc_grid_lo, -2)
  expect_equal(sdata$tcc_grid_step, 0.02)
  expect_equal(sdata$tcc_grid_rows, 703)
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

# the kernels run in a fixed_param program over the stanvars bmm generates, as
# the other parity tests do (expose_functions() fails inside the suite)
tcc_stan_lp <- function(model, rows) {
  stanvars <- bmm:::.tcc_stanvars(model, table = "stan")
  r_table <- bmm:::.tcc_table(bmm:::.tcc_similarity(model), bmm:::.tcc_grid(model))
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
    "  for (i in 1:N) lp[i] = tcc_lpmf({y[i]} | rep_vector(0, 1), rep_vector(d[i], 1), tcc_tab,\n",
    "                                  tcc_grid_lo, tcc_grid_step);\n",
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
  list(table_diff = csv$table_diff, lp_runs = csv$lp_runs, lp_runs_r = csv$lp_runs_r,
       lp = as.numeric(csv[1, paste0("lp.", seq_len(nrow(rows)))]))
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
})
