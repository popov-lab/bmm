# tcc() with forced-choice option sets (#566, D4): the constructor, the
# option-set data contract shared with sdm(), dtcc() and rtcc() over offered
# options, the cell tables, the Stan likelihood against R, and the
# cross-constructor bridge to sdm()

# exponential stand-in for a measured colour similarity, f(0) = 1, f(180 deg) = 0
opt_similarity <- function(n_options = 360) {
  delta <- seq(0, 180, length.out = n_options / 2 + 1)
  (exp(-delta / 20) - exp(-9)) / (1 - exp(-9))
}

# log P of each offered option (offsets in steps) under normal noise, by
# adaptive quadrature over the offered options only, written without the
# package's window or its Gauss-Legendre rule
ref_set_logp <- function(d, offsets, similarity) {
  f <- similarity[abs(offsets) + 1]
  vapply(seq_along(f), function(own) {
    logint <- function(x) {
      dnorm(x - d * f[own], log = TRUE) +
        rowSums(pnorm(outer(x, d * f[-own], "-"), log.p = TRUE))
    }
    grid <- seq(-12, max(0, d) + 12, by = 0.05)
    peak <- max(logint(grid))
    res <- integrate(function(x) exp(logint(x) - peak), -Inf, Inf,
                     rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)
    peak + log(res$value)
  }, numeric(1))
}

opt_data <- function() {
  # two- and four-option trials on a 360-option wheel, in degrees from the target
  deg <- data.frame(
    err = c(0, 180, 0, 90, -90, 0, 30, 0),
    o1 = c(0, 0, 0, 0, 0, 0, 0, 0),
    o2 = c(180, 180, 90, 90, 90, 30, 30, 180),
    o3 = c(NA, NA, 180, 180, 180, -30, -30, NA),
    o4 = c(NA, NA, -90, -90, -90, NA, NA, NA),
    cond = c("a", "a", "a", "b", "b", "b", "b", "b")
  )
  deg[c("err", paste0("o", 1:4))] <- deg[c("err", paste0("o", 1:4))] * pi / 180
  deg
}

opt_model <- function(noise = "gaussian", ...) {
  tcc("err", similarity = opt_similarity(360), options = paste0("o", 1:4), noise = noise, ...)
}

opt_mock_fit <- function(model = opt_model(), formula = bmf(d ~ cond), data = opt_data(), ...) {
  suppressMessages(bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE, ...))
}

############################################################################# !
# CONSTRUCTOR AND DATA                                                   ####
############################################################################# !

test_that("tcc() with options or n_alternatives joins the option-set contract after tcc", {
  f <- opt_similarity(360)
  expect_identical(class(tcc("err", f)), c("bmmodel", "tcc"))
  expect_null(tcc("err", f)$other_vars$options)
  expect_null(tcc("err", f)$other_vars$n_alternatives)

  sets <- tcc("err", f, options = c("o1", "o2"), n_alternatives = "m")
  expect_identical(class(sets), c("bmmodel", "tcc", "option_set"))
  expect_identical(sets$other_vars$options, c("o1", "o2"))
  expect_identical(sets$other_vars$n_alternatives, "m")
  short <- tcc("err", f, n_alternatives = 4, noise = "gumbel")
  expect_identical(class(short), c("bmmodel", "tcc", "option_set"))
  expect_identical(short$other_vars$n_alternatives, 4L)
  expect_named(short$parameters, "d")
})

test_that("tcc() checks the option arguments as sdm() does", {
  f <- opt_similarity(360)
  expect_error(tcc("err", f, n_alternatives = 7), "multiple of 7")
  expect_error(tcc("err", f, n_alternatives = 1), "whole number of at least 2")
  expect_error(tcc("err", f, n_alternatives = 2.5), "whole number of at least 2")
  expect_error(tcc("err", f, n_alternatives = "m"), "needs the `options` columns")
  expect_error(tcc("err", f, options = c("a", "b"), n_alternatives = 4), "not both")
  expect_error(tcc("err", f, options = 1:2), "column names")
})

test_that("check_data() builds the option columns and keeps the tcc bin without a reserved-column warning", {
  model <- opt_model()
  expect_no_warning(checked <- check_data(model, opt_data(), bmf(d ~ 1)))
  expect_identical(checked$tcc_bin, c(0L, 180L, 0L, 90L, -90L, 0L, 30L, 0L))
  expect_identical(checked$OptIdx1, rep(0L, 8))
  expect_identical(checked$OptIdx3[1], 360L)
  expect_identical(checked$opt_set_id, c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 1L))
  expect_setequal(built_data_columns(model), c("tcc_bin", paste0("OptIdx", 1:4), "opt_set_id"))

  dat <- opt_data()
  dat$o1[3] <- dat$o2[3]
  expect_error(check_data(model, dat, bmf(d ~ 1)), "target \\(an offset of 0\\) is not among the options in row 3")
  dat <- opt_data()
  dat$err[4] <- pi / 4
  expect_error(check_data(model, dat, bmf(d ~ 1)), "not one of the options in row 4")
})

test_that("n_alternatives = m builds the same option columns as the equivalent option columns", {
  f <- opt_similarity(360)
  dat <- data.frame(err = c(0, 90, 180, -90) * pi / 180)
  short <- check_data(tcc("err", f, n_alternatives = 4), dat, bmf(d ~ 1))
  dat[paste0("o", 1:4)] <- matrix(c(0, 90, 180, -90) * pi / 180, 4, 4, byrow = TRUE)
  long <- check_data(tcc("err", f, options = paste0("o", 1:4)), dat, bmf(d ~ 1))
  built <- c("tcc_bin", paste0("OptIdx", 1:4), "opt_set_id")
  expect_identical(short[built], long[built])
})

test_that("the stored frame gets the user's response and option columns back", {
  model <- opt_model()
  checked <- check_data(model, opt_data(), bmf(d ~ 1))
  kept <- checked[c("cond", built_data_columns(model))]
  reverted <- revert_check_data(model, kept)
  expect_equal(reverted$err, opt_data()$err, tolerance = 1e-12)
  expect_equal(reverted$o3, opt_data()$o3, tolerance = 1e-12)
  expect_false(any(built_data_columns(model) %in% names(reverted)))
})

############################################################################# !
# dtcc() AND rtcc() OVER OFFERED OPTIONS                                 ####
############################################################################# !

test_that("dtcc() over an offered set sums to one for both noises and is 0 off the set", {
  f <- opt_similarity(360)
  sets <- list(c(0, 180), c(0, 90, 180, -90), c(0, 1, 2, 3), c(0, 30, -30, 170),
               seq(0, 359, by = 22.5)[c(TRUE, FALSE)], c(0, 45, 90, 135, 180, -135, -90, -45))
  for (noise in c("gaussian", "gumbel")) {
    for (set in sets) {
      x <- set * pi / 180
      for (d in c(0.01, 1.5, 6, 25)) {
        expect_equal(sum(dtcc(x, d, f, options = x, noise = noise)), 1, tolerance = 1e-9,
                     label = paste(noise, d, paste(set, collapse = ",")))
      }
      expect_identical(dtcc(5 * pi / 180, 2, f, options = x, noise = noise), 0)
      expect_identical(dtcc(5 * pi / 180, 2, f, options = x, noise = noise, log = TRUE), -Inf)
    }
  }
})

test_that("dtcc() under normal noise matches adaptive quadrature over the offered options", {
  f <- opt_similarity(360)
  sharp <- tcc_similarity("laplace", kappa = 1 / deg2rad(1))
  cases <- list(
    list(set = c(0, 180), sim = f), list(set = c(0, 90, 180, -90), sim = f),
    list(set = c(0, 1, -1, 2, 7), sim = f), list(set = c(0, 3, 177, -100), sim = sharp),
    list(set = seq(0, 337.5, by = 22.5), sim = tcc_similarity("vonmises", 4, n_options = 720))
  )
  for (case in cases) {
    n_options <- 2 * (length(case$sim) - 1)
    steps <- case$set * n_options / 360
    x <- 2 * pi * steps / n_options
    for (d in c(0.05, 1, 3.3, 12)) {
      exact <- ref_set_logp(d, .bin_options(steps, n_options, "bins"), case$sim)
      got <- dtcc(x, d, case$sim, n_options = n_options, options = x, log = TRUE)
      expect_lt(max(abs(got - exact)), 1e-8)
    }
  }
})

test_that("the two-option cell with the foil opposite reduces to the 2-AFC d'", {
  # f(180 deg) = 0: the target signal is N(d, 1), the foil's N(0, 1), so
  # P(correct) = pnorm(d / sqrt(2)) and sqrt(2) qnorm(P(correct)) = d
  f <- opt_similarity(360)
  x <- c(0, pi)
  for (d in c(0.1, 0.8, 2, 4.5)) {
    p_correct <- dtcc(0, d, f, options = x)
    expect_equal(p_correct, pnorm(d / sqrt(2)), tolerance = 1e-10)
    expect_equal(sqrt(2) * qnorm(p_correct), d, tolerance = 1e-8)
  }
})

test_that("dtcc() under Gumbel noise is the softmax of d f over the offered options", {
  f <- opt_similarity(360)
  set <- c(0, 40, -40, 120)
  x <- set * pi / 180
  for (d in c(0.3, 2, 9)) {
    a <- d * f[abs(set) + 1]
    expect_equal(dtcc(x, d, f, options = x, noise = "gumbel"), exp(a) / sum(exp(a)), tolerance = 1e-12)
  }
})

test_that("tcc(options, noise = 'gumbel') with a von Mises similarity is sdm() over the same options", {
  # S7: c = d sqrt(2 pi / kappa) / (1 - exp(-2 kappa)) makes the SDM activation
  # c* exp(kappa (cos x - 1)) equal d f(x) plus a constant, and a softmax over
  # the same options ignores the constant
  set <- c(0, 45, 90, 180, -60)
  x <- set * pi / 180
  for (kappa in c(0.5, 4, 30)) {
    f <- tcc_similarity("vonmises", kappa = kappa)
    for (d in c(0.5, 3, 12)) {
      c_sdm <- d * sqrt(2 * pi / kappa) / -expm1(-2 * kappa)
      act <- c_sdm * sqrt(kappa / (2 * pi)) * exp(kappa * (cos(x) - 1))
      expect_equal(dtcc(x, d, f, options = x, noise = "gumbel", log = TRUE),
                   act - matrixStats::logSumExp(act), tolerance = 1e-10,
                   label = paste("kappa", kappa, "d", d))
      # sdm(options = ) reports a density over the option width; it exists
      # once the option sets of sdm() (#558) are on this branch
      if ("options" %in% names(formals(dsdm))) {
        expect_lt(max(abs(dtcc(x, d, f, options = x, noise = "gumbel", log = TRUE) -
                            dsdm(x, 0, c_sdm, kappa, n_options = 360, options = x, log = TRUE) -
                            log(2 * pi / 360))), 1e-8)
      }
    }
  }
})

test_that("rtcc() draws only offered options, with the frequencies dtcc() gives", {
  withr::local_seed(566)
  f <- opt_similarity(360)
  x <- c(0, 30, -30, 90, 180) * pi / 180
  for (noise in c("gaussian", "gumbel")) {
    draws <- rtcc(20000, d = 1.5, similarity = f, options = x, noise = noise)
    expect_true(all(round(draws * 180 / pi) %in% round(x * 180 / pi)))
    observed <- table(factor(round(draws * 180 / pi), levels = round(x * 180 / pi)))
    expect_gt(stats::chisq.test(observed, p = dtcc(x, 1.5, f, options = x, noise = noise))$p.value, 0.001)
  }
})

test_that("dtcc() and rtcc() check the options as dsdm() does", {
  f <- opt_similarity(360)
  expect_error(dtcc(0, 1, f, options = c(0, NA, pi)), "missing or empty")
  expect_error(rtcc(5, 1, f, options = numeric(0)), "missing or empty")
  x <- c(0, pi / 2)
  expect_identical(dtcc(numeric(0), 1, f, options = x), numeric(0))
  expect_identical(dtcc(0, numeric(0), f, options = x), numeric(0))
  expect_equal(dtcc(c(NA, 0, 0), c(1, NA, 1), f, options = x), c(NA, NA, dtcc(0, 1, f, options = x)))
  # a repeated option, or pi and -pi, is one option
  expect_equal(dtcc(x, 1, f, options = c(x, 0, pi / 2 - 2 * pi)), dtcc(x, 1, f, options = x))
  expect_error(dtcc(0, 1, f, options = "a"), "must be numeric, in radians")
  expect_error(rtcc(5, 1, f, options = "a"), "must be numeric, in radians")
})

############################################################################# !
# CELL TABLES AND THE GENERATED STAN CODE                                ####
############################################################################# !

test_that("the cells are runs of equal predictors times option sets, with each row's cell", {
  fit <- opt_mock_fit()
  sdata <- brms::standata(fit)
  # cond a: rows 1-2 set 1, row 3 set 2; cond b: rows 4-5 set 2, rows 6-7 set 3, row 8 set 1
  expect_equal(as.integer(sdata$tcc_row_cell), c(1L, 1L, 2L, 3L, 3L, 4L, 4L, 5L))
  expect_equal(as.integer(sdata$tcc_cell_set), c(1L, 2L, 2L, 3L, 1L))
  expect_equal(as.integer(sdata$tcc_set_size), c(2L, 4L, 3L))
  expect_equal(sdata$G_tcc_cells, 5)
  # sorted signed offsets per set, padded with 0
  expect_equal(unname(sdata$tcc_set_off[2, ]), c(-90L, 0L, 90L, 180L))
  expect_equal(unname(sdata$tcc_set_off[3, ]), c(-30L, 0L, 30L, 0L))
  expect_equal(sdata$tcc_noise, 1)

  code <- brms::stancode(fit)
  expect_match(code, "target += tcc_options_lpmf(Y | mu, d);", fixed = TRUE)
  expect_match(code, "target += tcc_options_cells(d, 1, N, ", fixed = TRUE)
  expect_match(code, "tcc_cell_n", fixed = TRUE)
  expect_false(grepl("tcc_tab", code, fixed = TRUE))
  expect_null(sdata$tcc_grid_rows)

  code <- brms::stancode(opt_mock_fit(threads = brms::threading(2)))
  expect_match(code, "tcc_options_cells(d, start, end, ", fixed = TRUE)
  expect_equal(brms::standata(opt_mock_fit(opt_model("gumbel")))$tcc_noise, 2)
})

test_that("adjacent rows with different predictor values never share a cell", {
  # interaction() pastes labels with ".", so (1, 5.2) and (1.5, 2) are both
  # "1.5.2", as are the levels ("a.b", "c") and ("a", "b.c"); one cell across
  # them would give every row the d of the first
  dat <- opt_data()[rep(1, 4), ]
  dat$x1 <- c(1, 1, 1.5, 1.5)
  dat$x2 <- c(2, 5.2, 2, 5.2)
  sdata <- brms::standata(opt_mock_fit(formula = bmf(d ~ x1 + x2), data = dat))
  expect_equal(as.integer(sdata$tcc_row_cell), 1:4)
  dat <- opt_data()[rep(1, 2), ]
  dat$a <- c("a.b", "a")
  dat$b <- c("c", "b.c")
  sdata <- brms::standata(opt_mock_fit(formula = bmf(d ~ a + b), data = dat))
  expect_equal(as.integer(sdata$tcc_row_cell), 1:2)
})

test_that("one cell and one set keep their tables arrays, which brms would declare as scalars", {
  model <- tcc("err", opt_similarity(360), n_alternatives = 4)
  dat <- data.frame(err = c(0, 90, 180, -90, 0) * pi / 180)
  code <- brms::stancode(opt_mock_fit(model, bmf(d ~ 1), dat))
  expect_match(code, "array[G_tcc_cells] int tcc_cell_set;", fixed = TRUE)
  expect_match(code, "array[S_tcc_sets] int tcc_set_size;", fixed = TRUE)
  expect_match(code, "array[N] int tcc_row_cell;", fixed = TRUE)
  code <- brms::stancode(opt_mock_fit(model, bmf(d ~ 1), dat[1, , drop = FALSE]))
  expect_match(code, "array[N] int tcc_row_cell;", fixed = TRUE)
})

test_that("rows with a missing response or predictor are left out of the cells as brms drops them", {
  dat <- opt_data()
  dat$err[2] <- NA
  dat$cond[6] <- NA
  # brms warns that it drops them, once per pass over the data
  fit <- suppressWarnings(opt_mock_fit(data = dat))
  sdata <- brms::standata(fit)
  expect_equal(sdata$N, 6)
  # rows 2 and 6 dropped: cond a rows 1, 3; cond b rows 4-5, 7, 8
  expect_equal(as.integer(sdata$tcc_row_cell), c(1L, 2L, 3L, 3L, 4L, 5L))
})

test_that("each row's own cell is reported, since a fit then sums once per row", {
  dat <- opt_data()
  dat$x <- seq_len(nrow(dat))
  expect_message(bmm(bmf(d ~ x), dat, opt_model(), backend = "mock", mock_fit = 1, rename = FALSE),
                 "cells of equal d and offered options")
  dat <- do.call(rbind, rep(list(opt_data()), 3))
  expect_no_message(bmm(bmf(d ~ cond), dat[order(dat$cond), ], opt_model(), backend = "mock",
                        mock_fit = 1, rename = FALSE))
})

test_that("the cells are reported once they outnumber half the rows, not only when each row is one", {
  # the rows' sets are 1 1 2 2 2 3 3 1, so cond a-a-a-b-b-b-b-b gives 5 cells
  # of 8 rows and a single b at the end gives 4, exactly half
  dat <- opt_data()
  expect_equal(brms::standata(opt_mock_fit(data = dat))$G_tcc_cells, 5)
  expect_message(bmm(bmf(d ~ cond), dat, opt_model(), backend = "mock", mock_fit = 1, rename = FALSE),
                 "The 8 rows fall into 5 cells")
  dat$cond <- c(rep("a", 7), "b")
  expect_equal(brms::standata(opt_mock_fit(data = dat))$G_tcc_cells, 4)
  expect_no_message(bmm(bmf(d ~ cond), dat, opt_model(), backend = "mock", mock_fit = 1,
                        rename = FALSE))
})

test_that("an options fit reads no d table and is not checked against one", {
  prep <- epred_prep(list(d = cbind(c(1, 45, 300), 2)))
  expect_silent(.tcc_check_grid(opt_model(), prep))
  expect_silent(postprocess_brm(opt_model(), opt_mock_fit()))
})

############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT, PP_CHECK                                   ####
############################################################################# !

test_that("log_lik_tcc_options() is dtcc() over each row's own options, for both noises", {
  d <- c(0.4, 2.2, 9)
  for (noise in c("gaussian", "gumbel")) {
    fit <- opt_mock_fit(opt_model(noise), formula = bmf(d ~ 1))
    sdata <- brms::standata(fit)
    prep <- epred_prep(list(d = matrix(d, 3, sdata$N)), sdata)
    checked <- check_data(opt_model(noise), opt_data(), bmf(d ~ 1))
    for (i in seq_len(prep$nobs)) {
      x <- unlist(checked[i, paste0("o", 1:4)])
      exact <- dtcc(checked$err[i], d, opt_similarity(360), options = x[!is.na(x)], noise = noise,
                    log = TRUE)
      expect_lt(max(abs(log_lik_tcc_options(i, prep) - exact)), 1e-12)
    }
  }
  family <- opt_mock_fit()$formula$family
  expect_identical(family$log_lik, log_lik_tcc_options)
  expect_identical(family$posterior_predict, posterior_predict_tcc_options)
  expect_error(family$posterior_epred(NULL), "not defined for the tcc model")
})

test_that("log_lik_tcc_options() reads each row's own d, and gives a response off the set -Inf", {
  # rows 1, 2 and 8 offer one set: consecutive rows with other d must not
  # reuse the probabilities of the row before
  fit <- opt_mock_fit()
  sdata <- brms::standata(fit)
  d <- matrix(rep(0.3 * seq_len(sdata$N), each = 2), 2, sdata$N) * c(1, 2)
  prep <- epred_prep(list(d = d), sdata)
  checked <- check_data(opt_model(), opt_data(), bmf(d ~ cond))
  for (i in seq_len(prep$nobs)) {
    x <- unlist(checked[i, paste0("o", 1:4)])
    exact <- dtcc(checked$err[i], d[, i], opt_similarity(360), options = x[!is.na(x)], log = TRUE)
    expect_lt(max(abs(log_lik_tcc_options(i, prep) - exact)), 1e-12)
  }
  # newdata whose response the trial did not offer
  prep$data$Y[1] <- 90L
  expect_identical(log_lik_tcc_options(1, prep), rep(-Inf, 2))
})

test_that("posterior_predict_tcc_options() and pp_simulate() draw from each row's options", {
  withr::local_seed(567)
  model <- opt_model()
  fit <- opt_mock_fit(model, formula = bmf(d ~ 1))
  sdata <- brms::standata(fit)
  prep <- epred_prep(list(d = matrix(1.5, 4000, sdata$N)), sdata)
  checked <- check_data(model, opt_data(), bmf(d ~ 1))
  sims <- pp_simulate(model, prep)$error
  expect_equal(dim(sims), c(4000, sdata$N))
  for (i in seq_len(sdata$N)) {
    offered <- unlist(checked[i, paste0("OptIdx", 1:4)])
    offered <- offered[offered != 360]
    one <- posterior_predict_tcc_options(i, prep)
    expect_true(all(one %in% offered))
    expect_true(all(sims[, i] %in% offered))
  }
  # row 3 offers 0, 90, 180, -90: the target is drawn about as often as dtcc() says
  p <- dtcc(0, 1.5, opt_similarity(360), options = c(0, 90, 180, -90) * pi / 180)
  expect_equal(mean(posterior_predict_tcc_options(3, prep) == 0), p, tolerance = 0.05)
})

############################################################################# !
# STAN AGAINST R                                                         ####
############################################################################# !

test_that("the Stan option-set likelihood matches dtcc() per set and per cell, for both noises", {
  skip_on_cran()
  skip_unless_stan_tests()
  f <- opt_similarity(360)
  sharp <- tcc_similarity("laplace", kappa = 1 / deg2rad(1))
  funs <- read_lines2(system.file("stan_chunks", "tcc_options_functions.stan", package = "bmm"))
  # the 16 options of Robinson et al. (2023) at 22.5 degrees, binned onto the 360 grid
  sets <- list(c(0, 180), c(-90, 0, 90, 180), c(-100, 0, 3, 177), seq(-135, 180, by = 45),
               sort(.bin_options(seq(0, 337.5, by = 22.5) * pi / 180, 360, "radians")),
               c(-2, -1, 0, 1, 2, 7, 120))
  d_values <- c(0.011, 1.234, 5.987, 25, 70)
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n  int D;\n  vector[D] sim;\n  vector[D] sharp;\n  int Q;\n  vector[Q] gl_nodes;\n",
    "  vector[Q] gl_logw;\n  int M;\n  int S;\n  array[S] int size;\n  array[S, M] int off;\n",
    "  int ND;\n  vector[ND] dv;\n}\n",
    "generated quantities {\n",
    "  array[2, 2, S, ND] vector[M] lp;\n",
    "  for (noise in 1:2) for (k in 1:2) for (s in 1:S) for (j in 1:ND) {\n",
    "    vector[size[s]] fs;\n",
    "    for (i in 1:size[s]) fs[i] = k == 1 ? sim[abs(off[s, i]) + 1] : sharp[abs(off[s, i]) + 1];\n",
    "    lp[noise, k, s, j] = rep_vector(0, M);\n",
    "    lp[noise, k, s, j, 1:size[s]] = tcc_set_logp(dv[j], fs, noise, gl_nodes, gl_logw);\n",
    "  }\n",
    "}\n"
  )
  rule <- .gl_rule(.tcc_gl_nodes)
  off <- matrix(0L, length(sets), max(lengths(sets)))
  for (s in seq_along(sets)) off[s, seq_along(sets[[s]])] <- as.integer(sets[[s]])
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = list(D = 181, sim = f, sharp = sharp, Q = length(rule$nodes), gl_nodes = rule$nodes,
                gl_logw = log(rule$weights), M = ncol(off), S = nrow(off), size = lengths(sets),
                off = off, ND = length(d_values), dv = d_values),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
    show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  for (noise in 1:2) for (k in 1:2) for (s in seq_along(sets)) for (j in seq_along(d_values)) {
    x <- 2 * pi * sets[[s]] / 360
    r <- dtcc(x, d_values[j], list(f, sharp)[[k]], options = x, noise = c("gaussian", "gumbel")[noise],
              log = TRUE)
    stan <- as.numeric(csv[1, sprintf("lp.%d.%d.%d.%d.%d", noise, k, s, j, seq_along(x))])
    expect_lt(max(abs(stan - r)), 1e-8, label = paste(noise, k, s, d_values[j]))
  }
})

test_that("the cells counted in transformed data give the summed log_lik, whole and in slices", {
  skip_on_cran()
  skip_unless_stan_tests()
  withr::local_seed(568)
  # 60 trials in three conditions, sorted, each offering one of four sets of 2
  # to 5 options, so that cells hold many rows with different responses
  n <- 60
  dat <- data.frame(cond = sort(sample(c("a", "b", "c"), n, TRUE)))
  pool <- list(c(0, 180, NA, NA, NA), c(0, 90, 180, -90, NA), c(0, 10, -10, NA, NA),
               c(0, 45, -120, 150, 170))
  offs <- do.call(rbind, pool[sample(4, n, TRUE)])
  dat[paste0("o", 1:5)] <- offs * pi / 180
  dat$err <- vapply(seq_len(n), function(i) sample(stats::na.omit(offs[i, ]), 1), numeric(1)) * pi / 180
  for (noise in c("gaussian", "gumbel")) {
    model <- tcc("err", opt_similarity(360), options = paste0("o", 1:5), noise = noise)
    fit <- suppressMessages(bmm(bmf(d ~ cond), dat, model, backend = "mock", mock_fit = 1,
                                rename = FALSE))
    sdata <- brms::standata(fit)
    expect_gt(max(table(sdata$tcc_row_cell)), 3)
    # the split after K cuts the largest cell: its first row in one slice, the rest in the other
    biggest <- which(sdata$tcc_row_cell == which.max(tabulate(sdata$tcc_row_cell)))
    K <- biggest[1]
    expect_gt(max(biggest), K)
    d <- c(a = 0.7, b = 2.1, c = 6)[fit$data$cond]
    # every random-foil row is its own cell, which configure reports
    stanvars <- suppressMessages(.tcc_options_configure(model, fit$data, bmf(d ~ cond)))$stanvars
    block <- function(b) paste(vapply(stanvars, function(s) if (s$block == b) s$scode else "", ""),
                               collapse = "\n")
    cells_args <- paste("tcc_noise, tcc_cell_first, tcc_cell_set, tcc_cell_n, tcc_set_size,",
                        "tcc_set_off, tcc_similarity, tcc_gl_nodes, tcc_gl_logw")
    program <- paste0(
      "functions {\n", block("functions"), "\n}\n",
      "data {\n  int<lower=1> N;\n  array[N] int Y;\n  vector[N] d;\n  int K;\n  int H;\n", block("data"),
      "\n}\n",
      "transformed data {\n", block("tdata"), "\n}\n",
      "generated quantities {\n",
      "  real lp = tcc_options_cells(d, 1, N, ", cells_args, ");\n",
      "  real lp_slices = tcc_options_cells(d[1:K], 1, K, ", cells_args, ")\n",
      "    + tcc_options_cells(d[(K + 1):N], K + 1, N, ", cells_args, ");\n",
      "  real lp_halves = tcc_options_cells(d[1:H], 1, H, ", cells_args, ")\n",
      "    + tcc_options_cells(d[(H + 1):N], H + 1, N, ", cells_args, ");\n",
      "}\n"
    )
    sv_data <- lapply(stanvars, function(s) s$sdata)
    sv_data <- sv_data[!vapply(sv_data, is.null, TRUE)]
    out <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
      data = c(sv_data, list(N = sdata$N, Y = sdata$Y, d = unname(d), K = K, H = sdata$N %/% 2)),
      fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
      show_messages = FALSE, sig_figs = 17
    )
    csv <- utils::read.csv(out$output_files()[1], comment.char = "#", check.names = FALSE)
    prep <- epred_prep(list(d = matrix(d, 1)), sdata)
    total <- sum(vapply(seq_len(sdata$N), log_lik_tcc_options, numeric(1), prep = prep))
    expect_equal(csv$lp, total, tolerance = 1e-10, label = noise)
    expect_equal(csv$lp_slices, total, tolerance = 1e-10, label = noise)
    # halves: the second slice's cells must read d from their own offset
    expect_equal(csv$lp_halves, total, tolerance = 1e-10, label = noise)
  }
})
