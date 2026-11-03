# tcc() with correlated and motor noise (#566): the constructor, the Monte
# Carlo table the likelihood reads (drift-conditioned estimator), its
# certification and seed, the motor kernel, and dtcc() / rtcc()

# the exponential stand-in similarity of test-model_tcc.R, on n options
corr_similarity <- function(n_options = 360) {
  delta <- seq(0, 180, length.out = n_options / 2 + 1)
  (exp(-delta / 20) - exp(-9)) / (1 - exp(-9))
}

# an exponential correlation by distance, positive definite on the circle
corr_correlation <- function(n_options = 360, scale = 10) {
  exp(-seq(0, 180, length.out = n_options / 2 + 1) / scale)
}

# the d grid of the core for d_max (log link)
corr_grid <- function(d_max) {
  step <- 0.02
  list(lo = 0, step = step, rows = as.integer(ceiling(round((d_max + step) / step, 6)) + 3))
}

corr_grid_d <- function(grid) grid$lo + (seq_len(grid$rows) - 2) * grid$step

# P of the chosen option at each distance 0, ..., n / 2 by brute force: the
# dense matrix C from the correlation, its symmetric square root by eigen()
# (what MASS::mvrnorm() uses), the argmax of d f + z, motor noise added by
# rounding a normal and wrapped; no FFT, no envelopes, no table
corr_brute_force <- function(d, similarity, correlation, n_draws, motor_steps = 0) {
  n <- 2 * (length(similarity) - 1)
  dist <- abs(outer(0:(n - 1), 0:(n - 1), "-"))
  dist <- pmin(dist, n - dist)
  e <- eigen(matrix(correlation[dist + 1], n), symmetric = TRUE)
  root <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  f <- similarity[pmin(0:(n - 1), n - 0:(n - 1)) + 1]
  z <- root %*% matrix(stats::rnorm(n * n_draws), n)
  x <- max.col(t(z + d * f), ties.method = "first") - 1
  if (motor_steps > 0) x <- x + round(stats::rnorm(n_draws, 0, motor_steps))
  x <- x %% n
  counts <- tabulate(pmin(x, n - x) + 1, n / 2 + 1)
  counts / (n_draws * c(1, rep(2, n / 2 - 1), 1))
}

############################################################################# !
# CONSTRUCTOR                                                            ####
############################################################################# !

test_that("tcc() stores correlated and motor noise only when given, and then reads the table from R", {
  f <- corr_similarity()
  plain <- tcc("err", similarity = f)
  expect_equal(plain$other_vars$table, "stan")
  expect_false(any(c("correlation", "correlation_draws", "motor_sd") %in% names(plain$other_vars)))
  expect_identical(tcc("err", similarity = f, motor_sd = 0)$other_vars, plain$other_vars)

  corr <- tcc("err", similarity = f, correlation = corr_correlation())
  expect_equal(corr$other_vars$table, "r")
  expect_equal(corr$other_vars$correlation, corr_correlation())
  # draws come in groups of six rotations
  expect_equal(corr$other_vars$correlation_draws, 1000002)
  expect_equal(
    tcc("err", f, correlation = corr_correlation(), correlation_draws = 6001)$other_vars$correlation_draws,
    6006
  )
  expect_match(corr$parameters$d, "d' of Schurgin et al., 2020, under the noise correlation")

  motor <- tcc("err", similarity = f, motor_sd = 0.05)
  expect_equal(motor$other_vars$table, "r")
  expect_equal(motor$other_vars$motor_sd, 0.05)
  expect_null(motor$other_vars$correlation_draws)
  expect_output(print(corr), "Noise: +correlated by distance \\(likelihood table from 1,000,002")
  expect_output(print(motor), "Noise: +motor noise with SD 0.05 radians")
})

test_that("tcc() refuses correlated or motor noise it cannot fit", {
  f <- corr_similarity()
  rho <- corr_correlation()
  expect_error(tcc("err", f, correlation = rho, noise = "gumbel"), "requires noise = 'gaussian'")
  expect_error(tcc("err", f, motor_sd = 0.1, noise = "gumbel"), "`motor_sd` requires noise")
  expect_error(tcc("err", f, correlation = rho, n_alternatives = 4), "planned for bmm 1.6.0")
  expect_error(tcc("err", f, motor_sd = 0.1, n_alternatives = 4), "a forced choice has none")
  expect_error(tcc("err", f, correlation = rho[-1]), "must be a numeric vector of 181 values")
  expect_error(tcc("err", f, correlation = c(1, rho[-1] * 2)), "between -1 and 1")
  expect_error(tcc("err", f, correlation = c(0.9, rho[-1])), "must be 1 at distance 0")
  # a measured correlation, such as Schurgin et al.'s, is 1 at 0 only to
  # rounding, which can fall on either side of 1
  expect_silent(tcc("err", f, correlation = c(1 - 6e-15, rho[-1])))
  expect_silent(tcc("err", f, correlation = c(1 + 2e-16, rho[-1])))
  expect_error(tcc("err", f, correlation = c(1 + 1e-6, rho[-1])), "must be 1 at distance 0")
  expect_error(tcc("err", f, correlation = c(1, rep(-0.9, 180))), "negative eigenvalue")
  # rank 3: the similarity lives mostly on modes without noise
  expect_error(
    tcc("err", f, correlation = (1 + cos((0:180) * pi / 180)) / 2),
    "leaves no noise in directions that carry 0.5"
  )
  expect_error(tcc("err", f, correlation = rho, correlation_draws = 100), "at least 6000")
  expect_error(tcc("err", f, motor_sd = -1), "non-negative")
  expect_error(tcc("err", f, motor_sd = c(1, 2)), "non-negative")
})

############################################################################# !
# THE MONTE CARLO TABLE                                                  ####
############################################################################# !

test_that("the drift splits each draw into a part along the similarity and an independent rest", {
  # Delta = b'z has variance 1 / kappa^2 and is uncorrelated with z - Delta
  # f_R exactly when b'f_R = 1 and C b = f_R / kappa^2, on the dense C
  f <- corr_similarity(36)
  rho <- corr_correlation(36, scale = 30)
  setup <- bmm:::.tcc_corr_setup(rho, f)
  dist <- abs(outer(0:35, 0:35, "-"))
  C <- matrix(rho[pmin(dist, 36 - dist) + 1], 36)
  expect_equal(sum(setup$b * setup$f_R), 1, tolerance = 1e-12)
  expect_equal(drop(C %*% setup$b), setup$f_R / setup$kappa^2, tolerance = 1e-10)
  expect_equal(sum(setup$f_R), 0, tolerance = 1e-12)
  expect_equal(setup$dropped, 0, tolerance = 1e-12)
})

test_that("the envelope's winners are the argmax at every drift, and its segments tile the line", {
  skip_on_cran()
  withr::local_seed(6)
  f <- corr_similarity(36)
  setup <- bmm:::.tcc_corr_setup(corr_correlation(36, scale = 30), f)
  rho_perp <- bmm:::.tcc_corr_split(bmm:::.tcc_corr_draw(200, setup), setup, shift = 7)
  env <- bmm:::.tcc_corr_envelope(rho_perp, setup$f_R)
  for (i in 1:200) {
    seg <- which(env$draw == i)
    seg <- seg[order(env$lo[seg])]
    expect_equal(env$lo[seg[1]], -Inf)
    expect_equal(env$hi[seg[length(seg)]], Inf)
    expect_equal(env$lo[seg[-1]], env$hi[seg[-length(seg)]])
    s <- stats::runif(5, -8, 14)
    winner <- vapply(s, function(x) env$opt[seg][findInterval(x, env$lo[seg])], numeric(1))
    expect_equal(winner, max.col(outer(s, setup$f_R) + rep(rho_perp[i, ], each = 5), ties.method = "first"))
  }
})

test_that("a segment binned on the s grid gives its exact probability at every d", {
  sg <- bmm:::.tcc_corr_sgrid(c(0, 1, 2), kappa = 1.2)
  acc <- list(bins = numeric(4 * sg$nb), above = numeric(4))
  env <- list(draw = c(1, 1, 1), opt = c(2, 3, 4), lo = c(-Inf, 0.537, 1.913), hi = c(0.537, 1.913, Inf))
  acc <- bmm:::.tcc_corr_accumulate(acc, env, sg)
  d <- c(0, 0.5, 1.2, 2)
  P <- matrix(acc$bins, 4) %*% stats::pnorm(outer(sg$s, d, "-") * 1.2) + acc$above
  exact <- rbind(0, stats::pnorm(1.2 * (0.537 - d)),
                 stats::pnorm(1.2 * (1.913 - d)) - stats::pnorm(1.2 * (0.537 - d)),
                 1 - stats::pnorm(1.2 * (1.913 - d)))
  expect_lt(max(abs(P - exact)), 2e-5)
})

test_that("the certification's per-draw terms add up to the table's own shift", {
  skip_on_cran()
  # with every group kept, mean psi = -sum_k mult_k s_k(d) P_k(d) over the
  # table it certifies: an identity up to the binning, which checks the score
  # of each option, its motor convolution and the segments' bookkeeping
  f <- corr_similarity(36)
  grid <- corr_grid(4)
  for (motor in c(0, 1.5)) {
    built <- bmm:::.tcc_build_noise_table(corr_correlation(36, 30), f, 6000, motor, grid, seed = 7)
    L <- t(built$tab)
    at <- round((built$info$rho_d - grid$lo) / grid$step) + 2
    slope <- (L[, at + 1] - L[, at - 1]) / (2 * grid$step)
    terms <- c(1, rep(2, 17), 1) * exp(L[, at]) * slope
    # the table bins the segment ends on the s grid and psi does not, which
    # leaves about 1e-5 of the terms; a score one option off gives 0.18
    expect_lt(max(abs(built$info$shift + colSums(terms)) / colSums(abs(terms))), 1e-4)
  }
})

test_that("the correlated table matches a brute-force argmax of correlated signals, with and without motor noise", {
  skip_on_cran()
  withr::local_seed(566)
  # 36 options keep the dense square root and the argmax cheap
  f <- corr_similarity(36)
  rho <- corr_correlation(36, scale = 30)
  grid <- corr_grid(3)
  d_values <- corr_grid_d(grid)
  for (motor in c(0, 1.5)) {
    tab <- bmm:::.tcc_build_noise_table(rho, f, 6e4, motor, grid, seed = 1)$tab
    for (d in c(1, 2.5)) {
      p_table <- exp(tab[which.min(abs(d_values - d)), ])
      p_brute <- corr_brute_force(d, f, rho, 2e5, motor)
      mult <- c(1, rep(2, 17), 1)
      # the table's error is at most that of crude Monte Carlo with its draws
      se <- sqrt(p_table * (1 - p_table) / mult * (1 / 2e5 + 1 / 6e4))
      expect_lt(max(abs(p_table - p_brute) / se), 4.5)
    }
  }
})

test_that("with uncorrelated noise the correlated table is the independent-noise TCC", {
  skip_on_cran()
  # C = I is the estimator's worst case (kappa is large, so the drift
  # integral is narrow), and .tcc_logp() gives the exact answer. Ten batches
  # give a batch-means SE; the z scores are over-dispersed at this size, so
  # the bounds catch gross errors only (a 10 % error in kappa gives z in the
  # tens), on cells with an expected count of at least 100
  f <- corr_similarity(60)
  identity <- c(1, rep(0, 30))
  setup <- bmm:::.tcc_corr_setup(identity, f)
  d_values <- seq(0, 4, by = 0.5)
  batches <- vapply(1:10, function(i) {
    withr::with_seed(i, bmm:::.tcc_corr_build(setup, d_values, 6000)$P)
  }, matrix(0, 60, length(d_values)))
  exact <- exp(vapply(d_values, bmm:::.tcc_logp, numeric(31), similarity = f))
  exact <- exact[pmin(0:59, 60 - 0:59) + 1, ]
  mean_p <- apply(batches, c(1, 2), mean)
  se <- apply(batches, c(1, 2), stats::sd) / sqrt(10)
  use <- exact * 6e4 >= 100 & se > 0
  z <- (mean_p - exact)[use] / se[use]
  expect_gt(sum(use), 200)
  expect_lt(sqrt(mean(z^2)), 2)
  expect_lt(max(abs(z)), 8)
})

test_that("each column of the table is a pmf over the options, before and after pooling the mirrors", {
  f <- corr_similarity(60)
  setup <- bmm:::.tcc_corr_setup(corr_correlation(60, 20), f)
  P <- withr::with_seed(2, bmm:::.tcc_corr_build(setup, c(0, 0.7, 2, 5), 6000)$P)
  expect_lt(max(abs(colSums(P) - 1)), 1e-12)
  pooled <- exp(bmm:::.tcc_by_distance(log(P)))
  expect_lt(max(abs(colSums(c(1, rep(2, 29), 1) * pooled) - 1)), 1e-12)
  expect_equal(pooled[2, ], (P[2, ] + P[60, ]) / 2)
})

test_that("the motor kernel is a rounded normal, and the log convolution keeps tiny and zero cells", {
  kernel <- bmm:::.tcc_motor_kernel(2)
  expect_length(kernel, 25)
  expect_equal(sum(kernel), 1, tolerance = 1e-14)
  expect_equal(kernel, rev(kernel))
  expect_equal(kernel[13], pnorm(0.25) - pnorm(-0.25), tolerance = 1e-6)
  withr::local_seed(3)
  rounded <- tabulate(round(stats::rnorm(1e6, 0, 2)) + 13, 25)
  # the outer cells expect 0.005 draws, where a single one gives z = 15
  # (2 of 20 other seeds); they are compared only where 5 or more are expected
  keep <- 1e6 * kernel >= 5
  expect_lt(max(abs(rounded - 1e6 * kernel)[keep] / sqrt(1e6 * kernel[keep])), 4.5)

  # 40 options, two columns: one with a P of 1e-300, one with exact zeros
  P <- cbind(c(0.5, 0.3, rep(0.2 / 37, 37), 1e-300), c(1, rep(0, 39)))
  L <- bmm:::.tcc_convolve_log(log(P), c(0.25, 0.5, 0.25))
  direct <- 0.25 * P[c(40, 1:39), ] + 0.5 * P + 0.25 * P[c(2:40, 1), ]
  expect_equal(exp(L), direct, tolerance = 1e-14)
  expect_equal(L[39, 1], log(0.25 * 0.2 / 37 + 0.5 * 0.2 / 37 + 0.25 * 1e-300), tolerance = 1e-14)
  expect_equal(L[3:39, 2], rep(-Inf, 37))
})

test_that("the table's seed depends on its inputs only, and the user's RNG is left as it was", {
  skip_on_cran()
  f <- corr_similarity(60)
  model <- tcc("err",
    similarity = f, n_options = 60, correlation = corr_correlation(60, 20),
    correlation_draws = 6000, d_max = 2
  )
  grid <- bmm:::.tcc_grid(model)
  ns <- asNamespace("bmm")
  old <- ns$.tcc_cache$noise
  withr::defer(assign("noise", old, envir = ns$.tcc_cache))
  build <- function() {
    assign("noise", NULL, envir = ns$.tcc_cache)
    bmm:::.tcc_noise_table(model, grid)
  }
  withr::local_seed(42)
  before <- .Random.seed
  first <- build()
  expect_identical(.Random.seed, before)
  second <- withr::with_seed(1, build(), .rng_kind = "L'Ecuyer-CMRG")
  expect_identical(second, first)
  more <- model
  more$other_vars$correlation_draws <- 6006
  expect_false(identical(build_seed <- bmm:::.tcc_noise_table(more, grid)$info$seed, first$info$seed))
  # the cache answers a second call: it still holds the two tables built
  # since the last reset, and no third
  expect_length(ns$.tcc_cache$noise, 2)
  expect_identical(bmm:::.tcc_noise_table(more, grid)$info$seed, build_seed)
  expect_length(ns$.tcc_cache$noise, 2)
  # a model that differs only in d_max or in motor noise gets its own table,
  # never the cached one of another grid or kernel
  wider <- model
  wider$other_vars$d_max <- 3
  wide_grid <- bmm:::.tcc_grid(wider)
  expect_equal(nrow(bmm:::.tcc_noise_table(wider, wide_grid)$tab), wide_grid$rows)
  motor <- model
  motor$other_vars$motor_sd <- 0.1
  expect_gt(max(abs(bmm:::.tcc_noise_table(motor, grid)$tab - first$tab)), 0.01)
})

test_that("the certification groups each draw with its rotations across chunks", {
  # 600-envelope chunks keep the segments of 300 groups from three chunks;
  # within one group the segments of its six envelopes tile the drift axis
  # six times, so its probabilities sum to 6 at any d, and ids that restart
  # in each chunk would merge groups and sum to 12 or 18
  setup <- bmm:::.tcc_corr_setup(corr_correlation(36, 30), corr_similarity(36))
  build <- withr::with_seed(8, bmm:::.tcc_corr_build(setup, c(0, 1), 1800, keep_groups = 250, chunk = 600))
  kept <- build$kept
  expect_setequal(unique(kept$draw), 1:300)
  total <- rowsum(stats::pnorm(setup$kappa * (kept$hi - 1.3)) - stats::pnorm(setup$kappa * (kept$lo - 1.3)), kept$draw)
  expect_equal(range(total), c(6, 6), tolerance = 1e-12)
})

test_that("the certified number of trials predicts the spread of tables built with other seeds", {
  skip_on_cran()
  # 2 x 30 tables at 6,000 draws on 36 options take about 5 s. The shift each
  # table induces in an estimate of d, in posterior SDs per root trial, should
  # spread across seeds as the influence-function SD of the certification
  # predicts (the SD of 30 values is known to about 13 %), with and without
  # motor noise, which the certification convolves back through its kernel
  f <- corr_similarity(36)
  rho <- corr_correlation(36, scale = 30)
  grid <- corr_grid(4)
  for (motor in c(0, 1.5)) {
    builds <- lapply(1:30, function(seed) bmm:::.tcc_build_noise_table(rho, f, 6000, motor, grid, seed))
    info <- builds[[1]]$info
    at <- round((info$rho_d - grid$lo) / grid$step) + 2
    mult <- c(1, rep(2, 17), 1)
    logs <- lapply(builds, function(b) t(b$tab))
    ref <- log(Reduce(`+`, lapply(logs, exp)) / 30)
    slope <- function(L) (L[, at + 1] - L[, at - 1]) / (2 * grid$step)
    info_per_trial <- colSums(mult * exp(ref[, at]) * slope(ref)^2)
    delta <- t(vapply(logs, function(L) {
      colSums(mult * exp(ref[, at]) * (slope(L) - slope(ref)))
    }, numeric(length(at))))
    spread <- apply(delta, 2, stats::sd) / sqrt(info_per_trial)
    predicted <- rowMeans(vapply(builds, function(b) b$info$rho, numeric(length(at))))
    expect_gt(stats::median(spread / predicted), 0.7)
    expect_lt(stats::median(spread / predicted), 1.4)
    }
})

############################################################################# !
# THE PIPELINE                                                           ####
############################################################################# !

corr_mock_fit <- function(model, data = NULL) {
  data <- data %||% data.frame(err = c(0, 0.1, -0.2, 1, 3, -2))
  suppressMessages(bmm(bmf(d ~ 1), data, model, backend = "mock", mock_fit = 1, rename = FALSE))
}

test_that("a correlated table reaches the lpmf as data, and another table leaves the Stan program unchanged", {
  skip_on_cran()
  f <- corr_similarity()
  plain <- tcc("err", similarity = f, d_max = 2)
  plain$other_vars$table <- "r"
  a <- tcc("err", similarity = f, d_max = 2, correlation = corr_correlation(), correlation_draws = 6000)
  b <- tcc("err",
    similarity = f, d_max = 2, correlation = corr_correlation(scale = 5),
    correlation_draws = 12000, motor_sd = 0.03
  )
  fits <- lapply(list(plain, a, b), corr_mock_fit)
  codes <- lapply(fits, brms::stancode)
  # update() and a refit with another table recompile nothing
  expect_identical(codes[[2]], codes[[1]])
  expect_identical(codes[[3]], codes[[1]])
  tabs <- lapply(fits, function(fit) brms::standata(fit)$tcc_tab)
  expect_equal(dim(tabs[[2]]), c(104, 181))
  expect_true(all(is.finite(tabs[[2]])))
  expect_gt(max(abs(tabs[[2]] - tabs[[3]])), 0.1)
  # the target's log P rises with d, so the padding row above the grid is
  # capped at the row below the last node, as in the quadrature table
  last <- nrow(tabs[[2]])
  expect_gt(tabs[[2]][last - 1, 1], tabs[[2]][last - 2, 1])
  expect_equal(tabs[[2]][last, 1], tabs[[2]][last - 2, 1])
  expect_equal(
    unname(tabs[[2]]),
    unname(bmm:::.tcc_noise_table(a, bmm:::.tcc_grid(a))$tab)
  )
})

test_that("a correlated fit stores and prints its certification, and warns past it", {
  skip_on_cran()
  model <- tcc("err",
    similarity = corr_similarity(), d_max = 2, correlation = corr_correlation(),
    correlation_draws = 6000
  )
  info <- bmm:::.tcc_noise_table(model, bmm:::.tcc_grid(model))$info
  expect_lt(info$n_pool, 200)
  fit <- corr_mock_fit(model)
  stored <- fit$bmm$model$other_vars$table_info
  expect_equal(stored[c("draws", "seed", "n_pool")], info[c("draws", "seed", "n_pool")])
  expect_match(
    summary_notes(fit$bmm$model, NULL),
    glue::glue("Monte Carlo table from 6,000 draws \\(seed {info$seed}\\).*certified for {info$n_pool} trials")
  )
  many <- data.frame(err = rep(c(0, 0.1, -0.2, 1, 3, -2), 50))
  msg <- expect_warning(
    corr_mock_fit(model, many),
    glue::glue("certified for {info$n_pool} trials.*hold 300 trials.*correlation_draws = ")
  )
  # the suggested draws, scaled in proportion, certify the 300 trials with
  # 10 % to spare (the n_pool of another seed differs by a few %), and are
  # the smallest such number with two significant digits
  needed <- as.numeric(gsub(",", "", sub(".*correlation_draws = ([0-9,]+).*", "\\1", conditionMessage(msg))))
  expect_gte(needed * info$n_pool / 6000, 1.1 * 300)
  expect_lt((needed - 10^(floor(log10(needed)) - 1)) * info$n_pool / 6000, 1.1 * 300)
  expect_null(summary_notes(tcc("err", similarity = corr_similarity()), NULL))
  # round counts print in full, as the warning asks users to type them
  round_info <- list(draws = 3e6, seed = 1L, n_pool = 1e5, zeros = 0L, zeros_from = NA_real_)
  expect_match(summary_notes(list(other_vars = list(table_info = round_info)) |>
                               structure(class = c("bmmodel", "tcc")), NULL),
               "3,000,000 draws.*100,000 trials")
  expect_warning(bmm:::.tcc_check_certified(round_info, 2e5, 30), "from 3,000,000 draws")
  # 1e6 draws certified for 10,000 trials, 10,040 in the data: 1.1 x 1,004,000
  # is rounded up to 1,200,000; to the nearest, 1,100,000 would certify only
  # 11,000 of the 11,044 trials the margin asks for
  expect_warning(
    bmm:::.tcc_check_certified(list(draws = 1e6, n_pool = 1e4, zeros_from = NA_real_), 10040, 30),
    "correlation_draws = 1,200,000"
  )
})

test_that("a weak correlation leaves cells no draw reaches, which the fit reports", {
  skip_on_cran()
  # with C = I the drift integral is narrow and far responses at large d are
  # never drawn (spike v2: from d = 10 at 2e5 draws)
  model <- tcc("err",
    similarity = corr_similarity(), d_max = 12, correlation = c(1, rep(0, 180)),
    correlation_draws = 6000
  )
  built <- bmm:::.tcc_noise_table(model, bmm:::.tcc_grid(model))
  info <- built$info
  expect_gt(info$zeros, 0)
  expect_lt(info$zeros_from, 12)
  # a floor of half a draw keeps every row a pmf to within the floored mass
  expect_true(all(is.finite(built$tab)))
  mass <- drop(exp(built$tab[-nrow(built$tab), ]) %*% c(1, rep(2, 179), 1))
  expect_lt(max(abs(mass - 1)), 360 * 0.5 / 6000)
  expect_warning(
    bmm:::.tcc_check_certified(info, 10, 12),
    glue::glue("reached {info$zeros} cells.*first at d = {signif(info$zeros_from, 3)}")
  )
  expect_silent(bmm:::.tcc_check_certified(info, 10, info$zeros_from - 0.1))
})

test_that("motor noise is converted from the response's unit to option steps", {
  radians <- tcc("err", similarity = corr_similarity(), motor_sd = 2 * pi / 360 * 1.5)
  bins <- tcc("err", similarity = corr_similarity(180), n_options = 180, units = "bins", motor_sd = 1.5)
  expect_equal(bmm:::.tcc_motor_steps(radians$other_vars), 1.5)
  expect_equal(bmm:::.tcc_motor_steps(bins$other_vars), 1.5)
})

test_that("the motor-noise table is the exact pmf convolved with the rounding kernel", {
  skip_on_cran()
  f <- corr_similarity()
  model <- tcc("err", similarity = f, d_max = 2, motor_sd = deg2rad(2))
  tab <- bmm:::.tcc_noise_table(model, bmm:::.tcc_grid(model))$tab
  d_values <- corr_grid_d(bmm:::.tcc_grid(model))
  for (row in c(2, 51, 101)) {
    x <- 2 * pi * (0:180) / 360
    expect_equal(tab[row, ], dtcc(x, d_values[row], f, motor_sd = deg2rad(2), log = TRUE),
      tolerance = 1e-12
    )
  }
})

############################################################################# !
# DTCC AND RTCC                                                          ####
############################################################################# !

test_that("dtcc(motor_sd) is a pmf, the convolution of the noise-free pmf, and the pmf of rtcc(motor_sd)", {
  skip_on_cran()
  f <- corr_similarity()
  x <- 2 * pi * (-179:180) / 360
  p <- dtcc(x, 2, f, motor_sd = deg2rad(2))
  plain <- dtcc(x, 2, f)
  # the convolution keeps the mass the quadrature gives (1 within 1e-8)
  expect_equal(sum(p), sum(plain), tolerance = 1e-13)
  expect_equal(sum(p), 1, tolerance = 1e-8)
  expect_equal(dtcc(x, 2, f, motor_sd = 0), plain)
  kernel <- bmm:::.tcc_motor_kernel(2)
  conv <- vapply(seq_along(x), function(i) {
    sum(kernel * plain[((i - 1 + (-12:12)) %% 360) + 1])
  }, numeric(1))
  expect_equal(p, conv, tolerance = 1e-12)
  withr::local_seed(4)
  r <- rtcc(1e5, 2, f, motor_sd = deg2rad(2))
  observed <- tabulate(round(r * 360 / (2 * pi)) + 180, 360)
  expected <- p * 1e5
  keep <- expected > 5
  chi2 <- sum((observed[keep] - expected[keep])^2 / expected[keep])
  expect_gt(stats::pchisq(chi2, sum(keep) - 1, lower.tail = FALSE), 0.001)
  expect_error(dtcc(x, 2, f, motor_sd = 0.1, noise = "gumbel"), "requires noise = 'gaussian'")
})

test_that("rtcc(correlation) draws from the distribution of the correlated table", {
  skip_on_cran()
  withr::local_seed(5)
  f <- corr_similarity(36)
  rho <- corr_correlation(36, scale = 30)
  grid <- corr_grid(3)
  p <- exp(bmm:::.tcc_build_noise_table(rho, f, 1.2e5, 0, grid, seed = 3)$tab[102, ])
  r <- rtcc(1e5, 2, f, n_options = 36, correlation = rho)
  steps <- round(r * 36 / (2 * pi))
  expect_true(all(steps > -18 & steps <= 18))
  observed <- tabulate(abs(steps) + 1, 19)
  expected <- p * c(1, rep(2, 17), 1) * 1e5
  keep <- expected > 5
  chi2 <- sum((observed[keep] - expected[keep])^2 / expected[keep])
  expect_gt(stats::pchisq(chi2, sum(keep) - 1, lower.tail = FALSE), 0.001)
  expect_equal(is.na(rtcc(3, c(1, NA, 2), f, n_options = 36, correlation = rho)), c(FALSE, TRUE, FALSE))
  expect_error(
    rtcc(5, 2, f, n_options = 36, correlation = rho, options = c(0, pi)),
    "planned for bmm 1.6.0"
  )
})
