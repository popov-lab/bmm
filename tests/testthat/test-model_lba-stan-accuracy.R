.stan_lba_softplus <- function(x, scale = 1e-12) {
  scaled <- x / scale
  ifelse(scaled > 30, x, scale * log1p(exp(scaled)))
}

.stan_lba_log_positive <- function(x) {
  log(.stan_lba_softplus(x))
}

.stan_lba_log_diff_exp <- function(log_x, log_y) {
  if (log_x <= log_y) {
    return(-Inf)
  }
  log_x + log1p(-exp(log_y - log_x))
}

.stan_lba_log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

.stan_lba_normal_single_lpdf <- function(t, v, b, A, s) {
  hi <- b / t
  lo <- (b - A) / t
  z_hi <- (hi - v) / s
  z_lo <- (lo - v) / s
  Phi_hi <- stats::pnorm(z_hi)
  Phi_lo <- stats::pnorm(z_lo)
  phi_hi <- stats::dnorm(z_hi)
  phi_lo <- stats::dnorm(z_lo)
  M <- v * (Phi_hi - Phi_lo) + s * (phi_lo - phi_hi)
  .stan_lba_log_positive(M) - log(A)
}

.stan_lba_normal_single_lccdf <- function(t, v, b, A, s) {
  hi <- b / t
  lo <- (b - A) / t
  z_hi <- (hi - v) / s
  z_lo <- (lo - v) / s
  Phi_hi <- stats::pnorm(z_hi)
  Phi_lo <- stats::pnorm(z_lo)
  phi_hi <- stats::dnorm(z_hi)
  phi_lo <- stats::dnorm(z_lo)
  M <- v * (Phi_hi - Phi_lo) + s * (phi_lo - phi_hi)
  surv_num <- (b * Phi_hi) - ((b - A) * Phi_lo) - (t * M)
  .stan_lba_log_positive(surv_num) - log(A)
}

.stan_lba_gamma_single_lpdf <- function(t, v, b, A, s) {
  hi <- b / t
  lo <- (b - A) / t
  log_M <- log(v) - log(s) + .stan_lba_log_diff_exp(
    stats::pgamma(hi, shape = v + 1, rate = s, log.p = TRUE),
    stats::pgamma(lo, shape = v + 1, rate = s, log.p = TRUE)
  )
  log_M - log(A)
}

.stan_lba_gamma_single_lccdf <- function(t, v, b, A, s) {
  hi <- b / t
  lo <- (b - A) / t
  log_M <- log(v) - log(s) + .stan_lba_log_diff_exp(
    stats::pgamma(hi, shape = v + 1, rate = s, log.p = TRUE),
    stats::pgamma(lo, shape = v + 1, rate = s, log.p = TRUE)
  )
  log_u <- .stan_lba_log_diff_exp(
    log(b) + stats::pgamma(hi, shape = v, rate = s, log.p = TRUE),
    log(b - A) + stats::pgamma(lo, shape = v, rate = s, log.p = TRUE)
  )
  surv_num <- exp(log_u) - t * exp(log_M)
  .stan_lba_log_positive(surv_num) - log(A)
}

.stan_lba_lognormal_single_lpdf <- function(t, v, b, A, s) {
  hi <- b / t
  lo <- (b - A) / t
  s2 <- s^2
  log_M <- v + s2 / 2 + .stan_lba_log_diff_exp(
    stats::pnorm((log(hi) - v - s2) / s, log.p = TRUE),
    stats::pnorm((log(lo) - v - s2) / s, log.p = TRUE)
  )
  log_M - log(A)
}

.stan_lba_lognormal_single_lccdf <- function(t, v, b, A, s) {
  hi <- b / t
  lo <- (b - A) / t
  s2 <- s^2
  log_M <- v + s2 / 2 + .stan_lba_log_diff_exp(
    stats::pnorm((log(hi) - v - s2) / s, log.p = TRUE),
    stats::pnorm((log(lo) - v - s2) / s, log.p = TRUE)
  )
  log_u <- .stan_lba_log_diff_exp(
    log(b) + stats::pnorm((log(hi) - v) / s, log.p = TRUE),
    log(b - A) + stats::pnorm((log(lo) - v) / s, log.p = TRUE)
  )
  surv_num <- exp(log_u) - t * exp(log_M)
  .stan_lba_log_positive(surv_num) - log(A)
}

.stan_lba_frechet_log_cdf <- function(x, v, s) {
  -(x / s)^(-v)
}

.stan_lba_frechet_log_M <- function(t, v, b, A, s) {
  nodes <- c(
    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318,
    -0.7554044083550030, -0.6178762444026438, -0.4580167776572274,
    -0.2816035507792589, -0.0950125098376374,  0.0950125098376374,
     0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
     0.7554044083550030,  0.8656312023878318,  0.9445750230732326,
     0.9894009349916499
  )
  weights <- c(
    0.0271524594117541, 0.0622535239386479, 0.0951585116824928,
    0.1246289712555339, 0.1495959888165767, 0.1691565193950025,
    0.1826034150449236, 0.1894506104550685, 0.1894506104550685,
    0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
    0.1246289712555339, 0.0951585116824928, 0.0622535239386479,
    0.0271524594117541
  )
  lo <- (b - A) / t
  hi <- b / t
  mid <- 0.5 * (hi + lo)
  half_range <- 0.5 * (hi - lo)
  log_terms <- vapply(seq_along(nodes), function(j) {
    u <- mid + half_range * nodes[j]
    log_z <- log(u / s)
    log(weights[j]) + log(v) - (v * log_z) - exp(-v * log_z)
  }, numeric(1))
  log(half_range) + .stan_lba_log_sum_exp(log_terms)
}

.stan_lba_frechet_single_lpdf <- function(t, v, b, A, s) {
  .stan_lba_frechet_log_M(t, v, b, A, s) - log(A)
}

.stan_lba_frechet_single_lccdf <- function(t, v, b, A, s) {
  log_M <- .stan_lba_frechet_log_M(t, v, b, A, s)
  hi <- b / t
  lo <- (b - A) / t
  log_u <- .stan_lba_log_diff_exp(
    log(b) + .stan_lba_frechet_log_cdf(hi, v, s),
    log(b - A) + .stan_lba_frechet_log_cdf(lo, v, s)
  )
  surv_num <- exp(log_u) - t * exp(log_M)
  .stan_lba_log_positive(surv_num) - log(A)
}

.stan_lba_simple_loglik <- function(rt, response, driftc, drifte, gap, sp, s,
                                    distribution) {
  b <- gap + sp
  A <- sp
  single_lpdf <- switch(distribution,
    normal = .stan_lba_normal_single_lpdf,
    gamma = .stan_lba_gamma_single_lpdf,
    lognormal = .stan_lba_lognormal_single_lpdf,
    frechet = .stan_lba_frechet_single_lpdf
  )
  single_lccdf <- switch(distribution,
    normal = .stan_lba_normal_single_lccdf,
    gamma = .stan_lba_gamma_single_lccdf,
    lognormal = .stan_lba_lognormal_single_lccdf,
    frechet = .stan_lba_frechet_single_lccdf
  )

  if (response == 1L) {
    single_lpdf(rt, driftc, b, A, s) + single_lccdf(rt, drifte, b, A, s)
  } else {
    single_lpdf(rt, drifte, b, A, s) + single_lccdf(rt, driftc, b, A, s)
  }
}

.ref_lba_simple_loglik <- function(rt, response, driftc, drifte, gap, sp, s,
                                   distribution) {
  b <- gap + sp
  A <- sp
  if (response == 1L) {
    .dlba_single(rt, driftc, b, A, s, distribution, log = TRUE) +
      log1p(-.plba_single(rt, drifte, b, A, s, distribution))
  } else {
    .dlba_single(rt, drifte, b, A, s, distribution, log = TRUE) +
      log1p(-.plba_single(rt, driftc, b, A, s, distribution))
  }
}

.lba_accuracy_grid <- function(n, distribution) {
  gap <- exp(stats::rnorm(n, -0.5, 0.35))
  sp <- exp(stats::rnorm(n, -1, 0.35))
  rt <- stats::runif(n, 0.08, 1.5)
  response <- sample(c(1L, 2L), n, replace = TRUE, prob = c(0.75, 0.25))
  s <- rep(1, n)
  drift <- switch(distribution,
    normal = cbind(
      stats::rnorm(n, 3, 0.35),
      stats::rnorm(n, 1, 0.35)
    ),
    gamma = cbind(
      stats::runif(n, 1.5, 3.5),
      stats::runif(n, 0.8, 2.5)
    ),
    lognormal = cbind(
      stats::rnorm(n, 0.5, 0.15),
      stats::rnorm(n, 0.1, 0.15)
    ),
    frechet = cbind(
      stats::runif(n, 1.2, 2.8),
      stats::runif(n, 0.8, 2.0)
    )
  )

  data.frame(
    rt = rt,
    response = response,
    driftc = drift[, 1],
    drifte = drift[, 2],
    gap = gap,
    sp = sp,
    s = s
  )
}

.check_lba_accuracy <- function(distribution, grid_n) {
  dat <- .lba_accuracy_grid(grid_n, distribution)

  ref <- mapply(
    .ref_lba_simple_loglik,
    dat$rt, dat$response, dat$driftc, dat$drifte, dat$gap, dat$sp, dat$s,
    MoreArgs = list(distribution = distribution)
  )
  stan_mirror <- mapply(
    .stan_lba_simple_loglik,
    dat$rt, dat$response, dat$driftc, dat$drifte, dat$gap, dat$sp, dat$s,
    MoreArgs = list(distribution = distribution)
  )

  diff <- abs(stan_mirror - ref)
  diff <- diff[is.finite(diff)]

  list(
    ref = ref,
    stan_mirror = stan_mirror,
    median = stats::median(diff),
    p95 = unname(stats::quantile(diff, 0.95)),
    max = max(diff)
  )
}

test_that("normal LBA Stan algebra stays close to the exact reference", {
  set.seed(101)
  acc <- .check_lba_accuracy("normal", 300)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$median, 1e-3)
  expect_lt(acc$p95, 1e-2)
})

test_that("gamma LBA Stan algebra matches the exact reference", {
  set.seed(102)
  acc <- .check_lba_accuracy("gamma", 300)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$max, 1e-10)
})

test_that("lognormal LBA Stan algebra stays close to the exact reference", {
  set.seed(103)
  acc <- .check_lba_accuracy("lognormal", 300)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$median, 1e-3)
  expect_lt(acc$p95, 1e-2)
})

test_that("frechet LBA Stan algebra stays close to the exact reference", {
  set.seed(104)
  acc <- .check_lba_accuracy("frechet", 100)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$p95, 1e-2)
})
