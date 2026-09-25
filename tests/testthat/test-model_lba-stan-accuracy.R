# Transcriptions of the CURRENT inst/stan_chunks/lba_*_functions.stan algebra,
# written independently of R/distributions.R's .lba_*_lpdf()/.lba_*_lsurv()
# twins so that a bug shared by both (e.g. in the derivation itself) is not
# hidden by comparing a function against itself. .ref_lba_simple_loglik()
# below is the production reference (.dlba_single()/.plba_single()).

.stan_lba_log_clip <- function(x) {
  log(max(x, 1e-300))
}

.stan_lba_log_floor <- function() log(1e-300)

.stan_lba_log_diff_exp <- function(log_x, log_y) {
  if (log_x <= log_y) return(-Inf)
  log_x + log1p(-exp(log_y - log_x))
}

.stan_lba_log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

.stan_lba_log_sum_exp2 <- function(a, b) .stan_lba_log_sum_exp(c(a, b))

.stan_lba_log_Phi_diff <- function(lo, hi) {
  if (lo + hi >= 0) {
    return(.stan_lba_log_diff_exp(stats::pnorm(-lo, log.p = TRUE), stats::pnorm(-hi, log.p = TRUE)))
  }
  .stan_lba_log_diff_exp(stats::pnorm(hi, log.p = TRUE), stats::pnorm(lo, log.p = TRUE))
}

.stan_lba_log_phi_int_narrow <- function(z_lo, z_hi) {
  dz <- z_hi - z_lo
  z_m <- 0.5 * (z_lo + z_hi)
  log(dz) + stats::dnorm(z_m, log = TRUE) + log1p(dz^2 / 24 * (z_m^2 - 1))
}

.stan_lba_normal_g <- function(z) z * stats::pnorm(z) + stats::dnorm(z)
.stan_lba_normal_h <- function(z) z * stats::pnorm(-z) - stats::dnorm(z)

.stan_lba_normal_log_M <- function(t, v, b, A, s) {
  z_lo <- ((b - A) / t - v) / s
  z_hi <- (b / t - v) / s
  delta <- A / (t * s)
  if (delta < 1e-4) {
    z_m <- 0.5 * (z_lo + z_hi)
    u_m <- (b - 0.5 * A) / t
    bracket <- u_m + delta^2 / 24 * (u_m * (z_m^2 - 1) - 2 * s * z_m)
    if (bracket <= 0) return(.stan_lba_log_floor())
    return(log(delta) + stats::dnorm(z_m, log = TRUE) + log(bracket))
  }
  log_dPhi <- .stan_lba_log_Phi_diff(z_lo, z_hi)
  sum_z <- z_lo + z_hi
  x <- 0.5 * delta * abs(sum_z)
  log_dphi <- if (x > 0) {
    log(s) + stats::dnorm(if (sum_z >= 0) z_lo else z_hi, log = TRUE) + log1p(-exp(-x))
  } else {
    -Inf
  }
  if (v > 0) {
    l1 <- log(v) + log_dPhi
    if (sum_z >= 0) return(.stan_lba_log_sum_exp2(l1, log_dphi))
    return(if (l1 > log_dphi) .stan_lba_log_diff_exp(l1, log_dphi) else .stan_lba_log_floor())
  }
  if (v < 0) {
    l1 <- log(-v) + log_dPhi
    return(if (log_dphi > l1) .stan_lba_log_diff_exp(log_dphi, l1) else .stan_lba_log_floor())
  }
  log_dphi
}

.stan_lba_normal_single_lpdf <- function(t, v, b, A, s) {
  .stan_lba_normal_log_M(t, v, b, A, s) - log(A) - stats::pnorm(v / s, log.p = TRUE)
}

.stan_lba_normal_single_lccdf <- function(t, v, b, A, s) {
  z_lo <- ((b - A) / t - v) / s
  z_hi <- (b / t - v) / s
  delta <- A / (t * s)
  log_denom <- stats::pnorm(v / s, log.p = TRUE)
  if (delta < 1e-4) {
    z_m <- 0.5 * (z_lo + z_hi)
    corr <- delta^2 / 24 * z_m * stats::dnorm(z_m)
    num <- (if (v >= 0) stats::pnorm(z_m) + expm1(log_denom) else exp(log_denom) - stats::pnorm(-z_m)) - corr
  } else if (v >= 0) {
    num <- (.stan_lba_normal_g(z_hi) - .stan_lba_normal_g(z_lo)) / delta + expm1(log_denom)
  } else {
    num <- exp(log_denom) - (.stan_lba_normal_h(z_hi) - .stan_lba_normal_h(z_lo)) / delta
  }
  .stan_lba_log_clip(num) - log_denom
}

.stan_lba_gamma_log_dF <- function(lo, hi, alpha, beta) {
  du <- hi - lo
  u_m <- 0.5 * (lo + hi)
  if (du * (abs(alpha - 1) / u_m + beta) < 1e-3) {
    r <- (alpha - 1) / u_m - beta
    return(log(du) + stats::dgamma(u_m, shape = alpha, rate = beta, log = TRUE) +
      log1p(du^2 / 24 * (r^2 - (alpha - 1) / u_m^2)))
  }
  if (lo * beta > alpha) {
    return(.stan_lba_log_diff_exp(
      stats::pgamma(lo, shape = alpha, rate = beta, lower.tail = FALSE, log.p = TRUE),
      stats::pgamma(hi, shape = alpha, rate = beta, lower.tail = FALSE, log.p = TRUE)
    ))
  }
  .stan_lba_log_diff_exp(
    stats::pgamma(hi, shape = alpha, rate = beta, log.p = TRUE),
    stats::pgamma(lo, shape = alpha, rate = beta, log.p = TRUE)
  )
}

.stan_lba_gamma_single_lpdf <- function(t, v, b, A, s) {
  log(v) - log(s) + .stan_lba_gamma_log_dF((b - A) / t, b / t, v + 1, s) - log(A)
}

.stan_lba_gamma_single_lccdf <- function(t, v, b, A, s) {
  lo <- (b - A) / t
  hi <- b / t
  log_tM <- log(t) + log(v) - log(s) + .stan_lba_gamma_log_dF(lo, hi, v + 1, s)
  log_u <- .stan_lba_log_sum_exp2(
    log(A) + stats::pgamma(lo, shape = v, rate = s, log.p = TRUE),
    log(b) + .stan_lba_gamma_log_dF(lo, hi, v, s)
  )
  (if (log_u > log_tM) .stan_lba_log_diff_exp(log_u, log_tM) else .stan_lba_log_floor()) - log(A)
}

.stan_lba_lognormal_log_dPhi <- function(z_lo, z_hi, dz) {
  if (dz < 1e-4) return(.stan_lba_log_phi_int_narrow(z_lo, z_hi))
  .stan_lba_log_Phi_diff(z_lo, z_hi)
}

.stan_lba_lognormal_single_lpdf <- function(t, v, b, A, s) {
  z_hi <- (log(b / t) - v) / s
  dz <- -log1p(-A / b) / s
  v + 0.5 * s^2 + .stan_lba_lognormal_log_dPhi(z_hi - dz - s, z_hi - s, dz) - log(A)
}

.stan_lba_lognormal_single_lccdf <- function(t, v, b, A, s) {
  z_hi <- (log(b / t) - v) / s
  dz <- -log1p(-A / b) / s
  z_lo <- z_hi - dz
  log_tM <- log(t) + v + 0.5 * s^2 + .stan_lba_lognormal_log_dPhi(z_lo - s, z_hi - s, dz)
  log_u <- .stan_lba_log_sum_exp2(
    log(A) + stats::pnorm(z_lo, log.p = TRUE),
    log(b) + .stan_lba_lognormal_log_dPhi(z_lo, z_hi, dz)
  )
  (if (log_u > log_tM) .stan_lba_log_diff_exp(log_u, log_tM) else .stan_lba_log_floor()) - log(A)
}

.stan_lba_frechet_log_F <- function(x, v, s) -(x / s)^(-v)

.stan_lba_frechet_log_dF <- function(lo, hi, v, s) {
  lF_lo <- .stan_lba_frechet_log_F(lo, v, s)
  lF_hi <- .stan_lba_frechet_log_F(hi, v, s)
  if (lF_lo > -log(2)) {
    return(.stan_lba_log_diff_exp(log1p(-exp(lF_lo)), log1p(-exp(lF_hi))))
  }
  .stan_lba_log_diff_exp(lF_hi, lF_lo)
}

.stan_lba_frechet_gl16_nodes <- c(
  -0.9894009349916499, -0.9445750230732326, -0.8656312023878318,
  -0.7554044083550030, -0.6178762444026438, -0.4580167776572274,
  -0.2816035507792589, -0.0950125098376374,  0.0950125098376374,
   0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
   0.7554044083550030,  0.8656312023878318,  0.9445750230732326,
   0.9894009349916499
)

.stan_lba_frechet_gl16_weights <- c(
  0.0271524594117541, 0.0622535239386479, 0.0951585116824928,
  0.1246289712555339, 0.1495959888165767, 0.1691565193950025,
  0.1826034150449236, 0.1894506104550685, 0.1894506104550685,
  0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
  0.1246289712555339, 0.0951585116824928, 0.0622535239386479,
  0.0271524594117541
)

.stan_lba_frechet_log_M <- function(t, v, b, A, s) {
  lo <- (b - A) / t
  hi <- b / t
  mid <- 0.5 * (hi + lo)
  half_range <- 0.5 * (hi - lo)
  log_terms <- vapply(seq_along(.stan_lba_frechet_gl16_nodes), function(j) {
    u <- mid + half_range * .stan_lba_frechet_gl16_nodes[j]
    log_z <- log(u / s)
    log_integrand <- log(v) - (v * log_z) - exp(-v * log_z)
    log(.stan_lba_frechet_gl16_weights[j]) + log_integrand
  }, numeric(1))
  log(half_range) + .stan_lba_log_sum_exp(log_terms)
}

.stan_lba_frechet_single_lpdf <- function(t, v, b, A, s) {
  .stan_lba_frechet_log_M(t, v, b, A, s) - log(A)
}

.stan_lba_frechet_single_lccdf <- function(t, v, b, A, s) {
  lo <- (b - A) / t
  hi <- b / t
  log_tM <- log(t) + .stan_lba_frechet_log_M(t, v, b, A, s)
  log_u <- .stan_lba_log_sum_exp2(
    log(A) + .stan_lba_frechet_log_F(lo, v, s),
    log(b) + .stan_lba_frechet_log_dF(lo, hi, v, s)
  )
  (if (log_u > log_tM) .stan_lba_log_diff_exp(log_u, log_tM) else .stan_lba_log_floor()) - log(A)
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
  withr::local_seed(101)
  acc <- .check_lba_accuracy("normal", 300)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$median, 1e-8)
  expect_lt(acc$p95, 1e-6)
})

test_that("gamma LBA Stan algebra matches the exact reference", {
  withr::local_seed(102)
  acc <- .check_lba_accuracy("gamma", 300)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$max, 1e-8)
})

test_that("lognormal LBA Stan algebra stays close to the exact reference", {
  withr::local_seed(103)
  acc <- .check_lba_accuracy("lognormal", 300)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$median, 1e-8)
  expect_lt(acc$p95, 1e-6)
})

test_that("frechet LBA Stan algebra stays close to the exact reference", {
  withr::local_seed(104)
  acc <- .check_lba_accuracy("frechet", 100)

  expect_true(all(is.finite(acc$stan_mirror)))
  expect_lt(acc$p95, 1e-4)
})
