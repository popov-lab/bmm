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

.stan_lba_simple_loglik <- function(rt, response, driftc, drifte, gap, sp, s,
                                    distribution) {
  b <- gap + sp
  A <- sp
  single_lpdf <- switch(distribution,
    normal = .stan_lba_normal_single_lpdf
  )
  single_lccdf <- switch(distribution,
    normal = .stan_lba_normal_single_lccdf
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
