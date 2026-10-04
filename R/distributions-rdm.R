############################################################################# !
# RACING DIFFUSION MODEL (RDM) DISTRIBUTION FUNCTIONS                   ####
############################################################################# !

#' @title Distribution functions for the Racing Diffusion Model (RDM)
#'
#' @description Density, random generation, CDF, and quantile functions for the
#'   Racing Diffusion Model (Tillman, Van Zandt, & Logan, 2020). The RDM is a
#'   multi-accumulator race model where each accumulator follows a Wald (inverse
#'   Gaussian) distribution. The first accumulator to finish determines the
#'   response and the RT.
#'
#' @name rdm_dist
#'
#' @param rt Numeric vector of response times (in seconds).
#' @param q Numeric vector of quantiles (response times in seconds).
#' @param p Numeric vector of probabilities.
#' @param n Number of observations to generate.
#' @param response Integer vector indicating which accumulator won (1-indexed,
#'   in `1:length(drift)`).
#' @param drift Numeric vector of drift rates, one per accumulator (all > 0).
#'   Its length sets the number of accumulators in the race.
#' @param gap Threshold gap (> 0). The distance between the maximum starting
#'   point and the decision threshold. The total threshold is computed as
#'   `b = gap + sp`, ensuring `b > sp` structurally.
#' @param ndt Non-decision time in seconds (>= 0).
#' @param s Diffusion constant (> 0), default = 1.
#' @param sp Maximum starting point (>= 0). Starting evidence is uniformly
#'   distributed on `[0, sp]`. Default = 0 (no starting point variability).
#' @param log Logical; if `TRUE`, values are returned on the log scale.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are P(X <= x).
#' @param log.p Logical; if `TRUE`, probabilities are given as log(p).
#'
#' @details `gap`, `ndt`, `s` and `sp` are recycled along `rt` (or `q`) in
#'   `drdm()` and `prdm()`, so trial-varying parameters can be passed as
#'   vectors. `rrdm()` and `qrdm()` take one value of each; `drift` is always
#'   one value per accumulator.
#'
#' @return
#'   - `drdm()` returns a numeric vector of (log-)densities.
#'   - `rrdm()` returns a data.frame with columns `rt` and `response`.
#'   - `prdm()` returns a numeric vector of (log-)probabilities.
#'   - `qrdm()` returns a numeric vector of quantiles (response times).
#'
#' @references
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling
#'   models without random between-trial variability: the racing diffusion
#'   model of speeded decision making. Psychonomic Bulletin & Review, 27,
#'   911-936.
#'
#' @keywords distribution
#'
#' @examples
#' dat <- rrdm(n = 1000, drift = c(3, 1.5), gap = 1, sp = 0, ndt = 0.2)
#' head(dat)
#' hist(dat$rt)
#'
#' # with starting point variability
#' dat2 <- rrdm(n = 1000, drift = c(3, 1.5), gap = 0.7, sp = 0.3, ndt = 0.2)
#' @export
drdm <- function(rt, response, drift, gap, ndt, s = 1, sp = 0,
                 log = FALSE) {
  validate_rdm_parameters(drift, gap, ndt, s, sp)
  K <- length(drift)
  stopif(
    any(is.na(response)) || any(response < 1 | response > K | response != round(response)),
    "response must index one of the {K} accumulators (integers in 1:{K})."
  )

  n <- max(length(rt), length(response), length(gap), length(ndt), length(s), length(sp))
  log_lik <- .rdm_race_lpdf(
    t = rep_len(rt, n) - rep_len(ndt, n),
    response = as.integer(rep_len(response, n)),
    drift = matrix(drift, n, K, byrow = TRUE),
    counts = matrix(1L, n, K),
    gap = rep_len(gap, n),
    A = rep_len(sp, n),
    s = rep_len(s, n)
  )
  if (log) log_lik else exp(log_lik)
}

#' @rdname rdm_dist
#' @export
rrdm <- function(n, drift, gap, ndt, s = 1, sp = 0) {
  validate_rdm_parameters(drift, gap, ndt, s, sp)
  stopif(
    any(lengths(list(gap, ndt, s, sp)) != 1),
    "rrdm() takes a single value for each of gap, ndt, s and sp."
  )
  K <- length(drift)
  race <- .rdm_race(
    drift = matrix(drift, n, K, byrow = TRUE),
    gap = rep(gap, n), A = rep(sp, n), s = rep(s, n),
    counts = matrix(1L, n, K)
  )
  data.frame(rt = race$rt + ndt, response = race$response)
}

#' @rdname rdm_dist
#' @export
prdm <- function(q, drift, gap, ndt, s = 1, sp = 0,
                 lower.tail = TRUE, log.p = FALSE) {
  validate_rdm_parameters(drift, gap, ndt, s, sp)
  n <- max(length(q), length(gap), length(ndt), length(s), length(sp))
  t <- rep_len(q, n) - rep_len(ndt, n)
  gap <- rep_len(gap, n)
  A <- rep_len(sp, n)
  s <- rep_len(s, n)

  log_surv <- numeric(n)
  for (j in seq_along(drift)) {
    log_surv <- log_surv +
      .pwald_full(t, drift = drift[j], bound = gap + A, A = A, s = s,
                  lower.tail = FALSE, log.p = TRUE)
  }
  log_p <- if (lower.tail) log1m_exp(log_surv) else log_surv
  if (log.p) log_p else exp(log_p)
}

#' @rdname rdm_dist
#' @export
qrdm <- function(p, drift, gap, ndt, s = 1, sp = 0,
                 lower.tail = TRUE, log.p = FALSE) {
  validate_rdm_parameters(drift, gap, ndt, s, sp)
  stopif(
    any(lengths(list(gap, ndt, s, sp)) != 1),
    "qrdm() takes a single value for each of gap, ndt, s and sp."
  )
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  cdf <- function(q) prdm(q, drift = drift, gap = gap, ndt = ndt, s = s, sp = sp)
  vapply(p, function(pi) {
    if (pi <= 0) return(ndt)
    if (pi >= 1) return(Inf)
    # the slowest accumulator's mean bounds the race from above only loosely,
    # so the bracket is widened until it contains the quantile
    upper <- ndt + max((gap + sp) / drift) * 5
    while (cdf(upper) < pi) upper <- ndt + (upper - ndt) * 2
    stats::uniroot(function(q) cdf(q) - pi, interval = c(ndt + 1e-10, upper),
                   tol = 1e-8)$root
  }, numeric(1))
}

validate_rdm_parameters <- function(drift, gap, ndt, s, sp) {
  stopif(any(drift <= 0), "drift rates must be positive.")
  stopif(any(gap <= 0), "gap (threshold gap) must be positive.")
  stopif(any(ndt < 0), "ndt (non-decision time) must be non-negative.")
  stopif(any(s <= 0), "s (diffusion constant) must be positive.")
  stopif(any(sp < 0), "sp (maximum starting point) must be non-negative.")
}

# Log-likelihood of one race per row: the winner's density times every other
# accumulator's survival, with n_j identical accumulators per category
# (log n_win for the winner's category, n_win - 1 survival copies for it).
# `drift` and `counts` are row-per-trial matrices, so per-draw parameters line
# up with the trials by row and every branch below subsets by the same mask.
.rdm_race_lpdf <- function(t, response, drift, counts, gap, A, s) {
  n <- length(t)
  K <- ncol(drift)
  win <- cbind(seq_len(n), response)
  bound <- gap + A
  log_lik <- log(counts[win]) +
    .dwald_full(t, drift = drift[win], bound = bound, A = A, s = s, log = TRUE)
  for (j in seq_len(K)) {
    # a zero-count or underflowed survival must stay out of the sum:
    # 0 * -Inf would poison the likelihood with NaN
    reps <- counts[, j] - (response == j)
    idx <- reps > 0 & t > 0
    if (!any(idx)) next
    log_lik[idx] <- log_lik[idx] + reps[idx] *
      .pwald_full(t[idx], drift = drift[idx, j], bound = bound[idx], A = A[idx],
                  s = s[idx], lower.tail = FALSE, log.p = TRUE)
  }
  log_lik[t <= 0] <- -Inf
  log_lik
}

# One race per row from the same simulator that pp_simulate() and
# posterior_predict() use: every accumulator draws its own uniform start point
# on [0, A] and a Wald finishing time from the remaining distance. `drift` and
# `counts` are row-per-trial matrices; a category with n_j accumulators draws
# n_j finishing times and keeps the fastest.
.rdm_race <- function(drift, gap, A, s, counts) {
  n <- nrow(drift)
  K <- ncol(drift)
  ft <- matrix(Inf, n, K)
  for (j in seq_len(K)) {
    for (k in seq_len(max(counts[, j]))) {
      active <- counts[, j] >= k
      m <- sum(active)
      start <- stats::runif(m, 0, A[active])
      ft[active, j] <- pmin(
        ft[active, j],
        .rwald_ig(m, drift = drift[active, j], bound = gap[active] + A[active] - start,
                  s = s[active])
      )
    }
  }
  list(
    rt = matrixStats::rowMins(ft),
    response = max.col(-ft, ties.method = "first")
  )
}

# Michael-Schucany-Haas algorithm for inverse Gaussian
.rwald_ig <- function(n, drift, bound, s) {
  mu_ig <- bound / drift
  lambda_ig <- (bound / s)^2
  y <- stats::rchisq(n, df = 1)
  x <- mu_ig + (mu_ig^2 * y -
    mu_ig * sqrt(4 * mu_ig * lambda_ig * y + mu_ig^2 * y^2)) /
    (2 * lambda_ig)
  u <- stats::runif(n)
  ifelse(u <= mu_ig / (mu_ig + x), x, mu_ig^2 / x)
}

# Wald with uniform start point on [0, A] and threshold `bound`, i.e. a distance
# to threshold uniform on [bound - A, bound] (Tillman et al., 2020, Eq. 5 and
# Appendix A, generalised to a diffusion scale s). Everything is computed in
# log space so that the tails stay finite where the raw forms underflow, and
# every argument is recycled to a common length so that per-draw parameter
# vectors line up with t. Keep these in step with rdm_functions.stan.
#
# The Stan twins (inst/stan_chunks/rdm_functions.stan) follow the same
# derivation but spell the normal tails for the gradient:
# - An upper normal tail is std_normal_lcdf(-z), never std_normal_lccdf(z).
#   The two are equal on Stan Math 5.4 (CmdStan 2.40), but on 5.3 (rstan /
#   StanHeaders 2.39) the lccdf is -Inf from z = 8.26 with an infinite partial,
#   which reached this model as a NaN trial for a loser with drift 5 from
#   t = 4 s. cogmod 0.3.2 (Makowski) exposed that flaw in bmm's cswald survivor
#   and spells its own RDM with the reflected lcdf; the derivation here is
#   bmm's own. A test pins that the chunk calls no lccdf.
# - log Phi(z) is rdm_log_Phi(). Stan's Phi() evaluates 0.5 (1 + erf(z / sqrt 2))
#   between -5 and 0, which cancels to a relative error of 1e-16 / Phi(z) (five
#   digits gone at z = -4.4), and the differences g(beta) - g(alpha) and
#   Phi(beta) - Phi(alpha) amplify that by their own cancellation (8.5e-7 in a
#   far-tail log survival, measured). erfc has no such cancellation, its value
#   is exact to 1e-16 down to its underflow at z = -37.5, and its derivative is
#   the exact normal density; std_normal_lcdf takes over below, at the cost of
#   its approximate (1e-5 relative) partial there.
# - Every log-difference goes through swald_log_diff_exp(): log_diff_exp(x, x)
#   has infinite partials that poison the gradient of the whole trial even in a
#   branch whose value is discarded.
# - rdm_log_g() sums y Phi(y) + phi(y) directly above y = -1 (g(-1) = 0.083, no
#   cancellation, exact gradient), takes the log-space difference between -12
#   and -1, and uses the series of .rdm_log_g() below -12. Twelve terms and the
#   log-space form agree to about 1e-12 at -12.
# - In rdm_log_surv() the E terms stay single exponents, so that they survive
#   where exp(2 u drift / s^2) and Phi(.) separately overflow and underflow.
#   D1, D2 and P are positive because g increases and the Wald reflection
#   identity signs the derivatives of the survivor's antiderivative pieces.
# Stan Math has no inverse Gaussian yet (stan-dev/math #3382 adds
# inv_gaussian_* in log space); nothing here depends on it.

# Below this ratio of start-point range to the diffusion spread s * sqrt(t)
# the difference quotients over the start point cancel; the plain Wald at the
# midpoint threshold is second-order accurate in A and takes over
.rdm_midpoint_ratio <- 1e-4

# log(y Phi(y) + phi(y)), the antiderivative of Phi, given log Phi(y). For y
# far below zero the two terms cancel to order y^-2, so g is written as
# phi(y) * (1 - |y| M(|y|)) with M the Mills ratio and the bracket expanded
# asymptotically (Abramowitz & Stegun 26.2.12): sum_j (-1)^(j+1) (2j-1)!! / y^2j.
# At the seam both forms are accurate to about 1e-12
.rdm_log_g <- function(y, log_Phi_y = stats::pnorm(y, log.p = TRUE)) {
  out <- numeric(length(y))
  series <- y < -12
  if (any(series)) {
    y2 <- y[series]^2
    term <- rep(1, length(y2))
    total <- 0
    for (j in seq_len(12)) {
      term <- term * (2 * j - 1) / y2
      total <- total + if (j %% 2 == 1) term else -term
    }
    out[series] <- stats::dnorm(y[series], log = TRUE) + log(total)
  }
  neg <- !series & y < 0
  out[neg] <- log_diff_exp_floored(
    stats::dnorm(y[neg], log = TRUE), log(-y[neg]) + log_Phi_y[neg]
  )
  pos <- y >= 0
  out[pos] <- log_sum_exp(log(y[pos]) + log_Phi_y[pos], stats::dnorm(y[pos], log = TRUE))
  out
}

# Tillman et al. (2020), Eq. 5: f = (1 / A) [drift (Phi(beta) - Phi(alpha)) +
# (s / sqrt(t)) (phi(alpha) - phi(beta))] with alpha, beta the standardised
# distances at the two ends of the start-point range. The first bracket is
# positive; the second takes the sign of |beta| - |alpha|
.dwald_full <- function(t, drift, bound, A, s, log = TRUE) {
  r <- recycle_args(t = t, drift = drift, bound = bound, A = A, s = s)
  t <- r$t; drift <- r$drift; bound <- r$bound; A <- r$A; s <- r$s
  gap <- bound - A
  st <- sqrt(pmax(t, 0))
  out <- rep(-Inf, length(t))
  ok <- t > 0

  small <- ok & A < .rdm_midpoint_ratio * s * st
  if (any(small)) {
    out[small] <- .dwald(t[small], drift = drift[small],
                         bound = gap[small] + A[small] / 2, s = s[small], log = TRUE)
  }

  gen <- ok & !small
  if (any(gen)) {
    t_g <- t[gen]; st_g <- st[gen]; drift_g <- drift[gen]; s_g <- s[gen]
    alpha <- (gap[gen] - drift_g * t_g) / (s_g * st_g)
    beta <- (bound[gen] - drift_g * t_g) / (s_g * st_g)
    l1 <- log(drift_g) + log_Phi_diff(alpha, beta)
    lpa <- stats::dnorm(alpha, log = TRUE)
    lpb <- stats::dnorm(beta, log = TRUE)
    second_positive <- abs(alpha) <= abs(beta)
    l2 <- log(s_g / st_g) +
      ifelse(second_positive, log_diff_exp_floored(lpa, lpb), log_diff_exp_floored(lpb, lpa))
    lnum <- ifelse(second_positive, log_sum_exp(l1, l2), log_diff_exp_floored(l1, l2))
    out[gen] <- lnum - log(A[gen])
  }

  if (log) out else exp(out)
}

# Tillman et al. (2020), Appendix A, assembled from the antiderivative of the
# plain Wald survival in the threshold. With G the antiderivative and
# g(y) = y Phi(y) + phi(y), q = s^2 / (2 drift), E(u) = exp(2 u drift / s^2)
# Phi(-(u + drift t) / (s sqrt(t))) and P = Phi(beta) - Phi(alpha):
#   S * A = D1 - D2,  D1 = s sqrt(t) (g(beta) - g(alpha)),
#                     D2 = q (P + E(bound) - E(gap)),
#   F * A = C1 + D2,  C1 = s sqrt(t) (g(-alpha) - g(-beta)).
# D1, D2, C1 and P are positive; E(bound) - E(gap) is signed and resolved
# explicitly. F never cancels; S cancels only in the far tail, by a factor of
# order drift t / gap, the same as the plain Wald survivor.
.pwald_full <- function(t, drift, bound, A, s, lower.tail = TRUE,
                        log.p = TRUE) {
  r <- recycle_args(t = t, drift = drift, bound = bound, A = A, s = s)
  t <- r$t; drift <- r$drift; bound <- r$bound; A <- r$A; s <- r$s
  gap <- bound - A
  st <- sqrt(pmax(t, 0))
  out <- rep(if (lower.tail) -Inf else 0, length(t))
  ok <- t > 0

  small <- ok & A < .rdm_midpoint_ratio * s * st
  if (any(small)) {
    out[small] <- .pwald(t[small], drift = drift[small],
                         bound = gap[small] + A[small] / 2, s = s[small],
                         lower.tail = lower.tail, log.p = TRUE)
  }

  gen <- ok & !small
  if (any(gen)) {
    t_g <- t[gen]; st_g <- st[gen]; drift_g <- drift[gen]; s_g <- s[gen]
    gap_g <- gap[gen]; bound_g <- bound[gen]
    denom <- s_g * st_g
    alpha <- (gap_g - drift_g * t_g) / denom
    beta <- (bound_g - drift_g * t_g) / denom
    log_q <- 2 * log(s_g) - log(2 * drift_g)
    lEb <- 2 * drift_g * bound_g / s_g^2 +
      stats::pnorm(-(bound_g + drift_g * t_g) / denom, log.p = TRUE)
    lEk <- 2 * drift_g * gap_g / s_g^2 +
      stats::pnorm(-(gap_g + drift_g * t_g) / denom, log.p = TRUE)
    lP <- log_Phi_diff(alpha, beta)
    lD2 <- log_q + ifelse(
      lEb >= lEk,
      log_sum_exp(lP, log_diff_exp_floored(lEb, lEk)),
      log_diff_exp_floored(log_sum_exp(lP, lEb), lEk)
    )
    if (lower.tail) {
      lC1 <- log(denom) + log_diff_exp_floored(.rdm_log_g(-alpha), .rdm_log_g(-beta))
      out[gen] <- pmin(log_sum_exp(lC1, lD2) - log(A[gen]), 0)
    } else {
      lD1 <- log(denom) + log_diff_exp_floored(.rdm_log_g(beta), .rdm_log_g(alpha))
      out[gen] <- pmin(log_diff_exp_floored(lD1, lD2) - log(A[gen]), 0)
    }
  }

  if (log.p) out else exp(out)
}

# log S(t) for one Wald accumulator with a start point uniform on [0, A], as a
# draws x length(t) matrix from draw-length parameter vectors
wald_log_surv <- function(t, drift, gap, A, s) {
  n_draws <- length(drift)
  matrix(
    .pwald_full(rep(t, each = n_draws), drift = drift, bound = gap + A, A = A,
                s = s, lower.tail = FALSE, log.p = TRUE),
    nrow = n_draws
  )
}
