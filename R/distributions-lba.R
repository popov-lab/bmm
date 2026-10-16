############################################################################# !
# LBA Distribution Functions                                             ####
############################################################################# !

#' @title LBA Distribution
#' @name lba_dist
#' @keywords distribution
#'
#' @description Density, distribution function, quantile function, and random
#'   generation for the Linear Ballistic Accumulator (LBA) model with normally
#'   distributed drift rates. `dlba()` is the joint density of the response
#'   time and the winning accumulator. `plba()` and `qlba()` describe the
#'   distribution of the race's first finishing time, regardless of which
#'   accumulator wins.
#'
#' @param rt Numeric vector of response times in seconds. Recycled with `gap`,
#'   `sp` and `ndt` to a common length.
#' @param n Integer. Number of samples to generate.
#' @param response Integer vector indicating the winning accumulator
#'   (1, 2, ..., K). May be recycled to the length of `rt`, and may contain
#'   `NA`; the corresponding density is `NA`.
#' @param drift Numeric vector of mean drift rates (one per accumulator). The
#'   mean may take any real value; the drawn drift rates are truncated at zero.
#' @param gap Numeric. Threshold gap (> 0). The distance between the maximum
#'   starting point and the decision threshold. The total threshold is computed
#'   as `b = gap + sp`, ensuring `b > sp` structurally. Recycled with `rt`,
#'   `sp` and `ndt` to a common length.
#' @param sp Numeric. Maximum starting point (`sp > 0`). Starting evidence is
#'   uniformly distributed on `[0, sp]`. Recycled with `rt`, `gap` and `ndt` to
#'   a common length.
#' @param ndt Numeric. Non-decision time in seconds (>= 0). Recycled with
#'   `rt`, `gap` and `sp` to a common length. `qlba()` requires a single value
#'   for `ndt`, `gap`, `sp` and `s`.
#' @param s Numeric. Drift rate SD (> 0, default = 1). Typically fixed to 1
#'   for identifiability.
#' @param distribution Character. The drift rate distribution. Currently
#'   `"normal"`.
#' @param log Logical. If `TRUE`, return log-densities. Default `FALSE`.
#' @param q Numeric vector of quantiles (response times).
#' @param p Numeric vector of probabilities.
#' @param lower.tail Logical. If `TRUE` (default), probabilities are
#'   `P(RT <= q)`.
#' @param log.p Logical. If `TRUE`, probabilities are on the log scale.
#'
#' @return
#'   - `dlba()` returns a numeric vector of (log-)densities.
#'   - `rlba()` returns a data.frame with columns `rt` and `response`.
#'   - `plba()` returns a numeric vector of (log-)probabilities.
#'   - `qlba()` returns a numeric vector of quantiles (response times).
#'
#' @references Brown, S. D., & Heathcote, A. (2008). The simplest complete
#'   model of choice response time: Linear ballistic accumulation. *Cognitive
#'   Psychology*, 57(3), 153-178.
#'
#' @examples
#' dat <- rlba(n = 1000, drift = c(3, 1.5), gap = 0.5, sp = 0.3, ndt = 0.3)
#' head(dat)
#' dlba(dat$rt[1:5], dat$response[1:5], drift = c(3, 1.5),
#'      gap = 0.5, sp = 0.3, ndt = 0.3)
#'
#' @export
dlba <- function(rt, response, drift, gap, sp, ndt, s = 1,
                 distribution = "normal",
                 log = FALSE) {
  distribution <- match.arg(distribution)
  validate_lba_parameters(drift, gap, sp, ndt, s, distribution)
  stopif(
    !is.numeric(rt) || !is.numeric(response),
    "rt and response must be numeric vectors."
  )
  stopif(
    !length(response) %in% c(1L, length(rt)),
    "response must have length 1 or the length of rt ({length(rt)}), \\
    not {length(response)}."
  )
  response <- rep_len(response, length(rt))
  stopif(
    !all(response[!is.na(response)] %in% seq_along(drift)),
    "response must index the accumulators, i.e. take values in 1:{length(drift)}."
  )
  b <- gap + sp
  A <- sp
  .dlba(rt, response, drift, b, A, ndt, s, distribution, log)
}

#' @rdname lba_dist
#' @export
rlba <- function(n, drift, gap, sp, ndt, s = 1,
                 distribution = "normal") {
  distribution <- match.arg(distribution)
  validate_lba_parameters(drift, gap, sp, ndt, s, distribution)
  b <- gap + sp
  A <- sp
  .rlba(n, drift, b, A, ndt, s, distribution)
}

#' @rdname lba_dist
#' @export
plba <- function(q, drift, gap, sp, ndt, s = 1,
                 distribution = "normal",
                 lower.tail = TRUE, log.p = FALSE) {
  distribution <- match.arg(distribution)
  validate_lba_parameters(drift, gap, sp, ndt, s, distribution)
  t <- q - ndt
  b <- gap + sp
  A <- sp
  # decision times t <= 0 never reach the kernels; those probabilities are 0 below
  t_kernel <- t
  t_kernel[!is.na(t) & t <= 0] <- NA_real_
  log_surv <- Reduce(`+`, lapply(drift, function(v) {
    .lba_lsurv_single(t_kernel, v, b, A, s, distribution)
  }))
  p_total <- -expm1(log_surv)
  p_total[!is.na(t) & t <= 0] <- 0

  if (!lower.tail) p_total <- 1 - p_total
  if (log.p) p_total <- log(p_total)
  p_total
}

#' @rdname lba_dist
#' @export
qlba <- function(p, drift, gap, sp, ndt, s = 1,
                 distribution = "normal",
                 lower.tail = TRUE, log.p = FALSE) {
  distribution <- match.arg(distribution)
  validate_lba_parameters(drift, gap, sp, ndt, s, distribution)

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  vapply(p, function(pi) {
    # p = 0 is ndt: the root search would extend into the flat region below ndt
    # and stop anywhere in it
    if (pi <= 0) {
      ndt
    } else {
      stats::uniroot(
        function(q) {
          plba(q, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s,
               distribution = distribution) - pi
        },
        interval = c(ndt + 1e-6, ndt + 20),
        extendInt = "upX",
        tol = 1e-8
      )$root
    }
  }, numeric(1))
}


.dlba <- function(rt, response, drift, b, A, ndt, s, distribution, log) {
  n <- length(rt)
  b <- rep_len(b, n)
  A <- rep_len(A, n)
  s <- rep_len(s, n)
  t <- rt - rep_len(ndt, n)
  # decision times t <= 0 never reach the kernels
  masked <- (!is.na(t) & t <= 0) | is.na(response)
  t_kernel <- t
  t_kernel[masked] <- NA_real_
  out <- .lba_lpdf_single(t_kernel, drift[response], b, A, s, distribution)
  out[is.na(response)] <- NA_real_
  for (j in seq_along(drift)) {
    i <- which(response != j)
    out[i] <- out[i] + .lba_lsurv_single(t_kernel[i], drift[j], b[i], A[i], s[i], distribution)
  }
  out[!is.na(t) & t <= 0] <- -Inf
  if (log) out else exp(out)
}

# Draw trial-to-trial drift rates for a single LBA accumulator. `mean` and `s`
# must be conformable (equal-length vectors or equal-shape matrices). The normal
# distribution resamples non-positive draws so drifts stay positive, matching the
# posdrift likelihood. This is the single source of drift sampling shared by
# rlba() and the posterior predict/epred methods.
.rlba_drift <- function(distribution, mean, s) {
  d <- switch(distribution,
    normal = {
      # brms passes a fixed s as one number, which s[neg] would turn into NA
      # and the redraw below into an endless loop
      s <- rep_len(s, length(mean))
      out <- stats::rnorm(length(mean), mean, s)
      neg <- which(out <= 0)
      while (length(neg) > 0) {
        out[neg] <- stats::rnorm(length(neg), mean[neg], s[neg])
        neg <- neg[out[neg] <= 0]
      }
      out
    }
  )
  if (is.matrix(mean)) dim(d) <- dim(mean)
  d
}

.rlba <- function(n, drift, b, A, ndt, s, distribution) {
  K <- length(drift)

  ft <- matrix(NA_real_, nrow = n, ncol = K)
  for (j in seq_len(K)) {
    start <- stats::runif(n, min = 0, max = A)
    d <- .rlba_drift(distribution, rep_len(drift[j], n), rep_len(s, n))
    ft[, j] <- (b - start) / d
  }

  data.frame(
    rt = apply(ft, 1, min) + ndt,
    response = apply(ft, 1, which.min)
  )
}


############################################################################# !
# Single-accumulator kernels, in log space                                ####
############################################################################# !
# Twins of inst/stan_chunks/lba_*_functions.stan. One accumulator with
# threshold b, start point k ~ U(0, A) and drift d finishes at (b - k) / d, so
# at decision time t the drifts that finish exactly then are u = (b - k) / t
# in (lo, hi) = ((b - A) / t, b / t). The density is M / A with
# M = int_{lo}^{hi} u f(u) du, and the survivor is (u_num - t M) / A with
# u_num = b F(hi) - (b - A) F(lo) = A F(lo) + b (F(hi) - F(lo)), the second
# form being a sum of positive terms and the one used on both sides.
#
# The differences of CDFs are taken from the tail where both terms are small,
# and where an interval is too narrow for a difference to carry digits the
# integral is replaced by the midpoint rule with its second-order term. No
# Stan branch evaluates an expression that can return -Inf with an infinite
# partial: Stan keeps the adjoint of a discarded branch, so log_diff_exp(x, x)
# or log(0) poison the gradient of the whole trial even when the value is not
# used. cogmod 0.3.2 (Makowski) exposed the Phi saturation and the -690 floor
# of the previous normal-drift kernel; the layouts were derived independently.
#
# The floor log(1e-300) is a last resort for a numerator that has rounded to
# zero or below. On the verification grid (t - ndt down to 1 ms, |z| to 40, A
# from 1e-8 to 2, gap from 1e-3 to 1.5) no branch reaches it; it is left for
# parameter values beyond that range (v / s below -37, or a decision time so
# long that the truncated survivor falls under 1e-300).
#
# Every kernel takes equal-length vectors (one call per observation over all
# posterior draws, the brms log_lik contract), selects its branches by mask,
# returns NA for an undefined decision time and floors a numerator that has
# rounded to zero at log(1e-300) as the Stan side does.

.lba_lpdf_single <- function(t, v, b, A, s, distribution) {
  args <- recycle_args(t = t, v = v, b = b, A = A, s = s)
  kernel <- switch(distribution,
    normal = .lba_normal_lpdf
  )
  do.call(kernel, unname(args))
}

.lba_lsurv_single <- function(t, v, b, A, s, distribution) {
  args <- recycle_args(t = t, v = v, b = b, A = A, s = s)
  kernel <- switch(distribution,
    normal = .lba_normal_lsurv
  )
  pmin(do.call(kernel, unname(args)), 0)
}

.dlba_single <- function(t, v, b, A, s, distribution, log = FALSE) {
  lp <- .lba_lpdf_single(t, v, b, A, s, distribution)
  if (log) lp else exp(lp)
}

.plba_single <- function(t, v, b, A, s, distribution) {
  pmax(-expm1(.lba_lsurv_single(t, v, b, A, s, distribution)), 0)
}

.lba_log_floor <- log(1e-300)

# int_{z_lo}^{z_hi} phi(z) dz for a narrow interval: midpoint rule with its
# second-order term (relative error O((dz max(|z|, 1))^4 / 1920))
.lba_log_phi_int_narrow <- function(z_lo, z_hi) {
  dz <- z_hi - z_lo
  z_m <- 0.5 * (z_lo + z_hi)
  log(dz) + stats::dnorm(z_m, log = TRUE) + log1p(dz^2 / 24 * (z_m^2 - 1))
}


# --- normal drift, truncated at zero -----------------------------------------

# Normal drift d ~ N(v, s^2) conditional on d > 0 (the posdrift convention of
# rtdists). In z-units, z = (u - v) / s, with delta = z_hi - z_lo = A / (t s):
#   M = int_{z_lo}^{z_hi} (v + s z) phi(z) dz = v [Phi(z_hi) - Phi(z_lo)] + s [phi(z_lo) - phi(z_hi)]
#   S - q = (1/delta) int_{z_lo}^{z_hi} [Phi(z) - Phi(-v/s)] dz,  q = Phi(-v/s)
# and the truncated density and survivor are M / (A Phi(v/s)) and
# (S - q) / Phi(v/s). Both Phi differences vanish with delta and both phi
# differences vanish with delta (z_lo + z_hi); the sum for M can also be a
# difference of two terms of the same size (v < 0, or v > 0 with the interval
# left of v). Every piece is therefore assembled in log space from a signed
# term, or by the midpoint rule below delta = 1e-4, where the direct form has
# lost four digits and the expansion has ~1e-12 left. In the midpoint branch
# u_m = v + s z_m is the midpoint of (lo, hi), written without the
# cancellation of v + s z_m. phi(near) - phi(far) = phi(near) (1 - exp(-x));
# the term is exactly zero when z_lo = -z_hi. With v < 0, z_lo > 0, so the
# phi term is positive and the larger one.
.lba_normal_log_M <- function(t, v, b, A, s) {
  z_lo <- ((b - A) / t - v) / s
  z_hi <- (b / t - v) / s
  delta <- A / (t * s)
  out <- rep(NA_real_, length(t))
  ok <- !is.na(t)

  i <- ok & delta < 1e-4
  if (any(i)) {
    z_m <- 0.5 * (z_lo[i] + z_hi[i])
    u_m <- (b[i] - 0.5 * A[i]) / t[i]
    bracket <- u_m + delta[i]^2 / 24 * (u_m * (z_m^2 - 1) - 2 * s[i] * z_m)
    out[i] <- ifelse(bracket > 0,
                     log(delta[i]) + stats::dnorm(z_m, log = TRUE) + log(pmax(bracket, 1e-300)),
                     .lba_log_floor)
  }

  i <- ok & delta >= 1e-4
  if (any(i)) {
    zl <- z_lo[i]
    zh <- z_hi[i]
    vv <- v[i]
    sum_z <- zl + zh
    x <- 0.5 * delta[i] * abs(sum_z)
    log_dPhi <- log_Phi_diff(zl, zh)
    log_dphi <- ifelse(x > 0,
                       log(s[i]) + stats::dnorm(ifelse(sum_z >= 0, zl, zh), log = TRUE) + log1m_exp(-x),
                       -Inf)
    l1 <- log(abs(vv)) + log_dPhi
    out[i] <- ifelse(vv > 0,
      ifelse(sum_z >= 0, log_sum_exp(l1, log_dphi), log_diff_exp_floored(l1, log_dphi, .lba_log_floor)),
      ifelse(vv < 0, log_diff_exp_floored(log_dphi, l1, .lba_log_floor), log_dphi)
    )
  }
  out
}

.lba_normal_lpdf <- function(t, v, b, A, s) {
  .lba_normal_log_M(t, v, b, A, s) - log(A) - stats::pnorm(v / s, log.p = TRUE)
}

# S - q is formed in probability space from the exact-gradient Phi(): from the
# lower tails through g(z) = z Phi(z) + phi(z) when v >= 0 (q < 1/2 is the
# small term), from the upper tails through h(z) = z Phi(-z) - phi(z) when
# v < 0 (1 - q and F are the small terms). The remaining cancellation is
# intrinsic to the closed form and bounded by (S - q) / S, which stays above
# 1e-5 for any decision time under 30 s with gap >= 1e-3.
.lba_normal_lsurv <- function(t, v, b, A, s) {
  z_lo <- ((b - A) / t - v) / s
  z_hi <- (b / t - v) / s
  delta <- A / (t * s)
  log_denom <- stats::pnorm(v / s, log.p = TRUE)
  g <- function(z) z * stats::pnorm(z) + stats::dnorm(z)
  h <- function(z) z * stats::pnorm(-z) - stats::dnorm(z)
  num <- rep(NA_real_, length(t))
  ok <- !is.na(t)

  i <- ok & delta < 1e-4
  if (any(i)) {
    z_m <- 0.5 * (z_lo[i] + z_hi[i])
    corr <- delta[i]^2 / 24 * z_m * stats::dnorm(z_m)
    num[i] <- ifelse(v[i] >= 0,
                     stats::pnorm(z_m) + expm1(log_denom[i]),
                     exp(log_denom[i]) - stats::pnorm(-z_m)) - corr
  }
  i <- ok & delta >= 1e-4 & v >= 0
  if (any(i)) {
    num[i] <- (g(z_hi[i]) - g(z_lo[i])) / delta[i] + expm1(log_denom[i])
  }
  i <- ok & delta >= 1e-4 & v < 0
  if (any(i)) {
    num[i] <- exp(log_denom[i]) - (h(z_hi[i]) - h(z_lo[i])) / delta[i]
  }
  log(pmax(num, 1e-300)) - log_denom
}


validate_lba_parameters <- function(drift, gap, sp, ndt, s, distribution) {
  stopif(any(gap <= 0), "gap (threshold gap) must be positive.")
  stopif(any(sp <= 0), "sp (maximum starting point) must be positive.")
  stopif(any(ndt < 0), "ndt (non-decision time) must be non-negative.")
  stopif(any(s <= 0), "s (scale parameter) must be positive.")
}
