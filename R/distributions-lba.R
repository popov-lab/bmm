############################################################################# !
# LBA Distribution Functions                                             ####
############################################################################# !

#' @title LBA Distribution
#' @name lba_dist
#' @keywords distribution
#'
#' @description Density, distribution function, quantile function, and random
#'   generation for the Linear Ballistic Accumulator (LBA) model with multiple
#'   drift rate distributions.
#'
#' @param rt Numeric vector of response times in seconds. Recycled with `gap`,
#'   `sp` and `ndt` to a common length.
#' @param n Integer. Number of samples to generate.
#' @param response Integer vector indicating the winning accumulator
#'   (1, 2, ..., K). May be recycled to the length of `rt`, and may contain
#'   `NA`; the corresponding density is `NA`.
#' @param drift Numeric vector of drift rate parameters (one per accumulator).
#'   Interpretation depends on `distribution`: mean drift for `"normal"`, shape
#'   for `"gamma"` and `"frechet"`, meanlog for `"lognormal"`. The lognormal
#'   meanlog may take any real value, unlike the gamma and Frechet shapes,
#'   which must be positive.
#' @param gap Numeric. Threshold gap (> 0). The distance between the maximum
#'   starting point and the decision threshold. The total threshold is computed
#'   as `b = gap + sp`, ensuring `b > sp` structurally. Recycled with `rt`,
#'   `sp` and `ndt` to a common length.
#' @param sp Numeric. Maximum starting point (>= 0 here; strictly positive,
#'   `sp > 0`, in the `lba()` model). Starting evidence is uniformly
#'   distributed on `[0, sp]`. These distribution functions additionally
#'   accept `sp = 0` as the no-starting-point-variability limit, which the
#'   fitted model does not. Recycled with `rt`, `gap` and `ndt` to a common
#'   length.
#' @param ndt Numeric. Non-decision time in seconds (>= 0). Recycled with
#'   `rt`, `gap` and `sp` to a common length.
#' @param s Numeric. Scale parameter (> 0, default = 1). Interpretation depends
#'   on `distribution`: drift SD for `"normal"`, rate for `"gamma"`, scale for
#'   `"frechet"`, sdlog for `"lognormal"`. Typically fixed to 1 for
#'   identifiability.
#' @param distribution Character. The drift rate distribution. One of
#'   `"normal"` (default), `"gamma"`, `"frechet"`, or `"lognormal"`.
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
                 distribution = c("normal", "gamma", "frechet", "lognormal"),
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
                 distribution = c("normal", "gamma", "frechet", "lognormal")) {
  distribution <- match.arg(distribution)
  validate_lba_parameters(drift, gap, sp, ndt, s, distribution)
  b <- gap + sp
  A <- sp
  .rlba(n, drift, b, A, ndt, s, distribution)
}

#' @rdname lba_dist
#' @export
plba <- function(q, drift, gap, sp, ndt, s = 1,
                 distribution = c("normal", "gamma", "frechet", "lognormal"),
                 lower.tail = TRUE, log.p = FALSE) {
  distribution <- match.arg(distribution)
  validate_lba_parameters(drift, gap, sp, ndt, s, distribution)
  t <- q - ndt
  b <- gap + sp
  A <- sp
  log_surv <- Reduce(`+`, lapply(drift, function(v) {
    .lba_lsurv_single(t, v, b, A, s, distribution)
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
                 distribution = c("normal", "gamma", "frechet", "lognormal"),
                 lower.tail = TRUE, log.p = FALSE) {
  distribution <- match.arg(distribution)
  validate_lba_parameters(drift, gap, sp, ndt, s, distribution)

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  vapply(p, function(pi) {
    stats::uniroot(
      function(q) {
        plba(q, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s,
             distribution = distribution) - pi
      },
      interval = c(ndt + 1e-6, ndt + 20),
      extendInt = "upX",
      tol = 1e-8
    )$root
  }, numeric(1))
}


.dlba <- function(rt, response, drift, b, A, ndt, s, distribution, log) {
  n <- length(rt)
  b <- rep_len(b, n)
  A <- rep_len(A, n)
  s <- rep_len(s, n)
  t <- rt - rep_len(ndt, n)
  out <- .lba_lpdf_single(t, drift[response], b, A, s, distribution)
  out[is.na(response)] <- NA_real_
  for (j in seq_along(drift)) {
    i <- which(response != j)
    out[i] <- out[i] + .lba_lsurv_single(t[i], drift[j], b[i], A[i], s[i], distribution)
  }
  out[!is.na(t) & t <= 0] <- -Inf
  if (log) out else exp(out)
}

# Draw trial-to-trial drift rates for a single LBA accumulator. `mean` and `s`
# must be conformable (equal-length vectors or equal-shape matrices). The normal
# distribution resamples non-positive draws so drifts stay positive, matching the
# posdrift likelihood; the other distributions are positive by construction. This
# is the single source of drift sampling shared by rlba() and the posterior
# predict/epred methods.
.rlba_drift <- function(distribution, mean, s) {
  d <- switch(distribution,
    normal = {
      out <- stats::rnorm(length(mean), mean, s)
      neg <- which(out <= 0)
      while (length(neg) > 0) {
        out[neg] <- stats::rnorm(length(neg), mean[neg], s[neg])
        neg <- neg[out[neg] <= 0]
      }
      out
    },
    gamma = stats::rgamma(length(mean), shape = mean, rate = s),
    frechet = .rfrechet(length(mean), shape = mean, scale = s),
    lognormal = stats::rlnorm(length(mean), meanlog = mean, sdlog = s)
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
    normal = .lba_normal_lpdf,
    gamma = .lba_gamma_lpdf,
    frechet = .lba_frechet_lpdf,
    lognormal = .lba_lognormal_lpdf
  )
  do.call(kernel, unname(args))
}

.lba_lsurv_single <- function(t, v, b, A, s, distribution) {
  args <- recycle_args(t = t, v = v, b = b, A = A, s = s)
  kernel <- switch(distribution,
    normal = .lba_normal_lsurv,
    gamma = .lba_gamma_lsurv,
    frechet = .lba_frechet_lsurv,
    lognormal = .lba_lognormal_lsurv
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


# --- gamma drift -------------------------------------------------------------

# Gamma drift d ~ Gamma(shape v, rate s). M = (v / s) [F(hi; v + 1) - F(lo; v + 1)].
#
# Two facts about Stan Math's gamma CDFs bind the Stan twin (stan-dev/math
# #3408, measured on 5.3.0 and 5.4.0): gamma_lcdf rounds to 0 once log Q < -37,
# so a difference of two lcdf values is -Inf for every fast response with
# b / t >= 10, and its shape partial is off by 2.5e-3 at shape 5, by 0.1 to
# 50 % for shape >= 10 in the lower tail, NaN once P saturates (shape >= 12,
# s b / t above ~1e3) and an exception above 5e4. gamma_lccdf is finite far
# beyond exp underflow. So the difference is taken from the upper tail
# (lccdf pair) beyond the mean, from the lower tail (lcdf pair) below it, and
# where the interval is narrow the midpoint rule through gamma_lpdf (exact
# partials) replaces both; the midpoint test uses the rule's expansion
# parameter, du |d log f / du| at the midpoint. The shape is v + 1 with v an
# autodiff parameter, and s b / t grows without bound as ndt approaches the
# fastest response: the gamma drift is the least robust of the four there.
.lba_gamma_log_dF <- function(lo, hi, alpha, beta) {
  du <- hi - lo
  u_m <- 0.5 * (lo + hi)
  out <- rep(NA_real_, length(lo))
  ok <- !is.na(lo)

  i <- ok & du * (abs(alpha - 1) / u_m + beta) < 1e-3
  if (any(i)) {
    r <- (alpha[i] - 1) / u_m[i] - beta[i]
    out[i] <- log(du[i]) + stats::dgamma(u_m[i], shape = alpha[i], rate = beta[i], log = TRUE) +
      log1p(du[i]^2 / 24 * (r^2 - (alpha[i] - 1) / u_m[i]^2))
  }
  i <- ok & is.na(out) & lo * beta > alpha
  if (any(i)) {
    out[i] <- log_diff_exp(
      stats::pgamma(lo[i], shape = alpha[i], rate = beta[i], lower.tail = FALSE, log.p = TRUE),
      stats::pgamma(hi[i], shape = alpha[i], rate = beta[i], lower.tail = FALSE, log.p = TRUE)
    )
  }
  i <- ok & is.na(out)
  if (any(i)) {
    out[i] <- log_diff_exp(
      stats::pgamma(hi[i], shape = alpha[i], rate = beta[i], log.p = TRUE),
      stats::pgamma(lo[i], shape = alpha[i], rate = beta[i], log.p = TRUE)
    )
  }
  out
}

.lba_gamma_lpdf <- function(t, v, b, A, s) {
  log(v) - log(s) + .lba_gamma_log_dF((b - A) / t, b / t, v + 1, s) - log(A)
}

# u_num - t M = t int_{lo}^{hi} F(u) du > 0; the difference loses digits only
# where the survivor itself is negligible
.lba_gamma_lsurv <- function(t, v, b, A, s) {
  lo <- (b - A) / t
  hi <- b / t
  log_tM <- log(t) + log(v) - log(s) + .lba_gamma_log_dF(lo, hi, v + 1, s)
  log_u <- log_sum_exp(
    log(A) + stats::pgamma(lo, shape = v, rate = s, log.p = TRUE),
    log(b) + .lba_gamma_log_dF(lo, hi, v, s)
  )
  log_diff_exp_floored(log_u, log_tM, .lba_log_floor) - log(A)
}


# --- lognormal drift ---------------------------------------------------------

# Lognormal drift d ~ LN(v, s), so z = (log u - v) / s and F(u) = Phi(z). The
# width of the interval in z-units is dz = log(b / (b - A)) / s, taken from
# log1m(A / b) rather than as a difference of two logs so that a tiny A keeps
# its digits. M = exp(v + s^2 / 2) [Phi(z_hi - s) - Phi(z_lo - s)].
.lba_lognormal_log_dPhi <- function(z_lo, z_hi, dz) {
  narrow <- !is.na(dz) & dz < 1e-4
  out <- log_Phi_diff(z_lo, z_hi)
  out[narrow] <- .lba_log_phi_int_narrow(z_lo[narrow], z_hi[narrow])
  out
}

.lba_lognormal_lpdf <- function(t, v, b, A, s) {
  z_hi <- (log(b / t) - v) / s
  dz <- -log1p(-A / b) / s
  v + 0.5 * s^2 + .lba_lognormal_log_dPhi(z_hi - dz - s, z_hi - s, dz) - log(A)
}

.lba_lognormal_lsurv <- function(t, v, b, A, s) {
  z_hi <- (log(b / t) - v) / s
  dz <- -log1p(-A / b) / s
  z_lo <- z_hi - dz
  log_tM <- log(t) + v + 0.5 * s^2 + .lba_lognormal_log_dPhi(z_lo - s, z_hi - s, dz)
  log_u <- log_sum_exp(
    log(A) + stats::pnorm(z_lo, log.p = TRUE),
    log(b) + .lba_lognormal_log_dPhi(z_lo, z_hi, dz)
  )
  log_diff_exp_floored(log_u, log_tM, .lba_log_floor) - log(A)
}


# --- Frechet drift -----------------------------------------------------------

.rfrechet <- function(n, shape, scale) {
  scale * (-log(stats::runif(n)))^(-1 / shape)
}

.lba_frechet_log_F <- function(x, shape, scale) {
  -(x / scale)^(-shape)
}

# Frechet drift d ~ Frechet(shape v, scale s), log F(u) = -(u / s)^-v. The
# survivor is spelled through this closed form and log1m_exp, from the upper
# tail once F(lo) > 1/2, never through Stan's frechet_lccdf, whose
# log1m(exp(.)) form is -Inf with an infinite partial once (s / u)^v < 5e-17
# (u / s = 113 at shape 8).
.lba_frechet_log_dF <- function(lo, hi, shape, scale) {
  lF_lo <- .lba_frechet_log_F(lo, shape, scale)
  lF_hi <- .lba_frechet_log_F(hi, shape, scale)
  upper <- !is.na(lF_lo) & lF_lo > -log(2)
  out <- log_diff_exp(lF_hi, lF_lo)
  out[upper] <- log_diff_exp(log1m_exp(lF_lo[upper]), log1m_exp(lF_hi[upper]))
  out
}

# The truncated first moment of a Frechet drift has no closed form: each
# element is integrated adaptively over the part of (lo, hi) where the
# integrand u f(u) is within exp(-80) of its peak (at u = scale, or the nearer
# end point), scaled by that peak. Without the window a fast or slow response
# concentrates the integrand in a sliver that integrate() cannot see. This is
# the independent reference for the 16-point rule of the Stan chunk.
#
# The Stan chunk integrates M by 16-point Gauss-Legendre in log space.
# Measured against this adaptive quadrature: within 1e-6 nats for
# A / gap <= 6, 6e-5 at A / gap = 10, 8.7e-3 at 40, 0.84 nats at gap = 1e-3
# with A = 2; the default priors put A / gap near 0.6. The survivor inherits
# that error multiplied by (u / s)^-2v, so for slow responses (b / t well
# below s) it is off by whole nats where its value is already below -100.
.lba_frechet_log_M <- function(t, shape, b, A, scale) {
  log_integrand <- function(u, shape, scale) {
    log(shape) - shape * log(u / scale) - (u / scale)^(-shape)
  }
  window_edge <- function(from, to, shape, scale, floor) {
    if (log_integrand(from, shape, scale) >= floor) return(from)
    stats::uniroot(function(u) log_integrand(u, shape, scale) - floor, c(from, to), tol = 1e-12)$root
  }
  out <- rep(NA_real_, length(t))
  ok <- !is.na(t)
  if (!any(ok)) return(out)
  out[ok] <- mapply(function(lo, hi, shape, scale) {
    peak <- min(max(scale, lo), hi)
    m <- log_integrand(peak, shape, scale)
    lo <- window_edge(lo, peak, shape, scale, m - 80)
    hi <- window_edge(hi, peak, shape, scale, m - 80)
    m + log(stats::integrate(
      function(u) exp(log_integrand(u, shape, scale) - m),
      lower = lo, upper = hi
    )$value)
  }, (b[ok] - A[ok]) / t[ok], b[ok] / t[ok], shape[ok], scale[ok])
  out
}

.lba_frechet_lpdf <- function(t, v, b, A, s) {
  .lba_frechet_log_M(t, v, b, A, s) - log(A)
}

.lba_frechet_lsurv <- function(t, v, b, A, s) {
  lo <- (b - A) / t
  hi <- b / t
  log_tM <- log(t) + .lba_frechet_log_M(t, v, b, A, s)
  log_u <- log_sum_exp(
    log(A) + .lba_frechet_log_F(lo, v, s),
    log(b) + .lba_frechet_log_dF(lo, hi, v, s)
  )
  log_diff_exp_floored(log_u, log_tM, .lba_log_floor) - log(A)
}


validate_lba_parameters <- function(drift, gap, sp, ndt, s, distribution) {
  stopif(any(gap <= 0), "gap (threshold gap) must be positive.")
  stopif(any(sp < 0), "sp (maximum starting point) must be non-negative.")
  stopif(any(ndt < 0), "ndt (non-decision time) must be non-negative.")
  stopif(any(s <= 0), "s (scale parameter) must be positive.")
  # the normal drift is a mean and the lognormal a meanlog, both unrestricted;
  # the gamma and Frechet drifts are shape parameters
  if (distribution %in% c("gamma", "frechet")) {
    stopif(any(drift <= 0),
           "drift must be positive for distribution '{distribution}'.")
  }
}
