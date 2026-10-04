############################################################################# !
# Log-Normal Race (LNR) distributions                                    ####
############################################################################# !

#' @title Distribution functions for the Log-Normal Race Model (LNR)
#'
#' @description Density, random generation, and cumulative distribution
#'   functions for the Log-Normal Race Model. The LNR (Rouder et al., 2015) is
#'   a K-accumulator race model where each accumulator's finishing time follows
#'   a lognormal distribution. The first accumulator to finish determines both
#'   the response and the response time.
#'
#' @name lnr_dist
#'
#' @param rt Numeric vector of response times in seconds. A response time at or
#'   below `ndt` lies outside the support, so the density there is 0 (`-Inf`
#'   with `log = TRUE`); `NA` propagates.
#' @param response Integer vector of responses (1:K, where K is the number of
#'   alternatives), of the same length as `rt` (or `q`).
#' @param m Numeric vector of meanlog parameters (one per accumulator).
#' @param s Numeric vector of sdlog parameters (one per accumulator, or a
#'   single value shared across all accumulators).
#' @param ndt Non-decision time in seconds, a single value or one per response
#'   time.
#' @param n Number of samples to generate.
#' @param log Logical; if `TRUE`, values are returned on the log scale.
#' @param q Numeric vector of quantiles (response times).
#' @param lower.tail Logical; if `TRUE` (default), probabilities are P(X <= x).
#' @param log.p Logical; if `TRUE`, probabilities are given as log(p).
#'
#' @param p Numeric vector of probabilities.
#'
#' @details
#' `plnr()` returns the marginal RT distribution function when `response` is
#' omitted, and the defective distribution function P(RT <= q, response = r)
#' when it is given; the defective values over all K responses sum to the
#' marginal one, and at large `q` each converges to the probability of that
#' response. The defective case has no closed form and is integrated
#' numerically, so it is slower than the marginal one. `qlnr()` inverts the
#' marginal distribution function only.
#'
#' @return
#'   - `dlnr()` returns a numeric vector of (log-)densities.
#'   - `rlnr()` returns a data.frame with columns `rt` and `response`.
#'   - `plnr()` returns a numeric vector of (log-)probabilities.
#'   - `qlnr()` returns a numeric vector of quantiles (response times).
#'
#' @references
#' Rouder, J. N., Province, J. M., Morey, R. D., Gomez, P., & Heathcote, A.
#'   (2015). The Lognormal Race: A Cognitive-Process Model of Choice and Latency
#'   with Desirable Psychometric Properties. Psychometrika, 80(2), 491-513.
#'
#' @keywords distribution
#'
#' @examples
#' dat <- rlnr(n = 1000, m = c(-1, 0), s = c(1, 1), ndt = 0.2)
#' head(dat)
#' hist(dat$rt)
#' @export
dlnr <- function(rt, response, m, s, ndt, log = FALSE) {
  validate_lnr_parameters(s, ndt, length(m), length(rt))
  stopif(
    !is.numeric(response) || anyNA(response),
    "response must contain integers in 1:{length(m)}."
  )
  stopif(
    any(response < 1) || any(response > length(m)) ||
      any(response != round(response)),
    "response must contain integers in 1:{length(m)}."
  )
  stopif(
    length(response) != length(rt),
    "response has {length(response)} entries but rt has {length(rt)}."
  )
  .dlnr(rt, response, m, s, ndt, log)
}

# A response time at or below the non-decision time is outside the support, so
# the density is 0 there, the way every stats::d*() reports an impossible value
.dlnr <- function(rt, response, m, s, ndt, log) {
  K <- length(m)
  if (length(s) == 1) s <- rep(s, K)
  t <- rt - ndt
  t[!is.na(t) & t <= 0] <- NA_real_

  log_lik <- stats::dlnorm(t, meanlog = m[response], sdlog = s[response],
                           log = TRUE)

  for (j in seq_len(K)) {
    is_loser <- (j != response)
    if (!any(is_loser)) next
    log_lik[is_loser] <- log_lik[is_loser] +
      stats::plnorm(t[is_loser], meanlog = m[j], sdlog = s[j],
                    lower.tail = FALSE, log.p = TRUE)
  }

  log_lik[is.na(log_lik) & !is.na(rt)] <- -Inf
  if (log) log_lik else exp(log_lik)
}

#' @rdname lnr_dist
#' @export
rlnr <- function(n, m, s, ndt) {
  validate_lnr_parameters(s, ndt, length(m), 1L)
  .rlnr(n, m, s, ndt)
}

.rlnr <- function(n, m, s, ndt) {
  if (n == 0) {
    return(data.frame(rt = numeric(0), response = integer(0)))
  }
  K <- length(m)
  if (length(s) == 1) s <- rep(s, K)
  race <- lnr_race(
    m = matrix(m, n, K, byrow = TRUE),
    s = matrix(s, n, K, byrow = TRUE),
    counts = matrix(1L, n, K)
  )
  data.frame(rt = race$rt + ndt, response = race$response)
}

#' @rdname lnr_dist
#' @export
plnr <- function(q, response, m, s, ndt, lower.tail = TRUE, log.p = FALSE) {
  validate_lnr_parameters(s, ndt, length(m), length(q))
  K <- length(m)
  if (length(s) == 1) s <- rep(s, K)
  t <- q - ndt

  if (missing(response)) {
    log_surv <- numeric(length(t))
    for (j in seq_len(K)) {
      log_surv <- log_surv + stats::plnorm(t, meanlog = m[j], sdlog = s[j],
                                            lower.tail = FALSE, log.p = TRUE)
    }
    # the survivor is the quantity with a closed form, so the upper tail is
    # returned as it stands: routing it through the CDF and back loses the
    # whole value once exp(log_surv) rounds below the resolution of 1
    log_p <- if (lower.tail) log1m_exp(log_surv) else log_surv
  } else {
    stopif(
      !length(response) %in% c(1L, length(q)),
      "response has {length(response)} entries but q has {length(q)}."
    )
    stopif(
      any(response < 1) || any(response > K) || any(response != round(response)),
      "response must contain integers in 1:{K}."
    )
    log_p <- .plnr_defective(t, rep_len(response, length(t)), m, s, ndt)
    if (!lower.tail) log_p <- log1m_exp(log_p)
  }

  if (log.p) log_p else exp(log_p)
}

# P(RT <= q, response = r) has no closed form: the winner's density is weighted
# by the losers' survivors, which do not factor out of the integral. The
# integration runs in log decision time and is bracketed by the extreme
# quantiles of the accumulators themselves; over an open or very wide interval
# integrate() samples past the decade that holds the mass and returns zero.
.plnr_defective <- function(t, response, m, s, ndt) {
  u_lo <- log(stats::qlnorm(1e-14, min(m), max(s)))
  u_cap <- log(stats::qlnorm(1e-14, max(m), max(s), lower.tail = FALSE))
  vapply(seq_along(t), function(i) {
    if (is.na(t[i])) return(NA_real_)
    if (t[i] <= 0) return(-Inf)
    log(stats::integrate(
      function(u) .dlnr(exp(u) + ndt, response[i], m, s, ndt, log = FALSE) *
        exp(u),
      lower = u_lo, upper = min(log(t[i]), u_cap),
      rel.tol = .Machine$double.eps^0.5
    )$value)
  }, numeric(1))
}

#' @rdname lnr_dist
#' @export
qlnr <- function(p, m, s, ndt, lower.tail = TRUE, log.p = FALSE) {
  validate_lnr_parameters(s, ndt, length(m), 1L)
  K <- length(m)
  if (length(s) == 1) s <- rep(s, K)
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  vapply(p, function(pi) {
    if (pi <= 0) return(ndt)
    if (pi >= 1) return(Inf)
    cdf_fn <- function(q) {
      t <- q - ndt
      if (t <= 0) return(-pi)
      log_surv <- sum(stats::plnorm(t, meanlog = m, sdlog = s,
                                     lower.tail = FALSE, log.p = TRUE))
      (1 - exp(log_surv)) - pi
    }
    upper <- ndt + stats::qlnorm(0.999, meanlog = max(m), sdlog = max(s))
    while (cdf_fn(upper) < 0) {
      upper <- ndt + 2 * (upper - ndt)
    }
    stats::uniroot(cdf_fn, interval = c(ndt + 1e-10, upper),
                   tol = 1e-8)$root
  }, numeric(1))
}

validate_lnr_parameters <- function(s, ndt, n_acc, n_obs) {
  stopif(
    any(!is.finite(s)) || any(s <= 0),
    "s (sdlog) must be finite and positive."
  )
  stopif(
    any(!is.finite(ndt)) || any(ndt < 0),
    "ndt (non-decision time) must be finite and non-negative."
  )
  # silent recycling of a mis-sized s pairs meanlogs with the wrong sdlogs
  stopif(
    !length(s) %in% c(1L, n_acc),
    "s has {length(s)} entries but the race has {n_acc} accumulators; \\
    pass one sdlog or one per accumulator."
  )
  stopif(
    !length(ndt) %in% c(1L, n_obs),
    "ndt has {length(ndt)} entries but there are {n_obs} response times."
  )
}

# log S(t) for a lognormal, as a draws x length(t) matrix from draw-length
# `meanlog`/`sdlog` vectors
lnorm_log_surv <- function(t, meanlog, sdlog) {
  matrix(
    stats::plnorm(rep(t, each = length(meanlog)), meanlog = meanlog,
                  sdlog = sdlog, lower.tail = FALSE, log.p = TRUE),
    nrow = length(meanlog)
  )
}

# The minimum of n i.i.d. finishing times has survivor S(t)^n, so one uniform
# inverts a whole group of identical accumulators; drawing the n copies and
# taking their minimum spends n times as many variates on the same
# distribution. `m`, `s` and `counts` are row-per-trial matrices, so per-trial
# parameters and per-trial accumulator counts need no loop over trials.
lnr_race <- function(m, s, counts) {
  ft <- matrix(Inf, nrow(m), ncol(m))
  for (j in seq_len(ncol(m))) {
    active <- counts[, j] > 0
    if (!any(active)) next
    surv <- stats::runif(sum(active))^(1 / counts[active, j])
    ft[active, j] <- stats::qlnorm(surv, meanlog = m[active, j],
                                   sdlog = s[active, j], lower.tail = FALSE)
  }
  list(
    rt = matrixStats::rowMins(ft),
    response = max.col(-ft, ties.method = "first")
  )
}
