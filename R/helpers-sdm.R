# The cells of sdm() (#558). Every density is exp(c* s(x - mu)) / Z(c*, kappa),
# with s(0) = 1 and c* the peak of the activation. The scaling decides what
# `c` measures: "density" divides it by the integral of s over the circle (the
# normal approximation sqrt(2 pi / kappa) for von Mises, which is bmm's
# historical sqrtexp form, and exactly 2 (1 - exp(-kappa pi)) / kappa for
# Laplace), "peak" takes it as the activation of the target, and "range" as
# the contrast between target and opposite, c* (1 - s(pi)), whose additive
# constant cancels against Z. Stan reads the cell as integer codes in this
# order (inst/stan_chunks/sdm_kernel_funs.stan)
.sdm_activations <- c("vonmises", "laplace")
.sdm_scalings <- c("density", "peak", "range")

# the scaling a cell gets when sdm() or the d/p/q/r functions are not given one;
# Laplace has no density scaling in Oberauer (2023), so its natural unit is the peak
.sdm_default_scaling <- function(activation) {
  if (activation == "laplace") "peak" else "density"
}

.sdm_is_default_cell <- function(activation, scaling) {
  activation == "vonmises" && scaling == "density"
}

# log(1 - exp(-a)) for a >= 0
.log1mexp_neg <- function(a) {
  ifelse(a < log(2), log(-expm1(-a)), log1p(-exp(-a)))
}

.sdm_log_cstar <- function(log_c, kappa, activation, scaling) {
  switch(scaling,
    density = if (activation == "vonmises") {
      log_c + 0.5 * log(kappa / (2 * pi))
    } else {
      # kappa / (1 - exp(-kappa pi)) tends to 1 / pi as kappa goes to 0;
      # ifelse() takes the length of its condition, so both are recycled first
      n <- max(length(log_c), length(kappa))
      log_c <- rep_len(log_c, n)
      kappa <- rep_len(kappa, n)
      ifelse(kappa == 0, log_c - log(2 * pi), log_c + log(kappa) - log(2) - .log1mexp_neg(kappa * pi))
    },
    peak = log_c,
    range = log_c - .log1mexp_neg(if (activation == "vonmises") 2 * kappa else kappa * pi)
  )
}

# log s(delta); cos(delta) - 1 is formed as -2 sin(delta / 2)^2, which keeps
# its digits near the target
.sdm_log_s <- function(delta, kappa, activation) {
  if (activation == "vonmises") {
    -2 * kappa * sin(delta / 2)^2
  } else {
    d <- abs(delta) %% (2 * pi)
    -kappa * pmin(d, 2 * pi - d)
  }
}

# log Z = log 2 + c* + log of the integral of exp(c* expm1(log s(u))) over
# [0, pi]: the integrand over its peak, so that exp() neither overflows at
# large c* (dsdm() before #558 stopped once log Z passed 709) nor loses digits
# near the peak. The rules are the Stan normalisers of sdm_simple_funs.stan
# and sdm_kernel_funs.stan, vectorised over cells; one integrate() call per
# cell would make log_lik(), and so loo(), several times slower than before
.sdm_log_z <- function(log_cstar, kappa, activation) {
  n <- max(length(log_cstar), length(kappa))
  log_cstar <- rep_len(log_cstar, n)
  kappa <- rep_len(kappa, n)
  out <- rep(NA_real_, n)
  ok <- !is.na(log_cstar) & !is.na(kappa)
  if (!any(ok)) {
    return(out)
  }
  out[ok] <- log(2) + exp(log_cstar[ok]) + .sdm_log_int(log_cstar[ok], kappa[ok], activation)
  out
}

.sdm_log_int <- function(log_cstar, kappa, activation) {
  peak <- exp(log_cstar)
  if (activation == "vonmises") {
    width <- pmin(pi, sqrt(2 * pmax(log_cstar, 1) / kappa))
    sd_peak <- exp(-0.5 * (log_cstar + log(kappa)))
    cuts <- cbind(0, pi, width * 0.5, width, width * 2, width * 5, sd_peak, sd_peak * 3, sd_peak * 9)
  } else {
    v0 <- pmax(log_cstar, 1)
    cuts <- cbind(0, pi, cbind(3 / peak, 30 / peak, v0, v0 + 6, v0 + 30) / kappa)
  }
  cuts <- pmin(pmax(cuts, 0), pi)
  cuts <- matrix(t(apply(cuts, 1, sort)), nrow = length(peak))
  gl <- .gl_rule(16L)
  out <- rep(-Inf, length(peak))
  for (j in seq_len(ncol(cuts) - 1)) {
    half <- (cuts[, j + 1] - cuts[, j]) / 2
    u <- (cuts[, j] + cuts[, j + 1]) / 2 + outer(half, gl$nodes)
    d <- peak * expm1(.sdm_log_s(u, kappa, activation))
    # d falls with u, and the first node is the smallest
    d_max <- d[, 1]
    part <- log(half) + d_max + log(drop(exp(d - d_max) %*% gl$weights))
    used <- half > 0
    out[used] <- .log_add_exp(out[used], part[used])
  }
  if (activation == "vonmises") {
    # the fixed rule of the #615 hybrid where the spike is not sharp
    fixed <- kappa * (peak + 1) <= 4000
    wide <- .sdm_wide_rule()
    out[fixed] <- log(drop(exp(outer(peak[fixed], rep(1, length(wide$cosm1))) *
      expm1(outer(kappa[fixed], wide$cosm1))) %*% wide$weights))
  }
  out
}

.log_add_exp <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}
