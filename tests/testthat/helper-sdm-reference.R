# log Z of an sdm() cell from integrate(), independent of the rules in
# R/helpers-sdm.R and the Stan chunks (#558). The integrand is
# exp(c* expm1(log s(u))), the density over its peak with the subtraction
# done without cancellation, broken at multiples of the widths of the spike
# at u = 0. von Mises: where the activation falls to about 1, and the SD of
# exp(eta) near the peak (#615). Laplace, in v = kappa u: multiples of 1 / c*,
# over which exp(eta) falls by e near the top, and log c* plus a few units,
# past which the integrand over its floor is exp(exp(-(v - log c*))).
# A panel whose left end times its width is below 1e-20 of the running total
# is skipped: the integrand falls with u, and integrate() calls a long panel
# that underflows divergent
sdm_log_z_integrate <- function(log_cstar, kappa, activation) {
  peak <- exp(log_cstar)
  cuts <- if (activation == "vonmises") {
    width <- sqrt(2 * max(log_cstar, 1) / kappa)
    sd_peak <- 1 / sqrt(peak * kappa)
    c(width * c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16),
      sd_peak * c(0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64))
  } else {
    v0 <- max(log_cstar, 1)
    c(c(0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 128) / peak,
      v0 * c(0.125, 0.25, 0.5), v0 + c(0, 1, 2, 4, 6, 10, 15, 20, 30, 40)) / kappa
  }
  breaks <- sort(unique(c(0, pmin(cuts, pi), pi)))
  log_s <- if (activation == "vonmises") {
    function(u) -2 * kappa * sin(u / 2)^2
  } else {
    function(u) -kappa * u
  }
  f <- function(u) exp(peak * expm1(log_s(u)))
  total <- 0
  for (i in seq_len(length(breaks) - 1)) {
    if (f(breaks[i]) * (breaks[i + 1] - breaks[i]) < 1e-20 * total) next
    total <- total + stats::integrate(f, breaks[i], breaks[i + 1],
      rel.tol = 1e-13, abs.tol = .Machine$double.xmin, subdivisions = 5000L
    )$value
  }
  log(2) + peak + log(total)
}

# c*, the peak of the activation, written out per cell from the definitions
# in ?sdm rather than taken from bmm
sdm_cstar_by_hand <- function(c, kappa, activation, scaling) {
  integral_of_s <- if (activation == "vonmises") {
    sqrt(2 * pi / kappa)
  } else {
    2 * (1 - exp(-kappa * pi)) / kappa
  }
  s_opposite <- if (activation == "vonmises") exp(-2 * kappa) else exp(-kappa * pi)
  switch(scaling,
    density = c / integral_of_s,
    peak = c,
    range = c / (1 - s_opposite)
  )
}

# every non-default cell, and the default
sdm_cells <- expand.grid(
  activation = c("vonmises", "laplace"),
  scaling = c("density", "peak", "range"),
  stringsAsFactors = FALSE
)
