# Recovery plan for sdm() with the von Mises activation (activation = "vonmises")
# and its default scaling, "density" (Oberauer, 2023; bmm's sqrtexp form). The
# response error x has the density exp(c* s(x)) / Z on the circle, with the
# activation kernel s(x) = exp(kappa * (cos(x) - 1)) and the activation of the
# target c* = c * sqrt(kappa / (2 * pi)); Z normalises over the circle. Both
# formulas are those of ?sdm ("Activation functions", "What `c` means") and
# of the SDM article (the approximated Bessel normaliser, not the exact
# I0 form).
#
# How this differs from bmm: the draws invert the cumulative distribution of
# that density on a fine grid, computed here; bmm's rsdm(), dsdm(), the
# normaliser it uses in Stan (Gauss-Legendre panels) and the likelihood are not
# called. A bug in bmm's normaliser or in its scaling would therefore show up
# as a recovery bias instead of cancelling.

# the von Mises activation, 1 at the target
sdv_kernel <- function(x, kappa) exp(kappa * (cos(x) - 1))

# the activation of the target under scaling = "density" (sqrtexp form)
sdv_cstar <- function(c, kappa) c * sqrt(kappa / (2 * pi))

# Inverse-CDF draws from the SDM response density on a fine grid. The density
# is exp(cstar * s(x)) / Z on the circle (x = response error minus the location
# mu), with s the activation kernel, equal to 1 at the target; Z is the
# integral of the numerator over the circle, taken here from the same grid
# (trapezoid rule), so no closed-form normaliser is used. The grid is dense at
# the target, where a large c or kappa puts all the mass: x = pi * t^2 on both
# sides, t from 0 to 1 in 6000 steps, so the spacing is 1e-3 at the edges and
# 1e-4 within 0.05 rad of the target. A uniform draw is inverted by linear
# interpolation of the cumulative sum.
sdv_grid <- local({
  t <- seq(0, 1, length.out = 6001L)
  c(-rev(pi * t[-1L]^2), pi * t^2)
})

# n draws for one set of parameters (cstar, kappa scalars; mu in radians)
sdv_draw <- function(n, cstar, kappa, mu = 0) {
  x <- sdv_grid
  logf <- cstar * sdv_kernel(x, kappa)
  f <- exp(logf - max(logf))
  cdf <- c(0, cumsum(diff(x) * (f[-1L] + f[-length(f)]) / 2))
  y <- stats::approx(cdf / cdf[length(cdf)], x, xout = stats::runif(n), ties = "ordered", rule = 2)$y
  (y + mu + pi) %% (2 * pi) - pi
}

# one response error per trial, for trial-wise c and kappa (a subject has one
# value of each, so the draws are made by distinct pair)
sdv_respond <- function(c, kappa, mu = 0) {
  y <- numeric(length(c))
  pair <- paste(c, kappa)
  for (p in unique(pair)) {
    rows <- which(pair == p)
    y[rows] <- sdv_draw(length(rows), sdv_cstar(c[rows[1L]], kappa[rows[1L]]), kappa[rows[1L]], mu)
  }
  y
}

# one subject's trials, for the recovery grid; c and kappa arrive on their
# natural scale, mu (fixed at 0) as given
sdv_generator <- function(pars, n_trials, model) {
  out <- data.frame(sdv_respond(rep(pars$c, n_trials), rep(pars$kappa, n_trials), pars$mu))
  names(out) <- model$resp_vars$resp_error
  out
}

# the calibration's data sets: bmmtools hands over one prior draw under the
# fit's names (c and kappa on the log scale) and the design data; the
# responses are drawn again from the subject values of that draw
sdv_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    exp(unname(draws[[sprintf("b_%s_Intercept", par)]] +
      draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))]))
  }
  data$y <- sdv_respond(subject("c"), subject("kappa"))
  data
}

bmmtools::validation_plan(
  bmm::sdm(resp_error = "y", activation = "vonmises"),
  # log scale. kappa: log 2.6 (kappa 13.5) is 1.1 scale units above the
  # centre of the default prior (student_t(5, 1.75, 0.75)); c: log 1.4 (c 4.1,
  # c* 6.0 at this kappa) is 0.8 scale units below the centre of
  # student_t(5, 2, 0.75), so the target is clearly above the floor but the
  # response error is wider than at the prior centre
  pars = c(kappa = 2.6, c = 1.4),
  sds = c(kappa = 0.3, c = 0.3),
  generator = sdv_generator,
  # prior = NULL: bmm's default priors (the sdm cell's own main priors)
  sbc = list(
    prior = NULL, generator = sdv_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
