# Recovery plan for sdm() with the Laplace activation (activation = "laplace")
# and its default scaling, "peak" (bmm 1.5.0 sets the default per activation:
# "density" for the von Mises, "peak" for the Laplace). The response error x
# has the density exp(c* s(x)) / Z on the circle, with the kernel
# s(x) = exp(-kappa * |x|), |x| the angular distance, and the activation of the
# target c* = c under "peak" (c is the activation of the target); Z normalises
# over the circle. The formulas are those of ?sdm ("Activation functions",
# "What `c` means").
#
# How this differs from bmm: the draws invert the cumulative distribution of
# that density on a fine grid, computed here; bmm's rsdm(), dsdm(), the
# normaliser it uses in Stan and the likelihood are not called. A bug in bmm's
# normaliser or in its scaling would therefore show up as a recovery bias
# instead of cancelling. The kernel has a kink at the target, so the grid is
# dense there.

# the Laplace activation, 1 at the target
sdl_kernel <- function(x, kappa) exp(-kappa * abs(x))

# the activation of the target under scaling = "peak"
sdl_cstar <- function(c, kappa) c

# Inverse-CDF draws from the SDM response density on a fine grid. The density
# is exp(cstar * s(x)) / Z on the circle (x = response error minus the location
# mu), with s the activation kernel, equal to 1 at the target; Z is the
# integral of the numerator over the circle, taken here from the same grid
# (trapezoid rule), so no closed-form normaliser is used. The grid is dense at
# the target, where a large c or kappa puts all the mass: x = pi * t^2 on both
# sides, t from 0 to 1 in 6000 steps, so the spacing is 1e-3 at the edges and
# 1e-4 within 0.05 rad of the target. A uniform draw is inverted by linear
# interpolation of the cumulative sum.
sdl_grid <- local({
  t <- seq(0, 1, length.out = 6001L)
  c(-rev(pi * t[-1L]^2), pi * t^2)
})

# n draws for one set of parameters (cstar, kappa scalars; mu in radians)
sdl_draw <- function(n, cstar, kappa, mu = 0) {
  x <- sdl_grid
  logf <- cstar * sdl_kernel(x, kappa)
  f <- exp(logf - max(logf))
  cdf <- c(0, cumsum(diff(x) * (f[-1L] + f[-length(f)]) / 2))
  y <- stats::approx(cdf / cdf[length(cdf)], x, xout = stats::runif(n), ties = "ordered", rule = 2)$y
  (y + mu + pi) %% (2 * pi) - pi
}

# one response error per trial, for trial-wise c and kappa (a subject has one
# value of each, so the draws are made by distinct pair)
sdl_respond <- function(c, kappa, mu = 0) {
  y <- numeric(length(c))
  pair <- paste(c, kappa)
  for (p in unique(pair)) {
    rows <- which(pair == p)
    y[rows] <- sdl_draw(length(rows), sdl_cstar(c[rows[1L]], kappa[rows[1L]]), kappa[rows[1L]], mu)
  }
  y
}

# one subject's trials, for the recovery grid; c and kappa arrive on their
# natural scale, mu (fixed at 0) as given
sdl_generator <- function(pars, n_trials, model) {
  out <- data.frame(sdl_respond(rep(pars$c, n_trials), rep(pars$kappa, n_trials), pars$mu))
  names(out) <- model$resp_vars$resp_error
  out
}

# the calibration's data sets: bmmtools hands over one prior draw under the
# fit's names (c and kappa on the log scale) and the design data; the
# responses are drawn again from the subject values of that draw
sdl_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    exp(unname(draws[[sprintf("b_%s_Intercept", par)]] +
      draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))]))
  }
  data$y <- sdl_respond(subject("c"), subject("kappa"))
  data
}

bmmtools::validation_plan(
  bmm::sdm(resp_error = "y", activation = "laplace"),
  # log scale. kappa: log 0.2 (kappa 1.2) is 0.9 scale units above the centre
  # of the default prior (student_t(5, -0.5, 0.75)); c: log 1.6 (c 5) is 1.2
  # scale units below the centre of student_t(5, 2.5, 0.75), a weaker target
  # activation than the prior centre (c 12)
  pars = c(kappa = 0.2, c = 1.6),
  sds = c(kappa = 0.3, c = 0.3),
  generator = sdl_generator,
  sbc = list(
    prior = NULL, generator = sdl_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
