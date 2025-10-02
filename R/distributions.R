#' Rejection Sampling
#'
#' Performs rejection sampling to generate samples from a target distribution.
#'
#' @param n Integer. The number of samples to generate.
#' @param f Function. The target density function from which to sample.
#' @param max_f Numeric. The maximum value of the target density function `f`.
#' @param proposal_fun Function. A function that generates samples from the proposal distribution.
#' @param ... Additional arguments to be passed to the target density function `f`.
#'
#' @return A numeric vector of length `n` containing samples from the target distribution.
#' @export
#' @keywords distribution
#'
#' @examples
#' target_density <- function(x) brms::dvon_mises(x, mu = 0, kappa = 10)
#' proposal <- function(n) runif(n, min = -pi, max = pi)
#' samples <- rejection_sampling(10000, target_density, max_f = target_density(0), proposal)
#' hist(samples, freq = FALSE)
#' curve(target_density, col = "red", add = TRUE)
rejection_sampling <- function(n, f, max_f, proposal_fun, ...) {
  stopifnot(is.numeric(n), length(n) == 1, n > 0)
  stopifnot(is.numeric(max_f), length(max_f) == 1 | length(max_f) == n, max_f > 0)

  inner <- function(n, f, max_f, proposal_fun, ..., acc = c()) {
    if (length(acc) > n) {
      return(acc[seq_len(n)])
    }
    x <- proposal_fun(n)
    y <- stats::runif(n) * max_f
    accept <- y < f(x, ...)
    inner(n, f, max_f, proposal_fun, ..., acc = c(acc, x[accept]))
  }

  inner(n, f, max_f, proposal_fun, ...)
}

#' @title Distribution functions for the Signal Discrimination Model (SDM)
#'
#' @description Density, distribution function, and random generation for the
#'   Signal Discrimination Model (SDM) Distribution with location `mu`,
#'   memory strength `c`, and precision `kappa`. Currently only a
#'   single activation source is supported.
#'
#' @name SDMdist
#'
#' @param x Vector of quantiles
#' @param q Vector of quantiles
#' @param p Vector of probabilities
#' @param n Number of observations to sample
#' @param mu Vector of location values in radians
#' @param c Vector of memory strength values
#' @param kappa Vector of precision values
#' @param log Logical; if `TRUE`, values are returned on the log scale.
#' @param parametrization Character; either `"bessel"` or `"sqrtexp"`
#'   (default). See [the online article](https://venpopov.github.io/bmm/articles/bmm_sdm_simple.html) for details on the
#'   parameterization.
#' @param log.p Logical; if `TRUE`, probabilities are returned on the log
#'   scale.
#' @param lower.bound Numeric; Lower bound of integration for the cumulative
#'   distribution
#' @param lower.tail Logical; If `TRUE` (default), return P(X <= x). Else,
#'   return P(X > x)
#' @keywords distribution
#'
#' @references Oberauer, K. (2023). Measurement models for visual working
#'   memory - A factorial model comparison. Psychological Review, 130(3), 841–852
#'
#' @return `dsdm` gives the density, `psdm` gives the distribution
#'   function, `qsdm` gives the quantile function, `rsdm` generates
#'   random deviates, and `.dsdm_integrate` is a helper function for
#'   calculating the density of the SDM distribution.
#'
#' @details **Parametrization**
#'
#' See [the online article](https://venpopov.github.io/bmm/articles/bmm_sdm_simple.html) for details on the parameterization.
#' Oberauer (2023) introduced the SDM with the bessel parametrization. The
#' sqrtexp parametrization is the default in the `bmm` package for
#' numerical stability and efficiency. The two parametrizations are related by
#' the functions `c_bessel2sqrtexp()` and `c_sqrtexp2bessel()`.
#'
#' **The cumulative distribution function**
#'
#' Since responses are on the circle, the cumulative distribution function
#' requires you to choose a lower bound of integration. The default is
#' \eqn{-\pi}, as for the brms::pvon_mises() function but you can choose any
#' value in the argument `lower_bound` of `psdm`. Another useful
#' choice is the mean of the response distribution minus \eqn{\pi}, e.g.
#' `lower_bound = mu-pi`. This is the default in
#' `circular::pvonmises()`, and it ensures that 50% of the cumulative
#' probability mass is below the mean of the response distribution.
#'
#' @export
#'
#' @examples
#' # plot the density of the SDM distribution
#' x <- seq(-pi, pi, length.out = 10000)
#' plot(x, dsdm(x, 0, 2, 3),
#'   type = "l", xlim = c(-pi, pi), ylim = c(0, 1),
#'   xlab = "Angle error (radians)",
#'   ylab = "density",
#'   main = "SDM density"
#' )
#' lines(x, dsdm(x, 0, 9, 1), col = "red")
#' lines(x, dsdm(x, 0, 2, 8), col = "green")
#' legend("topright", c(
#'   "c=2, kappa=3.0, mu=0",
#'   "c=9, kappa=1.0, mu=0",
#'   "c=2, kappa=8, mu=1"
#' ),
#' col = c("black", "red", "green"), lty = 1, cex = 0.8
#' )
#'
#' # plot the cumulative distribution function of the SDM distribution
#' p <- psdm(x, mu = 0, c = 3.1, kappa = 5)
#' plot(x, p, type = "l")
#'
#' # generate random deviates from the SDM distribution and overlay the density
#' r <- rsdm(10000, mu = 0, c = 3.1, kappa = 5)
#' d <- dsdm(x, mu = 0, c = 3.1, kappa = 5)
#' hist(r, breaks = 60, freq = FALSE)
#' lines(x, d, type = "l", col = "red")
dsdm <- function(x, mu = 0, c = 3, kappa = 3.5, log = FALSE,
                 parametrization = "sqrtexp") {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(c < 0)), "c must be non-negative")

  .dsdm_numer <- switch(parametrization,
    bessel = .dsdm_numer_bessel,
    sqrtexp = .dsdm_numer_sqrtexp,
    stop("Parametrization must be one of 'bessel' or 'sqrtexp'")
  )

  lnumerator <- .dsdm_numer(x, mu, c, kappa, log = TRUE)

  denom <- if (any(length(mu) > 1, length(c) > 1, length(kappa) > 1)) {
    .dsdm_integrate_numer_v(.dsdm_numer, mu, c, kappa, lower = mu, upper = mu + pi)
  } else {
    .dsdm_integrate_numer(.dsdm_numer, mu, c, kappa, lower = mu, upper = mu + pi)
  }

  denom <- 2 * denom

  if (!log) {
    return(exp(lnumerator) / denom)
  }
  lnumerator - log(denom)
}

#' @rdname SDMdist
#' @export
psdm <- function(q, mu = 0, c = 3, kappa = 3.5, lower.tail = TRUE, log.p = FALSE,
                 lower.bound = -pi, parametrization = "sqrtexp") {
  # parts adapted from brms::pvon_mises
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(c < 0)), "c must be non-negative")

  pi <- base::pi
  pi2 <- 2 * pi
  q <- (q + pi) %% pi2
  mu <- (mu + pi) %% pi2
  lower.bound <- (lower.bound + pi) %% pi2

  .dsdm_integrate <- function(mu, c, kappa, lower, upper, parametrization) {
    stats::integrate(dsdm,
      lower = lower, upper = upper, mu, c, kappa,
      parametrization = parametrization
    )$value
  }

  .dsdm_integrate_v <- Vectorize(.dsdm_integrate)

  if (any(length(q) > 1, length(mu) > 1, length(c) > 1, length(kappa) > 1)) {
    out <- .dsdm_integrate_v(mu, c, kappa,
      lower = lower.bound, upper = q,
      parametrization = parametrization
    )
  } else {
    out <- .dsdm_integrate(mu, c, kappa,
      lower = lower.bound, upper = q,
      parametrization = parametrization
    )
  }

  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

#' @rdname SDMdist
#' @export
qsdm <- function(p, mu = 0, c = 3, kappa = 3.5, parametrization = "sqrtexp") {
  .NotYetImplemented()
}

#' @rdname SDMdist
#' @export
rsdm <- function(n, mu = 0, c = 3, kappa = 3.5, parametrization = "sqrtexp") {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(c < 0)), "c must be non-negative")
  stopif(length(n) > 1, "n must be a single integer")

  .dsdm_numer <- switch(parametrization,
    bessel = .dsdm_numer_bessel,
    sqrtexp = .dsdm_numer_sqrtexp,
    stop("Parametrization must be one of 'bessel' or 'sqrtexp'")
  )

  rejection_sampling(
    n = n,
    f = function(x) .dsdm_numer(x, mu, c, kappa),
    max_f = .dsdm_numer(0, 0, c, kappa),
    proposal_fun = function(n) stats::runif(n, -pi, pi)
  )
}

# helper functions for calculating the density of the SDM distribution
.dsdm_numer_bessel <- function(x, mu, c, kappa, log = FALSE) {
  be <- besselI(kappa, nu = 0, expon.scaled = TRUE)
  out <- c * exp(kappa * (cos(x - mu) - 1)) / (2 * pi * be)
  if (!log) {
    out <- exp(out)
  }
  out
}

.dsdm_numer_sqrtexp <- function(x, mu, c, kappa, log = FALSE) {
  out <- c * exp(kappa * (cos(x - mu) - 1)) * sqrt(kappa) / sqrt(2 * pi)
  if (!log) {
    out <- exp(out)
  }
  out
}

.dsdm_integrate_numer <- function(fun, mu, c, kappa, lower, upper) {
  stats::integrate(fun, lower = lower, upper = upper, mu, c, kappa)$value
}

.dsdm_integrate_numer_v <- Vectorize(.dsdm_integrate_numer,
  vectorize.args = c("mu", "c", "kappa", "lower", "upper")
)


#' @title Distribution functions for the two-parameter mixture model (mixture2p)
#'
#' @description Density, distribution, and random generation functions for the
#'   two-parameter mixture model with the location of `mu`, precision of memory
#'   representations `kappa` and probability of recalling items from memory
#'   `p_mem`.
#'
#' @name mixture2p_dist
#'
#' @param x Vector of observed responses
#' @param q Vector of quantiles
#' @param p Vector of probability
#' @param n Number of observations to generate data for
#' @param mu Vector of locations
#' @param kappa Vector of precision values
#' @param p_mem Vector of probabilities for memory recall
#' @param log Logical; if `TRUE`, values are returned on the log scale.
#'
#' @keywords distribution
#'
#' @references Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
#'   representations in visual working memory. Nature, 453.
#'
#' @return `dmixture2p` gives the density of the two-parameter mixture model,
#'   `pmixture2p` gives the cumulative distribution function of the
#'   two-parameter mixture model, `qmixture2p` gives the quantile function of
#'   the two-parameter mixture model, and `rmixture2p` gives the random
#'   generation function for the two-parameter mixture model.
#'
#' @export
#'
#' @examples
#' # generate random samples from the mixture2p model and overlay the density
#' r <- rmixture2p(10000, mu = 0, kappa = 4, p_mem = 0.8)
#' x <- seq(-pi, pi, length.out = 10000)
#' d <- dmixture2p(x, mu = 0, kappa = 4, p_mem = 0.8)
#' hist(r, breaks = 60, freq = FALSE)
#' lines(x, d, type = "l", col = "red")
#'
dmixture2p <- function(x, mu = 0, kappa = 5, p_mem = 0.6, log = FALSE) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(p_mem < 0)), "p_mem must be larger than zero.")
  stopif(isTRUE(any(p_mem > 1)), "p_mem must be smaller than one.")

  density <- matrix(data = NaN, nrow = length(x), ncol = 2)

  density[, 1] <- log(p_mem) + brms::dvon_mises(x = x, mu = mu, kappa = kappa, log = T)
  density[, 2] <- log(1 - p_mem) + brms::dvon_mises(x = x, mu = 0, kappa = 0, log = T)

  density <- matrixStats::rowLogSumExps(density)

  if (!log) {
    return(exp(density))
  }

  density
}

#' @rdname mixture2p_dist
#' @export
pmixture2p <- function(q, mu = 0, kappa = 7, p_mem = 0.8) {
  .NotYetImplemented()
}

#' @rdname mixture2p_dist
#' @export
qmixture2p <- function(p, mu = 0, kappa = 5, p_mem = 0.6) {
  .NotYetImplemented()
}

#' @rdname mixture2p_dist
#' @export
rmixture2p <- function(n, mu = 0, kappa = 5, p_mem = 0.6) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(p_mem < 0)), "p_mem must be larger than zero.")
  stopif(isTRUE(any(p_mem > 1)), "p_mem must be smaller than one.")

  rejection_sampling(
    n = n,
    f = function(x) dmixture2p(x, mu, kappa, p_mem),
    max_f = dmixture2p(0, 0, kappa, p_mem),
    proposal_fun = function(n) stats::runif(n, -pi, pi)
  )
}

#' @title Distribution functions for the three-parameter mixture model (mixture3p)
#'
#' @description Density, distribution, and random generation functions for the
#'   three-parameter mixture model with the location of `mu`, precision of
#'   memory representations `kappa`, probability of recalling items from memory
#'   `p_mem`, and probability of recalling non-targets `p_nt`.
#'
#' @name mixture3p_dist
#'
#' @param x Vector of observed responses
#' @param q Vector of quantiles
#' @param p Vector of probability
#' @param n Number of observations to generate data for
#' @param mu Vector of locations. First value represents the location of the
#'   target item and any additional values indicate the location of non-target
#'   items.
#' @param kappa Vector of precision values
#' @param p_mem Vector of probabilities for memory recall
#' @param p_nt Vector of probabilities for swap errors
#' @param log Logical; if `TRUE`, values are returned on the log scale.
#'
#' @keywords distribution
#'
#' @references Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The
#'   precision of visual working memory is set by allocation of a shared
#'   resource. Journal of Vision, 9(10), 7.
#'
#' @return `dmixture3p` gives the density of the three-parameter mixture model,
#'   `pmixture3p` gives the cumulative distribution function of the
#'   two-parameter mixture model, `qmixture3p` gives the quantile function of
#'   the two-parameter mixture model, and `rmixture3p` gives the random
#'   generation function for the two-parameter mixture model.
#'
#' @export
#'
#' @examples
#' # generate random samples from the mixture3p model and overlay the density
#' r <- rmixture3p(10000, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
#' x <- seq(-pi, pi, length.out = 10000)
#' d <- dmixture3p(x, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
#' hist(r, breaks = 60, freq = FALSE)
#' lines(x, d, type = "l", col = "red")
#'
dmixture3p <- function(x, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2, log = FALSE) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(p_mem < 0)), "p_mem must be larger than zero.")
  stopif(isTRUE(any(p_nt < 0)), "p_nt must be larger than zero.")
  stopif(isTRUE(any(p_mem + p_nt > 1)), "The sum of p_mem and p_nt must be smaller than one.")

  density <- matrix(data = NaN, nrow = length(x), ncol = length(mu) + 1)
  probs <- c(
    p_mem,
    rep(p_nt / (length(mu) - 1), each = length(mu) - 1),
    (1 - p_mem - p_nt)
  )

  for (i in 1:(length(mu))) {
    density[, i] <- log(probs[i]) +
      brms::dvon_mises(x = x, mu = mu[i], kappa = kappa, log = T)
  }

  density[, length(mu) + 1] <- log(probs[length(mu) + 1]) +
    stats::dunif(x = x, -pi, pi, log = T)

  density <- matrixStats::rowLogSumExps(density)

  if (!log) {
    return(exp(density))
  }

  density
}

#' @rdname mixture3p_dist
#' @export
pmixture3p <- function(q, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2) {
  .NotYetImplemented()
}

#' @rdname mixture3p_dist
#' @export
qmixture3p <- function(p, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2) {
  .NotYetImplemented()
}

#' @rdname mixture3p_dist
#' @export
rmixture3p <- function(n, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(p_mem < 0)), "p_mem must be larger than zero.")
  stopif(isTRUE(any(p_nt < 0)), "p_nt must be larger than zero.")
  stopif(isTRUE(any(p_mem + p_nt > 1)), "The sum of p_mem and p_nt must be smaller than one.")

  xm <- seq(-pi, pi, length.out = 361)
  max_y <- max(dmixture3p(xm, mu, kappa, p_mem, p_nt))

  rejection_sampling(
    n = n,
    f = function(x) dmixture3p(x, mu, kappa, p_mem, p_nt),
    max_f = max_y,
    proposal_fun = function(n) stats::runif(n, -pi, pi)
  )
}

#' @title Distribution functions for the Interference Measurement Model (IMM)
#'
#' @description Density, distribution, and random generation functions for the
#'   interference measurement model with the location of `mu`, strength of cue-
#'   dependent activation `c`, strength of cue-independent activation `a`, the
#'   generalization gradient `s`, and the precision of memory representations
#'   `kappa`.
#'
#' @name IMMdist
#'
#' @param x Vector of observed responses
#' @param q Vector of quantiles
#' @param p Vector of probability
#' @param n Number of observations to generate data for
#' @param mu Vector of locations
#' @param dist Vector of distances of the item locations to the cued location
#' @param kappa Vector of precision values
#' @param c Vector of strengths for cue-dependent activation
#' @param a Vector of strengths for cue-independent activation
#' @param s Vector of generalization gradients
#' @param b Vector of baseline activation
#' @param log Logical; if `TRUE`, values are returned on the log scale.
#'
#' @keywords distribution
#'
#' @references Oberauer, K., Stoneking, C., Wabersich, D., & Lin, H.-Y. (2017).
#'   Hierarchical Bayesian measurement models for continuous reproduction of
#'   visual features from working memory. Journal of Vision, 17(5), 11.
#'
#' @return `dimm` gives the density of the interference measurement model,
#'   `pimm` gives the cumulative distribution function of the interference
#'   measurement model, `qimm` gives the quantile function of the interference
#'   measurement model, and `rimm` gives the random generation function for the
#'   interference measurement model.
#'
#' @export
#'
#' @examples
#' # generate random samples from the imm and overlay the density
#' r <- rimm(10000,
#'   mu = c(0, 2, -1.5), dist = c(0, 0.5, 2),
#'   c = 5, a = 2, s = 2, b = 1, kappa = 4
#' )
#' x <- seq(-pi, pi, length.out = 10000)
#' d <- dimm(x,
#'   mu = c(0, 2, -1.5), dist = c(0, 0.5, 2),
#'   c = 5, a = 2, s = 2, b = 1, kappa = 4
#' )
#' hist(r, breaks = 60, freq = FALSE)
#' lines(x, d, type = "l", col = "red")
#'
dimm <- function(x, mu = c(0, 2, -1.5), dist = c(0, 0.5, 2),
                 c = 5, a = 2, b = 1, s = 2, kappa = 5, log = FALSE) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  len_mu <- length(mu)
  stopif(
    len_mu != length(dist),
    "The number of items does not match the distances provided from the cued location."
  )
  stopif(isTRUE(any(s < 0)), "s must be non-negative")
  stopif(isTRUE(any(dist < 0)), "all distances have to be positive.")

  # compute activation for all items
  weights <- rep(c, len_mu) * exp(-s * dist) + rep(a, len_mu)

  # add activation of background noise
  weights <- c(weights, b)

  # compute probability for responding stemming from each distribution
  probs <- weights / sum(weights)
  density <- matrix(data = NaN, nrow = length(x), ncol = len_mu + 1)

  for (i in seq_along(mu)) {
    density[, i] <- log(probs[i]) +
      brms::dvon_mises(x, mu = mu[i], kappa = kappa, log = T)
  }

  density[, len_mu + 1] <- log(probs[len_mu + 1]) +
    stats::dunif(x = x, -pi, pi, log = T)

  density <- matrixStats::rowLogSumExps(density)

  if (!log) {
    return(exp(density))
  }

  density
}

#' @rdname IMMdist
#' @export
pimm <- function(q, mu = c(0, 2, -1.5), dist = c(0, 0.5, 2),
                 c = 1, a = 0.2, b = 0, s = 2, kappa = 5) {
  .NotYetImplemented()
}

#' @rdname IMMdist
#' @export
qimm <- function(p, mu = c(0, 2, -1.5), dist = c(0, 0.5, 2),
                 c = 1, a = 0.2, b = 0, s = 2, kappa = 5) {
  .NotYetImplemented()
}

#' @rdname IMMdist
#' @export
rimm <- function(n, mu = c(0, 2, -1.5), dist = c(0, 0.5, 2),
                 c = 1, a = 0.2, b = 1, s = 2, kappa = 5) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(s < 0)), "s must be non-negative")
  stopif(isTRUE(any(dist < 0)), "all distances have to be positive.")
  stopif(
    length(mu) != length(dist),
    "The number of items does not match the distances provided from the cued location."
  )

  xm <- seq(-pi, pi, length.out = 361)
  max_y <- max(dimm(xm, mu, dist, c, a, b, s, kappa))

  rejection_sampling(
    n = n,
    f = function(x) dimm(x, mu, dist, c, a, b, s, kappa),
    max_f = max_y,
    proposal_fun = function(n) stats::runif(n, -pi, pi)
  )
}

#' @title Distribution functions for the Memory Measurement Model (M3)
#'
#' @description Density and random generation functions for the memory
#'   measurement model. Please note that these functions are currently not
#'   vectorized.
#'
#' @name m3dist
#'
#' @param x Integer vector of length `K` where K is the number of response categories
#'   and each value is the number of observed responses per category
#' @param n Integer. Number of observations to generate data for
#' @param size The total number of observations in all categories
#' @param pars A named vector of parameters of the memory measurement model
#' @param m3_model A `bmmodel` object specifying the m3 model that densities or
#'   random samples should be generated for
#' @param act_funs A `bmmformula` object specifying the activation functions for
#'   the different response categories for the "custom" version of the M3. The
#'   default will attempt to construct the standard activation functions for the
#'   "ss" and "cs" model version. For a custom m3 model you need to specify the
#'   act_funs argument manually
#' @param log Logical; if `TRUE` (default), values are returned on the log scale.
#' @param ... can be used to pass additional variables that are used in the
#'   activation functions, but not parameters of the model
#'
#' @keywords distribution
#'
#' @references Oberauer, K., & Lewandowsky, S. (2019). Simple measurement models
#'   for complex working-memory tasks. Psychological Review, 126(6), 880–932.
#'   https://doi.org/10.1037/rev0000159
#'
#' @return `dm3` gives the density of the memory measurement model, and `rm3`
#'   gives the random generation function for the memory measurement model.
#'
#' @examples
#'   model <- m3(
#'    resp_cats = c("corr", "other", "npl"),
#'    num_options = c(1, 4, 5),
#'    choice_rule = "simple",
#'    version = "ss"
#'  )
#'  dm3(x = c(20, 10, 10), pars = c(a = 1, b = 1, c = 2), m3_model = model)
#' @export
dm3 <- function(x, pars, m3_model, act_funs = construct_m3_act_funs(m3_model, warnings = FALSE),
                log = TRUE, ...) {
  probs <- .compute_m3_probability_vector(pars, m3_model, act_funs, ...)
  dmultinom(x, prob = probs, log = log)
}

#' @rdname m3dist
#' @export
rm3 <- function(n, size, pars, m3_model, act_funs = construct_m3_act_funs(m3_model, warnings = FALSE),
                ...) {
  probs <- .compute_m3_probability_vector(pars, m3_model, act_funs, ...)
  t(rmultinom(n, size = size, prob = probs))
}

.compute_m3_probability_vector <-
  function(pars, m3_model, act_funs = construct_m3_act_funs(m3_model, warnings = FALSE), ...) {
    pars <- c(pars, unlist(list(...)))
    stopif(
      is_try_error(try(act_funs, silent = TRUE)),
      'No activation functions for version "custom" provided.
      Please pass activation functions for the different response categories
      using the "act_funs" argument.'
    )
    stopif(
      !identical(sort(rhs_vars(act_funs)), sort(names(pars))),
      'The names or number of parameters used in the activation functions mismatch the names or number
      of parameters ("pars") and additional arguments (i.e. ...) passed to the function.'
    )

    acts <- sapply(act_funs, function(pform) eval(pform[[length(pform)]], envir = as.list(pars)))

    num_options <- m3_model$other_vars$num_options
    choice_rule <- tolower(m3_model$other_vars$choice_rule)
    if (choice_rule == "softmax") acts <- exp(acts)
    acts <- acts * num_options
    acts / sum(acts)
  }


#' @title Distribution function for the Diffusion Decision Model (`ddm`)
#'
#' @description
#'   Density and random generation function for the Diffusion Decision Model.
#'
#' @name ddm_dist
#'
#' @param rt Vector of response times for which the density should be returned
#' @param response Vector of responses for which the density should be returned
#' @param n Number of random samples to generate
#' @param drift Drift rates of the ddm
#' @param bound Boundary separation of the ddm
#' @param ndt Non-decision time of the ddm
#' @param zr relative starting point of the ddm
#' @param sdrift Trial-to-trial variability of the drift rate
#' @param sndt Trial-to-trial variability of the non-decision time
#' @param szr Trial-to-trial variability of the relative starting point
#' @param log Logical, indicating if log-densities should be returned (default = TRUE)
#'
#' @keywords distributions
#'
#' @export
dddm <- function(rt, response,drift, bound, ndt, zr = 0.5, sdrift = 0, sndt = 0, szr = 0, log = TRUE) {

  stopif(any(rt < 0),
         glue("Negative RTs are not allowed.\n",
              "Please check your rt variable and pass only positive RTs."))

  if(!is.character(response)){
    response <- ifelse(response == 1, "upper","lower")
  }

  stopif(any(!response %in% c("upper","lower")),
         glue("Invalid responses passed! Either pass a numeric vector with 0 as a lower bound response,",
              "and 1 as upper bound response, or a character vector with \"upper\" and \"lower\" for the respecitve,",
              "responses given."))

  stopif(length(rt) != length(response),
         glue("Different number of rts and responses passed to dddm. Please pass vectors of equal length."))

  out <- rtdists::ddiffusion(rt = rt, response = response,
                             a = bound, v = drift, t0 = ndt, z = zr * bound, sz = szr, sv = sdrift, st0 = sndt)

  if(log) out <- log(out)
  out
}

#' @name ddm_dist
#' @export
rddm <- function(n, drift, bound, ndt, zr = 0.5, sdrift = 0, sndt = 0, szr = 0) {
  max_len <- max(lengths(list(drift, bound, ndt, zr, sdrift, sndt, szr)))

  if (max_len > 1L) {
    if (!n %in% c(1, max_len)) {
      stop2("Can only sample exactly once for each condition.")
    }
    n <- max_len
  }

  sim_data <- rtdists::rdiffusion(n = n, a = bound, v = drift, t0 = ndt, z = zr * bound, sz = szr, sv = sdrift, st0 = sndt)
  sim_data$response <- ifelse(sim_data$response == "upper",1,0)
  sim_data
}
