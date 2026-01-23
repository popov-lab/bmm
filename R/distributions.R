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


#' @title Distribution functions for the EZ-Diffusion Model (ezdm)
#'
#' @description Density and random generation functions for the EZ-Diffusion
#'   Model. The model operates on aggregated data: mean
#'   reaction time, variance of reaction time, and number of responses to the
#'   upper boundary.
#'
#' @name ezdm_dist
#'
#' @param mean_rt Observed mean reaction time(s) in seconds. For version
#'   "3par", a numeric vector or single value. For version "4par", either a vector
#'   of length 2 (c(mean_rt_upper, mean_rt_lower)) for single observation, or a matrix with 2
#'   columns for multiple observations.
#' @param var_rt Observed variance of reaction times in seconds^2. For version
#'   "3par", a numeric vector or single value. For version "4par", either a vector
#'   of length 2 (c(var_rt_upper, var_rt_lower)) for single observation, or a matrix with 2
#'   columns for multiple observations.
#' @param n_upper Number of responses to the upper boundary
#' @param n_trials Total number of trials
#' @param drift Drift rate (positive, evidence accumulation rate).
#' @param bound Boundary separation (distance between decision thresholds).
#' @param ndt Non-decision time (seconds).
#' @param zr Relative starting point (0 to 1). Only used for version "4par".
#' @param s Diffusion constant (standard deviation of noise), default = 1.
#' @param version Character; either "3par" (default) or "4par"
#' @param n Number of samples to generate
#' @param log Logical; if `TRUE`, values are returned on the log scale.
#'
#' @keywords distribution
#'
#' @references
#' Wagenmakers, E.-J., Van Der Maas, H. L. J., & Grasman, R. P. P. P. (2007).
#'   An EZ-diffusion model for response time and accuracy. Psychonomic Bulletin
#'   & Review, 14(1), 3-22.
#'
#' Chávez De la Peña, A. F., & Vandekerckhove, J. (2025). An EZ Bayesian
#'   hierarchical drift diffusion model for response time and accuracy.
#'   Psychonomic Bulletin & Review.
#'
#' @return `dezdm` gives the log-density of the observed summary statistics
#'   under the EZDM, and `rezdm` generates random summary statistics from the
#'   implied sampling distributions.
#'
#' @export
#'
#' @examples
#' # 3-parameter version (single observation)
#' dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
#'       drift = 2, bound = 1.5, ndt = 0.3)
#'
#' # 3-parameter version (vectorized)
#' dezdm(mean_rt = c(0.5, 0.55), var_rt = c(0.02, 0.025),
#'       n_upper = c(80, 75), n_trials = c(100, 100),
#'       drift = 2, bound = 1.5, ndt = 0.3)
#'
#' # 4-parameter version (single observation)
#' dezdm(mean_rt = c(0.45, 0.55), var_rt = c(0.018, 0.025),
#'       n_upper = 80, n_trials = 100,
#'       drift = 2, bound = 1.5, ndt = 0.3, zr = 0.55, version = "4par")
#'
#' # generate random summary statistics
#' rezdm(n = 100, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3)
#' rezdm(n = 100, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3,
#'       zr = 0.55, version = "4par")
#'
dezdm <- function(mean_rt, var_rt, n_upper, n_trials,
                  drift, bound, ndt, zr = 0.5, s = 1,
                  version = "3par", log = TRUE) {
  # validate version
  stopif(
    not_in(version, c("3par", "4par")),
    "version must be either '3par' or '4par'"
  )

  # parameter validation
  stopif(isTRUE(any(bound <= 0)), "bound must be positive")
  stopif(isTRUE(any(ndt <= 0)), "ndt must be positive")
  stopif(isTRUE(any(s <= 0)), "s must be positive")
  stopif(isTRUE(any(n_trials <= 0)), "n_trials must be positive")
  stopif(isTRUE(any(n_upper < 0)), "n_upper cannot be negative")
  stopif(isTRUE(any(n_upper > n_trials)), "n_upper cannot exceed n_trials")

  if (version == "4par") {
    stopif(isTRUE(any(zr <= 0 | zr >= 1)), "zr must be between 0 and 1")
    if (!is.matrix(mean_rt)) {
      stopif(
        length(mean_rt) != 2,
        "For version '4par', mean_rt must be length 2 or a matrix with 2 cols"
      )
      stopif(
        length(var_rt) != 2,
        "For version '4par', var_rt must be length 2 or a matrix with 2 cols"
      )
    } else {
      stopif(ncol(mean_rt) != 2, "mean_rt matrix must have 2 columns")
      stopif(ncol(var_rt) != 2, "var_rt matrix must have 2 columns")
    }
    ll <- .dezdm_4par(mean_rt, var_rt, n_upper, n_trials,
                      drift, bound, ndt, zr, s)
  } else {
    ll <- .dezdm_3par(mean_rt, var_rt, n_upper, n_trials,
                      drift, bound, ndt, s)
  }

  if (!log) {
    return(exp(ll))
  }
  ll
}

#' @rdname ezdm_dist
#' @export
rezdm <- function(n, n_trials, drift, bound, ndt, zr = 0.5, s = 1,
                  version = "3par") {

  # validate version
  stopif(
    not_in(version, c("3par", "4par")),
    "version must be either '3par' or '4par'"
  )

  # parameter validation
  stopif(isTRUE(any(bound <= 0)), "bound must be positive")
  stopif(isTRUE(any(ndt <= 0)), "ndt must be positive")
  stopif(isTRUE(any(s <= 0)), "s must be positive")
  stopif(length(n) > 1, "n must be a single integer")

  if (version == "4par") {
    stopif(isTRUE(any(zr <= 0 | zr >= 1)), "zr must be between 0 and 1")
    .rezdm_4par(n, n_trials, drift, bound, ndt, zr, s)
  } else {
    .rezdm_3par(n, n_trials, drift, bound, ndt, s)
  }
}

# Internal: 3par density - vectorized
# Handles both: (1) vector observations with scalar parameters, and
#               (2) scalar observations with vector parameters (for log_lik)
.dezdm_3par <- function(mean_rt, var_rt, n_upper, n_trials,
                        drift, bound, ndt, s) {
  # determine common length from observations AND parameters
  # this allows log_lik to work (scalar obs, vector params)
  n <- max(length(mean_rt), length(var_rt), length(n_upper),
           length(n_trials), length(drift), length(bound),
           length(ndt), length(s))

  # recycle all inputs to common length
  mean_rt <- rep_len(mean_rt, n)
  var_rt <- rep_len(var_rt, n)
  n_upper <- rep_len(n_upper, n)
  n_trials <- rep_len(n_trials, n)
  ndt <- rep_len(ndt, n)
  drift <- rep_len(drift, n)
  bound <- rep_len(bound, n)
  s <- rep_len(s, n)

  # compute moments (already vectorized)
  moments <- .ezdm_moments_3par(drift, bound, s)
  p_c <- moments$pC
  mdt <- moments$MDT
  vrt <- moments$VRT

  # binomial for n_upper
  ll <- stats::dbinom(n_upper, size = n_trials, prob = p_c, log = TRUE)

  # normal for mean RT
  ll <- ll + stats::dnorm(mean_rt, mean = ndt + mdt,
                          sd = sqrt(vrt / n_trials), log = TRUE)

  # gamma for variance RT
  shape <- (n_trials - 1) / 2
  rate <- (n_trials - 1) / (2 * vrt)
  ll <- ll + stats::dgamma(var_rt, shape = shape, rate = rate, log = TRUE)

  ll
}

# Internal: 4par density - vectorized
# Handles both: (1) vector observations with scalar parameters, and
#               (2) scalar observations with vector parameters (for log_lik)
# mean_rt, var_rt: matrices with 2 columns (upper, lower) and n rows
#                  OR vectors of length 2 for single observation
# n_upper, n_trials: vectors of length n (or scalars for single obs)
# drift, bound, ndt, zr, s: scalars or vectors (recycled)
.dezdm_4par <- function(mean_rt, var_rt, n_upper, n_trials,
                        drift, bound, ndt, zr, s) {
  # handle single observation (vectors) vs multiple (matrices)
  if (is.matrix(mean_rt)) {
    n_obs <- nrow(mean_rt)
    mean_rt_upper <- mean_rt[, 1]
    mean_rt_lower <- mean_rt[, 2]
    var_rt_upper <- var_rt[, 1]
    var_rt_lower <- var_rt[, 2]
  } else {
    # single observation - length 2 vectors
    n_obs <- 1
    mean_rt_upper <- mean_rt[1]
    mean_rt_lower <- mean_rt[2]
    var_rt_upper <- var_rt[1]
    var_rt_lower <- var_rt[2]
  }

  # determine common length from observations AND parameters

  # this allows log_lik to work (scalar obs, vector params)
  n <- max(n_obs, length(n_upper), length(n_trials),
           length(drift), length(bound), length(ndt),
           length(zr), length(s))

  # recycle observation values to common length
  mean_rt_upper <- rep_len(mean_rt_upper, n)
  mean_rt_lower <- rep_len(mean_rt_lower, n)
  var_rt_upper <- rep_len(var_rt_upper, n)
  var_rt_lower <- rep_len(var_rt_lower, n)
  n_upper <- rep_len(n_upper, n)
  n_trials <- rep_len(n_trials, n)

  n_lower <- n_trials - n_upper
  moments <- .ezdm_moments_4par(drift, bound, zr, s)

  # recycle moments and ndt to common length
  pC <- rep_len(moments$pC, n)
  mdt_upper <- rep_len(moments$mdt_upper, n)
  mdt_lower <- rep_len(moments$mdt_lower, n)
  vrt_upper <- rep_len(moments$vrt_upper, n)
  vrt_lower <- rep_len(moments$vrt_lower, n)
  ndt <- rep_len(ndt, n)

  # binomial for n_upper
  ll <- stats::dbinom(n_upper, size = n_trials, prob = pC, log = TRUE)

  # upper boundary contributions (vectorized)
  upper_valid <- n_upper >= 2
  if (any(upper_valid)) {
    ll[upper_valid] <- ll[upper_valid] +
      stats::dnorm(mean_rt_upper[upper_valid],
                   mean = ndt[upper_valid] + mdt_upper[upper_valid],
                   sd = sqrt(vrt_upper[upper_valid] / n_upper[upper_valid]),
                   log = TRUE) +
      stats::dgamma(var_rt_upper[upper_valid],
                    shape = (n_upper[upper_valid] - 1) / 2,
                    rate = (n_upper[upper_valid] - 1) / (2 * vrt_upper[upper_valid]),
                    log = TRUE)
  }

  # lower boundary contributions (vectorized)
  lower_valid <- n_lower >= 2
  if (any(lower_valid)) {
    ll[lower_valid] <- ll[lower_valid] +
      stats::dnorm(mean_rt_lower[lower_valid],
                   mean = ndt[lower_valid] + mdt_lower[lower_valid],
                   sd = sqrt(vrt_lower[lower_valid] / n_lower[lower_valid]),
                   log = TRUE) +
      stats::dgamma(var_rt_lower[lower_valid],
                    shape = (n_lower[lower_valid] - 1) / 2,
                    rate = (n_lower[lower_valid] - 1) / (2 * vrt_lower[lower_valid]),
                    log = TRUE)
  }

  ll
}

# Internal: truncated normal sampling via rejection sampling
# Samples from N(mean, sd) truncated to [lower, Inf)
# @param n Number of samples
# @param mean Mean of the normal distribution (scalar or vector of length n)
# @param sd Standard deviation (scalar or vector of length n)
# @param lower Lower truncation bound (scalar)
# @param max_iter Maximum rejection sampling iterations (default 1000)
# @return Numeric vector of n samples >= lower
.rtruncnorm_lower <- function(n, mean, sd, lower, max_iter = 1000) {
  samples <- stats::rnorm(n, mean = mean, sd = sd)
  rejected <- samples < lower
  iter <- 0

  while (any(rejected) && iter < max_iter) {
    n_rejected <- sum(rejected)
    # Resample only rejected values, using corresponding mean/sd if vectorized
    if (length(mean) == 1) {
      samples[rejected] <- stats::rnorm(n_rejected, mean = mean, sd = sd)
    } else {
      samples[rejected] <- stats::rnorm(
        n_rejected,
        mean = mean[rejected],
        sd = sd[rejected]
      )
    }
    rejected <- samples < lower
    iter <- iter + 1
  }

  # Fallback: clamp any remaining rejected samples (should be extremely rare)
  if (any(rejected)) {
    samples[rejected] <- lower
  }

  samples
}

# Internal: 3par random generation
.rezdm_3par <- function(n, n_trials, drift, bound, ndt, s) {
  moments <- .ezdm_moments_3par(drift, bound, s)

  n_upper <- stats::rbinom(n, size = n_trials, prob = moments$pC)
  var_rt <- moments$VRT * stats::rchisq(n, df = n_trials - 1) / (n_trials - 1)

  # Use truncated normal to ensure mean_rt >= ndt
  mean_rt <- .rtruncnorm_lower(
    n = n,
    mean = ndt + moments$MDT,
    sd = sqrt(var_rt / n_trials),
    lower = ndt
  )

  data.frame(
    mean_rt = mean_rt,
    var_rt = var_rt,
    n_upper = n_upper,
    n_trials = rep(n_trials, n)
  )
}

# Internal: 4par random generation
.rezdm_4par <- function(n, n_trials, drift, bound, ndt, zr, s) {
  moments <- .ezdm_moments_4par(drift, bound, zr, s)

  n_upper <- stats::rbinom(n, size = n_trials, prob = moments$pC)
  n_lower <- n_trials - n_upper

  # pre-allocate
  mean_rt_upper <- var_rt_upper <- rep(NA_real_, n)
  mean_rt_lower <- var_rt_lower <- rep(NA_real_, n)

  # generate upper boundary statistics where n_upper >= 2
  idx_upper <- n_upper >= 2
  if (any(idx_upper)) {
    n_u <- n_upper[idx_upper]
    var_rt_upper[idx_upper] <- moments$vrt_upper *
      stats::rchisq(sum(idx_upper), df = n_u - 1) / (n_u - 1)
    # Use truncated normal to ensure mean_rt_upper >= ndt
    mean_rt_upper[idx_upper] <- .rtruncnorm_lower(
      n = sum(idx_upper),
      mean = ndt + moments$mdt_upper,
      sd = sqrt(var_rt_upper[idx_upper] / n_u),
      lower = ndt
    )
  }

  # generate lower boundary statistics where n_lower >= 2
  idx_lower <- n_lower >= 2
  if (any(idx_lower)) {
    n_l <- n_lower[idx_lower]
    var_rt_lower[idx_lower] <- moments$vrt_lower *
      stats::rchisq(sum(idx_lower), df = n_l - 1) / (n_l - 1)
    # Use truncated normal to ensure mean_rt_lower >= ndt
    mean_rt_lower[idx_lower] <- .rtruncnorm_lower(
      n = sum(idx_lower),
      mean = ndt + moments$mdt_lower,
      sd = sqrt(var_rt_lower[idx_lower] / n_l),
      lower = ndt
    )
  }

  data.frame(
    mean_rt_upper = mean_rt_upper,
    mean_rt_lower = mean_rt_lower,
    var_rt_upper = var_rt_upper,
    var_rt_lower = var_rt_lower,
    n_upper = n_upper,
    n_trials = rep(n_trials, n)
  )
}

# Internal: compute 3par moments (zr = 0.5) - vectorized
.ezdm_moments_3par <- function(drift, bound, s) {
  # pre-allocate based on longest input
  n <- max(length(drift), length(bound), length(s))

  # recycle to common length
  drift <- rep_len(drift, n)
  bound <- rep_len(bound, n)
  s <- rep_len(s, n)

  # initialize outputs
  pC <- rep(NA_real_, n)
  MDT <- rep(NA_real_, n)
  VRT <- rep(NA_real_, n)

  # identify near-zero drift cases
  zero_drift <- abs(drift) < 1e-6

  # zero-drift formulas
  if (any(zero_drift)) {
    pC[zero_drift] <- 0.5
    MDT[zero_drift] <- bound[zero_drift]^2 / (4 * s[zero_drift]^2)
    VRT[zero_drift] <- bound[zero_drift]^4 / (24 * s[zero_drift]^4)
  }

  # non-zero drift formulas
  if (any(!zero_drift)) {
    idx <- !zero_drift
    # Use signed drift for pC calculation
    y <- -(bound[idx] * drift[idx]) / s[idx]^2
    expy <- exp(y)
    pC[idx] <- 1 / (1 + expy)
    # Use soft absolute value: sqrt(drift^2 + tau^2) with tau = 0.01
    # This avoids extreme curvature while maintaining smoothness
    tau <- 0.01
    drift_abs <- sqrt(drift[idx]^2 + tau^2)
    y_abs <- -(bound[idx] * drift_abs) / s[idx]^2
    expy_abs <- exp(y_abs)
    MDT[idx] <- (bound[idx] / (2 * drift_abs)) * ((1 - expy_abs) / (1 + expy_abs))
    VRT[idx] <- ((bound[idx] * s[idx]^2) / (2 * drift_abs^3)) *
      (2 * y_abs * expy_abs - exp(2 * y_abs) + 1) / ((expy_abs + 1)^2)
  }

  list(pC = pC, MDT = MDT, VRT = VRT)
}

# Internal: compute 4par moments (Srivastava et al. formulas) - vectorized
.ezdm_moments_4par <- function(drift, bound, zr, s) {
  # helper functions
  coth <- function(x) cosh(x) / sinh(x)
  csch <- function(x) 1 / sinh(x)

  # pre-allocate based on longest input
  n <- max(length(drift), length(bound), length(zr), length(s))

  # recycle to common length
  drift <- rep_len(drift, n)
  bound <- rep_len(bound, n)
  zr <- rep_len(zr, n)
  s <- rep_len(s, n)

  # compute intermediate values
  z <- bound / 2
  x0 <- (zr * bound) - z

  # Use signed drift for pC calculation
  k_z_signed <- (drift * z) / s^2
  k_x_signed <- (drift * x0) / s^2
  
  # proportion correct (works for all cases)
  pC <- 1 - (exp(-2 * k_x_signed) - exp(-2 * k_z_signed)) / (exp(2 * k_z_signed) - exp(-2 * k_z_signed))
  
  # Use soft absolute value: sqrt(drift^2 + tau^2) with tau = 0.01
  # This provides smooth gradients without extreme curvature
  tau <- 0.01
  a <- sqrt(drift^2 + tau^2)
  k_z <- (a * z) / s^2
  k_x <- (a * x0) / s^2

  # initialize outputs
  mdt_upper <- rep(NA_real_, n)
  mdt_lower <- rep(NA_real_, n)
  vrt_upper <- rep(NA_real_, n)
  vrt_lower <- rep(NA_real_, n)

  # identify near-zero drift cases
  zero_drift <- abs(drift) < 1e-6

  # zero-drift formulas
  if (any(zero_drift)) {
    idx <- zero_drift
    mdt_upper[idx] <- (4 * z[idx]^2 - (z[idx] + x0[idx])^2) / (3 * s[idx]^2)
    mdt_lower[idx] <- (4 * z[idx]^2 - (z[idx] - x0[idx])^2) / (3 * s[idx]^2)
    vrt_upper[idx] <- (32 * z[idx]^4 - 2 * (z[idx] + x0[idx])^4) / (45 * s[idx]^4)
    vrt_lower[idx] <- (32 * z[idx]^4 - 2 * (z[idx] - x0[idx])^4) / (45 * s[idx]^4)
  }

  # non-zero drift formulas
  if (any(!zero_drift)) {
    idx <- !zero_drift
    mdt_upper[idx] <- (s[idx]^2 / a[idx]^2) *
      (2 * k_z[idx] * coth(2 * k_z[idx]) -
         (k_x[idx] + k_z[idx]) * coth(k_x[idx] + k_z[idx]))
    mdt_lower[idx] <- (s[idx]^2 / a[idx]^2) *
      (2 * k_z[idx] * coth(2 * k_z[idx]) -
         (-k_x[idx] + k_z[idx]) * coth(-k_x[idx] + k_z[idx]))

    vrt_upper[idx] <- (s[idx]^4 / a[idx]^4) *
      (4 * k_z[idx]^2 * csch(2 * k_z[idx])^2 +
         2 * k_z[idx] * coth(2 * k_z[idx]) -
         (k_x[idx] + k_z[idx])^2 * csch(k_x[idx] + k_z[idx])^2 -
         (k_x[idx] + k_z[idx]) * coth(k_x[idx] + k_z[idx]))
    vrt_lower[idx] <- (s[idx]^4 / a[idx]^4) *
      (4 * k_z[idx]^2 * csch(2 * k_z[idx])^2 +
         2 * k_z[idx] * coth(2 * k_z[idx]) -
         (-k_x[idx] + k_z[idx])^2 * csch(-k_x[idx] + k_z[idx])^2 -
         (-k_x[idx] + k_z[idx]) * coth(-k_x[idx] + k_z[idx]))
  }

  list(
    pC = pC,
    mdt_upper = mdt_upper,
    mdt_lower = mdt_lower,
    vrt_upper = vrt_upper,
    vrt_lower = vrt_lower
  )
}
