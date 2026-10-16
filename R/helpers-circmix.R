############################################################################# !
# SHARED NUMERICS FOR THE CIRCULAR MIXTURE MODELS                        ####
############################################################################# !
# R mirror of inst/stan_chunks/circmix_funs.stan, used by the log_lik and
# posterior_predict methods and by the d* functions in distributions.R. The two
# implementations must agree to within floating point;
# tests/testthat/test-helpers-circmix.R checks that against the compiled Stan
# functions rather than against a second R implementation.

# I1(kappa) / I0(kappa). besselI() underflows to zero above kappa ~ 1e5, while
# the asymptotic series is already accurate to 2e-13 at kappa = 1e3.
.circmix_A <- function(kappa) {
  out <- besselI(kappa, 1, expon.scaled = TRUE) / besselI(kappa, 0, expon.scaled = TRUE)
  large <- which(!is.finite(out) | kappa > 1e3)
  if (length(large)) {
    k <- kappa[large]
    out[large] <- 1 - 1 / (2 * k) - 1 / (8 * k^2) - 1 / (8 * k^3)
  }
  out
}

.circmix_log_besselI0 <- function(kappa) {
  out <- log(besselI(kappa, 0, expon.scaled = TRUE)) + kappa
  large <- which(!is.finite(out))
  if (length(large)) {
    k <- kappa[large]
    out[large] <- k - 0.5 * log(2 * pi * k) +
      log1p(1 / (8 * k) + 9 / (128 * k^2) + 225 / (3072 * k^3))
  }
  out
}

# Fisher information of a von Mises about its location
.circmix_J <- function(kappa) {
  kappa * .circmix_A(kappa)
}

.circmix_dlogkappa_dlogJ <- function(kappa) {
  a <- .circmix_A(kappa)
  a / (kappa * (1 - a^2))
}

############################################################################# !
# INVERTING J(kappa)                                                     ####
############################################################################# !

# Newton on log kappa, where f'(u) is the reciprocal of d log kappa / d log J.
# The step is clamped because the asymptotic starting values are poor near the
# crossover at J ~ 2, where an undamped step can overshoot into the range where
# besselI() has already underflowed.
.circmix_invert_J <- function(logJ, tol = 1e-13, maxit = 40L) {
  J <- exp(logJ)
  u <- ifelse(J > 2, log(J + 0.5), 0.5 * log(2 * J))
  for (i in seq_len(maxit)) {
    kappa <- exp(u)
    step <- (log(.circmix_J(kappa)) - logJ) * .circmix_dlogkappa_dlogJ(kappa)
    u <- u - pmax(pmin(step, 1), -1)
    if (max(abs(step)) < tol) break
  }
  u
}

.build_circmix_kappa_table <- function(n = 2001L, logJ_min = log(1e-5),
                                       logJ_max = log(1e5)) {
  logJ <- seq(logJ_min, logJ_max, length.out = n)
  logkappa <- .circmix_invert_J(logJ)
  list(
    logkappa = logkappa,
    dlogkappa = .circmix_dlogkappa_dlogJ(exp(logkappa)),
    logJ_min = logJ_min,
    dlogJ = (logJ_max - logJ_min) / (n - 1)
  )
}

# The table is a pure function of its defaults but costs about a second to
# build, and configure_model() is re-run by stancode(), standata() and
# default_prior() as well as by bmm(), so hold it in the closure.
.circmix_kappa_table <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- .build_circmix_kappa_table()
    }
    cache
  }
})

# Cubic Hermite interpolation of log kappa against log J, with the asymptotics
# outside the tabulated range: below it, inverting
# J = kappa^2 / 2 - kappa^4 / 16 gives sqrt(2 J) (1 + J / 8), whose relative
# error is O(J^2), so the value jumps by about 2e-12 at the edge rather than by
# the 1.25e-6 of sqrt(2 J) alone; above it, J -> kappa - 1/2. Max relative
# error 1.7e-11. Stan's circmix_kappa() receives the same table as data, so the
# R and Stan inverses cannot drift apart.
.circmix_kappa <- function(J, tab = .circmix_kappa_table()) {
  n <- length(tab$logkappa)
  t <- log(J)
  out <- J + 0.5
  low <- which(t <= tab$logJ_min)
  out[low] <- sqrt(2 * J[low]) * (1 + J[low] / 8)
  mid <- which(t > tab$logJ_min & t < tab$logJ_min + (n - 1) * tab$dlogJ)
  if (!length(mid)) {
    return(out)
  }
  pos <- (t[mid] - tab$logJ_min) / tab$dlogJ
  i <- floor(pos) + 1L
  s <- pos - (i - 1)
  s2 <- s * s
  s3 <- s2 * s
  out[mid] <- exp(
    (2 * s3 - 3 * s2 + 1) * tab$logkappa[i] +
      (-2 * s3 + 3 * s2) * tab$logkappa[i + 1] +
      tab$dlogJ * ((s3 - 2 * s2 + s) * tab$dlogkappa[i] +
        (s3 - s2) * tab$dlogkappa[i + 1])
  )
  out
}

############################################################################# !
# MIXTURE DENSITIES                                                      ####
############################################################################# !
# cosd and logw are n x K matrices holding cos(y - mu_k) and the normalised log
# weight of each memory component (exp(logw) and exp(logw_guess) sum to one);
# logw_guess and kappa are length n. cosd is
# supplied by the caller because it does not depend on kappa, which is what lets
# the whole variable-precision grid re-use one set of cosines.

.circmix_ld <- function(cosd, logw, logw_guess, kappa) {
  memory <- matrixStats::rowLogSumExps(logw + kappa * cosd) -
    .circmix_log_besselI0(kappa)
  matrixStats::rowLogSumExps(cbind(memory, logw_guess)) - log(2 * pi)
}

.circmix_het_ld <- function(cosd, logw, logw_guess, kappa) {
  log_i0 <- matrix(.circmix_log_besselI0(as.vector(kappa)), nrow = nrow(cosd))
  memory <- matrixStats::rowLogSumExps(logw + kappa * cosd - log_i0)
  matrixStats::rowLogSumExps(cbind(memory, logw_guess)) - log(2 * pi)
}

# Marginalises the Fisher information of the memory components over
# J ~ gamma(shape = J(kappa) / tau, scale = tau). tau = 0 is the point mass at
# J(kappa), so rows with tau = 0 fall back to the constant-precision density,
# which is why tau is the parameter variable_precision frees. The density is
# linear in the components, sum_k w_k I_k(y) + w_guess / (2 pi), and each I_k
# is integrated on its own grid (.circmix_vp_component()). A grid shared by the
# components and centred on the gamma prior misses the integrand of a response
# far from a component, which peaks many prior SDs below E[log J]; without a
# guess weight to mask it, that underestimated the log density by up to 233
# nats (kappa 500, shape 100, y = pi).
# Below a shape of 60 / nodes the left tail of log J is too heavy for the grid:
# with 81 nodes the log density holds to about 1e-5 from that threshold (0.74)
# up, to 2e-6 from shape 1, and to 1e-9 from shape 10.
.circmix_vp_ld <- function(cosd, logw, logw_guess, kappa, tau, nodes = 81L,
                           tab = .circmix_kappa_table()) {
  out <- .circmix_ld(cosd, logw, logw_guess, kappa)
  shape <- .circmix_J(kappa) / tau
  half_width <- 8 * sqrt(trigamma(shape))
  vp <- which(tau > 0 & is.finite(half_width) & half_width >= 1e-6)
  if (!length(vp)) {
    return(out)
  }
  stopif(
    any(shape[vp] < 60 / nodes),
    "The variable-precision quadrature needs a gamma shape J(kappa)/tau of at \\
    least 60 / vp_nodes = {60 / nodes}, but the smallest value here is \\
    {min(shape[vp])}. Raise vp_nodes."
  )

  log_integral <- vapply(seq_len(ncol(cosd)), function(k) {
    .circmix_vp_component(cosd[vp, k], shape[vp], tau[vp], nodes, tab)
  }, numeric(length(vp)))
  memory <- matrixStats::rowLogSumExps(
    logw[vp, , drop = FALSE] + matrix(log_integral, nrow = length(vp))
  )
  out[vp] <- matrixStats::rowLogSumExps(cbind(memory, logw_guess[vp] - log(2 * pi)))
  out
}

# log of the integral over J ~ gamma(shape, scale = tau) of one von Mises
# component with cos(y - mu) = cos_d, on a composite Simpson grid in t = log J.
# The grid is centred on the mode of the integrand, found by a fixed number of
# damped Newton steps from the large-kappa mode (shape + 1/2) / (1/tau + 1 -
# cos_d); a fixed count keeps the grid a smooth function of the parameters, and
# four steps reach the mode to 1e-4 of the window scale from shape 0.5 to 1e6
# and kappa 0.01 to 1e4 (the start alone is within 0.4 of it there, which the
# window absorbs, so the steps guard parameters outside that range). The grid
# spans -12 to +6 times the integrand's curvature scale, wider below because
# log J is left-skewed, and never less than sd(log J) = sqrt(trigamma(shape)).
# The weights are not renormalised: the grid follows the integrand, not the
# prior, so their sum is not an approximation of the prior's mass.
.circmix_vp_component <- function(cos_d, shape, tau, nodes, tab) {
  t <- log((shape + 0.5) / (1 / tau + 1 - cos_d))
  for (i in 1:4) {
    slopes <- .circmix_vp_slopes(t, cos_d, shape, tau, tab)
    step <- ifelse(slopes$d2 < 0, -slopes$d1 / slopes$d2, ifelse(slopes$d1 > 0, 1, -1))
    t <- t + pmin(pmax(step, -1), 1)
  }
  d2 <- .circmix_vp_slopes(t, cos_d, shape, tau, tab)$d2
  scale <- pmax(1 / sqrt(ifelse(d2 < 0, -d2, Inf)), sqrt(trigamma(shape)))

  n <- length(t)
  grid <- t + outer(scale, seq(-12, 6, length.out = nodes))
  log_integrand <- matrix(.circmix_vp_log_integrand(
    as.vector(grid), rep(cos_d, nodes), rep(shape, nodes), rep(tau, nodes), tab
  ), nrow = n)
  simpson <- log(c(1, rep(c(4, 2), length.out = nodes - 2), 1))
  matrixStats::rowLogSumExps(log_integrand + rep(simpson, each = n)) +
    log(6 * scale / (nodes - 1))
}

.circmix_vp_log_integrand <- function(t, cos_d, shape, tau, tab) {
  J <- exp(t)
  kappa <- .circmix_kappa(J, tab)
  t + stats::dgamma(J, shape = shape, scale = tau, log = TRUE) +
    kappa * cos_d - .circmix_log_besselI0(kappa) - log(2 * pi)
}

# First and second derivatives of .circmix_vp_log_integrand() in t, from
# d kappa / d log J = A / (1 - A^2) and dA / d kappa = 1 - A / kappa - A^2.
.circmix_vp_slopes <- function(t, cos_d, shape, tau, tab) {
  J <- exp(t)
  a <- .circmix_A_terms(.circmix_kappa(J, tab))
  one_minus_A2 <- a$one_minus * (1 + a$A)
  dkappa <- a$A / one_minus_A2
  list(
    d1 = shape - J / tau + (cos_d - a$A) * dkappa,
    d2 = -J / tau - a$slope * dkappa^2 +
      (cos_d - a$A) * a$slope * (1 + a$A^2) * dkappa / one_minus_A2^2
  )
}

# A, 1 - A and dA / d kappa. Above kappa = 1e3, where .circmix_A() switches to
# its series, the last two come from that series too, because computing them
# from A cancels.
.circmix_A_terms <- function(kappa) {
  A <- .circmix_A(kappa)
  one_minus <- 1 - A
  slope <- 1 - A / kappa - A^2
  large <- which(kappa > 1e3)
  if (length(large)) {
    k <- kappa[large]
    one_minus[large] <- 1 / (2 * k) + 1 / (8 * k^2) + 1 / (8 * k^3)
    slope[large] <- 1 / (2 * k^2) + 1 / (4 * k^3) + 3 / (8 * k^4)
  }
  nlist(A, one_minus, slope)
}

# Draws from the mixture by picking a component and then sampling from it, which
# is exact and vectorised over rows, rather than by rejection sampling the
# marginal density. Under variable precision the per-draw concentration is drawn
# from the gamma on J first, so this samples the hierarchy the likelihood
# integrates out and is therefore an independent check on it.
.rcircmix <- function(mu, logw, logw_guess, kappa, tau) {
  n <- nrow(mu)
  n_components <- ncol(mu)
  component <- .circmix_sample_component(exp(cbind(logw, logw_guess)))

  out <- numeric(n)
  guessed <- which(component > n_components)
  out[guessed] <- stats::runif(length(guessed), -pi, pi)

  remembered <- which(component <= n_components)
  if (length(remembered) == 0) {
    return(out)
  }

  drawn_kappa <- kappa[remembered]
  variable <- which(tau[remembered] > 0)
  if (length(variable)) {
    mean_kappa <- kappa[remembered][variable]
    scale <- tau[remembered][variable]
    drawn_kappa[variable] <- .circmix_kappa(stats::rgamma(
      length(variable),
      shape = .circmix_J(mean_kappa) / scale, scale = scale
    ))
  }
  out[remembered] <- brms::rvon_mises(
    length(remembered),
    mu[cbind(remembered, component[remembered])],
    drawn_kappa
  )
  out
}

.circmix_sample_component <- function(weights) {
  cumulative <- matrixStats::rowCumsums(weights)
  total <- cumulative[, ncol(cumulative)]
  rowSums(cumulative < stats::runif(nrow(weights)) * total) + 1L
}

# An item receives floor(K / set_size) slots with probability 1 - extra and one
# more with probability extra, which is continuous in K across the crossings.
.circmix_slots <- function(K, set_size) {
  q <- K / set_size
  slots <- floor(q)
  nlist(slots, extra = q - slots)
}

# Slot averaging (Zhang & Luck, 2008): the reported item holds floor(K /
# set_size) or one more slot, and averaging independent samples adds their
# Fisher information, so a j-slot item has J_j = j * J(kappa) with kappa the
# precision of a single slot. logw holds the item weights given that the item is
# held; an item holding no slot is guessed. Under variable precision the sum of
# j draws from gamma(J(kappa)/tau, tau) is gamma(j J(kappa)/tau, tau), which is
# what passing kappa_j with the same tau gives, so the two mechanisms compose
# without a second quadrature over slot counts.
.circmix_slot_averaging_ld <- function(cosd, logw, K, set_size, kappa, tau, nodes) {
  n <- nrow(cosd)
  allocation <- .circmix_slots(K, set_size)
  no_guessing <- rep(-Inf, n)
  single_slot <- .circmix_J(kappa)

  # At a whole K / set_size the extra slot has weight zero, so it is skipped.
  # Stan cannot skip it: log(0) would give K a NaN gradient, and dropping the
  # branch would give the wrong slope, so it returns the same value in a form
  # whose gradient is the right-hand slope.
  high <- rep(-Inf, n)
  partial <- which(allocation$extra > 0)
  if (length(partial)) {
    high[partial] <- log(allocation$extra[partial]) + .circmix_vp_ld(
      cosd[partial, , drop = FALSE], logw[partial, , drop = FALSE],
      no_guessing[partial],
      .circmix_kappa((allocation$slots[partial] + 1) * single_slot[partial]),
      tau[partial], nodes
    )
  }

  low <- log1p(-allocation$extra) - log(2 * pi)
  held <- which(allocation$slots >= 0.5)
  if (length(held)) {
    low[held] <- log1p(-allocation$extra[held]) + .circmix_vp_ld(
      cosd[held, , drop = FALSE], logw[held, , drop = FALSE], no_guessing[held],
      .circmix_kappa(allocation$slots[held] * single_slot[held]),
      tau[held], nodes
    )
  }
  matrixStats::rowLogSumExps(cbind(low, high))
}

# Draws the number of slots the reported item holds, then samples at that
# precision: the mechanism .circmix_slot_averaging_ld() marginalises, which is
# what makes it an independent check on that density.
.rcircmix_slot_averaging <- function(mu, logw, K, set_size, kappa, tau) {
  n <- nrow(mu)
  allocation <- .circmix_slots(K, set_size)
  slots <- allocation$slots + stats::rbinom(n, 1, allocation$extra)

  out <- stats::runif(n, -pi, pi)
  held <- which(slots >= 0.5)
  if (length(held)) {
    out[held] <- .rcircmix(
      mu[held, , drop = FALSE], logw[held, , drop = FALSE], rep(-Inf, length(held)),
      .circmix_kappa(slots[held] * .circmix_J(kappa[held])), tau[held]
    )
  }
  out
}

############################################################################# !
# MODEL SPECIFICATION                                                    ####
############################################################################# !

.circmix_vp_citation <- glue(
  "van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. \\
  (2012). Variability in encoding precision accounts for visual \\
  short-term memory limitations. Proceedings of the National Academy \\
  of Sciences, 109(22), 8780-8785. https://doi.org/10.1073/pnas.1117465109"
)

# tau means the same in every model that offers variable precision, so it is
# defined once and inserted next to the kappa it modifies.
.circmix_add_variable_precision <- function(spec) {
  after_kappa <- which(names(spec$parameters) == "kappa")
  spec$parameters <- append(spec$parameters, list(tau = glue(
    "Scale of the trial-to-trial variability in memory precision. The Fisher \\
    information J of a memory representation is drawn from \\
    gamma(mean = J(kappa), scale = tau), so tau -> 0 recovers constant precision \\
    and larger values spread precision more widely across trials."
  )), after = after_kappa)
  spec$links <- append(spec$links, list(tau = "log"),
    after = which(names(spec$links) == "kappa")
  )
  spec$priors$tau <- list(
    main = "normal(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(2)"
  )
  spec$init_ranges$tau <- c(0.2, 1)
  spec
}

# The d* functions are called both with a vector of responses and scalar
# parameters, and from log_lik() with one response and a draw per parameter.
.circmix_recycle <- function(...) {
  args <- list(...)
  lapply(args, rep_len, length.out = max(lengths(args)))
}

# cos(y - mu) for the target and the active non-targets. nt is padded to
# max_set_size - 1, so only the first set_size - 1 entries are ever read and a
# trial never pays for components its set size does not have.
.circmix_cos <- function(x, mu, nt, set_size) {
  args <- .circmix_recycle(x = x, mu = mu)
  cosd <- matrix(cos(args$x - args$mu), nrow = length(args$x))
  for (j in seq_len(set_size - 1)) {
    cosd <- cbind(cosd, cos(args$x - nt[j]))
  }
  cosd
}

# tau is only a distributional parameter when the model estimates variable
# precision; otherwise the likelihood is the tau = 0 limit.
.circmix_prep_tau <- function(prep, i) {
  if (is.null(prep$dpars$tau)) {
    return(rep(0, prep$ndraws))
  }
  brms::get_dpar(prep, "tau", i = i)
}

# The node count travels on the family so that log_lik() uses the grid the Stan
# likelihood was fitted with.
.circmix_prep_nodes <- function(prep) {
  prep$family$vp_nodes %||% 81L
}

# brms wants bounds on the natural scale of each distributional parameter, which
# the declared link already determines. The circular location is left unbounded
# because tan_half maps onto the open circle, and an infinite end means no bound.
.circmix_bounds <- function(links) {
  ranges <- vapply(unlist(links), function(link) {
    if (identical(link, "tan_half")) c(NA, NA) else .link_ranges[[link]]
  }, numeric(2))
  ranges[is.infinite(ranges)] <- NA
  list(lb = unname(ranges[1, ]), ub = unname(ranges[2, ]))
}

# Called from the exported constructors, where bmm validates arguments.
.circmix_check_variable_precision <- function(variable_precision, vp_nodes) {
  stopif(
    !isTRUE(variable_precision) && !isFALSE(variable_precision),
    "The variable_precision argument must be either TRUE or FALSE."
  )
  stopif(
    !isTRUE(is.numeric(vp_nodes) && length(vp_nodes) == 1 &&
      vp_nodes >= 41 && vp_nodes %% 2 == 1),
    "vp_nodes must be a single, whole, odd number of at least 41, but is \\
    {collapse_comma(vp_nodes)}."
  )
}

# A posterior that the shape bound 60 / vp_nodes cuts off does not reach the
# reject message: the sampler only diverges at the bound. In a diverging
# mixture3p slot-averaging fit at 81 nodes, 34-50% of the draws sat within 10% of
# it; refitted at 121 nodes, at most 1% sat near the new bound. Draws stopped at
# the bound do not show how far below it the posterior reaches, so the node count
# suggested moves the bound to two thirds of its value instead of to the lowest
# draw.
#' @export
postprocess_brm.circular <- function(model, fit, ...) {
  fit <- NextMethod("postprocess_brm")
  if (isTRUE(model$variable_precision) && has_draws(fit) &&
    identical(fit$algorithm, "sampling")) {
    .circmix_warn_vp_bound(
      .circmix_vp_shape_min(fit),
      divergent = sum(brms::nuts_params(fit, pars = "divergent__")$Value),
      nodes = model$vp_nodes
    )
  }
  fit
}

.circmix_warn_vp_bound <- function(shape_min, divergent, nodes) {
  near <- mean(shape_min < 1.1 * 60 / nodes)
  suggested <- 2 * ceiling(0.75 * (nodes - 1)) + 1
  warnif(
    divergent > 0 && near >= 0.05,
    "{divergent} transitions after warmup were divergent, and in {round(100 * near)}% \\
    of the posterior draws the gamma shape J(kappa)/tau lies within 10% of the \\
    smallest value the variable-precision quadrature accepts, 60 / vp_nodes = \\
    {signif(60 / nodes, 3)}. That bound may cut off part of the posterior. Refit \\
    with vp_nodes = {suggested}, which lowers it to {signif(60 / suggested, 3)} \\
    and takes about 1.5 times as long."
  )
}

# The smallest J(kappa)/tau over the observations, per draw. kappa and tau can
# vary with predictors, so they are read per observation, in blocks that keep the
# draws x observations matrices small.
.circmix_vp_shape_min <- function(fit) {
  prep <- brms::prepare_predictions(fit)
  out <- rep(Inf, prep$ndraws)
  for (idx in split(seq_len(prep$nobs), ceiling(seq_len(prep$nobs) / 1000))) {
    shape <- .circmix_J(brms::get_dpar(prep, "kappa", i = idx)) /
      brms::get_dpar(prep, "tau", i = idx)
    out <- pmin(out, matrixStats::rowMins(matrix(shape, nrow = prep$ndraws)))
  }
  out
}

# The non-target locations reach the likelihood as vreal1..vrealK, padded to
# max_set_size - 1; the likelihood reads only as many as the trial's set size.
.circmix_prep_nt <- function(prep, i) {
  columns <- grep("^vreal[0-9]+$", names(prep$data), value = TRUE)
  columns <- columns[order(as.integer(sub("vreal", "", columns)))]
  vapply(columns, function(column) prep$data[[column]][i], numeric(1))
}

############################################################################# !
# STAN PLUMBING                                                          ####
############################################################################# !

# The distributional parameters are always the location, the concentration,
# optionally tau, and then whatever weights the version uses. core_dpars is
# carried on the family so that the Stan wrapper can pass a literal zero where a
# model without variable precision has no tau to pass. posterior_epred has no
# default because brms's prep path fails without one even for these circular
# models, whose expected response is refused; pass posterior_epred_undefined().
.circmix_custom_family <- function(model, family, weight_parameters,
                                   vint = FALSE, n_vreal = 0, log_lik,
                                   posterior_predict, posterior_epred) {
  variable_precision <- isTRUE(model$variable_precision)
  dpars <- c("mu", "kappa", if (variable_precision) "tau", weight_parameters)
  bounds <- .circmix_bounds(model$links[dpars])

  out <- brms::custom_family(
    family,
    dpars = dpars,
    links = unlist(model$links[dpars], use.names = FALSE),
    lb = bounds$lb, ub = bounds$ub,
    type = "real",
    vars = .circmix_family_vars(vint = vint, n_vreal = n_vreal),
    loop = TRUE,
    log_lik = log_lik,
    posterior_predict = posterior_predict,
    posterior_epred = posterior_epred
  )
  out$vp_nodes <- model$vp_nodes
  out$core_dpars <- c(
    "mu", "kappa", if (variable_precision) "tau" else "0.0", weight_parameters
  )
  out
}

.circmix_model_stanvars <- function(model, family, chunk, vint = NULL,
                                    vreal = list()) {
  .circmix_stanvars() +
    .circmix_chunk_stanvar(chunk) +
    brms::stanvar(
      scode = .circmix_stan_wrapper(
        family$name, family$dpars, family$core_dpars,
        vint = vint, vreal = vreal, nodes = model$vp_nodes
      ),
      block = "functions", name = paste0(family$name, "_lpdf")
    )
}

# Stan functions cannot reach the data block, so the kappa(J) table has to be
# handed to the likelihood through the family's `vars`. pll_args declares it for
# the threaded partial_log_lik, without which the model fails to compile under
# threading.
.circmix_stanvars <- function() {
  tab <- .circmix_kappa_table()
  .circmix_chunk_stanvar("circmix_funs.stan") +
    brms::stanvar(
      x = tab$logkappa, name = "circmix_logk",
      pll_args = "data vector circmix_logk"
    ) +
    brms::stanvar(
      x = tab$dlogkappa, name = "circmix_dlogk",
      pll_args = "data vector circmix_dlogk"
    ) +
    brms::stanvar(
      x = tab$logJ_min, name = "circmix_logJ_min",
      pll_args = "data real circmix_logJ_min"
    ) +
    brms::stanvar(
      x = tab$dlogJ, name = "circmix_dlogJ",
      pll_args = "data real circmix_dlogJ"
    )
}

.circmix_chunk_stanvar <- function(chunk) {
  brms::stanvar(
    scode = read_lines2(system.file("stan_chunks", chunk, package = "bmm")),
    block = "functions", name = sub("\\.stan$", "", chunk)
  )
}

.circmix_table_vars <- function() {
  c("circmix_logk", "circmix_dlogk", "circmix_logJ_min", "circmix_dlogJ")
}

# Only bare `name[n]` entries are rewritten to the global index by brms when
# threading is on. Anything more elaborate is passed through untouched and would
# pair each thread's response slice with the top of the data.
.circmix_family_vars <- function(vint = FALSE, n_vreal = 0) {
  c(
    if (vint) "vint1[n]",
    if (n_vreal > 0) paste0("vreal", seq_len(n_vreal), "[n]"),
    .circmix_table_vars()
  )
}

.circmix_aterm <- function(resp, vint = NULL, vreal = character(0)) {
  terms <- c(
    if (!is.null(vint)) glue("vint({vint})"),
    if (length(vreal)) glue("vreal({paste(vreal, collapse = ', ')})")
  )
  if (!length(terms)) {
    return(glue("{resp} ~ 1"))
  }
  glue("{resp} | {paste(terms, collapse = ' + ')} ~ 1")
}

# brms needs the likelihood to be a single `<family>_lpdf`, whose signature
# depends on max_set_size because every non-target needs its own vreal term, so
# this is the one piece of Stan that is generated. It only packs the scalar
# vreal arguments into the vectors `<family>_core` expects; the arithmetic stays
# in inst/stan_chunks. core_dpars may hold literals, which is how a
# constant-precision model passes 0 for the tau its family does not estimate.
.circmix_stan_wrapper <- function(family, dpars, core_dpars = dpars, vint = NULL,
                                  vreal = list(), nodes = 81L) {
  vreal_names <- lapply(names(vreal), function(g) paste0(g, seq_len(vreal[[g]])))
  vreal_args <- unlist(vreal_names)
  args <- c(
    "real y", paste("real", dpars),
    if (!is.null(vint)) paste("int", vint),
    if (length(vreal_args)) paste("real", vreal_args),
    "data vector circmix_logk", "data vector circmix_dlogk",
    "data real circmix_logJ_min", "data real circmix_dlogJ"
  )
  packed <- vapply(
    vreal_names,
    function(x) glue("to_vector({{{paste(x, collapse = ', ')}}})"),
    character(1)
  )
  call_args <- c(
    "y", core_dpars, vint, packed, as.character(as.integer(nodes)),
    .circmix_table_vars()
  )
  as.character(glue(
    "  real {family}_lpdf({paste(args, collapse = ', ')}) {{\n",
    "    return {family}_core({paste(call_args, collapse = ', ')});\n",
    "  }}\n"
  ))
}

############################################################################# !
# SET-SIZE-1 WARNINGS                                                    ####
############################################################################# !

# One wording for every model whose set-size-1 trials leave coefficients to
# their priors; only the reason the trials give differs between parameters
.warn_set_size_one <- function(coefs, pars, reason) {
  formulas <- if (length(pars) > 1) "formulas" else "formula"
  pars <- paste0("'", pars, "'", collapse = " and ")
  warnif(
    length(coefs) > 0,
    "Trials with a set size of 1 have no non-target, so they {reason}. The \\
    other trials do not determine {collapse_comma(coefs)} in the {formulas} for \\
    {pars}, so the posteriors of these coefficients depend on their priors. Give \\
    set size 1 a level of its own in a factor, so that only the set-size-1 \\
    coefficients depend on their priors."
  )
}

# The coefficients of the formulas for `pars` that the trials with more than one
# item leave undetermined, named as in the fit's draws. The trials with one item
# inform nothing about `pars` or, with `summed = TRUE`, only the sum of their
# effects. Group-level terms are judged by their columns pooled over groups,
# which is what their sd is estimated from, one grouping factor at a time.
.set_size_one_undetermined <- function(pars, formula, data, summed = FALSE) {
  forms <- lapply(pars, function(par) formula[[par]])
  if (any(vapply(forms, function(f) !is.null(f) && is_nl(f), logical(1)))) {
    return(character(0))
  }
  estimated <- vapply(forms, function(f) !is.null(f) && !is_constant(f), logical(1))
  pars <- pars[estimated]
  forms <- forms[estimated]
  # a variable missing from the data is left to brms, which names it
  vars <- unique(unlist(lapply(forms, function(f) all.vars(f[[3]]))))
  if (!length(pars) || !all(vars %in% names(data))) {
    return(character(0))
  }
  data <- data[stats::complete.cases(data[vars]), ]
  informed <- data[["ss_numeric"]] > 1
  bterms <- lapply(forms, function(f) brms::brmsterms(brms::bf(f))$dpars$mu)

  fixed <- lapply(bterms, function(bt) {
    get_model_matrix(bt$fe, data) %||% matrix(0, nrow(data), 0)
  })
  out <- .lost_coefficients(fixed, informed, summed, glue("b_{pars}_"))

  re <- lapply(bterms, function(bt) if (is.data.frame(bt$re)) bt$re)
  for (group in unique(unlist(lapply(re, function(r) r$group)))) {
    pooled <- lapply(re, function(r) {
      terms <- which(r$group == group)
      if (!length(terms)) {
        return(matrix(0, nrow(data), 0))
      }
      do.call(cbind, lapply(terms, function(i) get_model_matrix(r$form[[i]], data)))
    })
    out <- c(out, .lost_coefficients(pooled, informed, summed, glue("sd_{group}__{pars}_")))
  }
  as.character(out)
}

# Columns of the design matrices in `blocks`, one per parameter, involved in a
# direction that the informed rows (plus, when summed, the uninformed rows as
# one row across all blocks) leave free but every row informing every block
# separately would fix
.lost_coefficients <- function(blocks, informed, summed, prefixes) {
  blocks <- lapply(blocks, function(X) {
    norms <- sqrt(colSums(X^2))
    sweep(X, 2, ifelse(norms > 0, norms, 1), "/")
  })
  widths <- vapply(blocks, ncol, integer(1))
  if (!sum(widths)) {
    return(character(0))
  }
  full <- .block_diagonal(blocks)
  available <- rbind(
    .block_diagonal(lapply(blocks, function(X) X[informed, , drop = FALSE])),
    if (summed) do.call(cbind, lapply(blocks, function(X) X[!informed, , drop = FALSE]))
  )
  free <- .null_space(full)
  free_in_available <- .null_space(available)
  lost <- free_in_available - free %*% crossprod(free, free_in_available)
  names <- unlist(lapply(seq_along(blocks), function(i) {
    if (widths[i]) paste0(prefixes[i], colnames(blocks[[i]]))
  }))
  names[rowSums(abs(lost)) > 1e-8]
}

.block_diagonal <- function(blocks) {
  widths <- vapply(blocks, ncol, integer(1))
  heights <- vapply(blocks, nrow, integer(1))
  out <- matrix(0, sum(heights), sum(widths))
  rows <- cumsum(c(0, heights))
  columns <- cumsum(c(0, widths))
  for (i in seq_along(blocks)) {
    if (widths[i] && heights[i]) {
      out[rows[i] + seq_len(heights[i]), columns[i] + seq_len(widths[i])] <- blocks[[i]]
    }
  }
  out
}

.null_space <- function(M) {
  if (!nrow(M)) {
    return(diag(ncol(M)))
  }
  s <- svd(M, nu = 0, nv = ncol(M))
  d <- c(s$d, rep(0, ncol(M) - length(s$d)))
  s$v[, d <= 1e-8 * max(d), drop = FALSE]
}
