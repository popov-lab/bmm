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

# Cubic Hermite interpolation of log kappa against log J, with the exact
# asymptotics outside the tabulated range. Max relative error 1.7e-11.
.circmix_kappa <- function(J, tab = .circmix_kappa_table()) {
  n <- length(tab$logkappa)
  t <- log(J)
  out <- J + 0.5
  low <- which(t <= tab$logJ_min)
  out[low] <- sqrt(2 * J[low])
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
# weight of each memory component; logw_guess and kappa are length n. cosd is
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
# J ~ gamma(shape = J(kappa) / tau, scale = tau) on a composite Simpson grid in
# log J. tau = 0 is the point mass at J(kappa), so rows with tau = 0 fall back
# to the constant-precision density.
.circmix_vp_ld <- function(cosd, logw, logw_guess, kappa, tau, nodes = 41L,
                           tab = .circmix_kappa_table()) {
  out <- .circmix_ld(cosd, logw, logw_guess, kappa)
  shape <- .circmix_J(kappa) / tau
  half_width <- 8 * sqrt(trigamma(shape))
  vp <- which(tau > 0 & is.finite(half_width) & half_width >= 1e-6)
  if (!length(vp)) {
    return(out)
  }
  stopif(
    any(shape[vp] < 40 / nodes),
    "The variable-precision quadrature holds the log density to about 1e-5 \\
    while the gamma shape J(kappa)/tau stays above {40 / nodes}, but the \\
    smallest value here is {min(shape[vp])}. Raise vp_nodes."
  )

  centre <- digamma(shape[vp]) + log(tau[vp])
  offsets <- seq(-1, 1, length.out = nodes)
  simpson <- log(c(1, rep(c(4, 2), length.out = nodes - 2), 1))
  lp <- vapply(seq_len(nodes), function(i) {
    t <- centre + half_width[vp] * offsets[i]
    J <- exp(t)
    simpson[i] + t +
      stats::dgamma(J, shape = shape[vp], scale = tau[vp], log = TRUE) +
      .circmix_ld(
        cosd[vp, , drop = FALSE], logw[vp, , drop = FALSE], logw_guess[vp],
        .circmix_kappa(J, tab)
      )
  }, numeric(length(vp)))

  step <- 2 * half_width[vp] / (nodes - 1)
  out[vp] <- matrixStats::rowLogSumExps(matrix(lp, nrow = length(vp))) +
    log(step / 3)
  out
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

  high <- log(allocation$extra) + .circmix_vp_ld(
    cosd, logw, no_guessing,
    .circmix_kappa((allocation$slots + 1) * single_slot), tau, nodes
  )

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
  spec$priors$tau <- list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
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
  prep$family$vp_nodes %||% 41L
}

# brms wants bounds on the natural scale of each distributional parameter, which
# the declared link already determines.
.circmix_bounds <- function(links) {
  lower <- c(tan_half = NA, log = 0, softplus = 0, logit = 0, identity = NA)
  upper <- c(tan_half = NA, log = NA, softplus = NA, logit = 1, identity = NA)
  list(
    lb = unname(lower[unlist(links)]),
    ub = unname(upper[unlist(links)])
  )
}

# Called from the exported constructors, where bmm validates arguments.
.circmix_check_variable_precision <- function(variable_precision, vp_nodes) {
  stopif(
    !isTRUE(variable_precision) && !isFALSE(variable_precision),
    "The variable_precision argument must be either TRUE or FALSE."
  )
  stopif(
    vp_nodes < 41 || vp_nodes %% 2 == 0,
    "vp_nodes must be an odd number of at least 41, but is {vp_nodes}."
  )
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
# model without variable precision has no tau to pass.
.circmix_custom_family <- function(model, family, weight_parameters,
                                   vint = FALSE, n_vreal = 0, log_lik,
                                   posterior_predict) {
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
    posterior_predict = posterior_predict
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
                                  vreal = list(), nodes = 41L) {
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
