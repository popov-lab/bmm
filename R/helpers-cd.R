############################################################################# !
# CHANGE DETECTION ON THE CIRCULAR MIXTURES                              ####
############################################################################# !
# R twins of inst/stan_chunks/circmix_cd_funs.stan and cd_funs.stan, used by
# the log_lik and posterior_predict methods and by the d*_cd functions of the
# change-detection models. Every function takes one row per observation (or
# per posterior draw), as the circmix densities in helpers-circmix.R do.
#
# The decision rule is Lin and Oberauer's (2022, Eq. 8): the observer retrieves
# a feature x from the retrieval mixture of the continuous-reproduction model
# and answers "change" when
#   log[(1 / 2pi) / (p_s vM(x | probe, kappa) + (1 - p_s) / 2pi)] > criterion.
# The ratio grows with the distance between x and the probe, so the "same"
# region is an arc around the probe whose half-width hw solves
#   vM(hw | 0, kappa) = (e^-criterion - 1 + p_s) / (2 pi p_s)
# (their Appendix B), i.e. cos(hw) = (log I0(kappa) + offset) / kappa with
# offset = log((e^-criterion - 1 + p_s) / p_s). At criterion = 0 the offset is
# 0 whatever p_s, which is why p_s drops out of the unbiased model. P("same") is
# the retrieval mass inside the arc: each von Mises component contributes its
# arc integral, guessing contributes hw / pi. Non-targets enter only through
# the retrieval mixture (their Eq. B.20), which is what produces the intrusion
# cost. p_s is passed explicitly rather than read from the target weight,
# because under slot averaging the observer's prior of storage is not the
# target's retrieval weight.
#
# When precision varies (variable precision, slot averaging) the observer may
# know the precision of the current trial ("rich": one boundary per precision
# state) or only its distribution ("limited": one boundary from the marginal
# target density over the states, their Eq. B.11 and footnote 1). Lin and
# Oberauer's data favoured the limited observer for both mechanisms (their
# Table 2), so the models default to it.

# Gauss-Legendre rule on [-1, 1] by Golub-Welsch. 32 nodes integrate a von
# Mises over the arc to about 1e-10 against a 2e6-point reference rule
# (measured on the August change-detection branch). The rule reaches Stan as
# data, so the R and Stan arcs use identical nodes.
.cd_gl_rule <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      n <- 32L
      k <- seq_len(n - 1)
      b <- k / sqrt(4 * k^2 - 1)
      jacobi <- matrix(0, n, n)
      jacobi[cbind(k, k + 1)] <- b
      jacobi[cbind(k + 1, k)] <- b
      eig <- eigen(jacobi, symmetric = TRUE)
      ord <- order(eig$values)
      cache <<- list(x = eig$values[ord], w = 2 * eig$vectors[1, ord]^2)
    }
    cache
  }
})

# +Inf means never "same" and -Inf always "same". With p_s = 0 the likelihood
# ratio is identically 1, so the sign of the criterion alone decides; with
# e^-criterion - 1 + p_s <= 0 the decision density exceeds the threshold
# everywhere.
.circmix_cd_offset <- function(criterion, p_s) {
  args <- .circmix_recycle(criterion = criterion, p_s = p_s)
  criterion <- args$criterion
  p_s <- args$p_s
  scaled <- expm1(-criterion) + p_s
  out <- log(pmax(scaled, 0)) - log(p_s)
  out[scaled <= 0] <- -Inf
  out[p_s <= 0] <- ifelse(criterion[p_s <= 0] < 0, Inf, -Inf)
  out[criterion == 0] <- 0
  out
}

.circmix_cd_crit_angle <- function(kappa, offset) {
  acos(pmin(pmax((.circmix_log_besselI0(kappa) + offset) / kappa, -1), 1))
}

# Mass of each von Mises component inside the arc; d is n x K (probe minus each
# component location), hw and kappa are length n. A single 32-node rule holds
# the mass to about 1e-14 while the arc is narrower than 6 / sqrt(kappa), which
# covers every unbiased boundary (hw sqrt(kappa) is 1.3 at kappa = 1 and 3.2 at
# kappa = 5000) but not a wide arc around a sharp peak: at hw sqrt(kappa) = 20
# the plain rule errs by 1e-2. Wider arcs are split at the peak and at
# +-2, 6, 18 / sqrt(kappa) around it (and its images 2 pi away), which held the
# error below 1.3e-12 for kappa up to 2e4 against adaptive integration.
.circmix_cd_arc_mass <- function(d, hw, kappa, gl = .cd_gl_rule()) {
  hw <- rep_len(hw, nrow(d))
  kappa <- rep_len(kappa, nrow(d))
  log_norm <- .circmix_log_besselI0(kappa) + log(2 * pi)
  centre <- atan2(sin(d), cos(d))
  mass <- .circmix_cd_gl_vm(centre - hw, centre + hw, kappa, log_norm, gl)
  wide <- which(hw > 6 / sqrt(kappa))
  multiples <- c(-18, -6, -2, 0, 2, 6, 18)
  for (i in wide) {
    breaks <- as.vector(outer(multiples / sqrt(kappa[i]), 2 * pi * (-1:1), "+"))
    for (k in seq_len(ncol(d))) {
      a <- centre[i, k] - hw[i]
      b <- centre[i, k] + hw[i]
      cuts <- sort(c(a, breaks[breaks > a & breaks < b], b))
      n <- length(cuts)
      mass[i, k] <- sum(.circmix_cd_gl_vm(
        cuts[-n], cuts[-1], kappa[i], log_norm[i], gl
      ))
    }
  }
  mass
}

# Gauss-Legendre mass of vM(0, kappa) on [a, b]; a and b may be matrices with
# one row per kappa
.circmix_cd_gl_vm <- function(a, b, kappa, log_norm, gl = .cd_gl_rule()) {
  half <- (b - a) / 2
  mid <- (a + b) / 2
  total <- 0
  for (i in seq_along(gl$x)) {
    total <- total + gl$w[i] * exp(kappa * cos(mid + half * gl$x[i]) - log_norm)
  }
  half * total
}

.circmix_cd_psame_at <- function(d, logw, logw_guess, kappa, hw) {
  hw <- rep_len(hw, nrow(d))
  out <- exp(logw_guess) * hw / pi +
    rowSums(exp(logw) * .circmix_cd_arc_mass(d, hw, kappa))
  out[hw <= 0] <- 0
  out[hw >= pi] <- 1
  out
}

.circmix_cd_psame <- function(d, logw, logw_guess, kappa, criterion, p_s) {
  hw <- .circmix_cd_crit_angle(kappa, .circmix_cd_offset(criterion, p_s))
  .circmix_cd_psame_at(d, logw, logw_guess, kappa, hw)
}

# Knowledge-limited boundary: the root on (0, pi) of
#   f(h) = log sum_n exp(lw_n) exp(kappa_n cos h) / I0(kappa_n) - offset,
# the log of 2 pi times the marginal target density over the precision states.
# f decreases strictly on (0, pi), with f'(h) = -sin(h) sum_n softmax_n kappa_n,
# so Newton converges from the closed form at the mean kappa; a step leaving the
# bracket is replaced by bisection. The last step is a Newton step, which in
# Stan carries the implicit-function gradient. kappa and lw are n x S.
.circmix_cd_limited_hw <- function(kappa, lw, offset, start) {
  n <- nrow(kappa)
  base <- lw - .circmix_log_besselI0(kappa)
  out <- rep(NA_real_, n)
  out[is.infinite(offset)] <- ifelse(offset[is.infinite(offset)] > 0, 0, pi)
  todo <- which(is.na(out))
  out[todo][matrixStats::rowLogSumExps(base[todo, , drop = FALSE] +
    kappa[todo, , drop = FALSE]) <= offset[todo]] <- 0
  todo <- which(is.na(out))
  out[todo][matrixStats::rowLogSumExps(base[todo, , drop = FALSE] -
    kappa[todo, , drop = FALSE]) >= offset[todo]] <- pi
  todo <- which(is.na(out))

  lo <- rep(0, length(todo))
  hi <- rep(pi, length(todo))
  h <- start[todo]
  h[!(h > lo & h < hi)] <- pi / 2
  for (iter in seq_len(60)) {
    if (!length(todo)) break
    b <- base[todo, , drop = FALSE]
    k <- kappa[todo, , drop = FALSE]
    lp <- b + k * cos(h)
    lse <- matrixStats::rowLogSumExps(lp)
    f <- lse - offset[todo]
    lo[f > 0] <- h[f > 0]
    hi[f <= 0] <- h[f <= 0]
    h_new <- h + f / (sin(h) * rowSums(exp(lp - lse) * k))
    outside <- !(h_new > lo & h_new < hi)
    h_new[outside] <- (lo[outside] + hi[outside]) / 2
    done <- abs(h_new - h) < 1e-12
    out[todo[done]] <- h_new[done]
    keep <- !done
    todo <- todo[keep]
    h <- h_new[keep]
    lo <- lo[keep]
    hi <- hi[keep]
  }
  out[todo] <- h
  out
}

# P("same") over precision states: kappa and lw are n x S, lw holding the
# normalised log weight of each state (-Inf for padding). A row with a single
# state takes the closed-form boundary under either knowledge, because the two
# observers coincide when there is nothing to know.
.circmix_cd_mix_psame <- function(d, logw, logw_guess, kappa, lw, criterion,
                                  p_s, rich) {
  n <- nrow(kappa)
  offset <- .circmix_cd_offset(rep_len(criterion, n), rep_len(p_s, n))
  hw <- .circmix_cd_crit_angle(kappa, offset)
  shared <- which(!rich & rowSums(is.finite(lw)) > 1)
  if (length(shared)) {
    start <- .circmix_cd_crit_angle(
      exp(matrixStats::rowLogSumExps(lw[shared, , drop = FALSE] +
        log(kappa[shared, , drop = FALSE]))),
      offset[shared]
    )
    hw[shared, ] <- .circmix_cd_limited_hw(
      kappa[shared, , drop = FALSE], lw[shared, , drop = FALSE],
      offset[shared], start
    )
  }
  hw <- matrix(hw, nrow = n)
  ps <- vapply(seq_len(ncol(kappa)), function(s) {
    .circmix_cd_psame_at(d, logw, logw_guess, kappa[, s], hw[, s])
  }, numeric(n))
  rowSums(exp(lw) * matrix(ps, nrow = n))
}

# Precision states [kappa, normalised log weight] of a memory with mean kappa:
# the J grid of .circmix_vp_ld() under variable precision, else kappa alone.
# Rows with constant precision are padded to the grid width with weight zero,
# so all rows share one matrix.
.circmix_cd_vp_states <- function(kappa, tau, nodes, tab = .circmix_kappa_table()) {
  n <- length(kappa)
  vp <- .circmix_vp_rows(kappa, tau, nodes)
  if (!length(vp$rows)) {
    return(list(kappa = matrix(kappa), lw = matrix(0, nrow = n)))
  }
  states <- list(
    kappa = matrix(kappa, nrow = n, ncol = nodes),
    lw = cbind(0, matrix(-Inf, nrow = n, ncol = nodes - 1))
  )
  grid <- .circmix_vp_grid(vp$shape, tau[vp$rows], nodes)
  states$kappa[vp$rows, ] <- .circmix_kappa(exp(grid$t), tab)
  states$lw[vp$rows, ] <- grid$lw - matrixStats::rowLogSumExps(grid$lw)
  states
}

.circmix_cd_vp_psame <- function(d, logw, logw_guess, kappa, tau, criterion,
                                 p_s, rich, nodes = 41L) {
  states <- .circmix_cd_vp_states(kappa, tau, nodes)
  .circmix_cd_mix_psame(
    d, logw, logw_guess, states$kappa, states$lw, criterion, p_s, rich
  )
}

# branch from .circmix_slot_averaging_branches(). With a slot held, the
# observer knows the item is stored (p_s = 1) and the precision states are the
# two slot counts, each with its own J grid. With none held, retrieval is a
# guess with probability 1 - extra and the item at kappa_hi otherwise; the
# observer's prior of storage is extra, and both cases are judged with the
# precision the item would have at kappa_hi.
.circmix_cd_slot_averaging_psame <- function(d, logw, branch, tau, criterion,
                                             rich, nodes = 41L) {
  n <- nrow(d)
  criterion <- rep_len(criterion, n)
  out <- numeric(n)

  empty <- which(!branch$held)
  if (length(empty)) {
    out[empty] <- .circmix_cd_vp_psame(
      d[empty, , drop = FALSE],
      branch$log_w_hi[empty] + logw[empty, , drop = FALSE],
      branch$log_w_lo[empty], branch$kappa_hi[empty], tau[empty],
      criterion[empty], exp(branch$log_w_hi[empty]), rich, nodes
    )
  }

  held <- which(branch$held)
  if (length(held)) {
    lo <- .circmix_cd_vp_states(branch$kappa_lo[held], tau[held], nodes)
    hi <- .circmix_cd_vp_states(branch$kappa_hi[held], tau[held], nodes)
    out[held] <- .circmix_cd_mix_psame(
      d[held, , drop = FALSE], logw[held, , drop = FALSE], rep(-Inf, length(held)),
      cbind(lo$kappa, hi$kappa),
      cbind(branch$log_w_lo[held] + lo$lw, branch$log_w_hi[held] + hi$lw),
      criterion[held], 1, rich
    )
  }
  out
}

# y = 1 is "change"; the clamp matches machine_precision() in Stan
.cd_bernoulli_ld <- function(y, p_same) {
  p <- pmin(pmax(p_same, .Machine$double.eps), 1 - .Machine$double.eps)
  ifelse(y == 1, log1p(-p), log(p))
}

############################################################################# !
# MODEL SPECIFICATION AND STAN PLUMBING                                  ####
############################################################################# !

# criterion means the same in every change-detection model. It is fixed to 0 by
# default, which is the unbiased observer and makes p_s cancel; a constant(0)
# intercept reaches Stan as exactly 0, so circmix_cd_offset() tests
# criterion == 0 rather than taking a flag. Placed after the weight parameters,
# because it is the last distributional parameter of the family.
.cd_add_criterion <- function(spec) {
  spec$parameters$criterion <- glue(
    "Decision criterion. A 'change' response is given when the log-likelihood \\
    ratio of a change against no change exceeds the criterion, so larger \\
    values make 'change' responses less likely. Fixed to 0 (an unbiased \\
    observer) by default; include a formula for criterion to estimate it."
  )
  spec$links$criterion <- "identity"
  spec$fixed_parameters$criterion <- 0
  spec$priors$criterion <- list(
    main = "normal(0, 0.5)", effects = "normal(0, 0.5)", sd = "exponential(2)"
  )
  spec$init_ranges$criterion <- c(-0.2, 0.2)
  spec
}

.cd_gl_data <- function() {
  c(cd_gl_x = "vector", cd_gl_w = "vector")
}

# The Gauss-Legendre rule reaches the likelihood through the family's vars, like
# the kappa(J) table, and needs pll_args for the threaded partial_log_lik.
.cd_stanvars <- function() {
  gl <- .cd_gl_rule()
  .circmix_chunk_stanvar("cd_funs.stan") +
    .circmix_chunk_stanvar("circmix_cd_funs.stan") +
    brms::stanvar(x = gl$x, name = "cd_gl_x", pll_args = "data vector cd_gl_x") +
    brms::stanvar(x = gl$w, name = "cd_gl_w", pll_args = "data vector cd_gl_w")
}

# The family of a change-detection model: the circmix family with criterion as
# its last distributional parameter, an integer response, and the
# Gauss-Legendre rule as data. knowledge travels on the family so that
# log_lik() uses the rule the model was fitted with.
.cd_custom_family <- function(model, family, weight_parameters, vint = FALSE,
                              n_vreal = 1, log_lik, posterior_predict,
                              posterior_epred) {
  out <- .circmix_custom_family(
    model, family, weight_parameters,
    vint = vint, n_vreal = n_vreal, log_lik = log_lik,
    posterior_predict = posterior_predict, posterior_epred = posterior_epred,
    extra_dpars = "criterion", type = "int",
    data = c(.circmix_table_data(), .cd_gl_data())
  )
  out$knowledge <- model$knowledge
  out
}

# vreal must list the probe first, as a group of one, because bmf2bf() puts
# probe_centered first in vreal() and the wrapper passes it as a bare real.
.cd_model_stanvars <- function(model, family, chunks, vint = NULL,
                               vreal = list(probe = 1)) {
  .circmix_model_stanvars(
    model, family, chunks,
    vint = vint, vreal = vreal, vreal_scalars = "probe",
    literals = as.character(as.integer(identical(model$knowledge, "rich"))),
    data = c(.circmix_table_data(), .cd_gl_data()),
    extra_stanvars = .cd_stanvars()
  )
}

.cd_prep_probe <- function(prep, i) {
  prep$data$vreal1[i]
}

.cd_prep_rich <- function(prep) {
  identical(prep$family$knowledge, "rich")
}

# A logical response is accepted and coerced, as for the diffusion models.
.cd_check_binary_response <- function(data, resp) {
  values <- data[[resp]]
  stopif(anyNA(values), "The response variable '{resp}' contains missing values.")
  if (is.logical(values)) {
    warning2(
      "The response variable '{resp}' is logical and will be coded as 1 for \\
      TRUE ('change') and 0 for FALSE ('same')."
    )
    data[[resp]] <- as.integer(values)
    return(data)
  }
  stopif(
    !is.numeric(values) || !all(values %in% c(0, 1)),
    "The response variable '{resp}' must be coded as 0 ('same') and 1 \\
    ('change'), or be logical."
  )
  data
}
