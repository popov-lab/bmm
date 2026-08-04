# shared constants and helpers for change-detection models
.cd_sharpness <- 5

# Number of Gauss-Legendre nodes used to integrate the retrieval distribution
# over the "same" arc. 32 nodes match a 2e6-point reference rule to ~1e-10.
.cd_n_nodes <- 32

# Gauss-Legendre nodes and weights on [-1, 1] via Golub-Welsch: the nodes are
# the eigenvalues of the Jacobi matrix of the Legendre recurrence.
.cd_gauss_legendre <- function(n = .cd_n_nodes) {
  k <- seq_len(n - 1)
  b <- k / sqrt(4 * k^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(k, k + 1)] <- b
  jacobi[cbind(k + 1, k)] <- b
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  nlist(nodes = eig$values[ord], weights = 2 * eig$vectors[1, ord]^2)
}

.cd_log_i0 <- function(kappa) {
  log(besselI(kappa, nu = 0, expon.scaled = TRUE)) + kappa
}

# Half-width of the arc around the probe within which the Bayesian decision rule
# of Lin & Oberauer (2022) yields a "same" response. Their Eq. 8 sets
# LLR = log[(1/2pi) / (p_mem * vonMises(x | probe, kappa) + (1 - p_mem)/2pi)],
# and the response is "change" when LLR > criterion. Because the von Mises
# density decreases monotonically in |x - probe|, that region is always an arc.
# Their Appendix B shows p_mem drops out of the boundary when criterion = 0.
.cd_crit_angle <- function(kappa, criterion = 0, p_mem = 1) {
  offset <- if (all(criterion == 0)) {
    0
  } else {
    log(pmax((exp(-criterion) - 1 + p_mem) / p_mem, 0))
  }
  acos(pmin(pmax((.cd_log_i0(kappa) + offset) / kappa, -1), 1))
}

# Mass of a von Mises(mu, kappa) inside the arc centred on `centre` with
# half-width `half_width`, by Gauss-Legendre quadrature. All arguments are
# recycled to a common length n; `centre` and `half_width` are length n.
.cd_vm_arc_mass <- function(centre, half_width, kappa, mu = 0, gl = .cd_gauss_legendre()) {
  x <- outer(half_width, gl$nodes) + centre - mu
  dens <- exp(kappa * cos(x) - .cd_log_i0(kappa)) / (2 * pi)
  half_width * drop(dens %*% gl$weights)
}

.log1mexp <- function(x) {
  stopif(
    any(x > 0, na.rm = TRUE),
    "Argument 'x' must be less than or equal to 0 in .log1mexp()."
  )

  ifelse(x < log(2), log1p(-exp(x)), log(-expm1(x)))
}

.safe_log <- function(x) {
  out <- rep(-Inf, length(x))
  valid <- !is.na(x) & x > 0
  out[valid] <- log(x[valid])
  out
}

.vm_log_norm <- function(kappa) {
  log(2 * pi) + log(besselI(kappa, nu = 0, expon.scaled = TRUE)) + kappa
}

.cd_quadrature <- function(n_quad = 101) {
  x_grid <- seq(-pi, pi, length.out = n_quad)
  nlist(x_grid, dx = x_grid[2] - x_grid[1])
}

.recycle_cd_args <- function(...) {
  args <- list(...)
  if (length(args) == 0) {
    return(list(n = 0))
  }

  lens <- vapply(args, length, integer(1))
  n <- max(lens)
  stopif(
    any(!(lens %in% c(1L, n))),
    "All observation-level arguments must either have length 1 or the same length."
  )

  recycled <- lapply(args, rep_len, length.out = n)
  recycled$n <- n
  recycled
}

.as_cd_matrix <- function(x, n, name, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(matrix(numeric(0), nrow = n, ncol = 0))
    }
    stop2("Argument '{name}' is required.")
  }

  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }

  if (is.matrix(x)) {
    x <- unname(x)
    stopif(
      nrow(x) != n && nrow(x) != 1,
      "Argument '{name}' must have either 1 row or {n} rows."
    )
    if (nrow(x) == 1 && n > 1) {
      x <- x[rep(1, n), , drop = FALSE]
    }
    storage.mode(x) <- "double"
    return(x)
  }

  stopif(
    !is.atomic(x),
    "Argument '{name}' must be an atomic vector, matrix, or data.frame."
  )

  x <- as.numeric(x)
  stopif(length(x) == 0, "Argument '{name}' cannot be empty.")
  matrix(rep(x, each = n), nrow = n, byrow = FALSE)
}

.prepare_cd_nt_inputs <- function(nt_features, lure_idx = NULL, nt_distances = NULL, n) {
  nt_features <- .as_cd_matrix(nt_features, n, "nt_features")
  lure_idx <- if (is.null(lure_idx)) {
    matrix(1, nrow = nrow(nt_features), ncol = ncol(nt_features))
  } else {
    .as_cd_matrix(lure_idx, n, "lure_idx")
  }

  stopif(
    !all(dim(lure_idx) == dim(nt_features)),
    "Arguments 'nt_features' and 'lure_idx' must have matching dimensions."
  )

  lure_idx <- ifelse(lure_idx > 0, 1, 0)
  storage.mode(lure_idx) <- "double"

  out <- nlist(nt_features, lure_idx)

  if (!is.null(nt_distances)) {
    nt_distances <- .as_cd_matrix(nt_distances, n, "nt_distances")
    stopif(
      !all(dim(nt_distances) == dim(nt_features)),
      "Arguments 'nt_features' and 'nt_distances' must have matching dimensions."
    )
    out$nt_distances <- nt_distances
  }

  out
}

.resolve_mixture3p_cd_probs <- function(p_target = NULL, p_nontarget = NULL,
                                        thetat = NULL, thetant = NULL) {
  if (!is.null(p_target) || !is.null(p_nontarget)) {
    stopif(
      !is.null(thetat) || !is.null(thetant),
      "Please specify either probabilities ('p_target', 'p_nontarget') or legacy weights ('thetat', 'thetant'), not both."
    )
    p_target <- p_target %||% 0.6
    p_nontarget <- p_nontarget %||% 0.3
  } else if (!is.null(thetat) || !is.null(thetant)) {
    stopif(
      is.null(thetat) || is.null(thetant),
      "Both 'thetat' and 'thetant' must be provided together."
    )
    weights <- .recycle_cd_args(thetat = thetat, thetant = thetant)
    p_target <- numeric(weights$n)
    p_nontarget <- numeric(weights$n)
    for (i in seq_len(weights$n)) {
      probs <- softmax(c(weights$thetat[i], weights$thetant[i], 0))
      p_target[i] <- probs[1]
      p_nontarget[i] <- probs[2]
    }
  } else {
    p_target <- 0.6
    p_nontarget <- 0.3
  }

  stopif(isTRUE(any(p_target < 0 | p_target > 1, na.rm = TRUE)), "p_target must be in [0, 1].")
  stopif(
    isTRUE(any(p_nontarget < 0 | p_nontarget > 1, na.rm = TRUE)),
    "p_nontarget must be in [0, 1]."
  )
  stopif(
    isTRUE(any(p_target + p_nontarget > 1, na.rm = TRUE)),
    "The sum of p_target and p_nontarget must be less than or equal to 1."
  )

  nlist(p_target, p_nontarget)
}

.mean_cd_chain_time <- function(fit) {
  chain_times <- try(fit_info(fit, "time_mean"), silent = TRUE)
  if (is_try_error(chain_times)) {
    chain_times <- try(fit_info(fit, "time"), silent = TRUE)
  }
  if (is_try_error(chain_times)) {
    return(NA_real_)
  }

  times_num <- suppressWarnings(as.numeric(unlist(chain_times)))
  times_num <- times_num[is.finite(times_num)]
  if (length(times_num) == 0) {
    return(NA_real_)
  }

  mean(times_num)
}

#' Extract convergence diagnostics for change-detection fits
#'
#' @param fit A fitted `bmm`/`brmsfit` object.
#'
#' @return A one-row data.frame with backend, divergence, tree-depth, and
#'   convergence summaries.
#'
#' @keywords extract_info
#' @export
cd_fit_diagnostics <- function(fit) {
  stopif(
    !inherits(fit, "brmsfit"),
    "Argument 'fit' must inherit from 'brmsfit'."
  )

  backend <- class(fit$fit)[1]

  if (methods::is(fit$fit, "stanfit")) {
    sampler_params <- rstan::get_sampler_params(fit$fit, inc_warmup = FALSE)
    divergences <- sum(vapply(
      sampler_params,
      function(x) sum(x[, "divergent__"]),
      numeric(1)
    ))
    max_treedepth <- max(vapply(
      sampler_params,
      function(x) max(x[, "treedepth__"]),
      numeric(1)
    ))
  } else if (inherits(fit$fit, "CmdStanMCMC")) {
    sampler_diag <- fit$fit$sampler_diagnostics()
    divergences <- sum(sampler_diag[, , "divergent__"])
    max_treedepth <- max(sampler_diag[, , "treedepth__"])
  } else {
    stop2("Unsupported backend class: {backend}")
  }

  rhats <- suppressWarnings(brms::rhat(fit))
  neff <- suppressWarnings(brms::neff_ratio(fit))
  mean_chain_time <- .mean_cd_chain_time(fit)

  data.frame(
    backend = backend,
    divergences = divergences,
    max_treedepth = max_treedepth,
    max_rhat = suppressWarnings(max(rhats, na.rm = TRUE)),
    min_neff_ratio = suppressWarnings(min(neff, na.rm = TRUE)),
    mean_chain_time = mean_chain_time
  )
}
