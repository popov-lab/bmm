#' Population summaries of the model parameters, per posterior draw
#'
#' @description
#' Integrates each parameter over the distribution of the group-level effects
#' that `re_formula` selects, separately for every posterior draw and grid cell.
#' A parameter that enters `brms` through a linear formula, with plain
#' group-level terms, and whose native value is the elementwise inverse link of
#' its linear predictor is integrated over a normal distribution with standard
#' deviation `sqrt(z' Sigma z)`: in closed form or by Gauss-Hermite quadrature.
#' Every other parameter is integrated by Monte Carlo over simulated new levels
#' of the grouping factors.
#'
#' @param x A bmmfit object
#' @param model_pars Character vector of all parameters of the model
#' @param pars Character vector of the parameters the user asked for
#' @param newdata The prediction grid, one row per cell
#' @param grid_vars Variables spanning the grid
#' @param re_formula `NULL`, `NA` or a formula, as in [native_parameters()]
#' @param draw_ids Sorted indices of the posterior draws to use
#' @param population_summary `"mean"` or `"median"`
#' @param scale `"native"` or `"sampling"`
#' @param ndraws_population Number of simulated levels per draw on the Monte
#'   Carlo path
#' @param dots Further arguments for [brms::posterior_linpred()]
#' @return A list with `values`, a list of draws x cell matrices; `parameter`
#'   and `statistic`, one element per matrix; and `method`, a data frame with
#'   one row per parameter and statistic describing how it was computed.
#'
#' @keywords internal
#' @noRd
.np_population <- function(x, model_pars, pars, newdata, grid_vars, re_formula,
                           draw_ids, population_summary, scale,
                           ndraws_population, dots) {
  native <- scale == "native"
  eta <- .np_linpred(x, model_pars, pars, newdata, NA, draw_ids, dots)
  links <- .np_population_links(x, names(eta), native)

  if (.np_no_group_effects(x, re_formula)) {
    message2(
      "There are no group-level effects to integrate over, so the population \\
      mean and median equal the parameters at zero group-level effects, and \\
      their spread is zero."
    )
    values <- if (native) .np_transform_checked(x, eta, newdata) else eta
    stats <- lapply(names(values)[names(values) %in% pars], function(par) {
      .np_population_degenerate(values[[par]], links[[par]], population_summary)
    })
    names(stats) <- names(values)[names(values) %in% pars]
    return(.np_population_result(stats, .np_population_method_none(stats, links)))
  }

  prep <- .np_re_design(x, newdata, re_formula, draw_ids[1])
  analytic <- .np_analytic_pars(x, eta, prep, newdata, native, links)
  monte_carlo <- names(eta)[not_in(names(eta), analytic)]
  analytic <- analytic[analytic %in% pars]
  monte_carlo <- monte_carlo[monte_carlo %in% pars]

  stats <- list()
  method <- list()
  for (par in analytic) {
    latent <- .np_latent_sd(x, prep, par, draw_ids)
    stats[[par]] <- .np_population_normal(
      eta[[par]], latent$sd, links[[par]], population_summary
    )
    method[[par]] <- .np_population_method_analytic(
      par, names(stats[[par]]), links[[par]], latent
    )
  }

  if (length(monte_carlo) > 0) {
    mc_stats <- .np_population_mc(
      x, names(eta), monte_carlo, newdata, grid_vars, re_formula, draw_ids,
      population_summary, native, links, ndraws_population, dots
    )
    groups <- collapse_comma(.np_prep_groups(prep))
    for (par in monte_carlo) {
      stats[[par]] <- mc_stats[[par]]
      method[[par]] <- data.frame(
        parameter = par,
        statistic = names(mc_stats[[par]]),
        link = links[[par]],
        method = "Monte Carlo",
        n = ndraws_population,
        groups = groups,
        stringsAsFactors = FALSE
      )
    }
    circular <- native & unlist(links[monte_carlo]) == "tan_half"
    .np_check_mc_error(mc_stats[!circular], ndraws_population)
  }

  order <- names(eta)[names(eta) %in% names(stats)]
  .np_population_result(stats[order], do.call(rbind, unname(method[order])))
}


.np_population_result <- function(stats, method) {
  row.names(method) <- NULL
  list(
    values = unlist(unname(stats), recursive = FALSE),
    parameter = rep(names(stats), lengths(stats)),
    statistic = unlist(lapply(stats, names), use.names = FALSE),
    method = method
  )
}


#' The link each parameter is summarised under
#'
#' @description
#' On the sampling scale every parameter is its own normal linear predictor, so
#' the summaries are those of the identity link.
#'
#' @keywords internal
#' @noRd
.np_population_links <- function(x, pars, native) {
  links <- lapply(pars, function(par) {
    if (native) x$bmm$model$links[[par]] %||% "identity" else "identity"
  })
  stats::setNames(links, pars)
}


.np_no_group_effects <- function(x, re_formula) {
  nrow(x$ranef) == 0 || (!is.null(re_formula) && !is_formula(re_formula))
}


#' Statistics of a distribution with no spread
#'
#' @keywords internal
#' @noRd
.np_population_degenerate <- function(value, link, population_summary) {
  zero <- value
  zero[] <- 0
  if (population_summary == "mean") {
    return(list(mean = value, sd = zero))
  }
  if (link == "tan_half") {
    return(list(median = value))
  }
  list(median = value, q25 = value, q75 = value)
}


.np_population_method_none <- function(stats, links) {
  rows <- lapply(names(stats), function(par) {
    data.frame(
      parameter = par,
      statistic = names(stats[[par]]),
      link = links[[par]],
      method = "no group-level effects",
      n = NA_integer_,
      groups = "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


.np_population_method_analytic <- function(par, statistics, link, latent) {
  if (length(latent$groups) == 0) {
    method <- rep("no group-level effects", length(statistics))
    n <- rep(NA_integer_, length(statistics))
  } else {
    method <- vapply(statistics, .np_analytic_method, character(1), link = link)
    n <- vapply(statistics, .np_analytic_nodes, integer(1), link = link)
  }
  data.frame(
    parameter = par,
    statistic = statistics,
    link = link,
    method = unname(method),
    n = unname(n),
    groups = collapse_comma(latent$groups),
    stringsAsFactors = FALSE
  )
}


.np_analytic_method <- function(statistic, link) {
  if (link == "tan_half") {
    return(if (statistic == "median") "root search" else "Gauss-Hermite")
  }
  closed <- statistic %in% c("median", "q25", "q75") ||
    link %in% c("identity", "log") ||
    (link == "probit" && statistic == "mean")
  if (closed) "closed form" else "Gauss-Hermite"
}


.np_analytic_nodes <- function(statistic, link) {
  if (.np_analytic_method(statistic, link) != "Gauss-Hermite") {
    return(NA_integer_)
  }
  if (link == "tan_half") 200L else 100L
}


#' Group-level design of the prediction grid
#'
#' @description
#' [brms::prepare_predictions()] applies `re_formula` the way every `brms`
#' prediction function does: it replaces the group-level terms of each linear
#' formula with those of `re_formula`. The `Z` blocks it returns therefore hold
#' exactly the terms the population summary integrates over. `Z` does not depend
#' on the level, so the grid is predicted for new levels, whose effects `brms`
#' samples; the seed is restored so the caller's random number stream is not
#' advanced.
#'
#' @keywords internal
#' @noRd
.np_re_design <- function(x, newdata, re_formula, draw_id) {
  withr::with_preserve_seed(brms::prepare_predictions(
    x,
    newdata = newdata, re_formula = re_formula, draw_ids = draw_id,
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
  ))
}


# A distributional parameter fixed to a constant is stored as that number, with
# no terms to predict from
.np_prep_term <- function(prep, par) {
  term <- prep$dpars[[par]] %||% prep$nlpars[[par]]
  if (is.list(term)) term
}


.np_prep_groups <- function(prep) {
  terms <- Filter(is.list, c(prep$dpars, prep$nlpars))
  unique(unlist(lapply(terms, function(term) names(term$re$Z))))
}


#' Parameters whose population summaries have a closed form or quadrature
#'
#' @description
#' Three conditions: the parameter enters `brms` through a linear formula, so it
#' is normal given the draw; its group-level terms are plain `(... | g)` terms,
#' whose effects are normal with the covariance matrix of the `sd_` and `cor_`
#' draws; and on the native scale, `native_transform()` of the draws at zero
#' group-level effects equals the inverse link of the parameter's own linear
#' predictor. The last condition excludes softmax groups and transformations
#' that depend on the design or on other parameters, which a model's
#' `native_transform()` method may define.
#'
#' @keywords internal
#' @noRd
.np_analytic_pars <- function(x, eta, prep, newdata, native, links) {
  bterms <- brms::brmsterms(x$formula)
  pars <- names(eta)
  linear <- vapply(pars, function(par) {
    inherits(bterms$dpars[[par]] %||% bterms$nlpars[[par]], "btl")
  }, logical(1))
  plain <- vapply(pars, .np_plain_re, logical(1), x = x, prep = prep)
  ok <- linear & plain
  if (!native) {
    return(pars[ok])
  }

  transformed <- .np_transform_checked(x, eta, newdata)
  elementwise <- vapply(pars, function(par) {
    if (links[[par]] == "softmax") {
      return(FALSE)
    }
    isTRUE(all.equal(
      unname(transformed[[par]]),
      unname(link_transform(eta[[par]], links[[par]], inverse = TRUE)),
      tolerance = 1e-12
    ))
  }, logical(1))
  pars[ok & elementwise]
}


.np_plain_re <- function(x, prep, par) {
  term <- .np_prep_term(prep, par)
  if (is.null(term)) {
    return(FALSE)
  }
  special <- c(term$re$Zsp, term$re$Zcs)
  if (any(!vapply(special, is.null, logical(1)))) {
    return(FALSE)
  }
  ranef <- .np_par_ranef(x, par)
  all(
    ranef$type == "", ranef$gtype == "", ranef$by == "", ranef$cov == "",
    ranef$dist == "gaussian"
  )
}


.np_par_ranef <- function(x, par) {
  ranef <- as.data.frame(x$ranef)
  dpar <- ifelse(ranef$dpar == "", "mu", ranef$dpar)
  ranef[ranef$nlpar == par | (ranef$nlpar == "" & dpar == par), , drop = FALSE]
}


#' Latent between-level SD of a parameter per draw and grid cell
#'
#' @description
#' Different `id` blocks are independent, so their contributions add up. Within
#' a block, the coefficients of the parameter are correlated when the term was
#' written with `|` and the block estimates correlations.
#'
#' @return A list with `sd`, a draws x cells matrix, and `groups`, the grouping
#'   factors that contribute to it.
#'
#' @keywords internal
#' @noRd
.np_latent_sd <- function(x, prep, par, draw_ids) {
  re <- .np_prep_term(prep, par)$re
  n_cells <- prep$nobs
  variance <- matrix(0, length(draw_ids), n_cells)
  groups <- character(0)
  ranef <- .np_par_ranef(x, par)

  for (group in names(re$Z)) {
    coefs <- sub("^.*,([^,]*)\\]$", "\\1", colnames(re$r[[group]]))
    design <- as.matrix(re$Z[[group]])
    z <- vapply(unique(coefs), function(coef) {
      rowSums(design[, coefs == coef, drop = FALSE])
    }, numeric(n_cells))
    z <- matrix(z, nrow = n_cells, dimnames = list(NULL, unique(coefs)))
    if (all(z == 0)) {
      next
    }
    covariance <- .np_re_covariance(x, ranef[ranef$group == group, , drop = FALSE],
                                    colnames(z), group, draw_ids)
    variance <- variance + .np_quadratic_form(z, covariance$sd, covariance$cor)
    groups <- c(groups, group)
  }

  list(sd = sqrt(variance), groups = unique(groups))
}


#' SD and correlation draws of the group-level coefficients of one parameter
#'
#' @description
#' The variable names follow `brms`: `sd_<group>__<prefix><coef>` and
#' `cor_<group>__<name1>__<name2>`, where the prefix names the response,
#' distributional parameter (other than `mu`) and non-linear parameter. A
#' correlation is drawn only within an `id` block that estimates it.
#'
#' @return A list with `sd`, a draws x coefficients matrix, and `cor`, a list of
#'   `list(j, k, rho)` entries for the correlated pairs.
#'
#' @keywords internal
#' @noRd
.np_re_covariance <- function(x, ranef, coefs, group, draw_ids) {
  rows <- ranef[match(coefs, ranef$coef), , drop = FALSE]
  stopif(
    anyNA(rows$coef),
    "Internal error: the group-level coefficients {collapse_comma(coefs)} of \\
    group '{group}' were not found in the model. Please report this at \\
    https://github.com/popov-lab/bmm/issues"
  )
  prefix <- mapply(.np_ranef_prefix, rows$resp, rows$dpar, rows$nlpar)
  labels <- paste0(prefix, rows$coef)
  sd_vars <- paste0("sd_", group, "__", labels)

  pairs <- list()
  for (j in seq_along(coefs)) {
    for (k in seq_along(coefs)) {
      if (k > j && rows$id[j] == rows$id[k] && rows$cor[j] && rows$cor[k]) {
        pairs[[length(pairs) + 1]] <- list(
          j = j, k = k,
          vars = paste0("cor_", group, "__", labels[c(j, k)], "__", labels[c(k, j)])
        )
      }
    }
  }

  available <- posterior::variables(x)
  cor_vars <- vapply(pairs, function(pair) {
    found <- pair$vars[pair$vars %in% available]
    if (length(found) == 0) NA_character_ else found[1]
  }, character(1))
  missing <- c(sd_vars[not_in(sd_vars, available)], unlist(lapply(pairs, `[[`, "vars"))[is.na(cor_vars)])
  stopif(
    length(missing) > 0,
    "Internal error: the posterior has no variable(s) {collapse_comma(missing)}. \\
    Please report this at https://github.com/popov-lab/bmm/issues"
  )

  draws <- posterior::merge_chains(
    posterior::as_draws_matrix(x, variable = c(sd_vars, cor_vars))
  )
  draws <- unclass(posterior::subset_draws(draws, draw = draw_ids))
  list(
    sd = draws[, sd_vars, drop = FALSE],
    cor = lapply(seq_along(pairs), function(i) {
      list(j = pairs[[i]]$j, k = pairs[[i]]$k, rho = draws[, cor_vars[i]])
    })
  )
}


.np_ranef_prefix <- function(resp, dpar, nlpar) {
  parts <- c(resp, if (dpar != "mu") dpar, nlpar)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) "" else paste0(paste(parts, collapse = "_"), "_")
}


#' z' Sigma z for every draw and grid cell
#'
#' @param z Cells x coefficients matrix of group-level design values
#' @param sds Draws x coefficients matrix of SDs
#' @param cors List of `list(j, k, rho)`, `rho` a vector of draws
#' @return Draws x cells matrix of latent variances
#'
#' @keywords internal
#' @noRd
.np_quadratic_form <- function(z, sds, cors) {
  out <- tcrossprod(sds^2, z^2)
  for (pair in cors) {
    j <- pair$j
    k <- pair$k
    out <- out + 2 * tcrossprod(pair$rho * sds[, j] * sds[, k], z[, j] * z[, k])
  }
  out
}


#' Population statistics of g(eta + s Z), Z standard normal
#'
#' @param eta,s Matrices of the latent mean and SD, of the same dimensions
#' @param link The parameter's link, `"identity"` on the sampling scale
#' @param population_summary `"mean"` or `"median"`
#' @return A named list of matrices: `mean` and `sd`, or `median`, `q25` and
#'   `q75`. For `"tan_half"`, the circular mean and circular SD, or the circular
#'   median alone.
#'
#' @keywords internal
#' @noRd
.np_population_normal <- function(eta, s, link, population_summary) {
  if (population_summary == "median") {
    return(.np_population_quartiles(eta, s, link))
  }
  .np_population_moments(eta, s, link)
}


.np_population_moments <- function(eta, s, link) {
  if (link == "tan_half") {
    return(.np_circular_moments(eta, s))
  }
  if (link == "inverse") {
    warning2(
      "The population mean of a parameter with an 'inverse' link does not \\
      exist, because its linear predictor can be zero. Returning NA."
    )
    na <- eta
    na[] <- NA_real_
    return(list(mean = na, sd = na))
  }

  inverse <- function(values) link_transform(values, link, inverse = TRUE)
  out <- switch(link,
    identity = list(mean = eta, sd = s),
    log = {
      mean <- exp(eta + s^2 / 2)
      list(mean = mean, sd = mean * sqrt(expm1(s^2)))
    },
    probit = {
      mean <- stats::pnorm(eta / sqrt(1 + s^2))
      list(mean = mean, sd = .np_gh_sd(eta, s, inverse, mean, .np_gh_nodes(100)))
    },
    {
      mean <- .np_gh_mean(eta, s, inverse, .np_gh_nodes(100))
      list(mean = mean, sd = .np_gh_sd(eta, s, inverse, mean, .np_gh_nodes(100)))
    }
  )

  zero <- s == 0
  out$mean[zero] <- inverse(eta[zero])
  out$sd[zero] <- 0
  out
}


.np_population_quartiles <- function(eta, s, link) {
  if (link == "tan_half") {
    return(list(median = .np_circular_median(eta, s)))
  }
  inverse <- function(values) link_transform(values, link, inverse = TRUE)
  lower <- inverse(eta - stats::qnorm(0.75) * s)
  upper <- inverse(eta + stats::qnorm(0.75) * s)
  # pmin/pmax keep a decreasing inverse link's quartiles in order
  list(median = inverse(eta), q25 = pmin(lower, upper), q75 = pmax(lower, upper))
}


#' Gauss-Hermite nodes and weights for the expectation over a standard normal
#'
#' @description
#' Golub-Welsch: the nodes are the eigenvalues of the Jacobi matrix of the
#' probabilists' Hermite polynomials, the weights the squared first components
#' of its eigenvectors. Exact for polynomials up to degree `2n - 1`.
#'
#' @keywords internal
#' @noRd
.np_gh_nodes <- function(n) {
  off_diagonal <- sqrt(seq_len(n - 1))
  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n - 1), 2:n)] <- off_diagonal
  jacobi[cbind(2:n, seq_len(n - 1))] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(x = decomposition$values, w = decomposition$vectors[1, ]^2)
}


.np_gh_mean <- function(eta, s, f, nodes) {
  out <- 0
  for (k in seq_along(nodes$x)) {
    out <- out + nodes$w[k] * f(eta + s * nodes$x[k])
  }
  out
}


# Centred on the mean rather than E[f^2] - E[f]^2, which cancels when the
# spread is small relative to the mean
.np_gh_sd <- function(eta, s, f, mean, nodes) {
  out <- 0
  for (k in seq_along(nodes$x)) {
    out <- out + nodes$w[k] * (f(eta + s * nodes$x[k]) - mean)^2
  }
  sqrt(out)
}


#' Circular mean and circular SD of 2 * atan(eta + s Z)
#'
#' @description
#' The definitions of the `circular` package: the mean direction
#' `atan2(E sin, E cos)` and the circular SD `sqrt(-2 log rho)`, with
#' `rho = |E exp(i theta)|`. With `theta = 2 atan(v)`, `cos theta` and
#' `sin theta` are rational functions of `v`, integrated by Gauss-Hermite
#' quadrature with 200 nodes.
#'
#' @keywords internal
#' @noRd
.np_circular_moments <- function(eta, s) {
  nodes <- .np_gh_nodes(200)
  cosine <- 0
  sine <- 0
  for (k in seq_along(nodes$x)) {
    v <- eta + s * nodes$x[k]
    cosine <- cosine + nodes$w[k] * (1 - v^2) / (1 + v^2)
    sine <- sine + nodes$w[k] * 2 * v / (1 + v^2)
  }
  out <- list(mean = atan2(sine, cosine), sd = sqrt(pmax(-log(cosine^2 + sine^2), 0)))
  zero <- s == 0
  out$mean[zero] <- 2 * atan(eta[zero])
  out$sd[zero] <- 0
  out
}


#' Circular median of 2 * atan(eta + s Z)
#'
#' @description
#' The direction `m` minimising the expected arc distance `E d(theta, m)`. Its stationary points are the diameters that split the mass in half,
#' `P(m < theta < m + pi) = 1/2`; of the two ends of such a diameter the minimum
#' is the end with the higher density. The diameters are found by a grid over
#' `(-pi, 0]` refined by bisection, from the closed-form distribution function
#' `Phi((tan(t / 2) - eta) / s)`. When several diameters qualify, the
#' distribution is multimodal and the one with the smallest expected arc
#' distance is returned.
#'
#' @param eta,s Matrices or vectors of the latent mean and SD
#' @return An object of the shape of `eta`, in radians
#'
#' @keywords internal
#' @noRd
.np_circular_median <- function(eta, s) {
  out <- 2 * atan(eta)
  todo <- which(s > 0)
  if (length(todo) == 0) {
    return(out)
  }
  location <- eta[todo]
  scale <- s[todo]
  cdf <- function(t, i) stats::pnorm((tan(t / 2) - location[i]) / scale[i])
  half <- function(m, i) cdf(m + pi, i) - cdf(m, i) - 0.5
  density <- function(t, i) {
    stats::dnorm((tan(t / 2) - location[i]) / scale[i]) * (1 + tan(t / 2)^2)
  }

  n <- length(todo)
  grid <- seq(-pi, 0, length.out = 129)
  positive <- vapply(grid, function(m) half(m, seq_len(n)) > 0, logical(n))
  positive <- matrix(positive, nrow = n)
  change <- which(
    positive[, -1, drop = FALSE] != positive[, -length(grid), drop = FALSE],
    arr.ind = TRUE
  )

  pair <- change[, 1]
  lower <- grid[change[, 2]]
  upper <- grid[change[, 2] + 1]
  lower_positive <- positive[change]
  for (iteration in seq_len(60)) {
    middle <- (lower + upper) / 2
    same <- (half(middle, pair) > 0) == lower_positive
    lower[same] <- middle[same]
    upper[!same] <- middle[!same]
  }

  # a symmetric bimodal distribution can leave h <= 0 on the whole grid, with
  # the diameters at its ends
  bracketless <- setdiff(seq_len(n), pair)
  root <- c((lower + upper) / 2, rep(0, length(bracketless)))
  pair <- c(pair, bracketless)
  candidate <- ifelse(density(root, pair) >= density(root + pi, pair), root, root + pi)

  multiple <- pair %in% pair[duplicated(pair)]
  if (any(multiple)) {
    distance <- rep(Inf, length(candidate))
    distance[multiple] <- .np_expected_arc_distance(candidate[multiple], pair[multiple], cdf)
    keep <- order(pair, distance)
    keep <- keep[!duplicated(pair[keep])]
    candidate <- candidate[keep]
    pair <- pair[keep]
  }

  out[todo[pair]] <- candidate
  out
}


#' E d(theta, m) = integral over r in (0, pi) of P(d(theta, m) > r)
#'
#' @keywords internal
#' @noRd
.np_expected_arc_distance <- function(m, pair, cdf) {
  unwrapped_cdf <- function(t, i) {
    turns <- floor((t + pi) / (2 * pi))
    turns + cdf(t - 2 * pi * turns, i)
  }
  radii <- (seq_len(512) - 0.5) * pi / 512
  outside <- vapply(radii, function(r) {
    1 - (unwrapped_cdf(m + r, pair) - unwrapped_cdf(m - r, pair))
  }, numeric(length(m)))
  rowSums(matrix(outside, nrow = length(m))) * pi / 512
}


#' Fisher's circular median of a sample
#'
#' @description
#' The arc-distance objective is piecewise linear between the sample points and
#' their antipodes, and its minima lie at sample points. Sorting the sample and
#' unrolling it once around the circle gives the objective at every sample
#' point from cumulative sums.
#'
#' @param theta Numeric vector of angles in radians
#' @return The sample point minimising the summed arc distance
#'
#' @keywords internal
#' @noRd
.np_circular_sample_median <- function(theta) {
  n <- length(theta)
  sorted <- sort(theta)
  unrolled <- c(sorted, sorted + 2 * pi)
  cumulative <- c(0, cumsum(unrolled))
  k <- seq_len(n)
  ahead <- findInterval(sorted + pi, unrolled)
  n_ahead <- ahead - k + 1
  distance_ahead <- cumulative[ahead + 1] - cumulative[k] - n_ahead * sorted
  distance_behind <- (n - n_ahead) * (sorted + 2 * pi) - (cumulative[k + n] - cumulative[ahead + 1])
  sorted[which.min(distance_ahead + distance_behind)]
}


#' Population statistics by Monte Carlo over simulated new levels
#'
#' @description
#' Every grid cell is repeated for `ndraws_population` new levels of each
#' grouping factor, with the same level names in every cell so that contrasts
#' between cells share their simulated levels. `brms` samples the effects of all
#' new levels once per call, for all parameters at once, so replaying the random
#' number state before each parameter's call gives every parameter the same
#' simulated levels, and with them the correlations between parameters. The
#' parameters are then transformed jointly by `native_transform()`. The draws
#' are processed in chunks to bound memory. The caller's random number stream is
#' used and left advanced.
#'
#' @return A named list, per parameter, of named lists of draws x cell matrices
#'
#' @keywords internal
#' @noRd
.np_population_mc <- function(x, linpred_pars, mc_pars, newdata, grid_vars,
                              re_formula, draw_ids, population_summary, native,
                              links, ndraws_population, dots) {
  n_cells <- nrow(newdata)
  levels <- paste0("bmm_population_", seq_len(ndraws_population))
  persons <- newdata[rep(seq_len(n_cells), each = ndraws_population), , drop = FALSE]
  for (var in setdiff(.np_group_vars(x), grid_vars)) {
    persons[[var]] <- rep(levels, times = n_cells)
  }
  row.names(persons) <- NULL

  chunk_size <- max(1L, floor(1e7 / (n_cells * ndraws_population * length(linpred_pars))))
  chunks <- split(draw_ids, ceiling(seq_along(draw_ids) / chunk_size))
  dots <- c(dots, list(allow_new_levels = TRUE, sample_new_levels = "gaussian"))

  per_chunk <- lapply(chunks, function(ids) {
    linpred <- .np_linpred(
      x, linpred_pars, linpred_pars, persons, re_formula, ids, dots,
      replay_seed = TRUE
    )
    if (native) {
      linpred <- .np_transform_checked(x, linpred, persons)
    }
    lapply(stats::setNames(mc_pars, mc_pars), function(par) {
      .np_mc_stats(
        linpred[[par]], n_cells, ndraws_population, population_summary,
        circular = native && links[[par]] == "tan_half"
      )
    })
  })

  lapply(stats::setNames(mc_pars, mc_pars), function(par) {
    statistics <- names(per_chunk[[1]][[par]])
    stats::setNames(lapply(statistics, function(statistic) {
      do.call(rbind, lapply(per_chunk, function(chunk) chunk[[par]][[statistic]]))
    }), statistics)
  })
}


.np_group_vars <- function(x) {
  unique(unlist(lapply(names(x$bmm$model$parameters), function(par) {
    lapply(.formula_nodes(x, par), function(node) .extract_re_grouping_vars(node$rhs))
  })))
}


#' Statistics across the simulated levels of each grid cell
#'
#' @param values Draws x (cells * levels) matrix, the levels of a cell adjacent
#' @return A named list of draws x cells matrices
#'
#' @keywords internal
#' @noRd
.np_mc_stats <- function(values, n_cells, n_levels, population_summary, circular) {
  per_cell <- lapply(seq_len(n_cells), function(cell) {
    .np_level_stats(
      values[, (cell - 1) * n_levels + seq_len(n_levels), drop = FALSE],
      population_summary, circular
    )
  })
  statistics <- names(per_cell[[1]])
  stats::setNames(lapply(statistics, function(statistic) {
    matrix(
      vapply(per_cell, `[[`, numeric(nrow(values)), statistic),
      nrow = nrow(values)
    )
  }), statistics)
}


.np_level_stats <- function(values, population_summary, circular) {
  if (population_summary == "mean" && circular) {
    cosine <- rowMeans(cos(values))
    sine <- rowMeans(sin(values))
    return(list(
      mean = atan2(sine, cosine),
      sd = sqrt(pmax(-log(cosine^2 + sine^2), 0))
    ))
  }
  if (population_summary == "mean") {
    return(list(mean = rowMeans(values), sd = matrixStats::rowSds(values)))
  }
  if (circular) {
    return(list(median = apply(values, 1, .np_circular_sample_median)))
  }
  quartiles <- matrix(
    matrixStats::rowQuantiles(values, probs = c(0.5, 0.25, 0.75)),
    nrow = nrow(values)
  )
  list(median = quartiles[, 1], q25 = quartiles[, 2], q75 = quartiles[, 3])
}


#' Warn when the Monte Carlo error widens the posterior noticeably
#'
#' @description
#' The posterior variance of a Monte Carlo estimate across draws is the
#' posterior variance of the estimand plus the Monte Carlo variance within a
#' draw: `sd^2 / n` for the mean and, under a normal approximation,
#' `(pi / 2) (IQR / 1.349)^2 / n` for the median. Their ratio is how much wider
#' the reported posterior SD is than it would be with exact integration.
#' Circular statistics are not checked.
#'
#' @keywords internal
#' @noRd
.np_check_mc_error <- function(mc_stats, ndraws_population) {
  inflation <- vapply(mc_stats, function(stats) {
    .np_mc_inflation(stats, ndraws_population)
  }, numeric(1))
  inflation <- inflation[!is.na(inflation)]
  if (length(inflation) == 0 || max(inflation) <= 1.02) {
    return(invisible())
  }
  worst <- names(inflation)[which.max(inflation)]
  stats <- mc_stats[[worst]]
  within <- .np_mc_within_variance(stats, ndraws_population)
  between <- pmax(.np_mc_between_variance(stats, within), 0)
  # 1.02^2 - 1: the Monte Carlo share of the variance that widens the SD by 2%
  needed <- ceiling(max(ndraws_population * within / (0.0404 * between), na.rm = TRUE))
  advice <- if (is.finite(needed)) {
    "Set 'ndraws_population' to at least {needed} to bring this below 2%."
  } else {
    "Increase 'ndraws_population'."
  }
  percent <- if (is.finite(max(inflation))) {
    glue::glue("by {round(100 * (max(inflation) - 1))}%")
  } else {
    "so much that it hides the posterior uncertainty"
  }
  warning2(
    "The Monte Carlo error of ndraws_population = {ndraws_population} widens the \\
    posterior SD of the population {names(stats)[1]} of '{worst}' {percent}. ",
    advice
  )
}


.np_mc_inflation <- function(stats, ndraws_population) {
  within <- .np_mc_within_variance(stats, ndraws_population)
  if (is.null(within)) {
    return(NA_real_)
  }
  between <- .np_mc_between_variance(stats, within)
  observed <- between + within
  inflation <- ifelse(observed == 0, 1, sqrt(observed / pmax(between, 0)))
  max(inflation, na.rm = TRUE)
}


.np_mc_within_variance <- function(stats, ndraws_population) {
  if (nrow(stats[[1]]) < 2) {
    return(NULL)
  }
  if (!is.null(stats$sd) && !is.null(stats$mean)) {
    return(colMeans(stats$sd^2) / ndraws_population)
  }
  if (!is.null(stats$q25)) {
    iqr <- stats$q75 - stats$q25
    return(pi / 2 * colMeans((iqr / 1.349)^2) / ndraws_population)
  }
  NULL
}


.np_mc_between_variance <- function(stats, within) {
  observed <- apply(stats[[1]], 2, stats::var)
  observed - within
}
