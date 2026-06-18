############################################################################# !
# HELPERS-ANALYSIS.R                                                     ####
# Post-fitting analysis functions for SDT bmmfit objects (ROC / AUC).   ####
############################################################################# !


############################################################################# !
# ROC CURVES                                                             ####
############################################################################# !

#' Model-implied ROC curve from a fitted SDT model
#'
#' Computes the receiver operating characteristic (ROC) curve implied by the
#' posterior distribution of a signal detection model fit with [bmm()]. ROC
#' curves require a response criterion, so they are defined for [sdt_binary()]
#' and [sdt_rating()] only (not the criterion-free [sdt_mafc()] and
#' [sdt_ranking()]).
#'
#' For **binary** models the ROC is traced analytically from the posterior of
#' `dprime` (and `sdratio` for unequal-variance SDT) over a grid of criterion
#' values. When the criterion varies across conditions but `dprime`/`sdratio`
#' do not (e.g. a base-rate manipulation, as in [broeder_schuetz_2009_e3]), the
#' several criteria are operating points on a single curve: `roc_sdt()` returns
#' one smooth curve and attaches the model-implied points (one per criterion
#' level) as a `points` attribute.
#'
#' For **rating** models the K-1 confidence thresholds define K-1 empirical ROC
#' points per posterior draw.
#'
#' @param fit A `bmmfit` object returned by [bmm()] from an SDT model.
#' @param conditions Optional data frame of predictor values for which to
#'   compute the ROC curve. Column names must match predictor variables used in
#'   the formula. If `NULL` (default), unique predictor combinations are derived
#'   from the data.
#' @param n_points Integer. Number of equally-spaced points on the smooth ROC
#'   curve for binary models (default 100). Ignored for rating models, which
#'   return K+1 points (K-1 threshold points plus the (0,0) and (1,1)
#'   endpoints).
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible band (default `c(0.025, 0.975)`).
#' @param criterion_points Optional. Control of the binary multi-criteria
#'   behaviour. `NULL` (default) auto-detects predictors that vary the criterion
#'   only (not `dprime`/`sdratio`) and treats their levels as operating points
#'   on one curve. Pass a character vector of column names to force that
#'   classification, or `FALSE` to disable it (one separate curve per predictor
#'   combination). Ignored for rating models.
#' @param ... Additional arguments passed to [brms::posterior_linpred()].
#'
#' @return A data frame of class `"bmm_sdt_roc"` with columns `FA`, `Hit`,
#'   `.draw`, and any condition columns. The object carries a `summary`
#'   attribute (`FA`, `Hit_mean`, `Hit_lower`, `Hit_upper`) and, for binary
#'   multi-criteria fits, a `points` attribute with the model-implied operating
#'   points.
#'
#' @seealso [auc_sdt()], [roc_observed()], [plot.bmm_sdt_roc()]
#' @export
roc_sdt <- function(fit, conditions = NULL, n_points = 100,
                    probs = c(0.025, 0.975), criterion_points = NULL, ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "roc_sdt() is only available for SDT models (sdt_binary, sdt_rating)")
  stopif(inherits(model, "sdt_mafc"),
         "ROC curves are not defined for the m-AFC SDT model: it has no response criterion.")
  stopif(inherits(model, "sdt_ranking"),
         "ROC curves are not defined for the ranking SDT model: it has no response criterion.")

  is_rating  <- inherits(model, "sdt_rating")
  conditions <- .sdt_resolve_conditions(fit, conditions)

  if (is_rating) {
    roc_data    <- .roc_sdt_rating(fit, model, conditions, ...)
    points      <- NULL
    roc_summary <- .roc_sdt_summary(roc_data, probs)
  } else {
    binary      <- .roc_sdt_binary(fit, model, conditions, n_points,
                                   criterion_points = criterion_points,
                                   probs = probs, ...)
    roc_data    <- binary$curve
    points      <- binary$points
    roc_summary <- binary$summary
  }

  structure(
    roc_data,
    class       = c("bmm_sdt_roc", "data.frame"),
    summary     = roc_summary,
    points      = points,
    probs       = probs,
    model_class = class(model),
    dist        = model$other_vars$dist,
    is_rating   = is_rating,
    conditions  = conditions
  )
}


# Unique predictor combinations from fit data. Excludes response columns, the
# SDT design variables added by bmm (stimulus, n_trials, the multinomial Y
# matrix, dist_type), random-effect grouping variables, and matrix columns.
.sdt_resolve_conditions <- function(fit, conditions) {
  if (!is.null(conditions)) return(as.data.frame(conditions))

  data  <- fit$data
  model <- fit$bmm$model

  exclude <- unique(c(
    unlist(model$resp_vars),
    model$other_vars$stimulus,
    model$other_vars$n_trials,
    "Y", "nTrials", "dist_type",
    names(brms::ranef(fit))
  ))
  pred_cols <- setdiff(names(data), exclude)

  is_matrix_col <- vapply(pred_cols, function(col) is.matrix(data[[col]]), logical(1))
  pred_cols <- pred_cols[!is_matrix_col]

  if (length(pred_cols) == 0L) {
    return(data[1L, character(0), drop = FALSE])
  }
  out <- unique(data[, pred_cols, drop = FALSE])
  rownames(out) <- NULL
  out
}


# Classify the formula predictors of a binary fit into criterion "operating
# points" (predict the criterion only) and "curves" (predict dprime/sdratio).
# Random-effect grouping factors are stripped so (1 | id) is never a dimension.
.sdt_criterion_point_dims <- function(fit) {
  uf <- fit$bmm$user_formula
  preds <- if (inherits(uf, "bmmformula")) rhs_vars(uf, collapse = FALSE) else list()
  re_vars <- tryCatch(names(brms::ranef(fit)), error = function(e) character(0))
  strip <- function(v) setdiff(v %||% character(0), re_vars)

  curve_vars <- union(strip(preds[["dprime"]]), strip(preds[["sdratio"]]))
  list(points = setdiff(strip(preds[["criterion"]]), curve_vars),
       curves = curve_vars)
}


# Resolve point/curve dimensions for a binary fit, honouring the user override.
.sdt_resolve_point_dims <- function(fit, conditions, criterion_points) {
  all_vars <- names(conditions)
  if (isFALSE(criterion_points)) {
    return(list(points = character(0), curves = all_vars))
  }
  if (!is.null(criterion_points)) {
    pts <- intersect(criterion_points, all_vars)
    return(list(points = pts, curves = setdiff(all_vars, pts)))
  }
  dims   <- .sdt_criterion_point_dims(fit)
  pts    <- intersect(dims$points, all_vars)
  curves <- setdiff(all_vars, pts)
  list(points = pts, curves = curves)
}


# 1-row-0-column placeholder when `vars` is empty, else the unique rows of the
# named subset. Keeps the per-condition loops running exactly once for
# intercept-only / population-level cases.
.sdt_unique_subset <- function(conditions, vars) {
  vars <- intersect(vars, names(conditions))
  if (length(vars) == 0L) {
    return(conditions[1L, character(0), drop = FALSE])
  }
  out <- unique(conditions[, vars, drop = FALSE])
  rownames(out) <- NULL
  out
}


# Whether sdratio is an estimated (not fixed) parameter. The model object is
# authoritative once update_model_fixed_parameters() has dropped sdratio from
# fixed_parameters for a `sdratio ~ ...` formula; the brms::variables() check is
# a fail-safe against a desynchronised model object.
.sdt_has_estimated_sdratio <- function(model, fit = NULL) {
  obj_free <- "sdratio" %in% names(model$parameters) &&
    !"sdratio" %in% names(model$fixed_parameters)
  if (is.null(fit)) return(obj_free)
  vars <- tryCatch(brms::variables(fit), error = function(e) character(0))
  obj_free || any(grepl("(^|_)b(sp)?_sdratio($|_)", vars))
}


# posterior_linpred for one parameter, with conditions overriding the first data
# row so every formula variable is present. dprime/criterion/sdratio are dpars
# for binary and nlpars for rating; pass the literal name to brms.
.sdt_linpred <- function(fit, param, conditions, is_rating = FALSE, ...) {
  newdata <- fit$data[rep(1L, max(1L, nrow(conditions))), , drop = FALSE]
  rownames(newdata) <- NULL
  if (ncol(conditions) > 0L) {
    cols <- intersect(names(conditions), names(newdata))
    newdata[cols] <- conditions[cols]
  }
  if (is_rating) {
    brms::posterior_linpred(fit, nlpar = param, newdata = newdata,
                            re_formula = NA, allow_new_levels = TRUE, ...)
  } else {
    brms::posterior_linpred(fit, dpar = param, newdata = newdata,
                            re_formula = NA, allow_new_levels = TRUE, ...)
  }
}


# Smooth analytical ROC for binary EV- and UV-SDT plus, when the criterion
# varies across conditions only, the model-implied operating points. The curve
# is traced over a shared criterion grid and summarised as the posterior-mean
# operating point (FA, Hit) at each node -- the same functional as the
# model-implied points -- so those points fall on the curve. Averaging Hit at a
# fixed FA instead would put them slightly off it on the concave ROC.
.roc_sdt_binary <- function(fit, model, conditions, n_points,
                            criterion_points = NULL,
                            probs = c(0.025, 0.975), ...) {
  dist        <- model$other_vars$dist
  has_sdratio <- .sdt_has_estimated_sdratio(model, fit)
  cdf         <- .SDT_DISTS[[dist]]$cdf
  qf          <- .SDT_DISTS[[dist]]$qf

  dims       <- .sdt_resolve_point_dims(fit, conditions, criterion_points)
  curve_cond <- .sdt_unique_subset(conditions, dims$curves)

  n_curve     <- max(1L, nrow(curve_cond))
  dprime_mat  <- .sdt_linpred(fit, "dprime", curve_cond, is_rating = FALSE, ...)
  n_draws     <- nrow(dprime_mat)
  sdratio_mat <- if (has_sdratio) {
    exp(.sdt_linpred(fit, "sdratio", curve_cond, is_rating = FALSE, ...))
  } else {
    matrix(1, nrow = n_draws, ncol = n_curve)
  }

  fa_grid        <- seq(0.001, 0.999, length.out = n_points)
  curve_has_cols <- ncol(curve_cond) > 0L
  curve_list     <- vector("list", n_curve)
  summary_list   <- vector("list", n_curve)
  for (c_i in seq_len(n_curve)) {
    dp <- dprime_mat[, c_i]
    sr <- sdratio_mat[, c_i]
    # Anchor the criterion grid to the mean dprime so FA is ~evenly spaced; the
    # noise scale is 1, so FA = cdf(-dprime/2 - criterion) inverts to this.
    c_grid  <- -mean(dp) / 2 - qf(fa_grid)
    fa_mat  <- cdf(outer(-dp / 2, c_grid, "-"))
    hit_mat <- cdf(sweep(outer(dp / 2, c_grid, "-"), 1L, sr, "/"))

    fa_full  <- cbind(0, fa_mat, 1)
    hit_full <- cbind(0, hit_mat, 1)
    n_fa     <- ncol(fa_full)
    df <- data.frame(
      FA    = as.vector(t(fa_full)),
      Hit   = as.vector(t(hit_full)),
      .draw = rep(seq_len(n_draws), each = n_fa)
    )
    summ <- data.frame(
      FA        = c(0, colMeans(fa_mat), 1),
      Hit_mean  = c(0, colMeans(hit_mat), 1),
      Hit_lower = c(0, unname(apply(hit_mat, 2L, stats::quantile, probs[1L])), 1),
      Hit_upper = c(0, unname(apply(hit_mat, 2L, stats::quantile, probs[2L])), 1)
    )
    if (curve_has_cols) {
      cond_row <- curve_cond[c_i, , drop = FALSE]
      df   <- cbind(df,   cond_row[rep(1L, nrow(df)), , drop = FALSE],   row.names = NULL)
      summ <- cbind(summ, cond_row[rep(1L, nrow(summ)), , drop = FALSE], row.names = NULL)
    }
    curve_list[[c_i]]   <- df
    summary_list[[c_i]] <- summ
  }

  points <- if (length(dims$points) > 0L) {
    .roc_sdt_binary_points(fit, conditions, dims, has_sdratio, cdf, probs, ...)
  }

  list(curve   = do.call(rbind, curve_list),
       summary = do.call(rbind, summary_list),
       points  = points)
}


# Model-implied operating points: one (FA, Hit) per criterion level, predicted
# from the shared dprime/sdratio. FA = P("old" | noise), Hit = P("old" | signal).
.roc_sdt_binary_points <- function(fit, conditions, dims, has_sdratio,
                                   cdf, probs, ...) {
  point_cond <- .sdt_unique_subset(conditions, union(dims$curves, dims$points))
  n_pt <- max(1L, nrow(point_cond))

  dprime_mat  <- .sdt_linpred(fit, "dprime",    point_cond, is_rating = FALSE, ...)
  crit_mat    <- .sdt_linpred(fit, "criterion", point_cond, is_rating = FALSE, ...)
  sdratio_mat <- if (has_sdratio) {
    exp(.sdt_linpred(fit, "sdratio", point_cond, is_rating = FALSE, ...))
  } else {
    matrix(1, nrow = nrow(dprime_mat), ncol = n_pt)
  }

  rows <- vector("list", n_pt)
  for (c_i in seq_len(n_pt)) {
    fa  <- cdf(.sdt_eta(dprime_mat[, c_i], crit_mat[, c_i], 0L, sdratio_mat[, c_i]))
    hit <- cdf(.sdt_eta(dprime_mat[, c_i], crit_mat[, c_i], 1L, sdratio_mat[, c_i]))
    row <- data.frame(
      FA_mean   = mean(fa),
      FA_lower  = unname(stats::quantile(fa, probs[1L])),
      FA_upper  = unname(stats::quantile(fa, probs[2L])),
      Hit_mean  = mean(hit),
      Hit_lower = unname(stats::quantile(hit, probs[1L])),
      Hit_upper = unname(stats::quantile(hit, probs[2L]))
    )
    if (ncol(point_cond) > 0L) {
      row <- cbind(point_cond[c_i, , drop = FALSE], row, row.names = NULL)
    }
    rows[[c_i]] <- row
  }
  do.call(rbind, rows)
}


# Empirical (K+1-point) ROC for rating SDT models, one set of points per draw.
.roc_sdt_rating <- function(fit, model, conditions, ...) {
  dist           <- model$other_vars$dist
  threshold_type <- model$other_vars$threshold_type
  n_ratings      <- model$other_vars$n_ratings
  has_sdratio    <- .sdt_has_estimated_sdratio(model, fit)

  n_cond     <- max(1L, nrow(conditions))
  dprime_mat <- .sdt_linpred(fit, "dprime",    conditions, is_rating = TRUE, ...)
  crit_mat   <- .sdt_linpred(fit, "criterion", conditions, is_rating = TRUE, ...)
  n_draws    <- nrow(dprime_mat)

  spacing_mat <- if ("spacing" %in% names(model$parameters)) {
    .sdt_linpred(fit, "spacing", conditions, is_rating = TRUE, ...)
  }
  delta_names <- grep("^delta", names(model$parameters), value = TRUE)
  delta_mats  <- if (length(delta_names)) {
    stats::setNames(
      lapply(delta_names, function(nm) {
        .sdt_linpred(fit, nm, conditions, is_rating = TRUE, ...)
      }), delta_names
    )
  }
  sdratio_mat <- if (has_sdratio) {
    exp(.sdt_linpred(fit, "sdratio", conditions, is_rating = TRUE, ...))
  }

  cond_has_cols <- ncol(conditions) > 0L
  result <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    roc_rows <- vector("list", n_draws)
    for (d_i in seq_len(n_draws)) {
      dp   <- dprime_mat[d_i, c_i]
      sr   <- if (!is.null(sdratio_mat)) sdratio_mat[d_i, c_i] else 1
      sp   <- if (!is.null(spacing_mat)) spacing_mat[d_i, c_i] else NULL
      deltas_d <- if (!is.null(delta_mats)) {
        vapply(delta_mats, function(m) m[d_i, c_i], numeric(1))
      }

      thresholds <- .sdt_make_thresholds(crit_mat[d_i, c_i], n_ratings,
                                         threshold_type, spacing = sp,
                                         deltas = deltas_d)
      pn <- .sdt_category_probs(thresholds, dp, 1,  0L, dist)
      ps <- .sdt_category_probs(thresholds, dp, sr, 1L, dist)

      roc_rows[[d_i]] <- data.frame(
        FA    = c(1, 1 - cumsum(pn)[seq_len(n_ratings - 1L)], 0),
        Hit   = c(1, 1 - cumsum(ps)[seq_len(n_ratings - 1L)], 0),
        .draw = d_i
      )
    }
    cond_df <- do.call(rbind, roc_rows)
    if (cond_has_cols) {
      cond_df <- cbind(cond_df, conditions[rep(c_i, nrow(cond_df)), , drop = FALSE],
                       row.names = NULL)
    }
    result[[c_i]] <- cond_df
  }
  do.call(rbind, result)
}


# Summarise ROC draws to a mean + quantile band at each (FA, condition).
.roc_sdt_summary <- function(roc_data, probs) {
  cond_cols <- setdiff(names(roc_data), c("FA", "Hit", ".draw"))
  roc_data$.fa_r <- round(roc_data$FA, 8L)
  group_cols <- c(".fa_r", cond_cols)
  groups <- unique(roc_data[, group_cols, drop = FALSE])

  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    mask <- rep(TRUE, nrow(roc_data))
    for (col in group_cols) mask <- mask & (roc_data[[col]] == groups[i, col])
    hit <- roc_data$Hit[mask]
    rows[[i]] <- data.frame(
      FA        = roc_data$FA[mask][1L],
      Hit_mean  = mean(hit),
      Hit_lower = unname(stats::quantile(hit, probs[1L])),
      Hit_upper = unname(stats::quantile(hit, probs[2L])),
      groups[i, cond_cols, drop = FALSE],
      row.names = NULL, check.names = FALSE
    )
  }
  do.call(rbind, rows)
}


#' @export
print.bmm_sdt_roc <- function(x, ...) {
  model_name <- utils::tail(attr(x, "model_class"), 1L) %||% "sdt"
  n_draws    <- length(unique(x$.draw))
  n_cond     <- max(1L, nrow(attr(x, "conditions")))

  cat("SDT ROC curve (", model_name, ", dist = ", attr(x, "dist"), ")\n", sep = "")
  cat("  ", n_draws, " posterior draws",
      if (n_cond > 1L) paste0(" x ", n_cond, " conditions") else "", "\n", sep = "")
  if (isTRUE(attr(x, "is_rating"))) {
    cat("  Rating model: ", nrow(x) / (n_draws * n_cond) - 1L,
        " ROC points per draw\n", sep = "")
  } else {
    cat("  Smooth curve: ", nrow(x) / (n_draws * n_cond),
        " FA points per draw\n", sep = "")
    if (!is.null(attr(x, "points"))) {
      cat("  ", nrow(attr(x, "points")), " model-implied criterion points",
          " (see attr(x, 'points'))\n", sep = "")
    }
  }
  cat("Use plot() to visualise, or attr(x, 'summary') for summaries.\n")
  invisible(x)
}


############################################################################# !
# OBSERVED ROC                                                           ####
############################################################################# !

#' Observed ROC points from a fitted SDT model
#'
#' Computes empirical (observed) ROC points from the response-count data used to
#' fit a [sdt_rating()] or [sdt_binary()] model. Rating models pool counts
#' across observations within each stimulus type (and optional condition) to
#' produce cumulative hit/false-alarm rates at each confidence threshold; binary
#' models produce one operating point per criterion level (the empirical
#' counterpart of the model-implied points from [roc_sdt()]).
#'
#' @param fit A `bmmfit` object from a rating or binary SDT model.
#' @param conditions Optional character vector of column names to condition on.
#'   For rating models, `NULL` pools all observations into one ROC. For binary
#'   models, `NULL` auto-detects the criterion-varying predictor(s).
#'
#' @return A data frame of class `"bmm_sdt_roc_observed"` with columns `FA`,
#'   `Hit`, and any condition columns. Rating models additionally include the
#'   (0,0) and (1,1) endpoints.
#'
#' @seealso [roc_sdt()], [plot.bmm_sdt_roc()]
#' @export
roc_observed <- function(fit, conditions = NULL) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "roc_observed() is only available for SDT models (sdt_binary, sdt_rating)")
  if (inherits(model, "sdt_rating")) return(.roc_observed_rating(fit, model, conditions))
  if (inherits(model, "sdt_binary")) return(.roc_observed_binary(fit, model, conditions))
  stop2("roc_observed() requires a binary or rating SDT model with a response criterion.")
}


.roc_observed_rating <- function(fit, model, conditions = NULL) {
  data      <- fit$data
  stim_var  <- model$other_vars$stimulus
  n_ratings <- model$other_vars$n_ratings
  y_mat     <- data$Y

  cond_cols   <- conditions %||% character(0)
  cond_levels <- if (length(cond_cols) > 0L) {
    unique(data[, cond_cols, drop = FALSE])
  } else {
    data[1L, character(0), drop = FALSE]
  }
  n_cond <- max(1L, nrow(cond_levels))

  result <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    mask <- if (length(cond_cols) > 0L) {
      Reduce("&", lapply(cond_cols, function(col) data[[col]] == cond_levels[[col]][c_i]))
    } else {
      rep(TRUE, nrow(data))
    }
    stim_col <- data[[stim_var]]
    pn <- colSums(y_mat[mask & stim_col == 0, , drop = FALSE])
    ps <- colSums(y_mat[mask & stim_col == 1, , drop = FALSE])
    pn <- pn / sum(pn)
    ps <- ps / sum(ps)

    roc_df <- data.frame(
      FA  = c(1, 1 - cumsum(pn)[seq_len(n_ratings - 1L)], 0),
      Hit = c(1, 1 - cumsum(ps)[seq_len(n_ratings - 1L)], 0)
    )
    if (length(cond_cols) > 0L) {
      roc_df <- cbind(roc_df, cond_levels[rep(c_i, nrow(roc_df)), , drop = FALSE],
                      row.names = NULL)
    }
    result[[c_i]] <- roc_df
  }
  structure(do.call(rbind, result),
            class = c("bmm_sdt_roc_observed", "data.frame"),
            n_ratings = n_ratings, model_type = "rating")
}


.roc_observed_binary <- function(fit, model, conditions = NULL) {
  data     <- fit$data
  resp_var <- model$resp_vars$response
  stim_var <- model$other_vars$stimulus
  nt_var   <- model$other_vars$n_trials

  cond_cols <- intersect(conditions %||% .sdt_criterion_point_dims(fit)$points,
                         names(data))
  cond_levels <- if (length(cond_cols) > 0L) {
    unique(data[, cond_cols, drop = FALSE])
  } else {
    data[1L, character(0), drop = FALSE]
  }
  n_cond <- max(1L, nrow(cond_levels))

  rows <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    mask <- if (length(cond_cols) > 0L) {
      Reduce("&", lapply(cond_cols, function(col) data[[col]] == cond_levels[[col]][c_i]))
    } else {
      rep(TRUE, nrow(data))
    }
    stim_col <- data[[stim_var]]
    noise  <- mask & stim_col == 0
    signal <- mask & stim_col == 1
    row <- data.frame(
      FA  = sum(data[[resp_var]][noise])  / sum(data[[nt_var]][noise]),
      Hit = sum(data[[resp_var]][signal]) / sum(data[[nt_var]][signal])
    )
    if (length(cond_cols) > 0L) {
      row <- cbind(cond_levels[c_i, , drop = FALSE], row, row.names = NULL)
    }
    rows[[c_i]] <- row
  }
  structure(do.call(rbind, rows),
            class = c("bmm_sdt_roc_observed", "data.frame"),
            model_type = "binary")
}


############################################################################# !
# LATENT DECISION VARIABLE                                               ####
############################################################################# !

#' Model-implied latent decision-variable distributions
#'
#' Computes the noise and signal latent decision-variable densities implied by a
#' fitted [sdt_binary()] or [sdt_rating()] model, together with the response
#' criterion (binary) or the K-1 confidence thresholds (rating). This is the
#' canonical signal detection picture: two evidence distributions separated by
#' `dprime` and cut by one or more decision boundaries.
#'
#' Densities are evaluated at the posterior-mean parameters of each condition,
#' on the centred evidence axis the model uses internally (noise at `-dprime/2`
#' with unit SD, signal at `+dprime/2` with SD `exp(sdratio)`). The decision
#' boundaries additionally carry a posterior credible band on their location.
#' For the symmetric distributions (`"normal"`, `"logistic"`) the area beyond a
#' boundary equals the corresponding hit/false-alarm rate; for the Gumbel
#' distributions the densities and boundaries are faithful but that area
#' identity does not hold (the model applies the CDF to a signed distance, not
#' to a survivor probability).
#'
#' Like ROC curves, latent distributions require a response criterion and are
#' not defined for [sdt_mafc()] or [sdt_ranking()].
#'
#' @inheritParams roc_sdt
#' @param n_grid Integer. Number of points on the evidence-axis grid at which
#'   each density is evaluated (default 200).
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible band on the criterion/threshold locations (default
#'   `c(0.025, 0.975)`).
#'
#' @return A data frame of class `"bmm_sdt_latent"` with columns `x` (the
#'   evidence axis), `density`, `distribution` (`"noise"` or `"signal"`), and any
#'   condition columns. A `lines` attribute holds the decision-boundary
#'   positions (columns `position`, `lower`, `upper`, `marker`, plus condition
#'   columns).
#'
#' @seealso [roc_sdt()], [plot.bmm_sdt_latent()]
#' @export
latent_sdt <- function(fit, conditions = NULL, n_grid = 200,
                       probs = c(0.025, 0.975), ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "latent_sdt() is only available for SDT models (sdt_binary, sdt_rating)")
  stopif(inherits(model, "sdt_mafc"),
         "Latent distributions are not defined for the m-AFC SDT model: it has no response criterion.")
  stopif(inherits(model, "sdt_ranking"),
         "Latent distributions are not defined for the ranking SDT model: it has no response criterion.")

  dist        <- model$other_vars$dist
  pdf         <- .SDT_DISTS[[dist]]$pdf
  is_rating   <- inherits(model, "sdt_rating")
  has_sdratio <- .sdt_has_estimated_sdratio(model, fit)
  conditions  <- .sdt_resolve_conditions(fit, conditions)
  n_cond      <- max(1L, nrow(conditions))

  dprime_mat  <- .sdt_linpred(fit, "dprime",    conditions, is_rating, ...)
  crit_mat    <- .sdt_linpred(fit, "criterion", conditions, is_rating, ...)
  sdratio_mat <- if (has_sdratio) {
    exp(.sdt_linpred(fit, "sdratio", conditions, is_rating, ...))
  }

  if (is_rating) {
    n_ratings      <- model$other_vars$n_ratings
    threshold_type <- model$other_vars$threshold_type
    spacing_mat    <- if ("spacing" %in% names(model$parameters)) {
      .sdt_linpred(fit, "spacing", conditions, is_rating, ...)
    }
    delta_names <- grep("^delta", names(model$parameters), value = TRUE)
    delta_mats  <- if (length(delta_names)) {
      stats::setNames(lapply(delta_names, function(nm) {
        .sdt_linpred(fit, nm, conditions, is_rating, ...)
      }), delta_names)
    }
  }

  cond_has_cols <- ncol(conditions) > 0L
  dens_list  <- vector("list", n_cond)
  lines_list <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    dp <- mean(dprime_mat[, c_i])
    s  <- if (!is.null(sdratio_mat)) mean(sdratio_mat[, c_i]) else 1

    if (is_rating) {
      thr_draws <- t(vapply(seq_len(nrow(crit_mat)), function(d_i) {
        sp       <- if (!is.null(spacing_mat)) spacing_mat[d_i, c_i] else NULL
        deltas_d <- if (!is.null(delta_mats)) {
          vapply(delta_mats, function(m) m[d_i, c_i], numeric(1))
        }
        .sdt_make_thresholds(crit_mat[d_i, c_i], n_ratings, threshold_type,
                             spacing = sp, deltas = deltas_d)
      }, numeric(n_ratings - 1L)))
      positions <- colMeans(thr_draws)
      lower     <- unname(apply(thr_draws, 2L, stats::quantile, probs[1L]))
      upper     <- unname(apply(thr_draws, 2L, stats::quantile, probs[2L]))
      marker    <- paste0("t", seq_len(n_ratings - 1L))
    } else {
      crit_draws <- crit_mat[, c_i]
      positions  <- mean(crit_draws)
      lower      <- unname(stats::quantile(crit_draws, probs[1L]))
      upper      <- unname(stats::quantile(crit_draws, probs[2L]))
      marker     <- "criterion"
    }

    pad <- 4 * max(1, s)
    x   <- seq(min(-dp / 2, positions) - pad, max(dp / 2, positions) + pad,
               length.out = n_grid)
    dens <- data.frame(
      x            = rep(x, 2L),
      density      = c(pdf(x + dp / 2), pdf((x - dp / 2) / s) / s),
      distribution = rep(c("noise", "signal"), each = n_grid)
    )
    lines <- data.frame(position = positions, lower = lower, upper = upper,
                        marker = marker)
    if (cond_has_cols) {
      cond_row <- conditions[c_i, , drop = FALSE]
      dens  <- cbind(dens,  cond_row[rep(1L, nrow(dens)),  , drop = FALSE], row.names = NULL)
      lines <- cbind(lines, cond_row[rep(1L, nrow(lines)), , drop = FALSE], row.names = NULL)
    }
    dens_list[[c_i]]  <- dens
    lines_list[[c_i]] <- lines
  }

  structure(
    do.call(rbind, dens_list),
    class       = c("bmm_sdt_latent", "data.frame"),
    lines       = do.call(rbind, lines_list),
    probs       = probs,
    model_class = class(model),
    dist        = dist,
    is_rating   = is_rating,
    conditions  = conditions
  )
}


#' @export
print.bmm_sdt_latent <- function(x, ...) {
  model_name <- utils::tail(attr(x, "model_class"), 1L) %||% "sdt"
  n_cond     <- max(1L, nrow(attr(x, "conditions")))
  n_grid     <- nrow(x) / (2L * n_cond)
  n_lines    <- nrow(attr(x, "lines")) / n_cond

  cat("SDT latent distributions (", model_name, ", dist = ", attr(x, "dist"),
      ")\n", sep = "")
  cat("  noise + signal densities on ", n_grid, "-point grids",
      if (n_cond > 1L) paste0(" x ", n_cond, " conditions") else "", "\n",
      sep = "")
  cat("  ", n_lines, if (isTRUE(attr(x, "is_rating"))) " thresholds" else " criterion",
      " per condition\n", sep = "")
  cat("Use plot() to visualise, or attr(x, 'lines') for boundary locations.\n")
  invisible(x)
}


############################################################################# !
# AUC                                                                    ####
############################################################################# !

#' Area under the ROC curve from a fitted SDT model
#'
#' Computes the posterior area under the ROC curve (AUC). For Gaussian and
#' Gumbel-min equal-variance binary SDT the AUC is available in closed form from
#' the `dprime` draws; otherwise it is obtained by trapezoidal integration of
#' the ROC points from [roc_sdt()]. The returned AUC is always the area under
#' the full curve (for binary multi-criteria fits this is one value per curve,
#' not the trapezoid of the discrete operating points).
#'
#' @inheritParams roc_sdt
#' @param probs Numeric vector of length 2. Quantiles for the credible interval
#'   (default `c(0.025, 0.975)`).
#'
#' @return A data frame of class `"bmm_sdt_auc"` with columns `AUC`, `.draw`,
#'   and any condition columns, plus a `summary` attribute (`AUC_mean`,
#'   `AUC_lower`, `AUC_upper`).
#'
#' @details Analytical formulas: normal EV-SDT \eqn{AUC = \Phi(d'/\sqrt{2})};
#'   Gumbel-min EV-SDT \eqn{AUC = \mathrm{logistic}(g')}.
#'
#' @seealso [roc_sdt()], [plot.bmm_sdt_auc()]
#' @export
auc_sdt <- function(fit, conditions = NULL, probs = c(0.025, 0.975),
                    criterion_points = NULL, ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "auc_sdt() is only available for SDT models (sdt_binary, sdt_rating)")
  stopif(inherits(model, "sdt_mafc"),
         "AUC is not defined for the m-AFC SDT model: it has no response criterion.")
  stopif(inherits(model, "sdt_ranking"),
         "AUC is not defined for the ranking SDT model: it has no response criterion.")

  dist        <- model$other_vars$dist
  is_rating   <- inherits(model, "sdt_rating")
  has_sdratio <- .sdt_has_estimated_sdratio(model, fit)
  conditions  <- .sdt_resolve_conditions(fit, conditions)

  if (!is_rating) {
    dims <- .sdt_resolve_point_dims(fit, conditions, criterion_points)
    conditions <- .sdt_unique_subset(conditions, dims$curves)
  }

  use_analytical <- !is_rating && !has_sdratio && dist %in% c("normal", "gumbel_min")

  if (use_analytical) {
    auc_fn <- if (dist == "normal") function(d) stats::pnorm(d / sqrt(2)) else stats::plogis
    dprime_mat <- .sdt_linpred(fit, "dprime", conditions, is_rating = FALSE, ...)
    n_draws <- nrow(dprime_mat)
    result <- vector("list", ncol(dprime_mat))
    for (c_i in seq_len(ncol(dprime_mat))) {
      df <- data.frame(AUC = auc_fn(dprime_mat[, c_i]), .draw = seq_len(n_draws))
      if (ncol(conditions) > 0L) {
        df <- cbind(df, conditions[rep(c_i, n_draws), , drop = FALSE], row.names = NULL)
      }
      result[[c_i]] <- df
    }
    auc_data <- do.call(rbind, result)
  } else {
    auc_data <- .auc_sdt_numerical(fit, conditions, probs, criterion_points, ...)
  }

  structure(
    auc_data,
    class       = c("bmm_sdt_auc", "data.frame"),
    summary     = .auc_sdt_summary(auc_data, probs),
    model_class = class(model),
    dist        = dist,
    conditions  = conditions
  )
}


# Trapezoidal AUC over the full ROC curve, per draw x condition.
.auc_sdt_numerical <- function(fit, conditions, probs, criterion_points, ...) {
  roc_obj   <- roc_sdt(fit, conditions = conditions, probs = probs,
                       criterion_points = criterion_points, ...)
  cond_cols <- setdiff(names(roc_obj), c("FA", "Hit", ".draw"))
  group_key <- if (length(cond_cols) > 0L) {
    do.call(paste, c(roc_obj[c(".draw", cond_cols)], sep = "\r"))
  } else {
    as.character(roc_obj$.draw)
  }

  unique_keys <- unique(group_key)
  rows <- vector("list", length(unique_keys))
  for (k_i in seq_along(unique_keys)) {
    sub <- roc_obj[group_key == unique_keys[k_i], , drop = FALSE]
    sub <- sub[order(sub$FA), ]
    trap <- utils::head(sub$Hit, -1L) + utils::tail(sub$Hit, -1L)
    auc_val <- sum(diff(sub$FA) * trap) / 2
    row <- data.frame(AUC = auc_val, .draw = sub$.draw[1L])
    if (length(cond_cols) > 0L) {
      row <- cbind(row, sub[1L, cond_cols, drop = FALSE], row.names = NULL)
    }
    rows[[k_i]] <- row
  }
  do.call(rbind, rows)
}


.auc_sdt_summary <- function(auc_data, probs) {
  cond_cols <- setdiff(names(auc_data), c("AUC", ".draw"))
  if (length(cond_cols) == 0L) {
    return(data.frame(
      AUC_mean  = mean(auc_data$AUC),
      AUC_lower = unname(stats::quantile(auc_data$AUC, probs[1L])),
      AUC_upper = unname(stats::quantile(auc_data$AUC, probs[2L]))
    ))
  }
  groups <- unique(auc_data[, cond_cols, drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    mask <- rep(TRUE, nrow(auc_data))
    for (col in cond_cols) mask <- mask & (auc_data[[col]] == groups[i, col])
    auc <- auc_data$AUC[mask]
    rows[[i]] <- data.frame(
      AUC_mean  = mean(auc),
      AUC_lower = unname(stats::quantile(auc, probs[1L])),
      AUC_upper = unname(stats::quantile(auc, probs[2L])),
      groups[i, cond_cols, drop = FALSE],
      row.names = NULL, check.names = FALSE
    )
  }
  do.call(rbind, rows)
}


#' @export
print.bmm_sdt_auc <- function(x, ...) {
  model_name <- utils::tail(attr(x, "model_class"), 1L) %||% "sdt"
  cat("SDT AUC (", model_name, ", dist = ", attr(x, "dist"), ")\n", sep = "")
  print(attr(x, "summary"), digits = 3, row.names = FALSE)
  invisible(x)
}
