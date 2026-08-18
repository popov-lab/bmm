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
#' curves require a response criterion, so they are defined for [sdt_yn()]
#' and [sdt_rating()] only (not the criterion-free [sdt_mafc()] and
#' [sdt_ranking()]).
#'
#' For **binary** models the ROC is traced analytically from the posterior of
#' `d` (and `sdratio` for unequal-variance SDT) over a grid of criterion
#' values. When the criterion varies across conditions but `d`/`sdratio`
#' do not (e.g. a base-rate manipulation, as in [broeder_schuetz_2009_e3]), the
#' several criteria are operating points on a single curve: `roc_sdt()` returns
#' one smooth curve and attaches the model-implied points (one per criterion
#' level) as a `points` attribute.
#'
#' For **rating** models the K-1 confidence thresholds define K-1 empirical ROC
#' points per posterior draw (returned as the data frame). The smooth
#' model-implied curve is traced over a virtual cut from the posterior of
#' `d` (and `sdratio`) and attached as the `summary` attribute, with the K-1
#' thresholds attached as the `points` attribute (labelled `c1`..`c(K-1)`) so
#' they fall on the curve.
#'
#' @param fit A `bmmfit` object returned by [bmm()] from an SDT model.
#' @param conditions Optional data frame of predictor values for which to
#'   compute the ROC curve. Column names must match predictor variables used in
#'   the formula. If `NULL` (default), unique predictor combinations are derived
#'   from the data.
#' @param n_points Integer. Number of equally-spaced points on the smooth
#'   model-implied ROC curve (default 100); used for both binary and rating
#'   models. The rating data frame itself still holds K+1 points per draw (K-1
#'   threshold points plus the (0,0) and (1,1) endpoints).
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible band (default `c(0.025, 0.975)`).
#' @param criterion_points Optional. Control of the binary multi-criteria
#'   behaviour. `NULL` (default) auto-detects predictors that vary the criterion
#'   only (not `d`/`sdratio`) and treats their levels as operating points
#'   on one curve. Pass a character vector of column names to force that
#'   classification, or `FALSE` to disable it (one separate curve per predictor
#'   combination). Ignored for rating models.
#' @param ... Additional arguments passed to [brms::posterior_linpred()].
#'
#' @return A data frame of class `"bmm_sdt_roc"` with columns `FA`, `Hit`,
#'   `.draw`, and any condition columns. The object carries a `summary`
#'   attribute (`FA`, `Hit_mean`, `Hit_lower`, `Hit_upper`) with the smooth
#'   model-implied curve, and a `points` attribute with the model-implied
#'   operating points: one per criterion level for binary multi-criteria fits,
#'   or the K-1 confidence thresholds (labelled `c1`..`c(K-1)`) for rating fits.
#'
#' @seealso [auc_sdt()], [roc_observed()], [plot.bmm_sdt_roc()]
#' @export
roc_sdt <- function(fit, conditions = NULL, n_points = 100,
                    probs = c(0.025, 0.975), criterion_points = NULL, ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "roc_sdt() is only available for SDT models")
  stopif(inherits(model, "sdt_mafc"),
         "ROC curves are not defined for the m-AFC SDT model: it has no response criterion.")
  stopif(inherits(model, "sdt_ranking"),
         "ROC curves are not defined for the ranking SDT model: it has no response criterion.")

  is_rating  <- inherits(model, "sdt_rating")
  conditions <- .sdt_resolve_conditions(fit, conditions)

  if (is_rating) {
    rating      <- .roc_sdt_rating(fit, model, conditions, n_points,
                                   probs = probs, ...)
    roc_data    <- rating$curve
    points      <- rating$points
    roc_summary <- rating$summary
  } else {
    binary      <- .roc_sdt_yn(fit, model, conditions, n_points,
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
    "Y", "nTrials", "dist_type", "m_afc", "max_rank",
    names(brms::ranef(fit))
  ))
  pred_cols <- setdiff(names(data), exclude)

  is_matrix_col <- vapply(pred_cols, function(col) is.matrix(data[[col]]), logical(1))
  .sdt_unique_subset(data, pred_cols[!is_matrix_col])
}


# Formula predictors per parameter, with random-effect grouping factors
# stripped so (1 | id) is never a dimension.
.sdt_stripped_preds <- function(fit) {
  uf <- fit$bmm$user_formula
  preds <- if (inherits(uf, "bmmformula")) rhs_vars(uf, collapse = FALSE) else list()
  re_vars <- tryCatch(names(brms::ranef(fit)), error = function(e) character(0))
  lapply(preds, function(v) setdiff(v %||% character(0), re_vars))
}


# Classify the formula predictors of a binary fit into criterion "operating
# points" (predict the criterion only) and "curves" (predict d/sdratio).
.sdt_criterion_point_dims <- function(fit) {
  preds <- .sdt_stripped_preds(fit)
  curve_vars <- union(preds[["d"]], preds[["sdratio"]]) %||% character(0)
  list(points = setdiff(preds[["criterion"]] %||% character(0), curve_vars),
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


# Recycle a 1-row condition frame onto every row of `df` (no-op for NULL or
# zero-column conditions).
.sdt_bind_cond <- function(df, cond_row) {
  if (is.null(cond_row) || ncol(cond_row) == 0L) return(df)
  cbind(df, cond_row[rep(1L, nrow(df)), , drop = FALSE], row.names = NULL)
}


# Posterior mean and interval per column of a draws matrix (or one draws
# vector), as a data frame with (<prefix>_)mean/lower/upper columns. Pins the
# credible-interval convention every summary in this file uses.
.sdt_summarise_draws <- function(draws, probs, prefix = NULL) {
  draws <- cbind(draws)
  out <- data.frame(
    mean  = colMeans(draws),
    lower = unname(apply(draws, 2L, stats::quantile, probs[1L])),
    upper = unname(apply(draws, 2L, stats::quantile, probs[2L]))
  )
  rownames(out) <- NULL
  if (!is.null(prefix)) names(out) <- paste0(prefix, "_", names(out))
  out
}


# Split rows of `data` by the unique combinations of `cond_cols`: the unique
# level frame plus one logical row mask per level. A zero-column request yields
# a single all-TRUE mask so callers loop exactly once.
.sdt_condition_masks <- function(data, cond_cols) {
  if (length(cond_cols) == 0L) {
    return(list(levels = data[1L, character(0), drop = FALSE],
                masks  = list(rep(TRUE, nrow(data)))))
  }
  levels <- unique(data[, cond_cols, drop = FALSE])
  rownames(levels) <- NULL
  masks <- lapply(seq_len(nrow(levels)), function(i) {
    Reduce(`&`, lapply(cond_cols, function(col) data[[col]] == levels[[col]][i]))
  })
  list(levels = levels, masks = masks)
}


# Classify the predictor columns present in `conditions` for latent_sdt(). The
# latent densities depend only on d/sdratio, so `density` dims (predictors
# of either) get distinct density panels (faceted), while `boundary` dims
# (predictors of the criterion/threshold parameters only) leave the densities
# unchanged and are collapsed into one panel with the boundaries overlaid.
# Columns predicting neither (e.g. a set-size column for mafc/ranking) are
# dropped. `collapse = FALSE` facets every dimension; `collapse = <cols>` forces
# the named columns to collapse instead of facet.
.sdt_latent_dims <- function(fit, conditions, collapse = NULL) {
  all_vars <- names(conditions)
  if (length(all_vars) == 0L) {
    return(list(density = character(0), boundary = character(0)))
  }
  if (isFALSE(collapse)) {
    return(list(density = all_vars, boundary = character(0)))
  }

  preds <- .sdt_stripped_preds(fit)

  density_par   <- intersect(c("d", "sdratio"), names(preds))
  density_vars  <- unique(unlist(preds[density_par]))
  boundary_vars <- setdiff(
    unique(unlist(preds[setdiff(names(preds), density_par)])),
    density_vars
  )

  if (is.character(collapse)) {
    forced <- intersect(collapse, all_vars)
    return(list(density  = setdiff(all_vars, forced),
                boundary = intersect(forced, boundary_vars)))
  }
  list(density  = intersect(all_vars, density_vars),
       boundary = intersect(all_vars, boundary_vars))
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
# row so every formula variable is present. The parameters are nlpars for the
# multinomial models (rating, ranking) and dpars for the custom-family models
# (binary, mafc); which applies is derived from the model class so callers
# never need to know.
.sdt_linpred <- function(fit, param, conditions, ...) {
  newdata <- fit$data[rep(1L, max(1L, nrow(conditions))), , drop = FALSE]
  rownames(newdata) <- NULL
  if (ncol(conditions) > 0L) {
    cols <- intersect(names(conditions), names(newdata))
    newdata[cols] <- conditions[cols]
  }
  if (inherits(fit$bmm$model, c("sdt_rating", "sdt_ranking"))) {
    brms::posterior_linpred(fit, nlpar = param, newdata = newdata,
                            re_formula = NA, allow_new_levels = TRUE, ...)
  } else {
    brms::posterior_linpred(fit, dpar = param, newdata = newdata,
                            re_formula = NA, allow_new_levels = TRUE, ...)
  }
}


# Posterior draws of the latent geometry per condition: `d` is d_a, `sdratio`
# the natural-scale SD ratio (1 everywhere when fixed), and `sep` their product
# with the root-mean-square scale -- the separation between the distributions
# in noise-SD units, which equals d_a itself under equal variance. Every
# analysis function builds on these three matrices (draws x conditions).
.sdt_latent_geometry <- function(fit, conditions, has_sdratio, ...) {
  d <- .sdt_linpred(fit, "d", conditions, ...)
  sdratio <- if (has_sdratio) {
    exp(.sdt_linpred(fit, "sdratio", conditions, ...))
  } else {
    matrix(1, nrow = nrow(d), ncol = ncol(d))
  }
  list(d = d, sdratio = sdratio, sep = d * .sdt_rms_scale(sdratio))
}


# The K-1 rating thresholds are not posterior parameters: they are rebuilt per
# draw from the criterion plus the parameterization parameters (spacing and/or
# delta*), whose mapping to ordered thresholds differs by threshold_type.
# Keeping that reconstruction here lets roc_sdt(), latent_sdt(), and
# sdt_thresholds() share it. .sdt_make_thresholds is vectorized over draws, so
# one call per condition returns the full n_draws x (K-1) matrix.
.sdt_rating_thresholds <- function(fit, model, conditions, ...) {
  n_ratings      <- model$other_vars$n_ratings
  threshold_type <- model$other_vars$threshold_type
  crit    <- .sdt_linpred(fit, "criterion", conditions, ...)
  spacing <- if ("spacing" %in% names(model$parameters)) {
    .sdt_linpred(fit, "spacing", conditions, ...)
  }
  delta_names <- grep("^delta", names(model$parameters), value = TRUE)
  delta_mats  <- lapply(stats::setNames(delta_names, delta_names),
                        function(nm) .sdt_linpred(fit, nm, conditions, ...))

  lapply(seq_len(ncol(crit)), function(l_i) {
    deltas <- if (length(delta_mats)) {
      do.call(cbind, lapply(delta_mats, function(m) m[, l_i]))
    }
    rbind(.sdt_make_thresholds(crit[, l_i], n_ratings, threshold_type,
                               spacing = if (!is.null(spacing)) spacing[, l_i],
                               deltas = deltas))
  })
}


# Posterior-mean ROC + quantile band at each swept-cut node, with the (0,0) and
# (1,1) endpoints appended and any condition columns recycled in. Shared by the
# binary and rating smooth implied curves so a model's K-1 thresholds (rating) or
# criterion operating points (binary) fall on the displayed curve.
.roc_summary_from_mats <- function(fa_mat, hit_mat, probs, cond_row = NULL) {
  hit <- .sdt_summarise_draws(hit_mat, probs, prefix = "Hit")
  summ <- data.frame(
    FA        = c(0, colMeans(fa_mat), 1),
    Hit_mean  = c(0, hit$Hit_mean, 1),
    Hit_lower = c(0, hit$Hit_lower, 1),
    Hit_upper = c(0, hit$Hit_upper, 1)
  )
  .sdt_bind_cond(summ, cond_row)
}


# Smooth analytical ROC for binary EV- and UV-SDT plus, when the criterion
# varies across conditions only, the model-implied operating points. The curve
# is traced over a shared criterion grid and summarised as the posterior-mean
# operating point (FA, Hit) at each node -- the same functional as the
# model-implied points -- so those points fall on the curve. Averaging Hit at a
# fixed FA instead would put them slightly off it on the concave ROC.
.roc_sdt_yn <- function(fit, model, conditions, n_points,
                            criterion_points = NULL,
                            probs = c(0.025, 0.975), ...) {
  dist        <- model$other_vars$dist
  has_sdratio <- .sdt_has_estimated_sdratio(model, fit)
  cdf         <- .sdt_dists[[dist]]$cdf
  qf          <- .sdt_dists[[dist]]$qf

  dims       <- .sdt_resolve_point_dims(fit, conditions, criterion_points)
  curve_cond <- .sdt_unique_subset(conditions, dims$curves)

  n_curve <- max(1L, nrow(curve_cond))
  geom    <- .sdt_latent_geometry(fit, curve_cond, has_sdratio, ...)
  n_draws <- nrow(geom$d)

  fa_grid        <- seq(0.001, 0.999, length.out = n_points)
  curve_has_cols <- ncol(curve_cond) > 0L
  curve_list     <- vector("list", n_curve)
  summary_list   <- vector("list", n_curve)
  for (c_i in seq_len(n_curve)) {
    dp <- geom$sep[, c_i]
    sr <- geom$sdratio[, c_i]
    cond_row <- if (curve_has_cols) curve_cond[c_i, , drop = FALSE]
    # Anchor the criterion grid to the mean separation so FA is ~evenly spaced;
    # the noise scale is 1, so FA = 1 - cdf(sep/2 + criterion) inverts to this.
    c_grid  <- -mean(dp) / 2 + qf(1 - fa_grid)
    # P("old") is the survival function of the evidence distribution, matching
    # the likelihood; cdf(eta) would trace the mirror ROC for the asymmetric
    # extreme-value distributions.
    fa_mat  <- 1 - cdf(-outer(-dp / 2, c_grid, "-"))
    hit_mat <- 1 - cdf(-sweep(outer(dp / 2, c_grid, "-"), 1L, sr, "/"))

    fa_full  <- cbind(0, fa_mat, 1)
    hit_full <- cbind(0, hit_mat, 1)
    n_fa     <- ncol(fa_full)
    df <- data.frame(
      FA    = as.vector(t(fa_full)),
      Hit   = as.vector(t(hit_full)),
      .draw = rep(seq_len(n_draws), each = n_fa)
    )
    curve_list[[c_i]]   <- .sdt_bind_cond(df, cond_row)
    summary_list[[c_i]] <- .roc_summary_from_mats(fa_mat, hit_mat, probs, cond_row)
  }

  points <- if (length(dims$points) > 0L) {
    .roc_sdt_yn_points(fit, conditions, dims, has_sdratio, cdf, probs, ...)
  }

  list(curve   = do.call(rbind, curve_list),
       summary = do.call(rbind, summary_list),
       points  = points)
}


# Model-implied operating points: one (FA, Hit) per criterion level, predicted
# from the shared d/sdratio. FA = P("old" | noise), Hit = P("old" | signal).
.roc_sdt_yn_points <- function(fit, conditions, dims, has_sdratio,
                                   cdf, probs, ...) {
  point_cond <- .sdt_unique_subset(conditions, union(dims$curves, dims$points))
  n_pt <- max(1L, nrow(point_cond))

  geom     <- .sdt_latent_geometry(fit, point_cond, has_sdratio, ...)
  crit_mat <- .sdt_linpred(fit, "criterion", point_cond, ...)

  rows <- vector("list", n_pt)
  for (c_i in seq_len(n_pt)) {
    fa  <- 1 - cdf(-.sdt_eta(geom$d[, c_i], crit_mat[, c_i], 0L,
                             geom$sdratio[, c_i]))
    hit <- 1 - cdf(-.sdt_eta(geom$d[, c_i], crit_mat[, c_i], 1L,
                             geom$sdratio[, c_i]))
    row <- cbind(.sdt_summarise_draws(fa, probs, prefix = "FA"),
                 .sdt_summarise_draws(hit, probs, prefix = "Hit"))
    if (ncol(point_cond) > 0L) {
      row <- cbind(point_cond[c_i, , drop = FALSE], row, row.names = NULL)
    }
    rows[[c_i]] <- row
  }
  do.call(rbind, rows)
}


# ROC for rating SDT models. Returns three pieces (like .roc_sdt_yn): the
# discrete K+1-point ROC per draw (`curve`, also used by the numerical AUC), the
# smooth model-implied curve swept over a virtual cut (`summary`), and the K-1
# threshold operating points with a credible band (`points`, labelled c1..cK-1).
# The smooth curve uses the rating model's own probability map -- FA = 1 - cdf(t
# + sep/2), Hit = 1 - cdf((t - sep/2) / sdratio), where sep is the separation in
# noise-SD units -- the continuous envelope of the discrete points, so the
# thresholds fall on it for every distribution (the binary curve's criterion
# sign convention would mismatch for the asymmetric Gumbel distributions).
.roc_sdt_rating <- function(fit, model, conditions, n_points,
                            probs = c(0.025, 0.975), ...) {
  dist        <- model$other_vars$dist
  n_ratings   <- model$other_vars$n_ratings
  has_sdratio <- .sdt_has_estimated_sdratio(model, fit)
  cdf         <- .sdt_dists[[dist]]$cdf
  qf          <- .sdt_dists[[dist]]$qf
  K1          <- n_ratings - 1L

  n_cond   <- max(1L, nrow(conditions))
  geom     <- .sdt_latent_geometry(fit, conditions, has_sdratio, ...)
  n_draws  <- nrow(geom$d)
  thr_list <- .sdt_rating_thresholds(fit, model, conditions, ...)

  fa_grid       <- seq(0.001, 0.999, length.out = n_points)
  thr_levels    <- paste0("c", seq_len(K1))
  cond_has_cols <- ncol(conditions) > 0L
  curve_list   <- vector("list", n_cond)
  summary_list <- vector("list", n_cond)
  points_list  <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    d_vec   <- geom$d[, c_i]
    sep_vec <- geom$sep[, c_i]
    sr_vec  <- geom$sdratio[, c_i]

    # sdratio goes into both calls: .sdt_category_probs derives the noise/signal
    # scale from the stimulus flag, but it also needs sdratio to convert d_a into
    # the separation. Passing 1 for the noise row would silently drop that
    # conversion and place the noise density at -d_a/2 instead of -sep/2.
    pn <- .sdt_category_probs(thr_list[[c_i]], d_vec, sr_vec, 0L, dist)
    ps <- .sdt_category_probs(thr_list[[c_i]], d_vec, sr_vec, 1L, dist)
    fa_pts  <- 1 - matrixStats::rowCumsums(pn)[, seq_len(K1), drop = FALSE]
    hit_pts <- 1 - matrixStats::rowCumsums(ps)[, seq_len(K1), drop = FALSE]

    cond_row <- if (cond_has_cols) conditions[c_i, , drop = FALSE]
    cond_df <- data.frame(
      FA    = as.vector(t(cbind(1, fa_pts, 0))),
      Hit   = as.vector(t(cbind(1, hit_pts, 0))),
      .draw = rep(seq_len(n_draws), each = K1 + 2L)
    )

    t_grid  <- qf(1 - fa_grid) - mean(sep_vec) / 2
    fa_mat  <- 1 - cdf(outer(sep_vec / 2, t_grid, "+"))
    hit_mat <- 1 - cdf(sweep(outer(-sep_vec / 2, t_grid, "+"), 1L, sr_vec, "/"))

    pts <- data.frame(
      threshold = factor(thr_levels, levels = thr_levels),
      .sdt_summarise_draws(fa_pts, probs, prefix = "FA"),
      .sdt_summarise_draws(hit_pts, probs, prefix = "Hit")
    )
    curve_list[[c_i]]   <- .sdt_bind_cond(cond_df, cond_row)
    summary_list[[c_i]] <- .roc_summary_from_mats(fa_mat, hit_mat, probs, cond_row)
    points_list[[c_i]]  <- .sdt_bind_cond(pts, cond_row)
  }

  list(curve   = do.call(rbind, curve_list),
       summary = do.call(rbind, summary_list),
       points  = do.call(rbind, points_list))
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
    cat("  ", nrow(attr(x, "points")) / n_cond, " threshold operating points",
        " (see attr(x, 'points'))\n", sep = "")
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
#' fit a [sdt_rating()] or [sdt_yn()] model. Rating models pool counts
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
         "roc_observed() is only available for SDT models")
  if (inherits(model, "sdt_rating")) return(.roc_observed_rating(fit, model, conditions))
  if (inherits(model, "sdt_yn")) return(.roc_observed_yn(fit, model, conditions))
  stop2("roc_observed() requires a binary or rating SDT model with a response criterion.")
}


.roc_observed_rating <- function(fit, model, conditions = NULL) {
  data      <- fit$data
  stim_var  <- model$other_vars$stimulus
  n_ratings <- model$other_vars$n_ratings
  y_mat     <- data$Y

  split  <- .sdt_condition_masks(data, conditions %||% character(0))
  n_cond <- length(split$masks)

  stim_col <- data[[stim_var]]
  result <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    mask <- split$masks[[c_i]]
    pn <- colSums(y_mat[mask & stim_col == 0, , drop = FALSE])
    ps <- colSums(y_mat[mask & stim_col == 1, , drop = FALSE])
    pn <- pn / sum(pn)
    ps <- ps / sum(ps)

    roc_df <- data.frame(
      FA  = c(1, 1 - cumsum(pn)[seq_len(n_ratings - 1L)], 0),
      Hit = c(1, 1 - cumsum(ps)[seq_len(n_ratings - 1L)], 0)
    )
    result[[c_i]] <- .sdt_bind_cond(roc_df, split$levels[c_i, , drop = FALSE])
  }
  structure(do.call(rbind, result),
            class = c("bmm_sdt_roc_observed", "data.frame"),
            n_ratings = n_ratings, model_type = "rating")
}


.roc_observed_yn <- function(fit, model, conditions = NULL) {
  data     <- fit$data
  resp_var <- model$resp_vars$response
  stim_var <- model$other_vars$stimulus
  nt_var   <- model$other_vars$n_trials

  cond_cols <- intersect(conditions %||% .sdt_criterion_point_dims(fit)$points,
                         names(data))
  split  <- .sdt_condition_masks(data, cond_cols)
  n_cond <- length(split$masks)

  stim_col <- data[[stim_var]]
  rows <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    mask   <- split$masks[[c_i]]
    noise  <- mask & stim_col == 0
    signal <- mask & stim_col == 1
    row <- data.frame(
      FA  = sum(data[[resp_var]][noise])  / sum(data[[nt_var]][noise]),
      Hit = sum(data[[resp_var]][signal]) / sum(data[[nt_var]][signal])
    )
    if (length(cond_cols) > 0L) {
      row <- cbind(split$levels[c_i, , drop = FALSE], row, row.names = NULL)
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
#' fitted [sdt_yn()] or [sdt_rating()] model, together with the response
#' criterion (binary) or the K-1 confidence thresholds (rating). This is the
#' canonical signal detection picture: two evidence distributions separated by
#' the sensitivity and cut by one or more decision boundaries.
#'
#' Densities are evaluated at the posterior-mean parameters of each condition,
#' on the centred evidence axis the model uses internally: noise at `-sep/2`
#' with unit SD and signal at `+sep/2` with SD `exp(sdratio)`, where
#' `sep = d * sqrt((1 + exp(sdratio)^2) / 2)` is the separation in noise-SD
#' units. `d` is \eqn{d_a}, measured in root-mean-square SD units, so `sep`
#' equals `d` whenever `sdratio` is fixed at 0. The decision
#' boundaries additionally carry a posterior credible band on their location.
#' The area beyond a boundary equals the corresponding hit/false-alarm rate for
#' every noise distribution: the likelihood evaluates the survivor function of
#' the evidence distribution at the boundary, which is exactly that area.
#'
#' Available for all four SDT models. [sdt_yn()] draws the criterion and
#' [sdt_rating()] the K-1 confidence thresholds as boundary lines. [sdt_mafc()]
#' and [sdt_ranking()] have no response criterion -- the decision is a max/rank
#' rule over the `m` alternatives -- so only the densities are drawn (no boundary
#' lines); there the noise density represents each of the `m - 1` distractor
#' alternatives and the signal density the target.
#'
#' Because the densities depend only on `d`/`sdratio`, predictors that vary
#' the criterion/thresholds *only* (e.g. a base-rate manipulation, as in
#' [broeder_schuetz_2009_e3]) leave the densities unchanged: by default they are
#' collapsed into a single panel with their several boundaries overlaid and
#' colour-coded, rather than shown as repeated identical panels. Predictors that
#' vary `d`/`sdratio` produce distinct density panels (faceted). For
#' `sdt_mafc`/`sdt_ranking` the set size is likewise collapsed unless it predicts
#' `d`. Use `collapse` to override the auto-detection.
#'
#' @inheritParams roc_sdt
#' @param n_grid Integer. Number of points on the evidence-axis grid at which
#'   each density is evaluated (default 200).
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible band on the criterion/threshold locations (default
#'   `c(0.025, 0.975)`).
#' @param collapse Control of which predictor dimensions are collapsed into one
#'   panel. `NULL` (default) auto-detects: predictors of the criterion/thresholds
#'   only (and the set size) are collapsed; predictors of `d`/`sdratio` are
#'   faceted. Pass a character vector of column names to force those columns to
#'   collapse, or `FALSE` to facet every dimension (one panel per combination).
#' @param show_competitors Logical (default `FALSE`). For `sdt_mafc`/`sdt_ranking`
#'   only, additionally overlay the density of the maximum of the `m - 1`
#'   distractor samples (one curve per set size) -- the "effective competitor"
#'   the target must beat -- which shifts rightward as `m` grows and visualises
#'   why accuracy falls with set size. Ignored for `sdt_yn`/`sdt_rating`.
#'
#' @return A data frame of class `"bmm_sdt_latent"` with columns `x` (the
#'   evidence axis), `density`, `distribution` (`"noise"` or `"signal"`), and any
#'   faceting (density) condition columns. A `lines` attribute holds the
#'   decision-boundary positions (columns `position`, `lower`, `upper`, `marker`,
#'   `level`, plus condition columns), or `NULL` for the criterion-free
#'   `sdt_mafc`/`sdt_ranking` models. When `show_competitors = TRUE`, a
#'   `competitors` attribute holds the max-of-distractors densities.
#'
#' @seealso [roc_sdt()], [plot.bmm_sdt_latent()]
#' @export
latent_sdt <- function(fit, conditions = NULL, n_grid = 200,
                       probs = c(0.025, 0.975), collapse = NULL,
                       show_competitors = FALSE, ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "latent_sdt() is only available for SDT models")

  dist          <- model$other_vars$dist
  pdf           <- .sdt_dists[[dist]]$pdf
  cdf           <- .sdt_dists[[dist]]$cdf
  is_rating     <- inherits(model, "sdt_rating")
  # criterion exists only for the models with an explicit response boundary
  # (binary, rating)
  has_criterion <- "criterion" %in% names(model$parameters)
  has_sdratio   <- .sdt_has_estimated_sdratio(model, fit)
  competitors   <- isTRUE(show_competitors) &&
    inherits(model, c("sdt_mafc", "sdt_ranking"))

  conditions <- .sdt_resolve_conditions(fit, conditions)
  dims       <- .sdt_latent_dims(fit, conditions, collapse)
  panel_cond <- .sdt_unique_subset(conditions, dims$density)
  line_cond  <- .sdt_unique_subset(conditions, c(dims$density, dims$boundary))
  n_panel    <- max(1L, nrow(panel_cond))

  # The densities live on the noise-standardized axis, so the signal mean is
  # half the separation (geom$sep), not half of d_a. The two coincide under
  # equal variance.
  geom <- .sdt_latent_geometry(fit, panel_cond, has_sdratio, ...)

  # Boundary positions are evaluated over the density x boundary combinations so
  # that each density panel can carry all its criteria/thresholds.
  if (has_criterion) {
    crit_line <- .sdt_linpred(fit, "criterion", line_cond, ...)
    thr_list  <- if (is_rating) {
      .sdt_rating_thresholds(fit, model, line_cond, ...)
    }
    lines_all <- vector("list", nrow(line_cond))
    for (l_i in seq_len(nrow(line_cond))) {
      boundary_draws <- if (is_rating) thr_list[[l_i]] else crit_line[, l_i]
      lines <- stats::setNames(.sdt_summarise_draws(boundary_draws, probs),
                               c("position", "lower", "upper"))
      lines$marker <- if (is_rating) {
        paste0("t", seq_len(nrow(lines)))
      } else {
        "criterion"
      }
      lines_all[[l_i]] <- .sdt_bind_cond(lines, line_cond[l_i, , drop = FALSE])
    }
    lines_all <- do.call(rbind, lines_all)
    # Colour key: the boundary-dim value(s) when the criterion varies across a
    # predictor (e.g. base-rate conditions), else the marker so a rating model's
    # K-1 thresholds stay distinguishable.
    lines_all$level <- if (length(dims$boundary) > 0L) {
      factor(do.call(paste, c(lines_all[dims$boundary], sep = " / ")))
    } else {
      factor(lines_all$marker, levels = unique(lines_all$marker))
    }
  }

  m_col     <- if (competitors && is.character(model$other_vars$m)) model$other_vars$m
  set_sizes <- if (competitors) .sdt_competitor_sizes(model, fit$data)

  dens_list  <- vector("list", n_panel)
  lines_list <- vector("list", n_panel)
  comp_list  <- vector("list", n_panel)
  for (p_i in seq_len(n_panel)) {
    dp <- mean(geom$sep[, p_i])
    s  <- mean(geom$sdratio[, p_i])
    panel_row <- panel_cond[p_i, , drop = FALSE]

    # Each density panel carries only the boundaries sharing its density-dim
    # values, so restrict the pooled boundary frame to this panel's combination.
    lines <- if (has_criterion) {
      dcols <- intersect(dims$density, names(lines_all))
      if (length(dcols) == 0L) {
        lines_all
      } else {
        keep <- Reduce(`&`, lapply(dcols, function(col) lines_all[[col]] == panel_row[[col]]))
        lines_all[keep, , drop = FALSE]
      }
    }

    positions <- if (!is.null(lines)) lines$position else numeric(0)
    pad <- 4 * max(1, s)
    x   <- seq(min(-dp / 2, positions) - pad, max(dp / 2, positions) + pad,
               length.out = n_grid)
    dens <- data.frame(
      x            = rep(x, 2L),
      density      = c(pdf(x + dp / 2), pdf((x - dp / 2) / s) / s),
      distribution = rep(c("noise", "signal"), each = n_grid)
    )

    comp <- if (competitors) {
      ms <- if (!is.null(m_col) && m_col %in% dims$density) {
        as.integer(panel_cond[p_i, m_col])
      } else {
        set_sizes
      }
      do.call(rbind, lapply(ms, function(mv) {
        data.frame(
          x        = x,
          density  = (mv - 1) * cdf(x + dp / 2)^(mv - 2) * pdf(x + dp / 2),
          set_size = factor(mv, levels = set_sizes)
        )
      }))
    }

    dens_list[[p_i]]  <- .sdt_bind_cond(dens, panel_row)
    lines_list[[p_i]] <- lines
    comp_list[[p_i]]  <- if (!is.null(comp)) .sdt_bind_cond(comp, panel_row)
  }

  structure(
    do.call(rbind, dens_list),
    class        = c("bmm_sdt_latent", "data.frame"),
    lines        = do.call(rbind, lines_list),
    competitors  = if (competitors) do.call(rbind, comp_list),
    probs        = probs,
    model_class  = class(model),
    dist         = dist,
    is_rating    = is_rating,
    conditions   = panel_cond
  )
}


# Distinct set sizes for the competitor overlay. brms keeps only the derived
# per-row set-size column (m_afc / max_rank) in fit$data, not the user's original
# m column, so read that; fall back to a constant m for un-fitted/mock objects.
.sdt_competitor_sizes <- function(model, data) {
  col <- if (inherits(model, "sdt_ranking")) "max_rank" else "m_afc"
  m   <- model$other_vars$m
  if (!is.null(data[[col]])) {
    sort(unique(as.integer(data[[col]])))
  } else if (is.numeric(m)) {
    as.integer(m)
  } else if (m %in% names(data)) {
    sort(unique(as.integer(data[[m]])))
  } else {
    integer(0)
  }
}


#' @export
print.bmm_sdt_latent <- function(x, ...) {
  model_name <- utils::tail(attr(x, "model_class"), 1L) %||% "sdt"
  n_panel    <- max(1L, nrow(attr(x, "conditions")))
  n_grid     <- nrow(x) / (2L * n_panel)
  lines      <- attr(x, "lines")

  cat("SDT latent distributions (", model_name, ", dist = ", attr(x, "dist"),
      ")\n", sep = "")
  cat("  noise + signal densities on ", n_grid, "-point grids",
      if (n_panel > 1L) paste0(" x ", n_panel, " panels") else "", "\n",
      sep = "")
  if (is.null(lines)) {
    cat("  no response boundary (max/rank decision rule)\n")
  } else {
    cat("  ", nrow(lines) / n_panel,
        if (isTRUE(attr(x, "is_rating"))) " thresholds" else " criteria",
        " per panel\n", sep = "")
  }
  if (!is.null(attr(x, "competitors"))) {
    cat("  max-of-distractors densities overlaid\n")
  }
  cat("Use plot() to visualise, or attr(x, 'lines') for boundary locations.\n")
  invisible(x)
}


############################################################################# !
# THRESHOLDS                                                             ####
############################################################################# !

#' Latent decision thresholds from a fitted rating SDT model
#'
#' Extracts the \eqn{K-1} latent decision thresholds of a rating signal
#' detection model fit with [bmm()]. Thresholds are not posterior parameters:
#' they are reconstructed per draw from the `criterion` and the parameterization
#' parameters (`spacing` and/or `delta*`), whose mapping to ordered thresholds
#' depends on the model's `threshold_type`. This function returns those draws on
#' the latent decision-variable scale together with a posterior summary, so
#' threshold estimates are accessible without knowing the parameterization.
#'
#' @inheritParams roc_sdt
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible interval (default `c(0.025, 0.975)`).
#'
#' @return A data frame of class `"bmm_sdt_thresholds"` with columns `marker`
#'   (threshold label `t1`, `t2`, ...), `position` (latent location), `.draw`,
#'   and any condition columns. The object carries a `summary` attribute
#'   (`marker`, `position` posterior mean, `lower`, `upper`, plus condition
#'   columns). The `position`/`lower`/`upper` naming matches the `lines`
#'   attribute of [latent_sdt()], which visualises the same quantities.
#'
#' @seealso [latent_sdt()], [roc_sdt()]
#' @export
sdt_thresholds <- function(fit, conditions = NULL, probs = c(0.025, 0.975), ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt_rating"),
         "sdt_thresholds() is only available for rating SDT models (sdt_rating)")

  conditions <- .sdt_resolve_conditions(fit, conditions)
  cond_rows  <- .sdt_unique_subset(conditions, names(conditions))
  thr_list   <- .sdt_rating_thresholds(fit, model, cond_rows, ...)
  markers    <- paste0("t", seq_len(model$other_vars$n_ratings - 1L))

  draws_list   <- vector("list", length(thr_list))
  summary_list <- vector("list", length(thr_list))
  for (l_i in seq_along(thr_list)) {
    thr_draws <- thr_list[[l_i]]
    draws <- data.frame(
      marker   = rep(markers, each = nrow(thr_draws)),
      position = as.vector(thr_draws),
      .draw    = rep(seq_len(nrow(thr_draws)), times = ncol(thr_draws))
    )
    summ <- data.frame(
      marker = markers,
      stats::setNames(.sdt_summarise_draws(thr_draws, probs),
                      c("position", "lower", "upper"))
    )
    crow <- cond_rows[l_i, , drop = FALSE]
    draws_list[[l_i]]   <- .sdt_bind_cond(draws, crow)
    summary_list[[l_i]] <- .sdt_bind_cond(summ, crow)
  }

  structure(
    do.call(rbind, draws_list),
    class       = c("bmm_sdt_thresholds", "data.frame"),
    summary     = do.call(rbind, summary_list),
    probs       = probs,
    model_class = class(model),
    dist        = model$other_vars$dist,
    conditions  = cond_rows
  )
}


#' @export
print.bmm_sdt_thresholds <- function(x, ...) {
  model_name <- utils::tail(attr(x, "model_class"), 1L) %||% "sdt"
  cat("SDT decision thresholds (", model_name, ", dist = ", attr(x, "dist"),
      ")\n", sep = "")
  print(attr(x, "summary"), digits = 3, row.names = FALSE)
  invisible(x)
}


############################################################################# !
# AUC                                                                    ####
############################################################################# !

#' Area under the ROC curve from a fitted SDT model
#'
#' Computes the posterior area under the ROC curve (AUC). For Gaussian and
#' Gumbel-min equal-variance binary SDT the AUC is available in closed form from
#' the `d` draws; otherwise it is obtained by trapezoidal integration of
#' the ROC points from [roc_sdt()]. The returned AUC is always the area under
#' the full curve (for binary multi-criteria fits this is one value per curve,
#' not the trapezoid of the discrete operating points).
#'
#' The closed form is used only when `sdratio` is fixed, where `d` (which is
#' \eqn{d_a}) equals \eqn{d'}; every unequal-variance fit takes the numerical
#' route, so the AUC is invariant to the sensitivity parameterization.
#'
#' @inheritParams roc_sdt
#' @param probs Numeric vector of length 2. Quantiles for the credible interval
#'   (default `c(0.025, 0.975)`).
#'
#' @return A data frame of class `"bmm_sdt_auc"` with columns `AUC`, `.draw`,
#'   and any condition columns, plus a `summary` attribute (`AUC_mean`,
#'   `AUC_lower`, `AUC_upper`).
#'
#' @details Analytical formulas (equal variance, where \eqn{d_a = d'}): normal
#'   EV-SDT \eqn{AUC = \Phi(d'/\sqrt{2})}; Gumbel-min EV-SDT
#'   \eqn{AUC = \mathrm{logistic}(g')}.
#'
#' @seealso [roc_sdt()], [plot.bmm_sdt_auc()]
#' @export
auc_sdt <- function(fit, conditions = NULL, probs = c(0.025, 0.975),
                    criterion_points = NULL, ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "auc_sdt() is only available for SDT models")
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
    d_mat <- .sdt_linpred(fit, "d", conditions, ...)
    n_draws <- nrow(d_mat)
    result <- vector("list", ncol(d_mat))
    for (c_i in seq_len(ncol(d_mat))) {
      df <- data.frame(AUC = auc_fn(d_mat[, c_i]), .draw = seq_len(n_draws))
      result[[c_i]] <- .sdt_bind_cond(df, conditions[c_i, , drop = FALSE])
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
  split <- .sdt_condition_masks(auc_data, cond_cols)
  rows <- lapply(seq_along(split$masks), function(i) {
    summ <- .sdt_summarise_draws(auc_data$AUC[split$masks[[i]]], probs,
                                 prefix = "AUC")
    .sdt_bind_cond(summ, split$levels[i, , drop = FALSE])
  })
  do.call(rbind, rows)
}


#' @export
print.bmm_sdt_auc <- function(x, ...) {
  model_name <- utils::tail(attr(x, "model_class"), 1L) %||% "sdt"
  cat("SDT AUC (", model_name, ", dist = ", attr(x, "dist"), ")\n", sep = "")
  print(attr(x, "summary"), digits = 3, row.names = FALSE)
  invisible(x)
}


############################################################################# !
# SENSITIVITY                                                            ####
############################################################################# !

#' Sensitivity on the noise, signal, or root-mean-square scale
#'
#' Re-expresses the posterior sensitivity of a fitted SDT model against a
#' different reference standard deviation. bmm estimates \eqn{d_a}, which
#' measures the separation of the two evidence distributions in units of their
#' root-mean-square SD. Under unequal variance the same separation can also be
#' read against the noise SD (\eqn{d_N}, the classical \eqn{d'}) or against the
#' signal SD (\eqn{d_S}). This function returns any of the three as posterior
#' draws, so contrasts and intervals can be computed on whichever scale a
#' literature reports.
#'
#' @details
#' The model places the noise distribution at \eqn{-\delta/2} with SD
#' \eqn{\sigma_N = 1} and the signal distribution at \eqn{+\delta/2} with SD
#' \eqn{\sigma_S = \exp(\mathrm{sdratio})}, so \eqn{\delta} is the separation in
#' noise-SD units. Dividing that separation by each reference SD gives
#'
#' \deqn{d_N = \delta / \sigma_N, \quad d_S = \delta / \sigma_S, \quad
#'       d_a = \delta / \sqrt{(\sigma_N^2 + \sigma_S^2)/2}.}
#'
#' The estimated parameter is \eqn{d_a}, hence
#' \eqn{\delta = d_a \sqrt{(1 + \sigma_S^2)/2}} and
#'
#' \deqn{d_N = d_a \sqrt{(1 + \sigma_S^2)/2}, \qquad
#'       d_S = d_a \sqrt{(1 + \sigma_S^2)/2} \,/\, \sigma_S.}
#'
#' All three coincide when `sdratio` is 0 (equal variance), which is the case
#' for [sdt_mafc()] and for any fit that does not give `sdratio` a formula. The
#' conversion is applied draw by draw, so the returned intervals propagate the
#' joint posterior uncertainty in `d` and `sdratio` rather than combining
#' point estimates.
#'
#' Only \eqn{d_a} is invariant to which distribution is treated as the
#' reference. \eqn{d_N} and \eqn{d_S} are not comparable across conditions that
#' differ in `sdratio`: two conditions with identical discriminability can show
#' a large, confidently estimated \eqn{d_N} difference. Prefer \eqn{d_a} for
#' contrasts, and use \eqn{d_N}/\eqn{d_S} for comparison with published values.
#'
#' @inheritParams roc_sdt
#' @param measure Character vector naming the scales to return: `"da"`
#'   (root-mean-square SD, the estimated parameter), `"dn"` (noise SD), and/or
#'   `"ds"` (signal SD). Defaults to all three.
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible interval (default `c(0.025, 0.975)`).
#'
#' @return A data frame of class `"bmm_sdt_sensitivity"` with columns `measure`,
#'   `value`, `.draw`, and any condition columns. The object carries a `summary`
#'   attribute (`measure`, `mean`, `lower`, `upper`, plus condition columns).
#'
#' @seealso [auc_sdt()], [roc_sdt()], [latent_sdt()]
#' @export
sdt_sensitivity <- function(fit, measure = c("da", "dn", "ds"),
                            conditions = NULL, probs = c(0.025, 0.975), ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  measure <- match.arg(measure, several.ok = TRUE)
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "sdt_sensitivity() is only available for SDT models")
  stopif(!"d" %in% names(model$parameters),
         "sdt_sensitivity() requires a model with a sensitivity parameter `d`")

  conditions <- .sdt_resolve_conditions(fit, conditions)
  cond_rows  <- .sdt_unique_subset(conditions, names(conditions))

  geom <- .sdt_latent_geometry(fit, cond_rows,
                               .sdt_has_estimated_sdratio(model, fit), ...)
  scales <- list(da = geom$d, dn = geom$sep, ds = geom$sep / geom$sdratio)[measure]

  draws_list   <- vector("list", length(scales) * ncol(geom$d))
  summary_list <- vector("list", length(draws_list))
  i <- 0L
  for (nm in names(scales)) {
    for (c_i in seq_len(ncol(geom$d))) {
      i <- i + 1L
      crow <- cond_rows[c_i, , drop = FALSE]
      draws_list[[i]] <- .sdt_bind_cond(
        data.frame(measure = nm, value = scales[[nm]][, c_i],
                   .draw = seq_len(nrow(geom$d))),
        crow
      )
      summary_list[[i]] <- .sdt_bind_cond(
        data.frame(measure = nm,
                   .sdt_summarise_draws(scales[[nm]][, c_i], probs)),
        crow
      )
    }
  }

  structure(
    do.call(rbind, draws_list),
    class       = c("bmm_sdt_sensitivity", "data.frame"),
    summary     = do.call(rbind, summary_list),
    probs       = probs,
    model_class = class(model),
    dist        = model$other_vars$dist,
    conditions  = cond_rows
  )
}


#' @export
print.bmm_sdt_sensitivity <- function(x, ...) {
  model_name <- utils::tail(attr(x, "model_class"), 1L) %||% "sdt"
  cat("SDT sensitivity (", model_name, ", dist = ", attr(x, "dist"), ")\n",
      sep = "")
  cat("  da = RMS-SD units (estimated) | dn = noise-SD units (d') |",
      "ds = signal-SD units\n")
  print(attr(x, "summary"), digits = 3, row.names = FALSE)
  invisible(x)
}
