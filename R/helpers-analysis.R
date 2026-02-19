############################################################################# !
# HELPERS-ANALYSIS.R                                                     ####
# Post-fitting analysis functions for bmmfit objects.                   ####
############################################################################# !


############################################################################# !
# ROC CURVES                                                             ####
############################################################################# !

#' Model-implied ROC curve from a fitted SDT model
#'
#' Computes the receiver operating characteristic (ROC) curve implied by the
#' posterior distribution of an SDT model fit with [bmm()]. For binary SDT
#' models the ROC is traced analytically from the posterior of `dprime` (and
#' `sdratio` for UV-SDT) over a grid of criterion values. For confidence rating
#' models the K-1 thresholds define K-1 empirical ROC points per posterior draw.
#'
#' @param fit A `bmmfit` object returned by [bmm()] from an SDT model.
#' @param conditions Optional data frame of predictor values for which to
#'   compute the ROC curve. Each row defines one set of conditions. Column
#'   names must match the predictor variables used in the formula. If `NULL`
#'   (default), the function derives unique predictor combinations from the
#'   original data.
#' @param n_points Integer. Number of equally-spaced points on the smooth ROC
#'   curve for binary models (default 100). Ignored for rating models, which
#'   always return K+1 points (K-1 threshold points plus the two extreme
#'   endpoints (0,0) and (1,1)).
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible band (default `c(0.025, 0.975)` for a 95% CI).
#' @param ... Additional arguments passed to [brms::posterior_linpred()].
#'
#' @return A data frame of class `"bmm_sdt_roc"` with columns `FA`, `Hit`,
#'   `.draw`, and any condition columns from `conditions`. The object has a
#'   `summary` attribute with columns `FA`, `Hit_mean`, `Hit_lower`,
#'   `Hit_upper` (and condition columns) summarising over posterior draws.
#'
#' @details
#' ROC curves are meaningful only for SDT models with an explicit response
#' criterion. `roc_sdt()` raises an error for `version = "mafc"` and
#' `version = "ranking"`, which lack a criterion parameter.
#'
#' For Gumbel-min EV-SDT the ROC takes the power-ROC form
#' \eqn{Hit = FA^{\exp(-g')}}, which is available in closed form and can
#' be recognised as a testable signature of Gumbel-min noise.
#'
#' @seealso [auc_sdt()], [plot.bmm_sdt_roc()]
#' @export
roc_sdt <- function(fit, conditions = NULL, n_points = 100,
                    probs = c(0.025, 0.975), ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "roc_sdt() is only available for SDT models fit with sdt()")

  version <- model$version
  stopif(version %in% c("mafc", "ranking"),
         glue("ROC curves are not defined for version = '{version}'. ",
              "These models have no response criterion."))

  is_rating <- .is_sdt_rating(model)
  conditions <- .sdt_resolve_conditions(fit, conditions)

  roc_data <- if (is_rating) {
    .roc_sdt_rating(fit, model, conditions, ...)
  } else {
    .roc_sdt_binary(fit, model, conditions, n_points, ...)
  }

  summary_df <- .roc_sdt_summary(roc_data, conditions, probs)

  structure(
    roc_data,
    class      = c("bmm_sdt_roc", "data.frame"),
    summary    = summary_df,
    model_version = version,
    dist       = model$other_vars$dist,
    is_rating  = is_rating,
    conditions = conditions
  )
}


# Internal: unique predictor combinations from fit data.
# Excludes response columns, experimental-design variables added by bmm
# (stimulus, n_trials, dist_type), random-effect grouping variables, and
# matrix columns (e.g., the cbind() response for multinomial models).
.sdt_resolve_conditions <- function(fit, conditions) {
  if (!is.null(conditions)) return(as.data.frame(conditions))

  data  <- fit$data
  model <- fit$bmm$model

  # Variables added by bmm infrastructure that are not fixed-effect predictors
  exclude <- unique(c(
    unlist(model$resp_vars),
    model$other_vars$stimulus,
    model$other_vars$n_trials,
    "nTrials", "dist_type",
    names(brms::ranef(fit))   # random-effect grouping variables
  ))

  pred_cols <- setdiff(names(data), exclude)

  # Also drop matrix columns (e.g., cbind() combined response for rating models)
  is_matrix_col <- vapply(pred_cols, function(col) is.matrix(data[[col]]), logical(1))
  pred_cols <- pred_cols[!is_matrix_col]

  if (length(pred_cols) == 0L) {
    # Intercept-only: return a 1-row 0-column data frame as placeholder
    return(data.frame(.row = 1L)[1L, -1L, drop = FALSE])
  }
  unique(data[, pred_cols, drop = FALSE])
}


# Internal: call posterior_linpred for a single parameter, handling 0-column
# conditions. posterior_linpred() requires ALL formula variables in newdata, so
# we always base newdata on fit$data[1L, ] and override condition columns.
# Binary SDT models expose dprime/criterion as distributional parameters (dpar);
# rating models expose them as non-linear parameters (nlpar).
.sdt_linpred <- function(fit, param, conditions, is_rating = FALSE, ...) {
  n_rows <- max(1L, nrow(conditions))
  # Replicate first data row to get all required formula variables
  newdata <- fit$data[rep(1L, n_rows), , drop = FALSE]
  rownames(newdata) <- NULL

  # Override with user-supplied condition values
  if (ncol(conditions) > 0L) {
    for (col in names(conditions)) {
      if (col %in% names(newdata)) {
        newdata[[col]] <- conditions[[col]]
      }
    }
  }

  if (is_rating) {
    brms::posterior_linpred(fit, nlpar = param, newdata = newdata,
                            re_formula = NA, allow_new_levels = TRUE, ...)
  } else {
    brms::posterior_linpred(fit, dpar = param, newdata = newdata,
                            re_formula = NA, allow_new_levels = TRUE, ...)
  }
}


# Internal: smooth analytical ROC for binary EV- and UV-SDT
.roc_sdt_binary <- function(fit, model, conditions, n_points, ...) {
  dist      <- model$other_vars$dist
  variances <- model$other_vars$variances

  cdf  <- .SDT_DISTS[[dist]]$cdf
  qf   <- .SDT_DISTS[[dist]]$qf

  fa_grid  <- seq(0.001, 0.999, length.out = n_points)
  qf_vals  <- qf(fa_grid)              # quantile values at each FA point
  fa_full  <- c(0, fa_grid, 1)         # add trivial endpoints
  n_fa     <- length(fa_full)

  n_cond      <- max(1L, nrow(conditions))
  dprime_mat  <- .sdt_linpred(fit, "dprime", conditions, is_rating = FALSE, ...)
  n_draws     <- nrow(dprime_mat)

  sdratio_mat <- if (variances == "unequal") {
    exp(.sdt_linpred(fit, "sdratio", conditions, is_rating = FALSE, ...))
  } else {
    matrix(1, nrow = n_draws, ncol = n_cond)
  }

  cond_has_cols <- ncol(conditions) > 0L

  result <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    dp_vec <- dprime_mat[, c_i]   # [n_draws]
    sr_vec <- sdratio_mat[, c_i]  # [n_draws]

    # hit_interior[d, f] = cdf((dp[d] + qf_vals[f]) / sr[d])
    inner <- sweep(outer(dp_vec, qf_vals, "+"), 1L, sr_vec, "/")
    hit_interior <- cdf(inner)   # [n_draws × n_points]

    # Prepend Hit=0 at FA=0 and append Hit=1 at FA=1
    hit_full <- cbind(0, hit_interior, 1)  # [n_draws × n_fa]

    df <- data.frame(
      FA    = rep(fa_full, times = n_draws),
      Hit   = as.vector(t(hit_full)),
      .draw = rep(seq_len(n_draws), each = n_fa)
    )
    if (cond_has_cols) {
      df <- cbind(df, conditions[rep(c_i, nrow(df)), , drop = FALSE],
                  row.names = NULL)
    }
    result[[c_i]] <- df
  }
  do.call(rbind, result)
}


# Internal: empirical (K+1-point) ROC for rating SDT models
.roc_sdt_rating <- function(fit, model, conditions, ...) {
  dist           <- model$other_vars$dist
  variances      <- model$other_vars$variances
  version        <- model$version
  threshold_type <- model$other_vars$threshold_type
  n_ratings      <- model$other_vars$n_ratings

  n_cond      <- max(1L, nrow(conditions))
  dprime_mat  <- .sdt_linpred(fit, "dprime",    conditions, is_rating = TRUE, ...)
  crit_mat    <- .sdt_linpred(fit, "criterion", conditions, is_rating = TRUE, ...)
  n_draws     <- nrow(dprime_mat)

  # Threshold parameters
  if (threshold_type %in% c("equidistant", "parsimonious")) {
    spacing_mat <- .sdt_linpred(fit, "spacing", conditions, is_rating = TRUE, ...)
    delta_mats  <- NULL
  } else {
    delta_names <- grep("^delta", names(model$parameters), value = TRUE)
    delta_mats  <- lapply(delta_names, function(nm) {
      .sdt_linpred(fit, nm, conditions, is_rating = TRUE, ...)
    })
    names(delta_mats) <- delta_names
    spacing_mat <- NULL
  }

  sdratio_mat <- if (variances == "unequal") {
    exp(.sdt_linpred(fit, "sdratio", conditions, is_rating = TRUE, ...))
  } else NULL

  Ro_mat <- if (version == "dpsdt") {
    plogis(.sdt_linpred(fit, "Ro", conditions, is_rating = TRUE, ...))
  } else NULL

  metad_mat <- if (version == "metad") {
    .sdt_linpred(fit, "metad", conditions, is_rating = TRUE, ...)
  } else NULL

  n_roc <- n_ratings + 1L     # K-1 interior + 2 endpoints
  cond_has_cols <- ncol(conditions) > 0L

  result <- vector("list", n_cond)
  for (c_i in seq_len(n_cond)) {
    roc_rows <- vector("list", n_draws)
    for (d_i in seq_len(n_draws)) {
      dp   <- dprime_mat[d_i, c_i]
      crit <- crit_mat[d_i, c_i]
      sr   <- if (!is.null(sdratio_mat)) sdratio_mat[d_i, c_i] else 1

      # Reconstruct thresholds
      if (!is.null(spacing_mat)) {
        sp         <- spacing_mat[d_i, c_i]
        thresholds <- .sdt_make_thresholds(crit, n_ratings, threshold_type,
                                           spacing = sp)
      } else {
        deltas_d   <- vapply(delta_mats, function(m) m[d_i, c_i], numeric(1))
        thresholds <- .sdt_make_thresholds(crit, n_ratings, threshold_type,
                                           deltas = deltas_d)
      }

      # Category probability function, dispatched by model version
      cat_probs_fn <- switch(
        version,
        evsdt = function(stim) {
          shift <- dp / 2 * (2 * stim - 1)
          .sdt_category_probs(thresholds, shift, sr, stim, dist)
        },
        dpsdt = function(stim) {
          shift  <- dp / 2 * (2 * stim - 1)
          Ro_val <- Ro_mat[d_i, c_i]
          .sdt_dpsdt_category_probs(thresholds, shift, sr, stim, dist, Ro_val)
        },
        metad = function(stim) {
          metad_val <- metad_mat[d_i, c_i]
          .sdt_metad_category_probs(thresholds, dp, metad_val, stim, sr, dist)
        }
      )

      pn <- cat_probs_fn(0L)  # noise: length-K vector
      ps <- cat_probs_fn(1L)  # signal: length-K vector

      # Cumulative hit/FA: FA_k = P(resp > k | noise), Hit_k = P(resp > k | signal)
      # FA_1 = 1 (always "old"), FA_{K} = 0 (never "old")
      fa_pts <- c(1, 1 - cumsum(pn)[seq_len(n_ratings - 1L)], 0)
      hi_pts <- c(1, 1 - cumsum(ps)[seq_len(n_ratings - 1L)], 0)

      roc_rows[[d_i]] <- data.frame(FA = fa_pts, Hit = hi_pts, .draw = d_i)
    }

    cond_df <- do.call(rbind, roc_rows)
    if (cond_has_cols) {
      cond_df <- cbind(cond_df,
                       conditions[rep(c_i, nrow(cond_df)), , drop = FALSE],
                       row.names = NULL)
    }
    result[[c_i]] <- cond_df
  }
  do.call(rbind, result)
}


# Internal: summarise ROC draws → mean + quantile band at each (FA, condition)
.roc_sdt_summary <- function(roc_data, conditions, probs) {
  cond_cols <- setdiff(names(roc_data), c("FA", "Hit", ".draw"))

  # Round FA to avoid floating-point grouping mismatches
  roc_data$FA_r <- round(roc_data$FA, 8L)
  group_cols    <- c("FA_r", cond_cols)

  groups <- unique(roc_data[, group_cols, drop = FALSE])

  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    mask <- rep(TRUE, nrow(roc_data))
    for (col in group_cols) {
      mask <- mask & (roc_data[[col]] == groups[i, col])
    }
    hit_draws <- roc_data$Hit[mask]
    rows[[i]] <- data.frame(
      FA        = roc_data$FA[mask][1L],
      Hit_mean  = mean(hit_draws),
      Hit_lower = quantile(hit_draws, probs[1L]),
      Hit_upper = quantile(hit_draws, probs[2L]),
      groups[i, cond_cols, drop = FALSE],
      row.names = NULL,
      check.names = FALSE
    )
  }
  do.call(rbind, rows)
}


#' @export
print.bmm_sdt_roc <- function(x, ...) {
  version  <- attr(x, "model_version")
  dist     <- attr(x, "dist")
  is_rating <- attr(x, "is_rating")
  n_draws  <- length(unique(x$.draw))
  n_cond   <- nrow(attr(x, "conditions"))
  cond_cols <- setdiff(names(x), c("FA", "Hit", ".draw"))

  cat("SDT ROC curve (", version, ", dist = ", dist, ")\n", sep = "")
  cat("  ", n_draws, " posterior draws",
      if (n_cond > 1L) paste0(" x ", n_cond, " conditions") else "",
      "\n", sep = "")
  if (is_rating) {
    n_ratings <- nrow(x) / (n_draws * max(1L, n_cond))
    cat("  Rating model: ", n_ratings - 1L, " ROC points per draw\n", sep = "")
  } else {
    n_pts <- nrow(x) / (n_draws * max(1L, n_cond))
    cat("  Smooth curve: ", n_pts, " FA points per draw\n", sep = "")
  }
  cat("Use plot() to visualise, or attr(x, 'summary') for summaries.\n")
  invisible(x)
}


############################################################################# !
# AUC                                                                    ####
############################################################################# !

#' Area under the ROC curve from a fitted SDT model
#'
#' Computes the area under the ROC curve (AUC) from the posterior distribution
#' of an SDT model. For Gaussian and Gumbel-min EV-SDT the AUC is computed
#' analytically from `dprime` draws; for all other configurations it is
#' computed via trapezoidal integration on the ROC points from [roc_sdt()].
#'
#' @param fit A `bmmfit` object returned by [bmm()] from an SDT model.
#' @param conditions Optional data frame of predictor values (see [roc_sdt()]).
#' @param probs Numeric vector of length 2. Quantiles for the credible interval
#'   (default `c(0.025, 0.975)`).
#' @param ... Additional arguments passed to [brms::posterior_linpred()] (for
#'   analytical cases) or [roc_sdt()] (for numerical cases).
#'
#' @return A data frame of class `"bmm_sdt_auc"` with columns `AUC`, `.draw`,
#'   and any condition columns. The object has a `summary` attribute with
#'   columns `AUC_mean`, `AUC_lower`, `AUC_upper` (and condition columns).
#'
#' @details
#' Analytical formulas:
#' \itemize{
#'   \item Normal EV-SDT: \eqn{AUC = \Phi(d'/\sqrt{2})}
#'   \item Gumbel-min EV-SDT: \eqn{AUC = \text{logistic}(g') = 1/(1+\exp(-g'))}
#' }
#'
#' @seealso [roc_sdt()], [plot.bmm_sdt_auc()]
#' @export
auc_sdt <- function(fit, conditions = NULL, probs = c(0.025, 0.975), ...) {
  stopif(!inherits(fit, "bmmfit"),
         "fit must be a bmmfit object returned by bmm()")
  model <- fit$bmm$model
  stopif(!inherits(model, "sdt"),
         "auc_sdt() is only available for SDT models fit with sdt()")

  version   <- model$version
  stopif(version %in% c("mafc", "ranking"),
         glue("AUC is not defined for version = '{version}'."))

  dist      <- model$other_vars$dist
  variances <- model$other_vars$variances
  is_rating <- .is_sdt_rating(model)
  conditions <- .sdt_resolve_conditions(fit, conditions)

  # Analytical path: normal or gumbel_min EV-SDT, binary response
  use_analytical <- !is_rating && variances == "equal" &&
                    version == "evsdt" &&
                    dist %in% c("normal", "gumbel_min")

  if (use_analytical) {
    auc_fn <- if (dist == "normal") {
      function(d) pnorm(d / sqrt(2))
    } else {
      plogis   # logistic(g') for gumbel_min
    }

    dprime_mat <- .sdt_linpred(fit, "dprime", conditions, is_rating = FALSE, ...)
    n_draws    <- nrow(dprime_mat)
    n_cond     <- ncol(dprime_mat)

    result <- vector("list", n_cond)
    for (c_i in seq_len(n_cond)) {
      auc_vals <- auc_fn(dprime_mat[, c_i])
      df <- data.frame(AUC = auc_vals, .draw = seq_len(n_draws))
      if (ncol(conditions) > 0L) {
        df <- cbind(df, conditions[rep(c_i, n_draws), , drop = FALSE],
                    row.names = NULL)
      }
      result[[c_i]] <- df
    }
    auc_data <- do.call(rbind, result)

  } else {
    # Numerical path: trapezoidal integration on ROC points
    roc_obj  <- roc_sdt(fit, conditions = conditions, probs = probs, ...)
    cond_cols <- setdiff(names(roc_obj), c("FA", "Hit", ".draw"))

    # Group by draw × conditions, compute trapezoidal AUC
    if (length(cond_cols) > 0L) {
      group_key <- do.call(paste, c(roc_obj[c(".draw", cond_cols)], sep = "|||"))
    } else {
      group_key <- as.character(roc_obj$.draw)
    }

    unique_keys <- unique(group_key)
    auc_rows <- vector("list", length(unique_keys))

    for (k_i in seq_along(unique_keys)) {
      k   <- unique_keys[k_i]
      sub <- roc_obj[group_key == k, , drop = FALSE]
      sub <- sub[order(sub$FA), ]
      auc_val <- sum(diff(sub$FA) *
                     (utils::head(sub$Hit, -1L) + utils::tail(sub$Hit, -1L))) / 2
      row <- data.frame(AUC = auc_val, .draw = sub$.draw[1L])
      if (length(cond_cols) > 0L) {
        row <- cbind(row, sub[1L, cond_cols, drop = FALSE], row.names = NULL)
      }
      auc_rows[[k_i]] <- row
    }
    auc_data <- do.call(rbind, auc_rows)
  }

  summary_df <- .auc_sdt_summary(auc_data, conditions, probs)

  structure(
    auc_data,
    class         = c("bmm_sdt_auc", "data.frame"),
    summary       = summary_df,
    model_version = version,
    dist          = dist,
    conditions    = conditions
  )
}


# Internal: summarise AUC draws → mean + quantile band per condition
.auc_sdt_summary <- function(auc_data, conditions, probs) {
  cond_cols <- setdiff(names(auc_data), c("AUC", ".draw"))

  if (length(cond_cols) == 0L) {
    return(data.frame(
      AUC_mean  = mean(auc_data$AUC),
      AUC_lower = quantile(auc_data$AUC, probs[1L]),
      AUC_upper = quantile(auc_data$AUC, probs[2L])
    ))
  }

  groups <- unique(auc_data[, cond_cols, drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    mask <- rep(TRUE, nrow(auc_data))
    for (col in cond_cols) {
      mask <- mask & (auc_data[[col]] == groups[i, col])
    }
    auc_draws <- auc_data$AUC[mask]
    rows[[i]] <- data.frame(
      AUC_mean  = mean(auc_draws),
      AUC_lower = quantile(auc_draws, probs[1L]),
      AUC_upper = quantile(auc_draws, probs[2L]),
      groups[i, cond_cols, drop = FALSE],
      row.names = NULL,
      check.names = FALSE
    )
  }
  do.call(rbind, rows)
}


#' @export
print.bmm_sdt_auc <- function(x, ...) {
  version <- attr(x, "model_version")
  dist    <- attr(x, "dist")
  summ    <- attr(x, "summary")

  cat("SDT AUC (", version, ", dist = ", dist, ")\n", sep = "")
  print(summ, digits = 3, row.names = FALSE)
  invisible(x)
}

