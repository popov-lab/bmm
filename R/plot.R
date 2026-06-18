############################################################################# !
# PLOT.R                                                                  ####
# S3 plot() methods for bmm SDT analysis objects (ROC / AUC).            ####
############################################################################# !

utils::globalVariables(c(
  ".data", "FA", "Hit", "FA_mean", "FA_lower", "FA_upper",
  "Hit_mean", "Hit_lower", "Hit_upper", "AUC", "AUC_mean",
  "density", "distribution", "position", "lower", "upper"
))


############################################################################# !
# ROC                                                                     ####
############################################################################# !

# Per-threshold posterior summary for rating ROC plots. The rating roc_sdt()
# emits rows per (draw x condition) in the order (1,1), thresholds 1..K-1,
# (0,0), so the within-group position indexes the threshold; the two endpoints
# (first and last position) are dropped.
.roc_threshold_summary <- function(x, cond_cols, probs) {
  key <- if (length(cond_cols) > 0L) {
    do.call(paste, c(x[c(".draw", cond_cols)], sep = "\r"))
  } else {
    as.character(x$.draw)
  }
  x$.k <- stats::ave(seq_len(nrow(x)), key, FUN = seq_along)
  x$.n <- stats::ave(seq_len(nrow(x)), key, FUN = length)
  interior <- x[x$.k > 1L & x$.k < x$.n, , drop = FALSE]

  grp_cols <- c(".k", cond_cols)
  groups <- unique(interior[, grp_cols, drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    mask <- rep(TRUE, nrow(interior))
    for (col in grp_cols) mask <- mask & (interior[[col]] == groups[i, col])
    fa <- interior$FA[mask]
    hit <- interior$Hit[mask]
    rows[[i]] <- data.frame(
      FA_mean   = mean(fa),
      FA_lower  = unname(stats::quantile(fa, probs[1L])),
      FA_upper  = unname(stats::quantile(fa, probs[2L])),
      Hit_mean  = mean(hit),
      Hit_lower = unname(stats::quantile(hit, probs[1L])),
      Hit_upper = unname(stats::quantile(hit, probs[2L])),
      groups[i, cond_cols, drop = FALSE],
      row.names = NULL, check.names = FALSE
    )
  }
  do.call(rbind, rows)
}


#' Plot a model-implied SDT ROC curve
#'
#' For **rating** models, shows the posterior-mean hit vs. false-alarm rate at
#' each confidence threshold as points with crosshair error bars (uncertainty in
#' both directions), connected through the (0,0) and (1,1) endpoints. For
#' **binary** models, shows the smooth analytical ROC as a credible-band ribbon;
#' for binary multi-criteria fits the model-implied operating points are overlaid
#' with crosshair error bars. Pass `observed` to overlay empirical points from
#' [roc_observed()].
#'
#' With `scale = "z"` the rates are read on the distribution's quantile
#' (probit-style) axis: `z = qf(rate)`, where `qf` is the inverse CDF of the
#' fitted noise distribution. On this axis the binary model ROC is a straight
#' line with slope `1 / exp(sdratio)` and intercept `dprime / exp(sdratio)`, so
#' a slope below 1 is the unequal-variance signature (signal SD > noise SD), and
#' departures of the observed points from a straight line diagnose
#' misfit. This linearity holds for the symmetric distributions (`"normal"`,
#' `"logistic"`); for the Gumbel distributions the z-ROC is curved (their
#' natural linearising transform is the log-log power-ROC). The (0,0) and (1,1)
#' endpoints map to infinity and are dropped on the z scale.
#'
#' @param x A `"bmm_sdt_roc"` object from [roc_sdt()].
#' @param observed Optional `"bmm_sdt_roc_observed"` object from [roc_observed()]
#'   to overlay as empirical points.
#' @param condition_col Optional character. Condition column for colour/faceting.
#'   If `NULL` (default), auto-detected.
#' @param add_diagonal Logical. Draw the chance-level diagonal (default `TRUE`).
#' @param scale Either `"probability"` (default) for the usual hit vs. false-
#'   alarm-rate axes, or `"z"` for the z-transformed (quantile) axes.
#' @param ribbon_alpha Numeric. Transparency of the binary credible band
#'   (default `0.25`).
#' @param point_size Numeric. Size of operating-point markers (default `2.5`).
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @seealso [roc_sdt()], [roc_observed()], [auc_sdt()]
#' @export
plot.bmm_sdt_roc <- function(x, observed = NULL, condition_col = NULL,
                             add_diagonal = TRUE,
                             scale = c("probability", "z"),
                             ribbon_alpha = 0.25, point_size = 2.5, ...) {
  stopif(!requireNamespace("ggplot2", quietly = TRUE),
         "ggplot2 is required for plot.bmm_sdt_roc(). Please install it.")
  scale <- match.arg(scale)
  qf    <- if (scale == "z") .SDT_DISTS[[attr(x, "dist")]]$qf

  is_rating <- isTRUE(attr(x, "is_rating"))
  probs     <- attr(x, "probs") %||% c(0.025, 0.975)
  cond_cols <- setdiff(names(x), c("FA", "Hit", ".draw"))
  colour_col <- .roc_colour_col(condition_col, cond_cols)

  p <- ggplot2::ggplot()
  if (add_diagonal) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
                                  linetype = "dashed", colour = "grey50")
  }

  p <- if (is_rating) {
    .plot_roc_rating(p, x, cond_cols, colour_col, probs, point_size, qf)
  } else {
    .plot_roc_binary(p, x, colour_col, ribbon_alpha, point_size, qf)
  }

  if (!is.null(observed)) {
    obs <- .roc_observed_xy(observed)
    if (!is.null(qf)) {
      obs <- obs[obs$FA > 0 & obs$FA < 1 & obs$Hit > 0 & obs$Hit < 1, , drop = FALSE]
      obs$FA <- qf(obs$FA)
      obs$Hit <- qf(obs$Hit)
    }
    p <- p + ggplot2::geom_point(
      data = obs, ggplot2::aes(x = FA, y = Hit), shape = 1L,
      size = point_size, stroke = 0.9
    )
  }

  if (!is.null(colour_col) && length(cond_cols) > 1L) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", setdiff(cond_cols, colour_col)))
    )
  }

  coords <- if (is.null(qf)) {
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1))
  } else {
    ggplot2::coord_equal()
  }
  labs <- if (is.null(qf)) {
    ggplot2::labs(x = "False alarm rate", y = "Hit rate",
                  colour = colour_col, fill = colour_col)
  } else {
    ggplot2::labs(x = "z(False alarm rate)", y = "z(Hit rate)",
                  colour = colour_col, fill = colour_col)
  }
  p + coords + labs + ggplot2::theme_minimal()
}


.roc_colour_col <- function(condition_col, cond_cols) {
  if (!is.null(condition_col) && condition_col %in% cond_cols) {
    condition_col
  } else if (length(cond_cols) == 1L) {
    cond_cols
  }
}


.roc_observed_xy <- function(observed) {
  obs <- as.data.frame(observed)
  obs[, c("FA", "Hit", setdiff(names(obs), c("FA", "Hit"))), drop = FALSE]
}


.plot_roc_rating <- function(p, x, cond_cols, colour_col, probs, point_size,
                             qf = NULL) {
  summ <- .roc_threshold_summary(x, cond_cols, probs)
  if (!is.null(qf)) {
    summ[c("FA_mean", "FA_lower", "FA_upper",
           "Hit_mean", "Hit_lower", "Hit_upper")] <-
      lapply(summ[c("FA_mean", "FA_lower", "FA_upper",
                    "Hit_mean", "Hit_lower", "Hit_upper")], qf)
  }

  if (!is.null(qf)) {
    # Endpoints (0,0)/(1,1) are +/-Inf on the z scale: connect interior points.
    curve <- if (length(cond_cols) > 0L) {
      summ[do.call(order, summ[c(cond_cols, "FA_mean")]),
           c("FA_mean", "Hit_mean", cond_cols)]
    } else {
      summ[order(summ$FA_mean), c("FA_mean", "Hit_mean")]
    }
  } else if (length(cond_cols) > 0L) {
    conds <- unique(summ[, cond_cols, drop = FALSE])
    endpoints <- cbind(
      data.frame(FA_mean = rep(c(0, 1), each = nrow(conds)),
                 Hit_mean = rep(c(0, 1), each = nrow(conds))),
      rbind(conds, conds), row.names = NULL
    )
    curve <- rbind(summ[, c("FA_mean", "Hit_mean", cond_cols)], endpoints)
    curve <- curve[do.call(order, curve[c(cond_cols, "FA_mean")]), ]
  } else {
    curve <- rbind(summ[, c("FA_mean", "Hit_mean")],
                   data.frame(FA_mean = c(0, 1), Hit_mean = c(0, 1)))
    curve <- curve[order(curve$FA_mean), ]
  }

  aes_pt <- if (!is.null(colour_col)) {
    ggplot2::aes(x = FA_mean, y = Hit_mean, colour = .data[[colour_col]])
  } else {
    ggplot2::aes(x = FA_mean, y = Hit_mean)
  }
  aes_line <- aes_pt
  aes_h <- if (!is.null(colour_col)) {
    ggplot2::aes(y = Hit_mean, xmin = FA_lower, xmax = FA_upper,
                 colour = .data[[colour_col]])
  } else {
    ggplot2::aes(y = Hit_mean, xmin = FA_lower, xmax = FA_upper)
  }
  aes_v <- if (!is.null(colour_col)) {
    ggplot2::aes(x = FA_mean, ymin = Hit_lower, ymax = Hit_upper,
                 colour = .data[[colour_col]])
  } else {
    ggplot2::aes(x = FA_mean, ymin = Hit_lower, ymax = Hit_upper)
  }

  p +
    ggplot2::geom_line(data = curve, aes_line, linewidth = 0.7) +
    ggplot2::geom_errorbar(data = summ, aes_h, width = 0, linewidth = 0.5,
                           orientation = "y") +
    ggplot2::geom_errorbar(data = summ, aes_v, width = 0, linewidth = 0.5) +
    ggplot2::geom_point(data = summ, aes_pt, size = point_size)
}


.plot_roc_binary <- function(p, x, colour_col, ribbon_alpha, point_size,
                             qf = NULL) {
  summ   <- attr(x, "summary")
  points <- attr(x, "points")

  if (!is.null(qf)) {
    summ <- summ[summ$FA > 0 & summ$FA < 1, , drop = FALSE]
    summ$FA <- qf(summ$FA)
    summ[c("Hit_mean", "Hit_lower", "Hit_upper")] <-
      lapply(summ[c("Hit_mean", "Hit_lower", "Hit_upper")], qf)
    if (!is.null(points)) {
      points[c("FA_mean", "FA_lower", "FA_upper",
               "Hit_mean", "Hit_lower", "Hit_upper")] <-
        lapply(points[c("FA_mean", "FA_lower", "FA_upper",
                        "Hit_mean", "Hit_lower", "Hit_upper")], qf)
    }
  }

  aes_ribbon <- if (!is.null(colour_col)) {
    ggplot2::aes(x = FA, ymin = Hit_lower, ymax = Hit_upper,
                 fill = .data[[colour_col]])
  } else {
    ggplot2::aes(x = FA, ymin = Hit_lower, ymax = Hit_upper)
  }
  aes_line <- if (!is.null(colour_col)) {
    ggplot2::aes(x = FA, y = Hit_mean, colour = .data[[colour_col]])
  } else {
    ggplot2::aes(x = FA, y = Hit_mean)
  }

  p <- p +
    ggplot2::geom_ribbon(data = summ, aes_ribbon, alpha = ribbon_alpha) +
    ggplot2::geom_line(data = summ, aes_line, linewidth = 0.8)

  if (!is.null(points)) {
    pt_cols <- setdiff(names(points), c("FA_mean", "FA_lower", "FA_upper",
                                        "Hit_mean", "Hit_lower", "Hit_upper"))
    pt_colour <- if (length(pt_cols) == 1L) pt_cols
    aes_pt <- if (!is.null(pt_colour)) {
      ggplot2::aes(x = FA_mean, y = Hit_mean, colour = .data[[pt_colour]])
    } else {
      ggplot2::aes(x = FA_mean, y = Hit_mean)
    }
    eb_h <- ggplot2::aes(y = Hit_mean, xmin = FA_lower, xmax = FA_upper)
    eb_v <- ggplot2::aes(x = FA_mean, ymin = Hit_lower, ymax = Hit_upper)
    p <- p +
      ggplot2::geom_errorbar(data = points, mapping = eb_h, width = 0,
                             linewidth = 0.4, orientation = "y",
                             colour = "grey30") +
      ggplot2::geom_errorbar(data = points, mapping = eb_v, width = 0,
                             linewidth = 0.4, colour = "grey30") +
      ggplot2::geom_point(data = points, mapping = aes_pt, size = point_size)
  }
  p
}


############################################################################# !
# LATENT DISTRIBUTIONS                                                    ####
############################################################################# !

#' Plot the latent decision-variable distributions of an SDT model
#'
#' Draws the model-implied noise and signal evidence densities on the latent
#' decision axis, with the response criterion (binary) or the K-1 confidence
#' thresholds (rating) as dashed vertical lines, each carrying a shaded credible
#' band on its location. Distributions are distinguished by colour/fill;
#' conditions are faceted.
#'
#' @param x A `"bmm_sdt_latent"` object from [latent_sdt()].
#' @param condition_col Optional character. Condition column(s) for faceting. If
#'   `NULL` (default), all condition columns are used.
#' @param line_alpha Numeric. Opacity of the criterion/threshold lines (default
#'   `0.9`).
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @seealso [latent_sdt()], [roc_sdt()]
#' @export
plot.bmm_sdt_latent <- function(x, condition_col = NULL, line_alpha = 0.9, ...) {
  stopif(!requireNamespace("ggplot2", quietly = TRUE),
         "ggplot2 is required for plot.bmm_sdt_latent(). Please install it.")

  cond_cols <- setdiff(names(x), c("x", "density", "distribution"))
  lines     <- attr(x, "lines")

  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[["x"]], y = density)) +
    ggplot2::geom_rect(
      data = lines, inherit.aes = FALSE,
      ggplot2::aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
      fill = "grey70", alpha = 0.3
    ) +
    ggplot2::geom_vline(
      data = lines, inherit.aes = FALSE,
      ggplot2::aes(xintercept = position), linetype = "dashed",
      colour = "grey30", alpha = line_alpha
    ) +
    ggplot2::geom_area(ggplot2::aes(fill = distribution),
                       position = "identity", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(colour = distribution), linewidth = 0.7)

  if (length(cond_cols) > 0L) {
    facet <- condition_col %||% cond_cols
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", paste(facet, collapse = " + ")))
    )
  }

  p +
    ggplot2::labs(x = "Decision variable", y = "Density",
                  fill = "Distribution", colour = "Distribution") +
    ggplot2::theme_minimal()
}


############################################################################# !
# AUC                                                                     ####
############################################################################# !

#' Plot the posterior AUC distribution from a fitted SDT model
#'
#' @param x A `"bmm_sdt_auc"` object from [auc_sdt()].
#' @param condition_col Optional character. Condition column for colour. If
#'   `NULL` (default), auto-detected.
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @seealso [auc_sdt()], [roc_sdt()]
#' @export
plot.bmm_sdt_auc <- function(x, condition_col = NULL, ...) {
  stopif(!requireNamespace("ggplot2", quietly = TRUE),
         "ggplot2 is required for plot.bmm_sdt_auc(). Please install it.")

  cond_cols  <- setdiff(names(x), c("AUC", ".draw"))
  colour_col <- .roc_colour_col(condition_col, cond_cols)
  summ <- attr(x, "summary")

  p <- if (!is.null(colour_col)) {
    ggplot2::ggplot(x, ggplot2::aes(x = AUC, colour = .data[[colour_col]],
                                    fill = .data[[colour_col]]))
  } else {
    ggplot2::ggplot(x, ggplot2::aes(x = AUC))
  }

  p <- p +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50")

  p <- if (!is.null(colour_col)) {
    vline_aes <- ggplot2::aes(xintercept = AUC_mean, colour = .data[[colour_col]])
    p + ggplot2::geom_vline(data = summ, mapping = vline_aes, linewidth = 1)
  } else {
    p + ggplot2::geom_vline(xintercept = summ$AUC_mean[1L], linewidth = 1,
                            colour = "black")
  }

  p +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "AUC", y = "Density", colour = colour_col, fill = colour_col) +
    ggplot2::theme_minimal()
}
