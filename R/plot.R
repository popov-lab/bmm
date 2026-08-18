############################################################################# !
# PLOT.R                                                                  ####
# S3 plot() methods for bmm SDT analysis objects (ROC / AUC).            ####
############################################################################# !

utils::globalVariables(c(
  ".data", "FA", "Hit", "FA_mean", "FA_lower", "FA_upper",
  "Hit_mean", "Hit_lower", "Hit_upper", "AUC", "AUC_mean",
  "density", "distribution", "position", "lower", "upper",
  "level", "set_size"
))


############################################################################# !
# ROC                                                                     ####
############################################################################# !

#' Plot a model-implied SDT ROC curve
#'
#' Shows the smooth model-implied ROC as a credible-band ribbon, with the
#' model-implied operating points overlaid as colour-coded markers with
#' crosshair error bars (uncertainty in both directions). For **binary**
#' multi-criteria fits the points are the several criterion levels; for
#' **rating** models they are the K-1 confidence thresholds (labelled
#' `c1`..`c(K-1)`). Pass `observed` to overlay empirical points from
#' [roc_observed()].
#'
#' With `scale = "quantile"` (or its alias `scale = "z"`) the rates are read on
#' the distribution's quantile axis: `qf(rate)`, where `qf` is the inverse CDF of
#' the fitted noise distribution, and the axis label names that transform
#' (`z` for normal, `logit` for logistic, `loglog`/`cloglog` for the Gumbel
#' distributions). On this axis the binary model ROC is a straight line with
#' slope `1 / exp(sdratio)` and intercept `d * sqrt((1 + exp(sdratio)^2) / 2) /
#' exp(sdratio)` -- the separation in noise-SD units over the signal SD, since
#' `d` is \eqn{d_a} -- so a slope
#' below 1 is the unequal-variance signature (signal SD > noise SD), and
#' departures of the observed points from a straight line diagnose misfit. This
#' linearity holds for the symmetric distributions (`"normal"`, `"logistic"`);
#' for the Gumbel distributions the transformed ROC is curved (their natural
#' linearising transform is the log-log power-ROC). The (0,0) and (1,1)
#' endpoints map to infinity and are dropped on the transformed scale.
#'
#' @param x A `"bmm_sdt_roc"` object from [roc_sdt()].
#' @param observed Optional `"bmm_sdt_roc_observed"` object from [roc_observed()]
#'   to overlay as empirical points.
#' @param condition_col Optional character. Condition column for colour/faceting.
#'   If `NULL` (default), auto-detected.
#' @param add_diagonal Logical. Draw the chance-level diagonal (default `TRUE`).
#' @param scale Either `"probability"` (default) for the usual hit vs. false-
#'   alarm-rate axes, or `"quantile"` (alias `"z"`) for the quantile-transformed
#'   axes (the inverse CDF of the fitted noise distribution).
#' @param ribbon_alpha Numeric. Transparency of the credible band (default
#'   `0.25`).
#' @param point_size Numeric. Size of operating-point markers (default `2.5`).
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @seealso [roc_sdt()], [roc_observed()], [auc_sdt()]
#' @export
plot.bmm_sdt_roc <- function(x, observed = NULL, condition_col = NULL,
                             add_diagonal = TRUE,
                             scale = c("probability", "quantile", "z"),
                             ribbon_alpha = 0.25, point_size = 2.5, ...) {
  stopif(!requireNamespace("ggplot2", quietly = TRUE),
         "ggplot2 is required for plot.bmm_sdt_roc(). Please install it.")
  scale <- match.arg(scale)
  qf    <- if (scale != "probability") .sdt_dists[[attr(x, "dist")]]$qf

  cond_cols <- setdiff(names(x), c("FA", "Hit", ".draw"))
  colour_col <- .roc_colour_col(condition_col, cond_cols)

  p <- ggplot2::ggplot()
  if (add_diagonal) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
                                  linetype = "dashed", colour = "grey50")
  }

  p <- .plot_roc_curve(p, x, colour_col, ribbon_alpha, point_size, qf)

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
    zl <- .sdt_dists[[attr(x, "dist")]]$qf_label
    ggplot2::labs(x = paste0(zl, "(False alarm rate)"),
                  y = paste0(zl, "(Hit rate)"),
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


# Unified ROC drawing for binary and rating fits: both carry a `summary` (smooth
# implied curve) and `points` (operating points) attribute. The curve is drawn
# as a ribbon + line; the points get crosshair error bars and are coloured by
# `threshold` (rating: c1..cK-1) or, for binary, the single criterion-level
# column. Operating points fall on the curve because both are derived from the
# same probability map in roc_sdt().
.plot_roc_curve <- function(p, x, colour_col, ribbon_alpha, point_size,
                            qf = NULL) {
  summ      <- attr(x, "summary")
  points    <- attr(x, "points")
  cond_cols <- setdiff(names(x), c("FA", "Hit", ".draw"))

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
    pt_colour <- if ("threshold" %in% pt_cols) {
      "threshold"
    } else if (length(setdiff(pt_cols, cond_cols)) == 1L) {
      setdiff(pt_cols, cond_cols)
    }
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
#' decision axis. For [sdt_yn()]/[sdt_rating()] the response criterion or the
#' K-1 confidence thresholds are added as dashed vertical lines with a shaded
#' credible band; when several boundaries are shown together (e.g. base-rate
#' criteria or the rating thresholds) they are colour-coded. [sdt_mafc()] and
#' [sdt_ranking()] have no boundary; with `show_competitors = TRUE` in
#' [latent_sdt()] the max-of-distractors densities are overlaid as dashed lines,
#' one per set size. Distinct `d`/`sdratio` conditions are faceted.
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
  comp      <- attr(x, "competitors")
  n_levels  <- if (!is.null(lines)) length(unique(lines$level)) else 0L
  n_comp    <- if (!is.null(comp)) length(unique(comp$set_size)) else 0L
  # When several boundaries/competitors are colour-coded the colour scale is
  # spent on them, so distributions are shown by fill with a neutral outline;
  # otherwise the densities keep their coloured outline.
  colour_key <- n_levels > 1L || n_comp > 1L

  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[["x"]], y = density))

  if (!is.null(lines) && nrow(lines) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = lines, inherit.aes = FALSE,
      ggplot2::aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
      fill = "grey70", alpha = 0.3
    )
    p <- p + if (colour_key) {
      ggplot2::geom_vline(
        data = lines, inherit.aes = FALSE,
        ggplot2::aes(xintercept = position, colour = level),
        linetype = "dashed", alpha = line_alpha
      )
    } else {
      ggplot2::geom_vline(
        data = lines, inherit.aes = FALSE,
        ggplot2::aes(xintercept = position),
        linetype = "dashed", alpha = line_alpha, colour = "grey30"
      )
    }
  }

  p <- p + ggplot2::geom_area(ggplot2::aes(fill = distribution),
                              position = "identity", alpha = 0.3)
  p <- p + if (colour_key) {
    ggplot2::geom_line(ggplot2::aes(group = distribution),
                       linewidth = 0.7, colour = "grey25")
  } else {
    ggplot2::geom_line(ggplot2::aes(colour = distribution), linewidth = 0.7)
  }

  if (!is.null(comp) && nrow(comp) > 0L) {
    p <- p + if (colour_key) {
      ggplot2::geom_line(
        data = comp, inherit.aes = FALSE,
        ggplot2::aes(x = .data[["x"]], y = density, colour = set_size),
        linetype = "dashed", linewidth = 0.6
      )
    } else {
      ggplot2::geom_line(
        data = comp, inherit.aes = FALSE,
        ggplot2::aes(x = .data[["x"]], y = density),
        linetype = "dashed", linewidth = 0.6, colour = "grey40"
      )
    }
  }

  if (length(cond_cols) > 0L) {
    facet <- condition_col %||% cond_cols
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", paste(facet, collapse = " + ")))
    )
  }

  colour_lab <- if (n_comp > 1L) {
    "Set size (m)"
  } else if (colour_key) {
    if (isTRUE(attr(x, "is_rating"))) "Threshold" else "Criterion"
  } else {
    "Distribution"
  }
  p +
    ggplot2::labs(x = "Decision variable", y = "Density",
                  fill = "Distribution", colour = colour_lab) +
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
