############################################################################# !
# PLOT.R                                                                  ####
# S3 plot() methods for bmm model output objects.                        ####
############################################################################# !

# Suppress R CMD check notes for ggplot2 aes() bare column names and rlang .data
utils::globalVariables(c(".data", "FA", "FA_mean", "FA_lower", "FA_upper",
                          "Hit_mean", "Hit_lower", "Hit_upper",
                          "AUC", "AUC_mean", ".k", ".n"))


############################################################################# !
# ROC                                                                     ####
############################################################################# !

# Internal: per-threshold posterior summary for rating SDT ROC plots.
# Groups rows by position within each (draw × condition) group — the rating
# roc_sdt() always generates rows in the order: (1,1) endpoint, threshold 1
# through K-1, (0,0) endpoint — so row_number() is a reliable threshold index.
.roc_threshold_summary <- function(x, cond_cols, probs) {
  group_vars <- c(".draw", cond_cols)
  x |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::mutate(.k = dplyr::row_number(), .n = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(.k > 1L, .k < .n) |>
    dplyr::group_by(.k, dplyr::across(dplyr::all_of(cond_cols))) |>
    dplyr::summarise(
      FA_mean   = mean(FA),
      FA_lower  = quantile(FA,  probs[1L]),
      FA_upper  = quantile(FA,  probs[2L]),
      Hit_mean  = mean(Hit),
      Hit_lower = quantile(Hit, probs[1L]),
      Hit_upper = quantile(Hit, probs[2L]),
      .groups = "drop"
    )
}


#' Plot model-implied ROC curve
#'
#' For **rating SDT models**, shows the posterior mean hit vs. false-alarm rate
#' at each confidence threshold as points with crosshair error bars (uncertainty
#' in both the x/FA and y/Hit directions), connected by a line through the
#' posterior mean points.
#'
#' For **binary SDT models**, shows the smooth analytical ROC as a credible-band
#' ribbon around the posterior mean curve.
#'
#' @param x A `"bmm_sdt_roc"` object from [roc_sdt()].
#' @param condition_col Optional character. Name of a condition column to use
#'   for colour/faceting. If `NULL` (default), auto-detected from the object.
#' @param add_diagonal Logical. Whether to draw the chance-level diagonal
#'   (default `TRUE`).
#' @param ribbon_alpha Numeric. Transparency of the credible-band ribbon for
#'   binary models (default `0.25`).
#' @param probs Numeric vector of length 2. Lower and upper quantiles for
#'   credible intervals shown as error bars (rating) or ribbon (binary).
#'   Default `c(0.025, 0.975)`.
#' @param point_size Numeric. Size of threshold-position points for rating
#'   models (default `2.5`).
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @seealso [roc_sdt()], [auc_sdt()]
#' @export
plot.bmm_sdt_roc <- function(x, condition_col = NULL, add_diagonal = TRUE,
                              ribbon_alpha = 0.25, probs = c(0.025, 0.975),
                              point_size = 2.5, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot.bmm_sdt_roc(). Please install it.")
  }

  is_rating <- isTRUE(attr(x, "is_rating"))
  cond_cols <- setdiff(names(x), c("FA", "Hit", ".draw"))

  colour_col <- if (!is.null(condition_col) && condition_col %in% cond_cols) {
    condition_col
  } else if (length(cond_cols) == 1L) {
    cond_cols
  } else {
    NULL
  }

  p <- ggplot2::ggplot()

  if (add_diagonal) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
                                   linetype = "dashed", colour = "grey50")
  }

  if (is_rating) {
    threshold_summ <- .roc_threshold_summary(x, cond_cols, probs)

    # Connecting line through posterior mean (FA, Hit) at each threshold +
    # the trivial (0,0) and (1,1) endpoints
    endpoint_cols <- c("FA_mean", "Hit_mean", cond_cols)
    if (length(cond_cols) > 0L) {
      unique_conds <- unique(threshold_summ[, cond_cols, drop = FALSE])
      n_conds <- nrow(unique_conds)
      endpoints <- cbind(
        data.frame(FA_mean  = rep(c(0, 1), each = n_conds),
                   Hit_mean = rep(c(0, 1), each = n_conds)),
        rbind(unique_conds, unique_conds),
        row.names = NULL
      )
    } else {
      endpoints <- data.frame(FA_mean = c(0, 1), Hit_mean = c(0, 1))
    }
    curve_data <- dplyr::bind_rows(
      threshold_summ[, endpoint_cols, drop = FALSE],
      endpoints
    ) |>
      dplyr::arrange(dplyr::across(dplyr::all_of(c(cond_cols, "FA_mean"))))

    if (!is.null(colour_col)) {
      aes_line  <- ggplot2::aes(x = FA_mean, y = Hit_mean,
                                colour = .data[[colour_col]])
      aes_errh  <- ggplot2::aes(y = Hit_mean, xmin = FA_lower, xmax = FA_upper,
                                colour = .data[[colour_col]])
      aes_errv  <- ggplot2::aes(x = FA_mean, ymin = Hit_lower, ymax = Hit_upper,
                                colour = .data[[colour_col]])
      aes_point <- ggplot2::aes(x = FA_mean, y = Hit_mean,
                                colour = .data[[colour_col]])
    } else {
      aes_line  <- ggplot2::aes(x = FA_mean, y = Hit_mean)
      aes_errh  <- ggplot2::aes(y = Hit_mean, xmin = FA_lower, xmax = FA_upper)
      aes_errv  <- ggplot2::aes(x = FA_mean, ymin = Hit_lower, ymax = Hit_upper)
      aes_point <- ggplot2::aes(x = FA_mean, y = Hit_mean)
    }

    p <- p +
      ggplot2::geom_line(data = curve_data, aes_line, linewidth = 0.7) +
      ggplot2::geom_errorbar(data = threshold_summ, aes_errh,
                             width = 0, linewidth = 0.5,
                             orientation = "y") +
      ggplot2::geom_errorbar(data = threshold_summ, aes_errv,
                             width = 0, linewidth = 0.5) +
      ggplot2::geom_point(data = threshold_summ, aes_point, size = point_size)

  } else {
    summ <- attr(x, "summary")

    if (!is.null(colour_col)) {
      aes_ribbon <- ggplot2::aes(x = FA, ymin = Hit_lower, ymax = Hit_upper,
                                  fill = .data[[colour_col]])
      aes_line   <- ggplot2::aes(x = FA, y = Hit_mean,
                                  colour = .data[[colour_col]])
    } else {
      aes_ribbon <- ggplot2::aes(x = FA, ymin = Hit_lower, ymax = Hit_upper)
      aes_line   <- ggplot2::aes(x = FA, y = Hit_mean)
    }

    p <- p +
      ggplot2::geom_ribbon(data = summ, aes_ribbon, alpha = ribbon_alpha) +
      ggplot2::geom_line(data = summ, aes_line, linewidth = 0.8)
  }

  if (!is.null(colour_col) && length(cond_cols) > 1L) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", setdiff(cond_cols, colour_col)))
    )
  }

  p + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "False alarm rate", y = "Hit rate",
                  colour = colour_col, fill = colour_col) +
    ggplot2::theme_minimal()
}


############################################################################# !
# AUC                                                                     ####
############################################################################# !

#' Plot posterior AUC distribution from a fitted SDT model
#'
#' @param x A `"bmm_sdt_auc"` object from [auc_sdt()].
#' @param condition_col Optional character. Name of a condition column to use
#'   for colour. If `NULL` (default), auto-detected from the object.
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @seealso [auc_sdt()], [roc_sdt()]
#' @export
plot.bmm_sdt_auc <- function(x, condition_col = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot.bmm_sdt_auc(). Please install it.")
  }

  cond_cols <- setdiff(names(x), c("AUC", ".draw"))
  colour_col <- if (!is.null(condition_col) && condition_col %in% cond_cols) {
    condition_col
  } else if (length(cond_cols) == 1L) {
    cond_cols
  } else {
    NULL
  }

  summ <- attr(x, "summary")

  if (!is.null(colour_col)) {
    p <- ggplot2::ggplot(x, ggplot2::aes(x = AUC,
                                          colour = .data[[colour_col]],
                                          fill   = .data[[colour_col]]))
  } else {
    p <- ggplot2::ggplot(x, ggplot2::aes(x = AUC))
  }

  p <- p +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        colour = "grey50")

  # Add posterior mean markers from summary
  if (!is.null(colour_col)) {
    p <- p + ggplot2::geom_vline(
      data = summ,
      ggplot2::aes(xintercept = AUC_mean, colour = .data[[colour_col]]),
      linewidth = 1
    )
  } else {
    p <- p + ggplot2::geom_vline(xintercept = summ$AUC_mean[1L],
                                  linewidth = 1, colour = "black")
  }

  p +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "AUC", y = "Density",
                  colour = colour_col, fill = colour_col) +
    ggplot2::theme_minimal()
}
