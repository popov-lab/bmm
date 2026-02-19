############################################################################# !
# PLOT.R                                                                  ####
# S3 plot() methods for bmm model output objects.                        ####
############################################################################# !

# Suppress R CMD check notes for ggplot2 aes() bare column names and rlang .data
utils::globalVariables(c(".data", "FA", "Hit_mean", "Hit_lower", "Hit_upper",
                          "AUC", "AUC_mean"))


############################################################################# !
# ROC                                                                     ####
############################################################################# !

#' Plot model-implied ROC curve
#'
#' @param x A `"bmm_sdt_roc"` object from [roc_sdt()].
#' @param condition_col Optional character. Name of a condition column to use
#'   for colour/faceting. If `NULL` (default), auto-detected from the object.
#' @param add_diagonal Logical. Whether to draw the chance-level diagonal
#'   (default `TRUE`).
#' @param ribbon_alpha Numeric. Transparency of the credible-band ribbon
#'   (default `0.25`).
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @seealso [roc_sdt()], [auc_sdt()]
#' @export
plot.bmm_sdt_roc <- function(x, condition_col = NULL, add_diagonal = TRUE,
                              ribbon_alpha = 0.25, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot.bmm_sdt_roc(). Please install it.")
  }

  summ <- attr(x, "summary")
  cond_cols <- setdiff(names(summ), c("FA", "Hit_mean", "Hit_lower", "Hit_upper"))

  # Determine colour/facet column
  colour_col <- if (!is.null(condition_col) && condition_col %in% cond_cols) {
    condition_col
  } else if (length(cond_cols) == 1L) {
    cond_cols
  } else {
    NULL
  }

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = FA, y = Hit_mean))

  if (!is.null(colour_col)) {
    aes_ribbon <- ggplot2::aes(ymin = Hit_lower, ymax = Hit_upper,
                                fill = .data[[colour_col]])
    aes_line   <- ggplot2::aes(colour = .data[[colour_col]])
  } else {
    aes_ribbon <- ggplot2::aes(ymin = Hit_lower, ymax = Hit_upper)
    aes_line   <- ggplot2::aes()
  }

  p <- p +
    ggplot2::geom_ribbon(aes_ribbon, alpha = ribbon_alpha) +
    ggplot2::geom_line(aes_line, linewidth = 0.8)

  if (add_diagonal) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
                                   linetype = "dashed", colour = "grey50")
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
