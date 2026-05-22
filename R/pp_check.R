############################################################################# !
# PP_CHECK.R                                                              ####
# pp_check() S3 method for bmmfit objects.                               ####
############################################################################# !

#' Posterior predictive check for bmmfit objects
#'
#' For models where brms provides `pp_check` support, this method delegates to
#' [brms::pp_check()]. For models with multi-column responses (e.g., confidence
#' rating SDT), brms's `pp_check` is unavailable; this method dispatches to a
#' model-specific visualisation instead.
#'
#' For **rating SDT models**, the plot shows the observed proportion of
#' responses in each rating category (black points and dashed line) against the
#' posterior predictive mean and credible band (blue), faceted by stimulus type
#' and any experimental conditions.
#'
#' @param object A `bmmfit` object returned by [bmm()].
#' @param type Character. Type of pp_check. For non-rating models, passed to
#'   [brms::pp_check()]. Rating SDT models currently produce a response
#'   proportion profile regardless of the value supplied.
#' @param ndraws Integer. Number of posterior draws. Defaults to `100` for
#'   rating SDT models; otherwise passed to [brms::pp_check()].
#' @param ... For non-rating models: additional arguments forwarded to
#'   [brms::pp_check()]. For rating SDT models: arguments forwarded to
#'   [brms::posterior_predict()], plus `probs` (numeric vector of length 2,
#'   default `c(0.025, 0.975)`) to control the credible interval.
#' @return For rating SDT models, a `ggplot2` object. For other models, the
#'   result of [brms::pp_check()].
#' @seealso [roc_sdt()], [auc_sdt()], [brms::pp_check()]
#' @importFrom brms pp_check
#' @export
pp_check.bmmfit <- function(object, type = NULL, ndraws = NULL, ...) {
  model <- object$bmm$model
  is_rating_family <- inherits(model, "sdt") &&
    !inherits(model, "sdt_cdp") &&
    !inherits(model, "sdt_mafc") &&
    !inherits(model, "sdt_ranking") &&
    !is.null(model$other_vars$n_ratings)

  if (is_rating_family) {
    if (!is.null(type)) {
      warning2("Argument 'type' is ignored for rating SDT pp_check. ",
               "A response-proportion plot is always produced.")
    }
    if (is.null(ndraws)) ndraws <- 100L
    .pp_check_sdt_rating(object, ndraws = ndraws, ...)
  } else {
    NextMethod()
  }
}


# Internal: pp_check for rating SDT models.
# Aggregates observed and predicted response proportions by stimulus x
# condition and plots them as a faceted ribbon + observed line chart.
#
# Explicitly captures brms::pp_check arguments that do not apply here so they
# are not accidentally forwarded to brms::posterior_predict() and cause cryptic
# "unused argument" errors. type is warned about; layout args (group, x,
# prefix, draw_ids, nsamples, subset) are silently absorbed.
.pp_check_sdt_rating <- function(fit, ndraws = 100L,
                                 probs = c(0.025, 0.975),
                                 # brms pp_check args — intercepted, not forwarded
                                 type = NULL, group = NULL, x = NULL,
                                 prefix = c("ppc", "ppd"),
                                 draw_ids = NULL,
                                 nsamples = NULL, subset = NULL,
                                 ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2("ggplot2 is required for pp_check of rating SDT models. ",
          "Please install it.")
  }

  # nsamples is a deprecated brms alias for ndraws
  if (!is.null(nsamples)) {
    warning2("'nsamples' is deprecated; use 'ndraws' instead.")
    if (is.null(draw_ids)) ndraws <- nsamples
  }

  # draw_ids overrides ndraws for the posterior_predict call
  predict_args <- if (!is.null(draw_ids)) list(draw_ids = draw_ids) else list(ndraws = ndraws)

  model     <- fit$bmm$model
  stim_var  <- model$other_vars$stimulus
  n_ratings <- model$other_vars$n_ratings
  long_data <- fit$data
  has_long_counts <- all(c("count", "category") %in% names(long_data))

  # Posterior predictions: long-format custom families return a matrix of
  # category probabilities, while older wide-format rating fits return the
  # brms multinomial array [draw x obs x category].
  yrep <- do.call(brms::posterior_predict, c(list(fit), predict_args, list(...)))

  # Condition columns (stimulus and grouping vars are already excluded)
  conditions  <- .sdt_resolve_conditions(fit, NULL)
  cond_cols   <- names(conditions)

  # group adds an extra breakdown dimension (e.g., participant id, block).
  # .sdt_resolve_conditions() already excludes random-effect grouping vars, so
  # users can pass group = "id" to see per-subject proportion profiles.
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L) {
      stop2("'group' must be a single character string naming a column in the data.")
    }
    if (!group %in% names(fit$data)) {
      warning2("Ignoring 'group': column '", group, "' not found in model data.")
      group <- NULL
    } else if (!group %in% c(stim_var, cond_cols)) {
      cond_cols <- c(cond_cols, group)
    }
  }

  group_cols  <- c(stim_var, cond_cols)
  group_df    <- fit$data[, group_cols, drop = FALSE]
  unique_grps <- unique(group_df)
  rownames(unique_grps) <- NULL
  n_groups    <- nrow(unique_grps)

  obs_list  <- vector("list", n_groups)
  pred_list <- vector("list", n_groups)

  for (g in seq_len(n_groups)) {
    match_mask <- Reduce("&", lapply(group_cols, function(col) {
      group_df[[col]] == unique_grps[[col]][g]
    }))
    idx <- which(match_mask)

    if (has_long_counts) {
      total_cnt <- vapply(seq_len(n_ratings), function(cat_idx) {
        sum(long_data$count[idx][long_data$category[idx] == cat_idx])
      }, numeric(1))
    } else {
      resp_cols <- unlist(model$resp_vars)
      total_cnt <- colSums(as.matrix(long_data[idx, resp_cols, drop = FALSE]))
    }
    obs_props <- total_cnt / sum(total_cnt)

    obs_list[[g]] <- cbind(
      unique_grps[rep(g, n_ratings), , drop = FALSE],
      data.frame(category = seq_len(n_ratings), proportion = obs_props,
                 stringsAsFactors = FALSE)
    )

    if (has_long_counts) {
      yrep_sub <- yrep[, idx, drop = FALSE]
      weight_mat <- matrix(
        long_data$n_trials_total[idx],
        nrow = nrow(yrep_sub),
        ncol = length(idx),
        byrow = TRUE
      )
      yrep_weighted <- yrep_sub * weight_mat
      yrep_pool <- vapply(seq_len(n_ratings), function(cat_idx) {
        rowSums(yrep_weighted[, long_data$category[idx] == cat_idx, drop = FALSE])
      }, numeric(nrow(yrep_weighted)))
    } else {
      yrep_sub <- yrep[, idx, , drop = FALSE]
      yrep_pool <- apply(yrep_sub, c(1L, 3L), sum)
    }
    yrep_n <- rowSums(yrep_pool)
    yrep_props <- yrep_pool / yrep_n

    pred_list[[g]] <- cbind(
      unique_grps[rep(g, n_ratings), , drop = FALSE],
      data.frame(
        category   = seq_len(n_ratings),
        pred_mean  = colMeans(yrep_props),
        pred_lower = apply(yrep_props, 2L, stats::quantile, probs = probs[1L]),
        pred_upper = apply(yrep_props, 2L, stats::quantile, probs = probs[2L]),
        stringsAsFactors = FALSE
      )
    )
  }

  obs_df  <- do.call(rbind, obs_list)
  pred_df <- do.call(rbind, pred_list)
  plot_df <- merge(obs_df, pred_df, by = c(group_cols, "category"))
  rownames(plot_df) <- NULL

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["category"]])) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["pred_lower"]], ymax = .data[["pred_upper"]]),
      alpha = 0.30, fill = "#4393c3"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[["pred_mean"]]),
      colour = "#4393c3", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data[["proportion"]]),
      colour = "black", size = 2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[["proportion"]]),
      colour = "black", linewidth = 0.5, linetype = "dashed"
    )

  if (length(cond_cols) == 0L) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", stim_var)),
      labeller = ggplot2::label_both
    )
  } else {
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste(stim_var, "~",
                              paste(cond_cols, collapse = "+"))),
      labeller = ggplot2::label_both
    )
  }

  n_actual_draws <- dim(yrep)[1L]
  pct_lo <- round(probs[1L] * 100)
  pct_hi <- round(probs[2L] * 100)
  p +
    ggplot2::scale_x_continuous(breaks = seq_len(n_ratings)) +
    ggplot2::labs(
      x = "Rating category", y = "Proportion",
      subtitle = paste0(
        "Posterior predictive check (", n_actual_draws, " draws; band: ",
        pct_lo, "\u2013", pct_hi, "% CrI; points: observed)"
      )
    ) +
    ggplot2::theme_minimal()
}
