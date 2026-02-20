############################################################################# !
# PP_CHECK.R                                                              ####
# pp_check() S3 method for bmmfit objects.                               ####
############################################################################# !

#' Posterior predictive check for bmmfit objects
#'
#' For models where brms provides `pp_check` support, this method delegates to
#' [brms::pp_check()]. For models with multinomial families (e.g., the m3
#' model), brms's `pp_check` is unavailable; this method dispatches to a
#' model-specific visualisation instead.
#'
#' For **multinomial models**, the plot mirrors the bayesplot `ppc_bars` style:
#' observed proportions are shown as bars and posterior predictive medians with
#' credible intervals are shown as point-ranges, using the bayesplot default
#' colour scheme and theme.
#'
#' @param object A `bmmfit` object returned by [bmm()].
#' @param type Character. Type of pp_check (default `"dens_overlay"`). For
#'   non-multinomial models, passed to [brms::pp_check()]. When `group` is
#'   specified, the grouped variant (e.g., `"dens_overlay_grouped"`) is
#'   auto-selected if available. Multinomial models produce a response
#'   proportion profile regardless of the value supplied.
#' @param ndraws Integer. Number of posterior draws. Defaults to `100` for
#'   multinomial models; otherwise passed to [brms::pp_check()].
#' @param ... For non-multinomial models: additional arguments forwarded to
#'   [brms::pp_check()]. For multinomial models: `probs` (numeric vector of
#'   length 2, default `c(0.025, 0.975)`) to control the credible interval, and
#'   `group` (single character string) to facet by a predictor variable.
#'   Both model types accept `re_formula` (e.g., `re_formula = NA` to
#'   predict at the population level, excluding random effects).
#' @return For multinomial models, a `ggplot2` object. For other models, the
#'   result of [brms::pp_check()].
#' @seealso [brms::pp_check()]
#' @importFrom brms pp_check
#' @importFrom rlang .data
#' @export
pp_check.bmmfit <- function(object, type = "dens_overlay", ndraws = NULL, ...) {
  if (identical(family(object)$family, "multinomial")) {
    return(.pp_check_multinomial(object, type = type, ndraws = ndraws, ...))
  }
  dots <- list(...)
  if (!is.null(dots$group)) {
    type <- .auto_grouped_type(type)
  }
  NextMethod()
}


.auto_grouped_type <- function(type) {
  if (endsWith(type, "_grouped")) return(type)
  grouped <- paste0(type, "_grouped")
  ppc_fn <- paste0("ppc_", grouped)
  if (exists(ppc_fn, where = asNamespace("bayesplot"), mode = "function")) {
    return(grouped)
  }
  type
}


.pp_check_multinomial <- function(object, type = NULL, ndraws = NULL,
                                  probs = c(0.025, 0.975),
                                  group = NULL,
                                  # brms pp_check args — intercepted, not forwarded
                                  x = NULL, prefix = c("ppc", "ppd"),
                                  draw_ids = NULL,
                                  nsamples = NULL, subset = NULL,
                                  ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2("ggplot2 is required for pp_check of multinomial models.")
  }

  if (!is.null(type) && type != "dens_overlay") {
    warning2("Argument 'type' is ignored for multinomial pp_check. ",
             "A response-proportion plot is always produced.")
  }

  if (!is.null(nsamples)) {
    warning2("'nsamples' is deprecated; use 'ndraws' instead.")
    if (is.null(draw_ids)) ndraws <- nsamples
  }

  if (is.null(ndraws)) ndraws <- 100L

  predict_args <- if (!is.null(draw_ids)) {
    list(draw_ids = draw_ids)
  } else {
    list(ndraws = ndraws)
  }

  model <- object$bmm$model
  resp_cols <- unlist(model$resp_vars)
  yrep <- do.call(brms::posterior_predict,
                  c(list(object), predict_args, list(...)))

  y_mat <- object$data$Y
  cat_names <- colnames(y_mat) %||% resp_cols

  cond_cols <- .validate_group_arg(group, object)
  intercept_only <- length(cond_cols) == 0L
  grouping <- .build_grouping(object$data, cond_cols)

  plot_df <- .aggregate_pp_data(
    y_mat, yrep, cat_names, grouping, probs, intercept_only
  )

  .build_pp_plot(plot_df, cond_cols, intercept_only, yrep, probs)
}


.validate_group_arg <- function(group, fit) {
  if (is.null(group)) return(character(0))

  if (!is.character(group) || length(group) != 1L) {
    stop2("'group' must be a single character string naming a column.")
  }

  valid_cols <- .resolve_pp_conditions(fit)
  if (!group %in% valid_cols) {
    warning2("Ignoring 'group': column '{group}' is not a predictor ",
             "variable in the model data.")
    return(character(0))
  }
  group
}


.build_grouping <- function(data, cond_cols) {
  if (length(cond_cols) == 0L) {
    return(list(
      group_df = NULL,
      unique_grps = data.frame(.placeholder = 1L),
      cond_cols = character(0)
    ))
  }
  group_df <- data[, cond_cols, drop = FALSE]
  unique_grps <- unique(group_df)
  rownames(unique_grps) <- NULL
  list(group_df = group_df, unique_grps = unique_grps, cond_cols = cond_cols)
}


.aggregate_pp_data <- function(y_mat, yrep, cat_names, grouping,
                               probs, intercept_only) {
  unique_grps <- grouping$unique_grps
  group_df <- grouping$group_df
  cond_cols <- grouping$cond_cols
  n_groups <- nrow(unique_grps)
  n_cats <- length(cat_names)
  obs_list <- vector("list", n_groups)
  pred_list <- vector("list", n_groups)

  for (g in seq_len(n_groups)) {
    idx <- if (intercept_only) {
      seq_len(nrow(y_mat))
    } else {
      which(Reduce("&", lapply(cond_cols, function(col) {
        group_df[[col]] == unique_grps[[col]][g]
      })))
    }

    total_cnt <- colSums(y_mat[idx, , drop = FALSE])
    obs_props <- total_cnt / sum(total_cnt)

    obs_row <- data.frame(
      category = factor(cat_names, levels = cat_names),
      proportion = as.numeric(obs_props),
      stringsAsFactors = FALSE
    )

    yrep_pool <- apply(yrep[, idx, , drop = FALSE], c(1L, 3L), sum)
    yrep_props <- yrep_pool / rowSums(yrep_pool)

    pred_row <- data.frame(
      category    = factor(cat_names, levels = cat_names),
      pred_median = apply(yrep_props, 2L, stats::median),
      pred_lower  = apply(yrep_props, 2L, stats::quantile, probs = probs[1L]),
      pred_upper  = apply(yrep_props, 2L, stats::quantile, probs = probs[2L]),
      stringsAsFactors = FALSE
    )

    if (!intercept_only) {
      grp_row <- unique_grps[rep(g, n_cats), , drop = FALSE]
      obs_row <- cbind(grp_row, obs_row)
      pred_row <- cbind(grp_row, pred_row)
    }

    obs_list[[g]] <- obs_row
    pred_list[[g]] <- pred_row
  }

  obs_df <- do.call(rbind, obs_list)
  pred_df <- do.call(rbind, pred_list)
  merge_by <- if (intercept_only) "category" else c(cond_cols, "category")
  plot_df <- merge(obs_df, pred_df, by = merge_by)
  rownames(plot_df) <- NULL
  plot_df
}


.build_pp_plot <- function(plot_df, cond_cols, intercept_only, yrep, probs) {
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["category"]])) +
    ggplot2::geom_col(
      ggplot2::aes(y = .data[["proportion"]]),
      fill = "#d1e1ec", colour = "#b3cde0", width = 0.9
    ) +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        y    = .data[["pred_median"]],
        ymin = .data[["pred_lower"]],
        ymax = .data[["pred_upper"]]
      ),
      colour = "#03396c", size = 0.5, fatten = 2.5, linewidth = 1
    )

  if (!intercept_only) {
    facet_formula <- stats::as.formula(
      paste("~", paste(cond_cols, collapse = " + "))
    )
    p <- p + ggplot2::facet_wrap(facet_formula, labeller = ggplot2::label_both)
  }

  n_actual_draws <- dim(yrep)[1L]
  pct_lo <- round(probs[1L] * 100)
  pct_hi <- round(probs[2L] * 100)

  p +
    ggplot2::labs(
      x = "Response category", y = "Proportion",
      subtitle = paste0(
        "Posterior predictive check (", n_actual_draws, " draws; ",
        pct_lo, "\u2013", pct_hi, "% CrI)"
      )
    ) +
    bayesplot::theme_default()
}


.resolve_pp_conditions <- function(fit) {
  data <- fit$data
  model <- fit$bmm$model
  resp_cols <- unlist(model$resp_vars)
  re_vars <- names(brms::ranef(fit))
  exclude <- unique(c(resp_cols, "nTrials", re_vars))

  pred_cols <- setdiff(names(data), exclude)

  is_matrix <- vapply(pred_cols, function(col) is.matrix(data[[col]]), logical(1))
  pred_cols <- pred_cols[!is_matrix]

  # Drop m3 infrastructure columns (Idx_*, n_*)
  pred_cols <- pred_cols[!grepl("^(Idx_|n_)", pred_cols)]

  pred_cols
}
