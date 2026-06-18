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
#' @param group Character. Optional grouping variable for faceting. For
#'   non-multinomial models, passed to [brms::pp_check()]; when specified, the
#'   grouped variant of `type` (e.g., `"dens_overlay_grouped"`) is auto-selected
#'   if available. For multinomial models, facets by the named predictor. For
#'   the [sdt_rating()] model `group` defaults to the `stimulus` variable so the
#'   signal and noise rating distributions are checked separately; pass
#'   `group = NA` to pool them into a single profile. The [sdt_ranking()] model
#'   facets by rank position; for mixed set sizes, pass the set-size column as
#'   `group` so each facet has a homogeneous denominator. The custom-family
#'   models [sdt_binary()] and [sdt_mafc()] delegate to [brms::pp_check()]; for
#'   their aggregated counts an overlaid empirical CDF (`type = "ecdf_overlay"`,
#'   adding `group = "stimulus"` for [sdt_binary()]) is usually clearer than the
#'   default bar view.
#' @param ... Additional arguments forwarded to [brms::pp_check()] (non-multinomial)
#'   or to [brms::posterior_predict()] (multinomial). For multinomial models,
#'   `probs` (numeric vector of length 2, default `c(0.025, 0.975)`) controls the
#'   credible interval. Both model types accept `re_formula` (e.g.,
#'   `re_formula = NA` to predict at the population level, excluding random
#'   effects).
#' @return For multinomial models, a `ggplot2` object. For other models, the
#'   result of [brms::pp_check()].
#' @seealso [brms::pp_check()]
#' @aliases pp_check
#' @importFrom brms pp_check
#' @importFrom rlang .data
#' @export
pp_check.bmmfit <- function(object, type = "dens_overlay", ndraws = NULL,
                            group = NULL, ...) {
  if (identical(family(object)$family, "multinomial")) {
    group <- .pp_check_resolve_group(object, group)
    .pp_check_multinomial(object, type = type, ndraws = ndraws %||% 100L,
                          group = group, ...)
  } else {
    if (!is.null(group)) {
      type <- .auto_grouped_type(type)
    }
    NextMethod()
  }
}


# Resolve the multinomial pp_check grouping. sdt_rating defaults to faceting by
# stimulus (signal vs noise — the meaningful SDT check); group = NA is the
# explicit "pool everything" opt-out. Other multinomial models (e.g. m3) keep
# the user-supplied group (NULL = pool).
.pp_check_resolve_group <- function(object, group) {
  if (length(group) == 1L && is.na(group)) {
    return(NULL)
  }
  if (is.null(group) && inherits(object$bmm$model, "sdt_rating")) {
    return(object$bmm$model$other_vars$stimulus)
  }
  group
}


.auto_grouped_type <- function(type) {
  if (endsWith(type, "_grouped")) {
    type
  } else {
    grouped <- paste0(type, "_grouped")
    ppc_fn <- paste0("ppc_", grouped)
    if (exists(ppc_fn, where = asNamespace("bayesplot"), mode = "function")) {
      grouped
    } else {
      type
    }
  }
}


.pp_check_multinomial <- function(object, type = NULL, ndraws = 100L,
                                  probs = c(0.025, 0.975),
                                  group = NULL, draw_ids = NULL, ...) {
  if (!is.null(type) && type != "dens_overlay") {
    warning2("Argument 'type' is ignored for multinomial pp_check. ",
             "A response-proportion plot is always produced.")
  }

  model <- object$bmm$model
  resp_cols <- unlist(model$resp_vars)
  yrep <- brms::posterior_predict(object, draw_ids = draw_ids,
                                  ndraws = ndraws, ...)

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
  if (is.null(group)) {
    character(0)
  } else {
    stopif(
      !is.character(group) || length(group) != 1L,
      "'group' must be a single character string naming a column."
    )
    valid_cols <- .resolve_pp_conditions(fit)
    if (group %in% valid_cols) {
      group
    } else {
      warning2("Ignoring 'group': column '{group}' is not a predictor ",
               "variable in the model data.")
      character(0)
    }
  }
}


.build_grouping <- function(data, cond_cols) {
  if (length(cond_cols) == 0L) {
    list(
      group_df = NULL,
      unique_grps = data.frame(.placeholder = 1L),
      cond_cols = character(0)
    )
  } else {
    group_df <- data[, cond_cols, drop = FALSE]
    unique_grps <- unique(group_df)
    rownames(unique_grps) <- NULL
    list(group_df = group_df, unique_grps = unique_grps, cond_cols = cond_cols)
  }
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

  pred_cols <- pred_cols[!grepl("^(Idx_|n_)", pred_cols)]

  pred_cols
}
