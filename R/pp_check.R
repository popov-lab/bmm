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
#' Some models describe several observables jointly (e.g. [ddm()]: response
#' times *and* responses; [ezdm()]: mean RT, RT variance *and* accuracy), but
#' brms only ever checks the primary response. For these models the `resp_var`
#' argument selects which observable to check; [pp_check_vars()] lists the
#' available checks. All selected observables are drawn from **one** joint
#' posterior predictive simulation, so `resp_var = "all"` panels are mutually
#' consistent.
#'
#' @param object A `bmmfit` object returned by [bmm()].
#' @param type Character. Type of pp_check. When `NULL` (default), resolves to
#'   `"dens_overlay"`, or to the selected observable's default type when
#'   `resp_var` is specified. When `group` is specified, the grouped variant
#'   (e.g., `"dens_overlay_grouped"`) is auto-selected if available.
#'   Multinomial models produce a response proportion profile regardless of
#'   the value supplied.
#' @param ndraws Integer. Number of posterior draws. Defaults to `100` for
#'   multinomial models and `10` when `resp_var` is specified; otherwise
#'   passed to [brms::pp_check()].
#' @param group Character. Optional grouping variable for faceting. For
#'   non-multinomial models, passed to [brms::pp_check()]; when specified, the
#'   grouped variant of `type` (e.g., `"dens_overlay_grouped"`) is auto-selected
#'   if available. For multinomial models, facets by the named predictor.
#' @param resp_var Character. For models that declare several observables,
#'   the name of the observable to check, or `"all"` for a panel of all
#'   available checks built from one shared simulation. See [pp_check_vars()]
#'   for the options of a fitted model. The default `NULL` checks the primary
#'   response via [brms::pp_check()]. For the RT models, passing
#'   `negative_rt = TRUE` (a [brms::posterior_predict()] argument) is
#'   redirected to `resp_var = "signed_rt"`, so that observed and predicted
#'   response times are both signed by the response.
#' @param ... Additional arguments. Without `resp_var`, forwarded to
#'   [brms::pp_check()], or for multinomial models to
#'   [brms::posterior_predict()] (`probs`, a numeric vector of length 2 with
#'   default `c(0.025, 0.975)`, sets the credible interval). With `resp_var`,
#'   `draw_ids` and `re_formula` go to [brms::prepare_predictions()] and the
#'   rest to the `bayesplot::ppc_*` function. `re_formula = NA` predicts at
#'   the population level on every path.
#' @return For multinomial models or when `resp_var` is specified, a `ggplot2`
#'   object (a `bayesplot_grid` for `resp_var = "all"`). For other models, the
#'   result of [brms::pp_check()].
#' @seealso [brms::pp_check()], [pp_check_vars()]
#' @aliases pp_check
#' @importFrom brms pp_check
#' @importFrom rlang .data
#' @export
pp_check.bmmfit <- function(object, type = NULL, ndraws = NULL,
                            group = NULL, resp_var = NULL, ...) {
  dots <- list(...)
  spec <- pp_observables(object$bmm$model)

  if (isTRUE(dots$negative_rt)) {
    # brms would compare signed predicted RTs against the unsigned observed Y
    # -- a silent mismatch that looks like severe misfit (#401)
    stopif(!is.null(resp_var) && !identical(resp_var, "signed_rt"),
           "'negative_rt = TRUE' cannot be combined with resp_var = '{resp_var}'.")
    stopif(is.null(spec) || !"signed_rt" %in% names(spec$checks),
           "'negative_rt = TRUE' is not supported in pp_check() for this model.")
    message2("negative_rt = TRUE: checking the 'signed_rt' observable \\
              (response times signed by the response).")
    resp_var <- "signed_rt"
  }

  if (!is.null(resp_var)) {
    stopif(is.null(spec),
           "'resp_var' is not supported for the {object$bmm$model$name}: \\
            it declares no additional observables.")
    valid <- c(names(spec$checks), "all")
    stopif(!is.character(resp_var) || length(resp_var) != 1L ||
             !resp_var %in% valid,
           "'resp_var' must be one of {collapse_comma(valid)}. \\
            See pp_check_vars(fit).")
    stopif(!is.null(group) && (!is.character(group) || length(group) != 1L ||
             !group %in% names(object$data)),
           "'group' must name a column of the model data.")
    stopif(!is.null(dots$newdata),
           "'newdata' is not supported when 'resp_var' is specified.")
    dots$negative_rt <- NULL
    return(.pp_check_observable(object, spec, resp_var, type, ndraws, group,
                                dots))
  }

  if (identical(family(object)$family, "multinomial")) {
    return(.pp_check_multinomial(object, type = type, ndraws = ndraws %||% 100L,
                                 group = group, ...))
  }

  # a bare NextMethod() re-supplies only the formals present in the original
  # call, so the resolved type is passed explicitly (#401)
  type <- type %||% "dens_overlay"
  if (!is.null(group)) {
    type <- .auto_grouped_type(type)
  }
  NextMethod(type = type)
}


.ppc_fun <- function(type) {
  name <- paste0("ppc_", type)
  if (name %in% as.character(bayesplot::available_ppc(""))) {
    get(name, asNamespace("bayesplot"))
  }
}


.auto_grouped_type <- function(type) {
  grouped <- paste0(type, "_grouped")
  if (endsWith(type, "_grouped") || is.null(.ppc_fun(grouped))) type else grouped
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
