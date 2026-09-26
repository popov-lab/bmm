############################################################################# !
# PP_OBSERVABLES.R                                                        ####
# Model-declared observables for multi-observable pp_check() support.     ####
############################################################################# !

# brms::pp_check() only ever compares the brms `Y` slot; observables carried by
# addition terms (vreal, vint, trials, dec) are invisible to it. Models declare
# them here so that pp_check(fit, resp_var = ...) can build y and yrep from one
# brms::prepare_predictions() object and one joint simulation.

#' Declare the observables of a bmm model for posterior predictive checks
#'
#' `pp_observables()` returns the model's observable declaration used by
#' [pp_check.bmmfit()] when `resp_var` is specified, or `NULL` for models
#' that delegate fully to [brms::pp_check()]. `pp_simulate()` draws all
#' observables jointly from the posterior predictive distribution.
#'
#' A `pp_observables()` method returns `list(observed, checks)`:
#' * `observed`: named character vector mapping observable names to slots of
#'   the brms standata (`"Y"`, `"vreal1"`, `"vint1"`, `"trials"`, `"dec"`).
#'   The observable mapped to `"Y"` is the default check.
#' * `checks`: named list of entries built by the internal `.pp_observable()`
#'   constructor. Each `compute` closure receives a named list keyed by
#'   `names(observed)` and must be elementwise, so the identical closure
#'   produces `y` from length-N vectors and `yrep` from ndraws x N matrices.
#'
#' A `pp_simulate()` method returns a named list of ndraws x nobs matrices
#' drawn jointly, typically through the internal `.pp_simulate_joint()` helper
#' around the model's `r*()` function. Simulating observables independently
#' would break their joint distribution (e.g. rt and response under the DDM).
#' Names not in `observed` are ignored; declared observables that are not
#' simulated (design quantities such as trial counts) are filled in from the
#' data.
#'
#' Register exactly one method per model at the most general class level
#' where the declaration is identical across versions.
#'
#' These two generics are exported so that model methods defined outside bmm
#' can be registered against them, but they are an internal developer
#' interface documented for bmm's own model authors and carry no stability
#' guarantee across releases.
#'
#' @param model A `bmmodel` object.
#' @param prep A `brmsprep` object from [brms::prepare_predictions()].
#' @return `pp_observables()` returns `NULL` for a model that delegates fully
#'   to [brms::pp_check()], or a list with elements `observed` (a named
#'   character vector mapping observable names to brms standata slots) and
#'   `checks` (a named list of check definitions, each with a `compute`
#'   closure, a `label` and a default plot `type`: a bayesplot `ppc_*` type
#'   or bmm's `"bars_binned"`). `pp_simulate()`
#'   returns a named list of `ndraws` x `nobs` matrices, one per simulated
#'   observable.
#' @keywords internal developer
#' @export
pp_observables <- function(model) {
  UseMethod("pp_observables")
}

#' @export
pp_observables.default <- function(model) {
  NULL
}

#' @rdname pp_observables
#' @export
pp_simulate <- function(model, prep) {
  UseMethod("pp_simulate")
}

#' @export
pp_simulate.default <- function(model, prep) {
  stop2("Model '{model$name}' declares pp_observables() but has no \\
         pp_simulate() method.")
}

.pp_observable <- function(compute, label, type = "dens_overlay") {
  nlist(compute, label, type)
}

# get_dpar() returns a scalar for dpars that brms stores fixed and an
# ndraws x nobs matrix otherwise; rep_len() flattens both to one column-major
# vector. Any other length would recycle into the wrong draw-major layout and
# attribute every simulated value to the wrong observation, silently.
.pp_dpar_vector <- function(prep, name) {
  v <- as.vector(brms::get_dpar(prep, name))
  stopif(!length(v) %in% c(1L, prep$ndraws * prep$nobs),
         "Cannot map dpar '{name}' (length {length(v)}) onto \\
          {prep$ndraws} draws x {prep$nobs} observations.")
  rep_len(v, prep$ndraws * prep$nobs)
}

# One RNG call over all draws and observations, reshaped column-major into
# ndraws x nobs matrices -- the layout .pp_expand_data() gives the data
.pp_simulate_joint <- function(prep, rfun, dpars, ...) {
  args <- c(list(n = prep$ndraws * prep$nobs),
            lapply(stats::setNames(dpars, dpars), .pp_dpar_vector, prep = prep),
            list(...))
  lapply(do.call(rfun, args), matrix, nrow = prep$ndraws)
}

# as.vector() of the result is rep(x, each = ndraws)
.pp_expand_data <- function(x, ndraws) {
  matrix(x, nrow = ndraws, ncol = length(x), byrow = TRUE)
}

#' List the posterior predictive checks available for a model fit
#'
#' Lists the values `resp_var` accepts in [pp_check.bmmfit()] for a fitted
#' model whose likelihood involves several observables.
#'
#' @param fit A `bmmfit` object returned by [bmm()].
#' @return A `data.frame` with one row per available check (columns `resp_var`,
#'   `label`, `default_type`, `slot`, listing the brms standata slots the check
#'   reads, and `default`, flagging the observable that `pp_check()` plots when
#'   `resp_var` is not specified), or `NULL` invisibly for models without
#'   multi-observable support.
#' @seealso [pp_check.bmmfit()]
#' @keywords extract_info
#' @examples
#' \dontrun{
#' fit <- bmm(bmmformula(drift ~ condition), data, ddm(rt = "rt", response = "response"))
#' pp_check_vars(fit)
#' pp_check(fit, resp_var = "response")
#' }
#' @export
pp_check_vars <- function(fit) {
  stopif(!inherits(fit, "bmmfit"), "'fit' must be a bmmfit object.")
  spec <- pp_observables(fit$bmm$model)
  if (is.null(spec)) {
    message2("Model '{fit$bmm$model$name}' declares no additional observables; \\
              pp_check() delegates to brms::pp_check().")
    return(invisible(NULL))
  }
  data.frame(
    resp_var = names(spec$checks),
    label = vapply(spec$checks, `[[`, character(1), "label"),
    default_type = vapply(spec$checks, `[[`, character(1), "type"),
    slot = vapply(spec$checks, .pp_check_slots, character(1),
                  observed = spec$observed),
    default = names(spec$checks) == names(spec$observed)[spec$observed == "Y"],
    row.names = NULL
  )
}

# Derived checks read several observables, so the slots are recovered from the
# closure rather than declared twice. Anchoring on "$" and a trailing word
# boundary keeps mean_rt from matching d$mean_rt_upper.
.pp_check_slots <- function(check, observed) {
  code <- paste(deparse(body(check$compute)), collapse = " ")
  reads <- vapply(names(observed), function(nm) {
    grepl(paste0("\\$", nm, "\\b"), code)
  }, logical(1))
  paste(observed[reads], collapse = ", ")
}

.pp_check_observable <- function(object, spec, resp_var, type, ndraws, group,
                                 dots) {
  if (is.null(dots$draw_ids)) {
    ndraws <- ndraws %||% 10L
  }
  prep <- brms::prepare_predictions(object, ndraws = ndraws,
                                    draw_ids = dots$draw_ids,
                                    re_formula = dots$re_formula)

  observed <- lapply(spec$observed, function(slot) prep$data[[slot]])
  yrep_inputs <- lapply(observed, .pp_expand_data, ndraws = prep$ndraws)
  sims <- pp_simulate(object$bmm$model, prep)
  sims <- sims[intersect(names(sims), names(spec$observed))]
  yrep_inputs[names(sims)] <- sims

  all_checks <- identical(resp_var, "all")
  checks <- if (all_checks) spec$checks else spec$checks[resp_var]
  plot_dots <- dots[setdiff(names(dots), c("draw_ids", "re_formula"))]

  reduced <- .pp_reduce_na(checks, observed, yrep_inputs)
  group_vec <- if (!is.null(group)) object$data[[group]][reduced$keep]

  plots <- lapply(names(checks), function(nm) {
    .pp_build_ppc_plot(checks[[nm]], reduced$values[[nm]], type, group_vec,
                       plot_dots)
  })
  if (all_checks) {
    bayesplot::bayesplot_grid(plots = unname(plots))
  } else {
    plots[[1L]]
  }
}

# A simulated statistic can be undefined (rezdm() returns NA for a boundary's
# mean RT when fewer than 2 responses reach it). Reducing the OBSERVATION
# dimension would make the retained count decay as (1 - p)^ndraws, so asking
# for more draws would check less data; only observations whose observed value
# is undefined are dropped, and the remaining NAs are absorbed by dropping
# exchangeable draws. All panels of resp_var = "all" share one reduction so
# that they are computed on the same observations and draws.
.pp_reduce_na <- function(checks, observed, yrep_inputs) {
  label <- collapse_comma(vapply(checks, `[[`, character(1), "label"))
  values <- lapply(checks, function(check) {
    list(y = check$compute(observed), yrep = check$compute(yrep_inputs))
  })

  keep <- Reduce(`&`, lapply(values, function(v) !is.na(v$y)))
  stopif(!any(keep), "All observations of {label} are undefined in the data.")
  warnif(!all(keep),
         "Dropped {sum(!keep)} of {length(keep)} observations because the \\
          observed {label} is undefined (too few responses at a boundary).")
  values <- lapply(values, function(v) {
    list(y = v$y[keep], yrep = v$yrep[, keep, drop = FALSE])
  })

  keep_draws <- Reduce(`&`, lapply(values, function(v) {
    rowSums(is.na(v$yrep)) == 0L
  }))
  stopif(!any(keep_draws),
         "Every posterior draw of {label} contains an undefined observation; \\
          try a model with more trials per cell.")
  warnif(!all(keep_draws),
         "Dropped {sum(!keep_draws)} of {length(keep_draws)} posterior draws \\
          because {label} was undefined for some observations.")
  values <- lapply(values, function(v) {
    list(y = v$y, yrep = v$yrep[keep_draws, , drop = FALSE])
  })
  nlist(values, keep)
}

.pp_build_ppc_plot <- function(check, value, type, group_vec, plot_dots) {
  type <- type %||% check$type
  if (!is.null(group_vec)) {
    type <- .auto_grouped_type(type)
  }
  ppc_fun <- .ppc_fun(type)
  stopif(is.null(ppc_fun) || startsWith(type, "loo_"),
         "'{type}' is not a supported pp_check type for resp_var.")

  args <- c(list(y = value$y, yrep = value$yrep), plot_dots)
  if ("group" %in% names(formals(ppc_fun))) {
    stopif(is.null(group_vec), "Argument 'group' is required for type '{type}'.")
    args$group <- group_vec
  }
  do.call(ppc_fun, args) + ggplot2::labs(subtitle = check$label)
}

# bmm's own plot type (see .ppc_fun()): ppc_bars() for a continuous statistic.
# bayesplot requires whole numbers there, so the statistic is binned and
# bayesplot counts the bins; the bars are then drawn on the statistic's scale
# instead of at bin indices.
.ppc_bars_binned <- function(y, yrep, ..., breaks = NULL, prob = 0.9,
                             freq = TRUE) {
  .pp_binned_bars(y, yrep, group = NULL, breaks, prob, freq)
}

.ppc_bars_binned_grouped <- function(y, yrep, group, ..., breaks = NULL,
                                     facet_args = list(), prob = 0.9,
                                     freq = TRUE) {
  facet_args$facets <- "group"
  facet_args$scales <- facet_args$scales %||% "free"
  .pp_binned_bars(y, yrep, group, breaks, prob, freq) +
    do.call(ggplot2::facet_wrap, facet_args)
}

# The default breaks span y and yrep together, so predicted mass outside the
# observed range gets a bin. User breaks that miss a value would put it in
# bin 0 or K + 1 and silently drop it from the plot.
.pp_binned_bars <- function(y, yrep, group, breaks, prob, freq) {
  values <- c(y, yrep)
  breaks <- breaks %||% pretty(range(values), n = ceiling(log2(length(y)) + 1))
  stopif(min(values) < min(breaks) || max(values) > max(breaks),
         "'breaks' must cover the observed and predicted values \\
          ({signif(min(values), 3)} to {signif(max(values), 3)}).")
  bin <- function(x) findInterval(x, breaks, rightmost.closed = TRUE)

  data <- bayesplot::ppc_bars_data(bin(y), matrix(bin(yrep), nrow = nrow(yrep)),
                                   group = group, prob = prob, freq = freq)
  data$lower <- breaks[data$x]
  data$upper <- breaks[data$x + 1L]
  data$x <- (data$lower + data$upper) / 2

  scheme <- bayesplot::color_scheme_get()
  ggplot2::ggplot(data, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper, ymin = 0,
                   ymax = .data$y_obs, fill = "y"),
      colour = scheme$light_highlight
    ) +
    ggplot2::geom_pointrange(
      ggplot2::aes(y = .data$m, ymin = .data$l, ymax = .data$h,
                   colour = "yrep"),
      size = 0.5, linewidth = 1
    ) +
    ggplot2::scale_fill_manual(NULL, values = c(y = scheme$light),
                               labels = expression(italic(y)),
                               guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_colour_manual(NULL, values = c(yrep = scheme$dark),
                                 labels = expression(italic(y)[rep]),
                                 guide = ggplot2::guide_legend(order = 2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = if (freq) "Count" else "Proportion") +
    bayesplot::bayesplot_theme_get()
}

# response is 0/1, so the sign flip is (2 * response - 1): sign(response)
# would zero out every lower-boundary rt
.pp_spec_rt_response <- function() {
  list(
    observed = c(rt = "Y", response = "dec"),
    checks = list(
      rt = .pp_observable(function(d) d$rt, label = "Response time"),
      response = .pp_observable(function(d) d$response,
                                label = "Response (0 = lower, 1 = upper)",
                                type = "bars"),
      signed_rt = .pp_observable(function(d) d$rt * (2 * d$response - 1),
                                 label = "Signed response time (lower = negative)")
    )
  )
}
