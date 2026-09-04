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
#' @param model A `bmmodel` object.
#' @param prep A `brmsprep` object from [brms::prepare_predictions()].
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
# vector
.pp_dpar_vector <- function(prep, name) {
  rep_len(as.vector(brms::get_dpar(prep, name)), prep$ndraws * prep$nobs)
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
#'   `label`, `default_type`, and `default`, flagging the observable that
#'   `pp_check()` plots when `resp_var` is not specified), or `NULL` invisibly
#'   for models without multi-observable support.
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
    default = names(spec$checks) == names(spec$observed)[spec$observed == "Y"],
    row.names = NULL
  )
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
  if (all_checks && !is.null(type)) {
    warning2("'type' is ignored for resp_var = 'all'; \\
              each panel uses its default type.")
    type <- NULL
  }
  checks <- if (all_checks) spec$checks else spec$checks[resp_var]
  group_vec <- if (!is.null(group)) object$data[[group]]
  plot_dots <- dots[setdiff(names(dots), c("draw_ids", "re_formula"))]

  plots <- lapply(checks, function(check) {
    .pp_build_ppc_plot(check, observed, yrep_inputs, type, group_vec, plot_dots)
  })
  if (all_checks) {
    bayesplot::bayesplot_grid(plots = unname(plots))
  } else {
    plots[[1L]]
  }
}

.pp_build_ppc_plot <- function(check, observed, yrep_inputs, type, group_vec,
                               plot_dots) {
  y <- check$compute(observed)
  yrep <- check$compute(yrep_inputs)

  # a simulated statistic can be undefined for some draws (rezdm() returns NA
  # for a boundary's mean RT when fewer than 2 responses reach it); y and yrep
  # must be reduced together or the two halves would be misaligned
  keep <- !is.na(y) & colSums(is.na(yrep)) == 0L
  stopif(!any(keep),
         "All observations of '{check$label}' are undefined in the \\
          posterior predictive simulation.")
  warnif(!all(keep),
         "Dropped {sum(!keep)} of {length(keep)} observations from the \\
          '{check$label}' check because the statistic was undefined for some \\
          posterior draws (e.g. too few simulated responses at a boundary).")
  y <- y[keep]
  yrep <- yrep[, keep, drop = FALSE]

  type <- type %||% check$type
  if (!is.null(group_vec)) {
    type <- .auto_grouped_type(type)
  }
  ppc_fun <- .ppc_fun(type)
  stopif(is.null(ppc_fun) || startsWith(type, "loo_"),
         "'{type}' is not a supported pp_check type for resp_var.")

  args <- c(list(y = y, yrep = yrep), plot_dots)
  if ("group" %in% names(formals(ppc_fun))) {
    stopif(is.null(group_vec), "Argument 'group' is required for type '{type}'.")
    args$group <- group_vec[keep]
  }
  do.call(ppc_fun, args) + ggplot2::labs(subtitle = check$label)
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
