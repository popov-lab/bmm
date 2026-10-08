#' Conditional Effects for BMM Models
#'
#' @description
#' Compute conditional effects for parameters of a bmmfit object.
#' This method provides a more intuitive interface than directly calling
#' [brms::conditional_effects()] on bmmfit objects, by:
#' \itemize{
#'   \item Accepting model parameter names directly (e.g., `"kappa"`, `"thetat"`)
#'   \item Automatically determining whether parameters are distributional or non-linear
#'   \item Optionally applying inverse link transformations to show parameters on their natural scale
#' }
#'
#' @param x A bmmfit object (created by [bmm()])
#' @param par Character string. Name of the model parameter to compute effects for.
#'   This should be one of the parameter names from the original model specification
#'   (see `names(x$bmm$model$parameters)`). If `NULL` (the default), conditional
#'   effects are computed for all estimated (non-fixed) parameters.
#' @param scale Character. Scale on which to show the parameter:
#'   \describe{
#'     \item{`"native"` (default)}{Show on natural scale using inverse link transformation.
#'       For example, `kappa` with log link shown on exp scale, `thetat` with
#'       logit (mixture2p) or softmax (mixture3p) link shown on the probability scale.}
#'     \item{`"sampling"`}{Show on the sampling scale (as used during MCMC).
#'       For example, `kappa` with log link shown on log scale.}
#'   }
#' @param ... Additional arguments passed to [brms::conditional_effects()].
#'   Common arguments include:
#'   \itemize{
#'     \item `effects`: Character vector specifying which predictor effects to plot
#'     \item `conditions`: Named list for setting values of covariates
#'     \item `int_conditions`: Conditions for interactions
#'     \item `prob`: Probability mass to include in credible intervals (default 0.95)
#'     \item `spaghetti`: Logical, whether to add spaghetti lines
#'     \item `method`: Method for computing effects ("posterior_predict" or
#'       "posterior_epred"; "posterior_predict" is not available for models with
#'       a multinomial family, see Details)
#'   }
#'
#' @return A `brms_conditional_effects` object (from brms), which can be:
#' \itemize{
#'   \item Plotted directly using [plot()]
#'   \item Converted to a data frame for custom plotting
#'   \item Combined with other conditional effects plots
#' }
#'
#' @details
#' ## Parameter Types
#'
#' bmm models use two types of parameters internally:
#' \itemize{
#'   \item **Non-linear parameters (`nlpar`)**: Model parameters that enter the
#'     likelihood through a non-linear formula, like `kappa` and `thetat` of
#'     [mixture3p()] or `c` and `a` of [m3()]
#'   \item **Distributional parameters (`dpar`)**: Parameters of the response
#'     distribution itself, like `kappa` and `c` of [sdm()]
#' }
#'
#' Users should not need to know this distinction - `conditional_effects.bmmfit()`
#' automatically routes to the correct parameter type.
#'
#' ## Scale Transformations
#'
#' By default (`scale = "native"`), parameters are shown on their natural scale by
#' applying inverse link transformations:
#' \itemize{
#'   \item `log` link → exp transformation
#'   \item `logit` link → inverse logit (probability scale)
#'   \item `tan_half` link → 2*atan transformation (radians)
#'   \item `identity` link → no transformation
#' }
#'
#' Use `scale = "sampling"` to see parameters on the scale used during MCMC sampling.
#'
#' `estimate__`, `lower__`, `upper__` and the spaghetti lines follow `scale`.
#' `se__` does not: it is the spread of the draws on the scale `brms` or bmm
#' summarises them on, so for a link other than identity it can sit on a
#' different scale than `estimate__`. For example, `se__` of `kappa` is on the
#' sampling scale in [mixture3p()] and on the native scale in [sdm()], at either
#' `scale`. Use `lower__` and `upper__` to describe uncertainty on the requested
#' scale. With `robust = FALSE`, `estimate__` is the mean and `se__` the standard
#' deviation of the draws on that scale, and `estimate__` is then mapped to the
#' requested scale: at `scale = "native"` the estimate of a non-linear parameter
#' with a log link is the exponential of the mean log value, not the mean of the
#' exponentiated draws.
#'
#' The default of `robust` is `TRUE` (median and MAD), as in [brms::conditional_effects()],
#' with one exception: `thetat` and `thetant` of [mixture3p()] at
#' `scale = "native"` default to `robust = FALSE` (mean and standard deviation).
#'
#' ## Models with a multinomial family
#'
#' In models fitted with a multinomial response (`m3()`, `sdt_rating()`,
#' `sdt_cdp()`, `sdt_ranking()`), the parameters are latent quantities: they
#' have no predictive distribution, no response categories and no response
#' points. For these models:
#' \itemize{
#'   \item The effect is the posterior of the parameter itself, summarised by
#'     the median and MAD of the draws; use `robust = FALSE` for their mean and
#'     standard deviation.
#'   \item The grid of the plotted predictor and the values of all other
#'     predictors are set by the rules of [brms::conditional_effects()].
#'   \item Without `effects`, only the effects `brms` lists for the parameter's own
#'     formula are returned. A parameter without predictors has no effects.
#'   \item `categorical`, `ordinal`, `select_points` and `transform` are
#'     errors, and `method` must be `"posterior_epred"` or `"posterior_linpred"`
#'     (or an alias [brms::conditional_effects()] reads as one of them, such as
#'     `"fitted"`), which agree for a parameter.
#' }
#'
#' @seealso [brms::conditional_effects()] for the underlying brms function
#'
#' @aliases conditional_effects
#' @method conditional_effects bmmfit
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a mixture model with set size effect on kappa
#' fit <- bmm(
#'   formula = bmf(kappa ~ 0 + setsize, thetat ~ 1),
#'   data = zhang_luck_2008,
#'   model = mixture3p(
#'     resp_error = "response_error",
#'     nt_features = paste0("col_lure", 1:5),
#'     set_size = "setsize"
#'   )
#' )
#'
#' # Get conditional effects for kappa on natural scale (exp of log)
#' ce_kappa <- conditional_effects(fit, par = "kappa")
#' plot(ce_kappa)
#'
#' # Get conditional effects for kappa on log scale (sampling scale)
#' ce_kappa_log <- conditional_effects(fit, par = "kappa", scale = "sampling")
#' plot(ce_kappa_log)
#'
#' # Get effects for thetat (memory probability)
#' ce_thetat <- conditional_effects(fit, par = "thetat")
#' plot(ce_thetat)
#'
#' # Specify which effects to plot
#' ce_specific <- conditional_effects(fit, par = "kappa", effects = "setsize")
#'
#' # Combine with other brms options
#' ce_detailed <- conditional_effects(
#'   fit,
#'   par = "kappa",
#'   effects = "setsize",
#'   spaghetti = TRUE,
#'   ndraws = 100
#' )
#' }
conditional_effects.bmmfit <- function(x,
                                       par = NULL,
                                       scale = c("native", "sampling"),
                                       ...) {
  x <- restructure(x)
  scale <- match.arg(scale)
  stopif(
    !is.null(par) && (!is.character(par) || length(par) != 1),
    "Argument 'par' must be a single character string"
  )
  if (.has_category_dpars(x)) {
    .ce_check_category_family_args(x, ...)
  }

  if (is.null(par)) {
    .ce_all_parameters(x, scale, ...)
  } else {
    .ce_single_parameter(x, par, scale, ...)
  }
}


#' Refuse the arguments a parameter of a multinomial-family model cannot honour
#'
#' @description
#' Every parameter of these models is a non-linear parameter, so the check
#' applies to all of them. `method` is read as brms reads it: its aliases of
#' "posterior_epred" and "posterior_linpred" are accepted, the predictive ones
#' are not.
#'
#' The arguments are checked here, on the path of the exported method, because
#' brms cannot check them: the routes ask brms for a category dpar, which accepts
#' them, and evaluate the parameter themselves.
#'
#' @param x A bmmfit object
#' @param ... Arguments of [brms::conditional_effects()]
#'
#' @return `NULL`, or an error that names what to change
#'
#' @keywords internal
#' @noRd
.ce_check_category_family_args <- function(x, ...) {
  args <- .ce_match_args(list(...))
  remove <- c(
    if (isTRUE(args[["categorical"]])) "categorical",
    if (isTRUE(args[["ordinal"]])) "ordinal",
    if (isTRUE(args[["select_points"]] > 0)) "select_points",
    if (!is.null(args[["transform"]])) "transform"
  )
  method_ok <- length(args[["method"]]) <= 1 &&
    all(args[["method"]] %in% c("posterior_epred", "fitted", "pp_expect",
                                "posterior_linpred"))
  fix <- c(
    if (length(remove) > 0) paste0("Remove ", collapse_comma(remove), "."),
    if (!method_ok) paste0("Set 'method' to \"posterior_epred\" or ",
                           "\"posterior_linpred\", which agree for a parameter.")
  )
  stopif(
    length(fix) > 0,
    "The parameters of the '{intersect(class(x$bmm$model), model_names())[1]}' \\
    model are latent quantities without predictive distribution, response \\
    categories or response points. {paste(fix, collapse = ' ')}"
  )
}


.ce_all_parameters <- function(x, scale, ...) {
  model <- x$bmm$model
  estimated_pars <- setdiff(names(model$parameters),
                            names(model$fixed_parameters))
  all_effects <- list()
  for (p in estimated_pars) {
    ce <- .ce_single_parameter(x, par = p, scale = scale, ...)
    if (length(ce) > 0) {
      names(ce) <- paste0(p, ".", names(ce))
      all_effects <- c(all_effects, ce)
    }
  }
  class(all_effects) <- c("brms_conditional_effects")
  all_effects
}


.ce_single_parameter <- function(x, par, scale, ...) {
  par_info <- .get_parameter_info(x, par)

  if (par_info$softmax && scale == "native") {
    softmax_result <- .compute_softmax_conditional_effects(x, par, ...)
    if (!is.null(softmax_result)) {
      .filter_internal_effects(softmax_result, x)
    } else {
      warning2(
        "Parameter '{par}' uses softmax transformation.\n",
        "Native scale display not available for this model configuration.\n",
        "Showing on sampling scale instead."
      )
      .ce_compute_and_transform(x, par, par_info, "sampling", ...)
    }
  } else {
    .ce_compute_and_transform(x, par, par_info, scale, ...)
  }
}


.ce_compute_and_transform <- function(x, par, par_info, scale, ...) {
  ce_result <- if (par_info$type == "dpar") {
    .brms_conditional_effects(x, dpar = par_info$brms_name, ...)
  } else if (par_info$type == "nlpar" && .has_category_dpars(x)) {
    .ce_nlpar_category_family(x, par, ...)
  } else if (par_info$type == "nlpar") {
    .brms_conditional_effects(x, nlpar = par_info$brms_name, ...)
  } else {
    stop2("Internal error: parameter type must be 'dpar' or 'nlpar'")
  }

  # nlpars: brms returns on sampling scale; dpars: on native scale
  if (par_info$link != "identity") {
    if (par_info$type == "nlpar" && scale == "native") {
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = TRUE)
    } else if (par_info$type == "dpar" && scale == "sampling") {
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = FALSE)
    }
  }

  .filter_internal_effects(ce_result, x)
}


#' Call brms conditional_effects without infinite recursion
#'
#' @description
#' Strips the `"bmmfit"` class so that S3 dispatch reaches
#' `brms::conditional_effects.brmsfit()` instead of recursing back to
#' `conditional_effects.bmmfit()`.
#'
#' @param x A bmmfit object
#' @param ... Arguments forwarded to [brms::conditional_effects()]
#'
#' @return A `brms_conditional_effects` object
#'
#' @keywords internal
#' @noRd
.brms_conditional_effects <- function(x, ...) {
  class(x) <- class(x)[class(x) != "bmmfit"]
  conditional_effects(x, ...)
}


#' Filter internal variables from conditional_effects results
#'
#' @description
#' Removes conditional effects plots for internal model variables
#' (like LureIdx, Idx_*, inv_ss, etc.) that are created during data
#' preprocessing but are not part of the user's formula.
#'
#' @param ce_result A brms_conditional_effects object
#' @param bmmfit A bmmfit object
#'
#' @return Filtered conditional_effects object with only user-specified predictors
#'
#' @keywords internal
#' @noRd
.filter_internal_effects <- function(ce_result, bmmfit) {
  internal_patterns <- c(
    "^LureIdx",
    "^Idx_",
    "^inv_ss$",
    "^Item[0-9]+_",
    "^expS$"
  )
  
  model <- bmmfit$bmm$model
  if (!is.null(model$other_vars$nt_features)) {
    nt_features <- model$other_vars$nt_features
    escaped <- gsub("([][(){}^$*+?.|\\\\])", "\\\\\\1", nt_features)
    internal_patterns <- c(internal_patterns, paste0("^", escaped, "$"))
  }
  if (!is.null(model$other_vars$nt_distances)) {
    nt_distances <- model$other_vars$nt_distances
    escaped <- gsub("([][(){}^$*+?.|\\\\])", "\\\\\\1", nt_distances)
    internal_patterns <- c(internal_patterns, paste0("^", escaped, "$"))
  }
  
  effect_names <- names(ce_result)
  combined_pattern <- paste(internal_patterns, collapse = "|")
  keep_effects <- !vapply(effect_names, function(name) {
    vars <- strsplit(name, ":")[[1]]
    any(grepl(combined_pattern, vars))
  }, logical(1))

  if (any(keep_effects)) {
    ce_result <- ce_result[keep_effects]
    class(ce_result) <- c("brms_conditional_effects")
  }

  ce_result
}


#' Extract grouping variable names from random effects in a formula
#'
#' @description
#' Parses the RHS of a formula to identify the data columns that define the
#' grouping levels of random effects, e.g. to exclude them from conditional
#' effects. Handles all brms grouping specifications:
#' \itemize{
#'   \item Bare names: `(1 | id)`, `(1 || id)`
#'   \item Correlation IDs: `(1 |p| id)` — returns `id` only; `p` labels the
#'     correlation structure and is not a data column (see
#'     `.extract_re_cor_ids()`)
#'   \item `gr()`: `(1 | gr(id, by = exp))` — extracts `id`, not `exp`
#'   \item `mm()`: `(1 | mm(g1, g2))` — extracts all positional args
#'   \item Interaction and nesting: `(1 | id:group)`, `(1 | id/group)` — extracts
#'     both `id` and `group`, the columns brms combines into the grouping factor
#' }
#'
#' @param formula A formula object
#'
#' @return Character vector of grouping variable names to exclude
#'
#' @keywords internal
#' @noRd
.extract_re_grouping_vars <- function(formula) {
  .re_bar_parts(formula[[length(formula)]])$groups
}


#' Extract correlation IDs from random effects in a formula
#'
#' @description
#' Returns the labels `p` of terms written `(1 |p| id)`. They tie the random
#' effects of several parameters into one correlation matrix and are not
#' columns of the data.
#'
#' @param formula A formula object
#'
#' @return Character vector of correlation IDs
#'
#' @keywords internal
#' @noRd
.extract_re_cor_ids <- function(formula) {
  .re_bar_parts(formula[[length(formula)]])$cor_ids
}


# `1 |p| g` parses as `(1 | p) | g`, so a bar whose left side is itself a bar
# carries a correlation ID. `group_cols` is every column the grouping term
# names, including `by =` and `weights =` arguments of gr() and mm().
.re_bar_parts <- function(expr) {
  none <- list(groups = character(0), cor_ids = character(0), group_cols = character(0))
  if (!is.call(expr)) {
    return(none)
  }
  if (.is_bar_call(expr)) {
    return(list(
      groups = .re_group_names(expr[[3]]),
      cor_ids = if (.is_bar_call(expr[[2]])) deparse(expr[[2]][[3]]) else character(0),
      group_cols = all.vars(expr[[3]])
    ))
  }
  parts <- lapply(seq_along(expr)[-1], function(i) {
    if (rlang::is_missing(expr[[i]])) none else .re_bar_parts(expr[[i]])
  })
  fields <- names(none)
  stats::setNames(
    lapply(fields, function(field) unique(as.character(unlist(lapply(parts, `[[`, field))))),
    fields
  )
}


.is_bar_call <- function(expr) {
  is.call(expr) && is.name(expr[[1]]) && as.character(expr[[1]]) %in% c("|", "||")
}


# `(1 |p| g)` becomes `(1 | g)`: the data columns of a term exclude its
# correlation ID.
.drop_re_cor_ids <- function(expr) {
  if (!is.call(expr)) {
    return(expr)
  }
  if (.is_bar_call(expr) && .is_bar_call(expr[[2]])) {
    expr[[2]] <- expr[[2]][[2]]
  }
  for (i in seq_along(expr)[-1]) {
    if (!rlang::is_missing(expr[[i]])) {
      expr[[i]] <- .drop_re_cor_ids(expr[[i]])
    }
  }
  expr
}


.re_group_names <- function(group) {
  if (is.name(group)) {
    return(as.character(group))
  }
  fun <- if (is.call(group) && is.name(group[[1]])) as.character(group[[1]]) else ""
  args <- if (is.call(group)) as.list(group)[-1]
  unnamed <- args[!nzchar(names(args) %||% character(length(args)))]
  unique(unlist(
    if (identical(fun, "gr")) {
      # by = and cor = may precede the grouping variable, which brms takes
      # from `group =` or else the first unnamed argument
      .re_group_names(if ("group" %in% names(args)) args$group else unnamed[[1]])
    } else if (identical(fun, "mm")) {
      lapply(unnamed, .re_group_names)
    } else if (fun %in% c(":", "/")) {
      lapply(args, .re_group_names)
    } else {
      all.vars(group)
    }
  ))
}


#' Summarize posterior draws into conditional-effect statistics
#'
#' @description
#' Takes a draws matrix (n_draws x n_grid_points) and computes summary
#' statistics suitable for `brms_conditional_effects` data frames.
#'
#' @param draws Matrix. Posterior draws (rows = draws, columns = grid points).
#' @param prob Numeric. Probability mass for credible intervals (default 0.95).
#' @param robust Logical. If `TRUE`, use median/MAD instead of mean/SD.
#' @param probs Numeric vector of the two quantiles bounding the interval. If
#'   `NULL`, derived from `prob`.
#'
#' @return A list with elements `estimate`, `lower`, `upper`, `se` — each a
#'   numeric vector of length `ncol(draws)`.
#'
#' @keywords internal
#' @noRd
.ce_summarize_draws <- function(draws, prob = 0.95, robust = FALSE,
                                probs = NULL) {
  probs <- probs %||% c((1 - prob) / 2, 1 - (1 - prob) / 2)
  if (robust) {
    estimate <- apply(draws, 2, stats::median)
    se <- apply(draws, 2, stats::mad)
  } else {
    estimate <- colMeans(draws)
    se <- apply(draws, 2, stats::sd)
  }
  lower <- apply(draws, 2, stats::quantile, probs = probs[1])
  upper <- apply(draws, 2, stats::quantile, probs = probs[2])
  list(estimate = estimate, lower = lower, upper = upper, se = se)
}


#' Compute softmax transformation for multinomial parameters
#'
#' @description
#' For mixture models with multinomial logit (softmax), manually computes
#' the softmax transformation by extracting conditional effects for all
#' relevant nlpars and applying the softmax formula.
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name to return (thetat or thetant)
#' @param ... Arguments of [brms::conditional_effects()], see `.ce_split_args()`
#'
#' @return A brms_conditional_effects object with softmax-transformed values
#'
#' @keywords internal
#' @noRd
.compute_softmax_conditional_effects <- function(bmmfit, par, ...) {
  if (!"mixture3p" %in% class(bmmfit$bmm$model)) {
    NULL
  } else {
    .compute_softmax_ce_inner(bmmfit, par, ...)
  }
}


.compute_softmax_ce_inner <- function(bmmfit, par, ...) {
  args <- .ce_split_args(...)
  grid <- do.call(.brms_conditional_effects, c(list(bmmfit, nlpar = par), args$grid),
                  quote = TRUE)

  softmax <- function(linpred) {
    draws_t <- linpred("thetat")
    draws_nt <- linpred("thetant")
    # numerically stable softmax: subtract max before exponentiating
    shift <- pmax(draws_t, draws_nt, 0)
    exp_t <- exp(draws_t - shift)
    exp_nt <- exp(draws_nt - shift)
    denom <- exp_t + exp_nt + exp(-shift)
    if (par == "thetat") exp_t / denom else exp_nt / denom
  }
  structure(
    .ce_summarise_grid(bmmfit, grid, args, robust = FALSE, draws_of = softmax),
    class = class(grid)
  )
}


#' Does the fit's family have one distributional parameter per category?
#'
#' @description
#' These are the families for which [brms::conditional_effects()] refuses a
#' request for a non-linear parameter and asks for `categorical = TRUE`
#' (`conv_cats_dpars()` in brms 2.23.0): the categorical, multinomial and
#' simplex-valued families. bmm's **m3** and the **sdt_rating**, **sdt_cdp** and
#' **sdt_ranking** models use `brms::multinomial()`.
#'
#' @param x A bmmfit object
#'
#' @return Logical scalar
#'
#' @keywords internal
#' @noRd
.has_category_dpars <- function(x) {
  x$family$family %in% c("categorical", "dirichlet", "dirichlet_multinomial",
                         "dirichlet2", "logistic_normal", "multinomial")
}


#' Formals of brms's conditional_effects method
#'
#' @return Pairlist, as [formals()]
#'
#' @keywords internal
#' @noRd
.brms_ce_formals <- function() {
  formals(utils::getS3method("conditional_effects", "brmsfit"))
}


#' Argument defaults of brms's conditional_effects method
#'
#' @param arg Character string. Name of an argument of
#'   `conditional_effects.brmsfit()`.
#'
#' @return The default of `arg`, as brms declares it
#'
#' @keywords internal
#' @noRd
.brms_ce_default <- function(arg) {
  .brms_ce_formals()[[arg]]
}


#' Match the arguments of a call to conditional_effects() to brms's formals
#'
#' @description
#' The routes of a multinomial family and of the softmax take their arguments
#' from `...`. R matches the arguments of a call to the formals of a function
#' by exact name, then by partial name, then by position, and refuses a partial
#' name that fits several formals. brms does this in stages: first for the
#' formals of its conditional_effects method, then, among what remains, for
#' `ndraws` and `draw_ids` of the prediction method, and then for the deprecated
#' `nsamples` and `subset`. The routes need the same reading of a call, so R's
#' [match.call()] does it here, stage by stage, so that `su` is `surface` and
#' `n` is `ndraws` as in brms.
#'
#' The positions count from `effects`, the first formal of brms's method after
#' `x`: `par` and `scale` are the formals of bmm's own method.
#'
#' @param args List of the arguments in `...`, named or not
#'
#' @return `args` with the formal names of exact, partial and positional
#'   matches; a name that matches none is kept
#'
#' @keywords internal
#' @noRd
.ce_match_args <- function(args) {
  stages <- list(
    .brms_ce_formals()[-1],
    alist(ndraws = , draw_ids = , ... = ),
    alist(nsamples = , subset = , ... = )
  )
  matched <- list()
  for (formals_ in stages) {
    stub <- function() NULL
    formals(stub) <- formals_
    args <- as.list(match.call(stub, as.call(c(quote(stub), args))))[-1]
    own <- names(args) %in% setdiff(names(formals_), "...")
    matched <- c(matched, args[own])
    args <- args[!own]
  }
  c(matched, args)
}


#' Split the arguments of conditional_effects() between brms's grid and draws
#'
#' @description
#' The routes of a multinomial family and of the softmax first let brms build
#' the conditions grid, then evaluate the parameter on that grid with
#' [brms::posterior_linpred()]. An argument belongs to one of them:
#' \itemize{
#'   \item The arguments of [brms::conditional_effects()] (`conditions`,
#'     `resolution`, `spaghetti`, `prob`, ...) reach the grid call, where brms
#'     validates them. The grid does not depend on the draws, so brms builds it
#'     at one draw.
#'   \item `ndraws` and `draw_ids` (and their deprecated aliases `nsamples` and
#'     `subset`) are resolved once per call by `.ce_resolve_draw_ids()`, so that
#'     every evaluation of one parameter's call, e.g. `thetat` and `thetant`
#'     of a softmax, uses the same draws. A call for several parameters
#'     (`par = NULL`) makes this call once per parameter, and `ndraws` then
#'     draws new ids for each of them.
#'   \item The remaining arguments prepare the draws (`sample_new_levels`, ...)
#'     and only reach [brms::posterior_linpred()].
#' }
#' `re_formula` reaches all of them, with the default of
#' [brms::conditional_effects()] if the call has none.
#'
#' @param ... Arguments of [brms::conditional_effects()]
#'
#' @return List with the argument lists for the `grid` call and for the
#'   `linpred` calls, the `effects` of the call, its `draws` arguments and the
#'   `summary` arguments (`prob`, `probs`, `robust`) that were given
#'
#' @keywords internal
#' @noRd
.ce_split_args <- function(...) {
  args <- .ce_match_args(list(...))
  draw_args <- c("ndraws", "draw_ids", "nsamples", "subset")
  re_formula <- if ("re_formula" %in% names(args)) {
    args["re_formula"]
  } else {
    list(re_formula = .brms_ce_default("re_formula"))
  }
  rest <- args[!names(args) %in% c("re_formula", draw_args)]

  list(
    grid = c(rest, re_formula, list(draw_ids = 1)),
    linpred = c(
      rest[!names(rest) %in% names(.brms_ce_formals())],
      re_formula,
      list(allow_new_levels = TRUE)
    ),
    effects = args[["effects"]],
    draws = args[names(args) %in% draw_args],
    summary = list(
      prob = args[["prob"]] %||% .brms_ce_default("prob"),
      probs = args[["probs"]],
      robust = args[["robust"]]
    )
  )
}


#' Resolve and validate the draws a call selects
#'
#' @description
#' Reads `ndraws` and `draw_ids` in the sequence of brms's own check and with
#' its messages, and draws `ndraws` ids if no `draw_ids` are given. brms
#' reports a call without draws only when it summarises the draws, so the
#' message of `posterior_summary()` stands in for that step.
#'
#' The draws are validated here and not by brms because the routes of a
#' multinomial family and of the softmax resolve them themselves, once per call,
#' and hand brms only the ids; this runs on the path of the exported method.
#'
#' @param x A bmmfit object
#' @param ndraws,draw_ids As in [brms::posterior_linpred()]
#'
#' @return Integer vector of draw ids, or `NULL` if the call selects none
#'
#' @keywords internal
#' @noRd
.validate_draw_ids <- function(x, ndraws = NULL, draw_ids = NULL) {
  ndraws_total <- brms::ndraws(x)
  if (is.null(draw_ids) && !is.null(ndraws)) {
    ndraws <- as_one_integer(ndraws)
    stopif(
      ndraws < 1 || ndraws > ndraws_total,
      "Argument 'ndraws' should be between 1 and the maximum number of draws \\
      ({ndraws_total})."
    )
    draw_ids <- sample.int(ndraws_total, ndraws)
  }
  if (!is.null(draw_ids)) {
    draw_ids <- as.integer(draw_ids)
    stopif(length(draw_ids) == 0L, "No posterior draws supplied.")
    stopif(
      any(draw_ids < 1L) || any(draw_ids > ndraws_total),
      "Some 'draw_ids' indices are out of range."
    )
  }
  draw_ids
}


#' Resolve the draws of a call into the draw ids to evaluate
#'
#' @description
#' Reads `ndraws` and `draw_ids` as brms does: a deprecated alias that is given
#' replaces its argument with a warning, and `draw_ids` wins over `ndraws`.
#'
#' @param fit A bmmfit object
#' @param draws The `draws` of `.ce_split_args()`
#'
#' @return Integer vector of draw ids, or `NULL` for all draws
#'
#' @keywords internal
#' @noRd
.ce_resolve_draw_ids <- function(fit, draws) {
  aliases <- c(nsamples = "ndraws", subset = "draw_ids")
  for (alias in names(aliases)) {
    arg <- aliases[[alias]]
    if (!is.null(draws[[alias]])) {
      warning2("Argument '{alias}' is deprecated. Please use argument '{arg}' instead.")
      draws[[arg]] <- draws[[alias]]
    }
  }
  .validate_draw_ids(fit, draws[["ndraws"]], draws[["draw_ids"]])
}


#' Summarise a parameter on the grid of a conditional effect
#'
#' @description
#' Evaluates the parameter on each element of brms's conditions grid with
#' [brms::posterior_linpred()], on the draws the call selected, and replaces the
#' summary columns and spaghetti lines of the element by those of the draws.
#'
#' @param x A bmmfit object
#' @param grid A `brms_conditional_effects` object, the conditions grid
#' @param args The list of `.ce_split_args()`
#' @param robust Default of `robust` on this route
#' @param draws_of Function of `linpred`, a function of the name of an `nlpar`
#'   that returns its draws on the grid element, returning the draws of the
#'   parameter
#'
#' @return Named list with one element per element of `grid`
#'
#' @keywords internal
#' @noRd
.ce_summarise_grid <- function(x, grid, args, robust, draws_of) {
  linpred_args <- c(
    args$linpred, list(draw_ids = .ce_resolve_draw_ids(x, args$draws))
  )
  lapply(grid, function(ce) {
    newdata <- .ce_cond_data(ce)
    draws <- draws_of(function(nlpar) {
      do.call(brms::posterior_linpred, c(
        list(x, newdata = newdata, nlpar = nlpar), linpred_args
      ), quote = TRUE)
    })
    .ce_set_summary(ce, draws, args$summary, robust)
  })
}


#' Conditional effects of a non-linear parameter in a family with category dpars
#'
#' @description
#' brms computes the predictions of an `nlpar` for these families and only then
#' refuses them, because it expects the categories to be plotted. A non-linear
#' parameter does not depend on the categories, so brms is asked for one of the
#' category dpars instead, at a single draw, which gives the conditions grid
#' (see `.ce_split_args()`). The parameter is then evaluated on that grid and
#' summarised as brms would.
#'
#' Without `effects`, brms is asked only for those of its default effects that
#' it would list for the parameter's own formula. The category dpar is evaluated
#' on every element of the grid, and on the grid of another formula's effect it
#' can leave its domain: brms holds a numeric moderator at mean ± sd, which is
#' below 0 for a 0/1 indicator or a count that is mostly 0, and the dpar takes
#' its log. A parameter without such effects needs no grid.
#'
#' @param x A bmmfit object
#' @param par Character string. The bmm parameter name, which is also its brms
#'   name.
#' @param ... Arguments of [brms::conditional_effects()]
#'
#' @return A `brms_conditional_effects` object
#'
#' @keywords internal
#' @noRd
.ce_nlpar_category_family <- function(x, par, ...) {
  args <- .ce_split_args(...)
  if (is.null(args$effects)) {
    args$grid$effects <- vapply(
      .ce_select_own_effects(.ce_brms_default_effects(x, args$grid$re_formula),
                             .ce_own_effects(x, par)),
      paste, "", collapse = ":"
    )
  }
  # an empty `effects` of the user's is for brms to refuse
  grid <- if (is.null(args$effects) && length(args$grid$effects) == 0) {
    structure(list(), names = character(0), class = "brms_conditional_effects")
  } else {
    .ce_category_grid(x, args$grid)
  }

  structure(
    lapply(
      .ce_summarise_grid(x, grid, args, robust = TRUE,
                         draws_of = function(linpred) linpred(par)),
      function(ce) {
        attr(ce, "response") <- par
        # a parameter has no observed response, but plot(points = TRUE) needs a
        # frame with a numeric `resp__`
        attr(ce, "points") <- attr(ce, "points")[0, ]
        ce
      }
    ),
    class = "brms_conditional_effects"
  )
}


#' brms's conditions grid of the first category dpar
#'
#' @param x A bmmfit object
#' @param grid_args The `grid` arguments of `.ce_split_args()`
#'
#' @return A `brms_conditional_effects` object
#'
#' @keywords internal
#' @noRd
.ce_category_grid <- function(x, grid_args) {
  # the parameter does not depend on the trials that brms sets to 1 for the grid
  withCallingHandlers(
    do.call(.brms_conditional_effects, c(
      list(x, dpar = names(brms::brmsterms(x$formula)$dpars)[1]), grid_args
    ), quote = TRUE),
    message = function(m) {
      if (grepl("Setting all 'trials' variables", conditionMessage(m), fixed = TRUE)) {
        invokeRestart("muffleMessage")
      }
    }
  )
}


#' The effects brms plots without `effects`, in its order and orientation
#'
#' @description
#' The list `brms::conditional_effects()` builds when it is given no `effects`
#' (brms 2.23.0): the effects of every formula of the model, each pair with a
#' numeric variable before a factor or grouping variable, which brms does only
#' for this list. An `effects` argument is plotted in the order it is written,
#' so these are the strings that reproduce brms's default grid.
#'
#' @param x A bmmfit object
#' @param re_formula As in [brms::conditional_effects()]
#'
#' @return List of character vectors of one or two variable names
#'
#' @keywords internal
#' @noRd
.ce_brms_default_effects <- function(x, re_formula) {
  brms_internal <- function(name) utils::getFromNamespace(name, "brms")
  bterms <- brms::brmsterms(
    brms_internal("update_re_terms")(x$formula, re_formula = re_formula)
  )
  group_vars <- brms_internal("get_group_vars")(bterms)
  is_like_factor <- brms_internal("is_like_factor")
  lapply(
    brms_internal("get_all_effects")(bterms, rsv_vars = brms_internal("rsv_vars")(bterms)),
    function(effect) {
      effect[order(vapply(x$data[effect], is_like_factor, logical(1)) |
                     effect %in% group_vars)]
    }
  )
}


#' The conditions grid of a conditional effect, without its summary columns
#'
#' @param ce One element of a `brms_conditional_effects` object
#'
#' @return A data frame
#'
#' @keywords internal
#' @noRd
.ce_cond_data <- function(ce) {
  ce[setdiff(names(ce), c("estimate__", "se__", "lower__", "upper__"))]
}


#' Put the summary of the draws of a parameter into a conditional effect
#'
#' @param ce One element of a `brms_conditional_effects` object
#' @param draws Matrix of draws (rows = draws, columns = grid rows)
#' @param summary List with `prob`, `probs` and `robust`, see `.ce_split_args()`
#' @param robust Used if `summary$robust` is `NULL`
#'
#' @return `ce`, its summary columns and its spaghetti lines (if it has any)
#'   replaced by those of `draws`
#'
#' @keywords internal
#' @noRd
.ce_set_summary <- function(ce, draws, summary, robust) {
  summ <- .ce_summarize_draws(draws, prob = summary$prob,
                              robust = summary$robust %||% robust,
                              probs = summary$probs)
  ce$estimate__ <- summ$estimate
  ce$se__ <- summ$se
  ce$lower__ <- summ$lower
  ce$upper__ <- summ$upper
  if (!is.null(attr(ce, "spaghetti"))) {
    attr(ce, "spaghetti") <- .ce_spaghetti(.ce_cond_data(ce), draws,
                                           attr(ce, "effects"))
  }
  ce
}


#' Keep the effects of brms's default grid that are the parameter's own effects
#'
#' @description
#' An effect that two formulas list with its variables in another order
#' (`x:z` and `z:x`) is one effect; it is kept in the order the parameter's own
#' formula lists it, or else the first one. brms plots the first variable on the
#' x-axis and refuses an effect given twice.
#'
#' @param effects List of the effects of brms's default grid, see
#'   `.ce_brms_default_effects()`
#' @param own_effects List of the variable sets of the parameter's effects
#'
#' @return `effects` without the other elements
#'
#' @keywords internal
#' @noRd
.ce_select_own_effects <- function(effects, own_effects) {
  own <- vapply(effects, function(e) any(vapply(own_effects, setequal, logical(1), y = e)),
                logical(1))
  in_own_order <- vapply(effects, function(e) any(vapply(own_effects, identical, logical(1), y = e)),
                         logical(1))
  set <- vapply(effects, function(e) paste(sort(e), collapse = ":"), "")
  first <- vapply(unique(set[own]), function(s) {
    candidates <- which(own & set == s)
    c(candidates[in_own_order[candidates]], candidates)[1]
  }, integer(1))
  effects[sort(first)]
}


#' The effects brms lists for a parameter's own formula
#'
#' @description
#' The variable sets [brms::conditional_effects()] would list by default for a
#' model with only this parameter's formula, and the formulas of the parameters
#' it contains. For a linear formula that is one set per fixed-effects term,
#' and for the variables of a smooth, `gp()` and monotonic or measurement-error
#' term each of them and each pair of them, see `.ce_term_effects()`. For a
#' non-linear formula it is each of its data variables and each pair of them.
#'
#' @param x A bmmfit object
#' @param par Character string. The bmm parameter name.
#'
#' @return List of character vectors, the effects of one or two variables that
#'   [brms::conditional_effects()] can plot
#'
#' @keywords internal
#' @noRd
.ce_own_effects <- function(x, par) {
  effects <- unname(unlist(lapply(.formula_nodes(x, par), function(node) {
    if (length(node$sub_pars) == 0) {
      unlist(Map(.ce_term_effects, node$labels, node$term_vars), recursive = FALSE)
    } else {
      .ce_var_combs(unique(unlist(node$term_vars)))
    }
  }), recursive = FALSE))
  unique(effects[lengths(effects) <= 2])
}


#' Each of some variables and each pair of them
#'
#' @param vars Character vector of variable names
#'
#' @return List of character vectors
#'
#' @keywords internal
#' @noRd
.ce_var_combs <- function(vars) {
  c(as.list(vars), if (length(vars) > 1) utils::combn(vars, 2, simplify = FALSE))
}


#' The effects brms lists for one fixed-effects term
#'
#' @description
#' A plain term is one effect of all its variables. brms lists the variables of
#' the terms `s()`, `t2()`, `te()`, `ti()`, `gp()`, `mo()`, `me()` and `mi()`,
#' alone and in pairs, and takes them from the arguments that name variables:
#' `by` and those of the smooth, `gp()` or its first argument.
#'
#' @param label A term label of a formula
#' @param vars The variables of the term
#'
#' @return List of character vectors
#'
#' @keywords internal
#' @noRd
.ce_term_effects <- function(label, vars) {
  parts <- strsplit(label, ":", fixed = TRUE)[[1]]
  special <- grepl("^(s|t2|te|ti|gp|mo|me|mi)\\([^:]*\\)$", parts)
  if (!any(special)) {
    return(list(vars))
  }
  .ce_var_combs(unique(unlist(Map(function(part, is_special) {
    if (is_special) .ce_special_term_vars(str2lang(part)) else all.vars(str2lang(part))
  }, parts, special))))
}


#' The variables a special term names
#'
#' @param call The call of the term, e.g. `s(x, by = z)` or `mo(x)`
#'
#' @return Character vector of variable names: those of the unnamed arguments
#'   and of `by` for a smooth or `gp()`, the variables of `x` for `mo()`,
#'   `me()` and `mi()`
#'
#' @keywords internal
#' @noRd
.ce_special_term_vars <- function(call) {
  if (as.character(call[[1]]) %in% c("mo", "me", "mi")) {
    all.vars(match.call(getExportedValue("brms", as.character(call[[1]])), call)$x)
  } else {
    all.vars(call[(names(call) %||% character(length(call))) %in% c("", "by")])
  }
}


#' Spaghetti lines from the draws of a conditional effect
#'
#' @description
#' Builds the `spaghetti` attribute the way `brms` does: one line per draw,
#' and per level of the second effect when there is one.
#'
#' @param cond_data The conditions grid, without the summary columns
#' @param draws Matrix of draws (rows = draws, columns = grid rows)
#' @param effects Character vector of the effect names
#'
#' @return A data frame with the grid repeated per draw, plus `estimate__` and
#'   `sample__`
#'
#' @keywords internal
#' @noRd
.ce_spaghetti <- function(cond_data, draws, effects) {
  sample <- rep(seq_len(nrow(draws)), each = ncol(draws))
  if (length(effects) == 2L) {
    sample <- paste0(sample, "_", cond_data[[effects[2]]])
  }
  cbind(
    cond_data[rep(seq_len(nrow(cond_data)), times = nrow(draws)), , drop = FALSE],
    estimate__ = as.numeric(t(draws)),
    sample__ = factor(sample)
  )
}


#' Apply link transformation to conditional effects
#'
#' @description
#' Internal function that applies link transformation to a
#' conditional_effects object from brms. Can apply either forward
#' or inverse transformation to the estimate and credible interval bounds.
#'
#' @param ce_object A brmsfit_conditional_effects object from brms::conditional_effects()
#' @param link Character string. Link function name
#' @param inverse Logical. If TRUE, apply inverse link (sampling → native).
#'   If FALSE, apply forward link (native → sampling).
#'
#' @return Modified conditional_effects object with transformed values
#'
#' @keywords internal
#' @noRd
.apply_link_transform <- function(ce_object, link, inverse = TRUE) {
  if (link == "identity") {
    ce_object
  } else {
    .apply_link_transform_inner(ce_object, link, inverse)
  }
}


.apply_link_transform_inner <- function(ce_object, link, inverse) {
  transform_cols <- function(d) {
    cols <- intersect(c("estimate__", "lower__", "upper__"), names(d))
    d[cols] <- lapply(d[cols], link_transform, link = link, inverse = inverse)
    d
  }

  # spaghetti draws live in an attribute, so a column-wise transform misses them
  result <- lapply(ce_object, function(df) {
    df <- transform_cols(df)
    if (!is.null(attr(df, "spaghetti"))) {
      attr(df, "spaghetti") <- transform_cols(attr(df, "spaghetti"))
    }
    df
  })
  
  names(result) <- names(ce_object)
  class(result) <- class(ce_object)
  result
}
