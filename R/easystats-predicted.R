#' Predictions of a model parameter for easystats and marginaleffects
#'
#' @description
#' With these methods, \pkg{insight} and \pkg{marginaleffects}, and through
#' them `modelbased::estimate_means()`, `modelbased::estimate_contrasts()`,
#' `modelbased::estimate_relation()` and the \pkg{see} plots of their results,
#' predict a model parameter on its native scale, as [native_parameters()]
#' computes it: with the inverse link applied draw by draw and the mixture
#' weights softmaxed jointly, so that `predict = "thetat"` returns the
#' probability of a target response.
#'
#' @details
#' `predict` (or `dpar` for `get_predict()`, which is how \pkg{modelbased}
#' passes it on) names the model parameter. Without it, `get_predicted()`
#' predicts the first parameter the model estimates and says which. The
#' prediction types of \pkg{insight} (`"expectation"`, `"link"`,
#' `"prediction"`, ...) and of \pkg{marginaleffects} are left to the methods
#' for \pkg{brms} fits; the expected response is not defined for most bmm
#' models.
#'
#' Random effects enter the predictions only if every grouping factor is a
#' column of the data to predict for and holds only levels of the fit, as
#' \pkg{insight} does for \pkg{brms} fits, or new levels with
#' `allow_new_levels = TRUE`. Otherwise the predictions are for the population
#' level, with a warning if you asked for `include_random = TRUE`. A grouping
#' factor that is `NA`, as in the grids of \pkg{modelbased}, always means the
#' population level. For a mixed model, `modelbased::estimate_means()`
#' averages the predictions over the participants in its grid, so its means
#' differ from those of `native_parameters(re_formula = NA)`.
#'
#' The point prediction of `get_predicted()` is the `centrality_function` of
#' the draws, `base::mean` by default, as for \pkg{brms} fits.
#'
#' @param x,model A `bmmfit` object (created by [bmm()]).
#' @param data,newdata The data to predict for. `get_predicted()` also takes
#'   `newdata` in `...`, and predicts for the data of the fit if both are
#'   `NULL`.
#' @param predict The name of a model parameter, or one of the prediction
#'   types of [insight::get_predicted()]. A `type` in `...` replaces it.
#' @param iterations The number of posterior draws to predict from, all of
#'   them if `NULL`.
#' @param ci The probability mass of the credible intervals.
#' @param ci_method `"quantile"` (the default), `"eti"` or `"hdi"`.
#' @param include_random `FALSE` or `NA` for population-level predictions, a
#'   formula for the random effects to include, or `"default"` or `TRUE` for
#'   all of them if the data contain the grouping factors.
#' @param type The prediction type of \pkg{marginaleffects}, used only when
#'   `dpar` is not a model parameter.
#' @param ... For a model parameter: `draw_ids`, `allow_new_levels` and
#'   `sample_new_levels`, passed to [native_parameters()]; for
#'   `get_predicted()` also `newdata`, `type` and `centrality_function`, and
#'   for `get_predict()` `dpar`, `ndraws` and `re_formula` (or `re.form`). For
#'   any other prediction, passed to the methods for \pkg{brms} fits.
#'
#' @return `get_predicted()`: a `get_predicted` object, as described in
#'   [insight::get_predicted()]. `get_predict()`: a `data.table` of the
#'   estimates with the draws in the `posterior_draws` attribute, as
#'   \pkg{marginaleffects} expects.
#'
#' @seealso [native_parameters()], [easystats-bmmfit]
#'
#' @examples
#' \dontrun{
#' insight::get_predicted(fit, predict = "c")
#' modelbased::estimate_means(fit, by = "set_size", predict = "c")
#' modelbased::estimate_contrasts(fit, contrast = "set_size", predict = "c")
#' }
#'
#' @name predicted-bmmfit
NULL


#' @rdname predicted-bmmfit
#' @exportS3Method insight::get_predicted
get_predicted.bmmfit <- function(x, data = NULL, predict = NULL, iterations = NULL,
                                 ci = NULL, ci_method = NULL, include_random = "default", ...) {
  model <- x$bmm$model
  dots <- list(...)
  # NextMethod() passes `predict` and `type` on as given, so that insight
  # resolves them for brms fits; only the parameter path reads `par`
  par <- dots$type %||% predict
  if (is.null(par)) {
    par <- .estimated_parameters(x)[1]
    message2("Predicting the model parameter '{par}'; choose another with `predict`.")
  }
  if (!par[1] %in% names(model$parameters)) {
    if (par[1] %in% c("expectation", "expected", "response", "link",
                      "prediction", "predicted", "classification")) {
      return(NextMethod())
    }
    stop2("'{par[1]}' is not a parameter of the {model$name} model. \\
          Its parameters are {collapse_comma(names(model$parameters))}.")
  }
  if (length(par) > 1) {
    par <- par[1]
    warning2("`predict` takes one value; predicting '{par}'.")
  }
  data <- data %||% dots$newdata %||% get_data.bmmfit(x)
  draws <- .native_draws(
    x, newdata = data, par = par,
    re_formula = .easystats_re_formula(x, data, include_random, dots$allow_new_levels),
    ndraws = iterations, draw_ids = dots$draw_ids,
    allow_new_levels = isTRUE(dots$allow_new_levels),
    sample_new_levels = dots$sample_new_levels %||% "uncertainty"
  )
  predictions <- apply(draws, 1, dots$centrality_function %||% base::mean)
  attr(predictions, "iterations") <- as.data.frame(draws)
  ci_data <- insight::get_predicted_ci(x, predictions, ci = ci, ci_method = ci_method %||% "quantile")
  structure(
    predictions,
    ci_data = ci_data, data = data, ci = ci, predict = par,
    class = c("get_predicted", "numeric")
  )
}


#' @rdname predicted-bmmfit
#' @exportS3Method marginaleffects::get_predict
get_predict.bmmfit <- function(model, newdata, type = "response", ...) {
  dots <- list(...)
  if (!isTRUE(dots$dpar %in% names(model$bmm$model$parameters))) {
    return(NextMethod())
  }
  # NA asks for the population level and NULL for all random effects, so the
  # argument's presence decides, not its value
  re_formula <- if ("re_formula" %in% names(dots)) {
    dots$re_formula
  } else if ("re.form" %in% names(dots)) {
    dots$re.form
  } else {
    .easystats_re_formula(model, newdata, allow_new_levels = dots$allow_new_levels)
  }
  draws <- .native_draws(
    model, newdata, dots$dpar, re_formula,
    ndraws = dots$ndraws, draw_ids = dots$draw_ids,
    allow_new_levels = isTRUE(dots$allow_new_levels),
    sample_new_levels = dots$sample_new_levels %||% "uncertainty"
  )
  # marginaleffects' comparisons split the predictions by group with data.table
  out <- data.table::data.table(
    group = "main_marginaleffect", rowid = newdata$rowid %||% seq_len(nrow(newdata)),
    estimate = apply(draws, 1, stats::median)
  )
  data.table::setattr(out, "posterior_draws", draws)
  out
}


# The draws of one parameter on its native scale, one row per row of newdata
# and one column per draw, as insight and marginaleffects hold them
.native_draws <- function(x, newdata, par, re_formula, ndraws = NULL, draw_ids = NULL,
                          allow_new_levels = FALSE, sample_new_levels = "uncertainty") {
  draws <- t(.native_draw_matrices(
    x, newdata, par, re_formula, ndraws = ndraws, draw_ids = draw_ids,
    dots = nlist(allow_new_levels, sample_new_levels)
  )$draws[[par]])
  dimnames(draws) <- list(NULL, paste0("iter_", seq_len(ncol(draws))))
  draws
}

# insight includes the random effects of a brms fit only if the data name
# every grouping factor and hold only levels the fit knows. modelbased sets
# the grouping factors of a grid built from the fixed effects to NA, and that
# grid asks for the population level: with allow_new_levels it would
# otherwise predict a new participant
.easystats_re_formula <- function(x, data, include_random = "default",
                                  allow_new_levels = FALSE) {
  if (inherits(include_random, "formula")) {
    return(include_random)
  }
  if (isFALSE(include_random) || identical(include_random, NA)) {
    return(NA)
  }
  groups <- insight::find_random(x, split_nested = TRUE, flatten = TRUE)
  absent <- groups[not_in(groups, names(data))]
  if (length(absent) > 0) {
    warnif(
      isTRUE(include_random),
      "`include_random = TRUE`, but the data have no column {collapse_comma(absent)}; \\
      predicting the population level."
    )
    return(NA)
  }
  unknown <- groups[vapply(groups, function(g) {
    anyNA(data[[g]]) || !isTRUE(allow_new_levels) && !all(data[[g]] %in% x$data[[g]])
  }, logical(1))]
  if (length(unknown) > 0) {
    warnif(
      isTRUE(include_random),
      "`include_random = TRUE`, but {collapse_comma(unknown)} in the data hold \\
      missing values, or levels the fit does not know without \\
      `allow_new_levels = TRUE`; predicting the population level."
    )
    return(NA)
  }
  NULL
}

# The parameters the fit estimates: those not fixed to a constant, or fixed by
# the model but given a formula by the user
.estimated_parameters <- function(x) {
  model <- x$bmm$model
  pars <- setdiff(names(model$parameters), internal_parameters(model))
  fixed <- setdiff(names(model$fixed_parameters), names(x$bmm$user_formula))
  setdiff(pars, fixed)
}
