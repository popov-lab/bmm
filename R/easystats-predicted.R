#' Predictions of a model parameter for easystats and marginaleffects
#'
#' @description
#' With these methods, \pkg{insight} and \pkg{marginaleffects}, and through
#' them `modelbased::estimate_means()`, `modelbased::estimate_contrasts()`,
#' `modelbased::estimate_relation()` and the \pkg{see} plots of their results,
#' predict a model parameter on its native scale, as [native_parameters()]
#' computes it: with the inverse link applied draw by draw and the mixture
#' weights softmaxed jointly, so that `predict = "thetat"` returns the
#' probability of a memory response.
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
#' column of the data to predict for and holds only levels of the fit (or
#' `allow_new_levels = TRUE`); otherwise the predictions are for the
#' population level, as \pkg{insight} does for \pkg{brms} fits.
#'
#' The point prediction is the median of the draws.
#'
#' @param x,model A `bmmfit` object (created by [bmm()]).
#' @param data,newdata The data to predict for. `get_predicted()` predicts for
#'   the data of the fit if `data` is `NULL`.
#' @param predict The name of a model parameter, or one of the prediction
#'   types of [insight::get_predicted()].
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
#'   `sample_new_levels`, passed to [native_parameters()], and for
#'   `get_predict()` `dpar`, `ndraws` and `re_formula` (or `re.form`). For
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
  if (is.null(predict)) {
    predict <- .estimated_parameters(x)[1]
    message2("Predicting the model parameter '{predict}'; choose another with `predict`.")
  }
  if (length(predict) > 1) {
    predict <- predict[1]
    warning2("`predict` takes one value; predicting '{predict}'.")
  }
  if (!predict %in% names(model$parameters)) {
    if (predict %in% c("expectation", "expected", "response", "link",
                       "prediction", "predicted", "classification")) {
      return(NextMethod())
    }
    stop2("'{predict}' is not a parameter of the {model$name} model. \\
          Its parameters are {collapse_comma(names(model$parameters))}.")
  }
  data <- data %||% get_data.bmmfit(x)
  draws <- .native_draws(x, data, predict, .easystats_re_formula(x, data, include_random, ...),
                         ndraws = iterations, ...)
  predictions <- apply(draws, 1, stats::median)
  attr(predictions, "iterations") <- as.data.frame(draws)
  ci_data <- insight::get_predicted_ci(x, predictions, ci = ci, ci_method = ci_method %||% "quantile")
  structure(
    predictions,
    ci_data = ci_data, data = data, ci = ci, predict = predict,
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
    .easystats_re_formula(model, newdata, ...)
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
# and one column per draw. native_parameters() returns the draws of each grid
# cell together, the draws varying fastest.
.native_draws <- function(x, newdata, par, re_formula, ndraws = NULL, draw_ids = NULL,
                          allow_new_levels = FALSE, sample_new_levels = "uncertainty", ...) {
  np <- native_parameters(x, newdata = newdata, pars = par, re_formula = re_formula,
                          ndraws = ndraws, draw_ids = draw_ids,
                          allow_new_levels = allow_new_levels,
                          sample_new_levels = sample_new_levels)
  draws <- matrix(np$value, nrow = nrow(newdata), byrow = TRUE)
  colnames(draws) <- paste0("iter_", seq_len(ncol(draws)))
  draws
}

# insight includes the random effects of a brms fit only if the data name
# every grouping factor and hold only levels the fit knows, so a grid built
# from the fixed effects, whose grouping factors modelbased sets to NA,
# predicts the population level instead of failing or predicting the first
# participant
.easystats_re_formula <- function(x, data, include_random = "default",
                                  allow_new_levels = FALSE, ...) {
  if (inherits(include_random, "formula")) {
    return(include_random)
  }
  if (isFALSE(include_random) || identical(include_random, NA)) {
    return(NA)
  }
  groups <- insight::find_random(x, split_nested = TRUE, flatten = TRUE)
  known <- all(groups %in% names(data)) && (isTRUE(allow_new_levels) || all(vapply(
    groups, function(g) all(unique(data[[g]]) %in% x$data[[g]]), logical(1)
  )))
  if (known) NULL else NA
}

# The parameters the fit estimates: those not fixed to a constant, or fixed by
# the model but given a formula by the user
.estimated_parameters <- function(x) {
  model <- x$bmm$model
  pars <- setdiff(names(model$parameters), internal_parameters(model))
  fixed <- setdiff(names(model$fixed_parameters), names(x$bmm$user_formula))
  setdiff(pars, fixed)
}
