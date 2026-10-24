#' easystats support for bmmfit objects
#'
#' @description
#' These S3 methods let the \pkg{insight} package, and through it the easystats
#' packages (\pkg{parameters}, \pkg{bayestestR}, \pkg{modelbased},
#' \pkg{performance}, \pkg{report}), read a bmm fit as the model you specified
#' with [bmm()]: its parameters, their formulas, the response variables and the
#' data, rather than the \pkg{brms} model that bmm builds from them.
#'
#' @details
#' Each method reshapes what the fit already stores:
#' \itemize{
#'   \item `find_formula()` returns one formula per model parameter, with the
#'     population-level and group-level parts as separate elements
#'     (`<parameter>` and `<parameter>_random`), plus a `conditional` formula
#'     of the response variables on all population-level predictors.
#'   \item `find_predictors()` returns the predictors of these formulas.
#'   \item `find_response()` returns the response variables of the model.
#'   \item `model_info()` describes the model as not linear and as having no
#'     dispersion parameter; the family is the one of the \pkg{brms} model.
#'   \item `link_function()` and `link_inverse()` return the link function of
#'     the first parameter that has one and its inverse. bmm models have one
#'     link per parameter, so for anything else use `model$links`.
#' }
#'
#' @param x A `bmmfit` object (created by [bmm()]).
#' @param effects Which coefficients to return: `"fixed"`, `"random"`, `"all"`
#'   and for `find_parameters()` also `"full"`, `"random_variance"` and
#'   `"grouplevel"`, as for the methods of \pkg{insight} for \pkg{brms} fits.
#' @param component For `find_parameters()`, `"all"`, `"conditional"`,
#'   `"location"` or the name of a model parameter; for `find_predictors()`
#'   as in \pkg{insight}; ignored by `get_data()`.
#' @param flatten Should a character vector be returned instead of a list?
#' @param combine,verbose Ignored, for consistency with the generics.
#' @param ... Ignored, for consistency with the generics.
#'
#' @return See the \pkg{insight} generics, e.g. [insight::find_formula()] and
#'   [insight::find_parameters()].
#'
#' @seealso [insight::find_formula()], [insight::find_parameters()],
#'   [insight::get_priors()]
#'
#' @examples
#' \dontrun{
#' insight::find_formula(fit)
#' parameters::model_parameters(fit)
#' bayestestR::describe_posterior(fit)
#' modelbased::estimate_means(fit, by = "set_size", predict = "c")
#' }
#'
#' @name easystats-bmmfit
NULL


#' @rdname easystats-bmmfit
#' @exportS3Method insight::find_formula
find_formula.bmmfit <- function(x, verbose = TRUE, ...) {
  parts <- lapply(.insight_parameter_formulas(x), .split_bar_terms)
  fixed <- lapply(parts, `[[`, "fixed")
  random <- lapply(parts, `[[`, "random")
  names(random) <- paste0(names(random), "_random")
  random <- Filter(length, random)
  predictors <- setdiff(unique(unlist(lapply(fixed, all.vars))), names(fixed))
  conditional <- stats::as.formula(paste(
    paste(find_response.bmmfit(x), collapse = " + "), "~",
    paste(if (length(predictors)) predictors else "1", collapse = " + ")
  ))
  structure(c(list(conditional = conditional), fixed, random),
            class = c("insight_formula", "list"))
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::find_predictors
find_predictors.bmmfit <- function(x, effects = "fixed", component = "all",
                                   flatten = FALSE, verbose = TRUE, ...) {
  # insight adds every variable of the brms formulas of a non-linear model as
  # a `nonlinear` element, the columns bmm builds among them
  out <- NextMethod(flatten = FALSE)
  out$nonlinear <- NULL
  if (!length(out)) {
    return(NULL)
  }
  if (flatten) unique(unlist(out, use.names = FALSE)) else out
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::find_response
find_response.bmmfit <- function(x, combine = TRUE, ...) {
  unlist(x$bmm$model$resp_vars, use.names = FALSE)
}


#' @rdname easystats-bmmfit
#' @param response Ignored, for consistency with the generic.
#' @exportS3Method insight::model_info
model_info.bmmfit <- function(x, response = NULL, ...) {
  out <- NextMethod()
  out$is_linear <- FALSE
  out$is_dispersion <- FALSE
  out
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::link_function
link_function.bmmfit <- function(x, ...) {
  link <- x$bmm$model$links[[1]]
  function(mu) link_transform(mu, link)
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::link_inverse
link_inverse.bmmfit <- function(x, ...) {
  link <- x$bmm$model$links[[1]]
  function(eta) link_transform(eta, link, inverse = TRUE)
}


# The user's formula as bmm completed it before fitting, restricted to the
# model's declared parameters: m3's activation formulas and the technical
# parameters a family requires are not parameters of the model
.insight_parameter_formulas <- function(x) {
  model <- x$bmm$model
  formula <- suppressMessages(add_missing_parameters(model, x$bmm$user_formula))
  formula[intersect(names(formula), names(model$parameters))]
}

# insight expects the group-level terms of a formula apart from its
# population-level terms, as brms formulas are split by insight's brmsfit
# method; stats::terms() keeps each `(... | group)` term whole. insight reads
# only the grouping factors and slopes from these terms and cannot parse `||`,
# which says nothing about either.
.split_bar_terms <- function(formula) {
  terms <- stats::terms(formula)
  labels <- attr(terms, "term.labels")
  is_bar <- grepl("|", labels, fixed = TRUE)
  env <- environment(formula)
  random <- lapply(sub("||", "|", labels[is_bar], fixed = TRUE), function(label) {
    stats::as.formula(paste("~", label), env = env)
  })
  list(
    fixed = stats::as.formula(
      paste("~", paste(c(attr(terms, "intercept"), labels[!is_bar]), collapse = " + ")),
      env = env
    ),
    random = if (length(random) == 1) random[[1]] else random
  )
}
