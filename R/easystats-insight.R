#' easystats support for bmmfit objects
#'
#' @description
#' These S3 methods let the \pkg{insight} package, and through it the easystats
#' packages (\pkg{parameters}, \pkg{bayestestR}, \pkg{modelbased},
#' \pkg{performance}, \pkg{report}), read a bmm fit as the model you specified
#' with [bmm()]: its parameters, their formulas, the response variables and the
#' data, rather than the \pkg{brms} model that bmm builds from them. With them,
#' for example, [parameters::model_parameters()] and
#' [bayestestR::describe_posterior()] summarise each model parameter in its own
#' component, and `modelbased::estimate_means(fit, by = "set_size", predict =
#' "c")` returns the estimates of a model parameter on its native scale.
#'
#' @details
#' Each method reshapes what the fit already stores:
#' \itemize{
#'   \item `find_formula()` returns one formula per model parameter, with the
#'     population-level and group-level parts as separate elements
#'     (`<parameter>` and `<parameter>_random`), plus a `conditional` formula
#'     of the response variables on all population-level predictors.
#'   \item `find_predictors()` returns the predictors of these formulas:
#'     the population-level predictors for `effects = "fixed"`, the grouping
#'     factors for `"random"` and both for `"all"`. Every model parameter
#'     counts as both a location and a distributional parameter, so each
#'     `component` returns the same variables.
#'   \item `find_response()` returns the names of the response columns, as
#'     `get_data()` returns them.
#'   \item `get_data()` returns the data the model was fitted to as you passed
#'     them: without the columns bmm adds for its own use, and with the
#'     response categories of a model fitted to counts (e.g. [m3()]) as
#'     columns of their own. The model frame stored in the fit is always the
#'     source; `source` is ignored.
#'   \item `find_parameters()` returns the names of the posterior draws in one
#'     element per model parameter (`<parameter>` for the population-level
#'     coefficients, `<parameter>_random` for the group-level ones), without
#'     the parameters that the \pkg{brms} family needs but the model does not
#'     have, such as the fixed `mu` of the response-time models. `component`
#'     selects a model parameter by name; `"all"`, `"conditional"` and
#'     `"location"` select all of them. [insight::get_parameters()] and
#'     [insight::clean_parameters()] build on it.
#'   \item `get_priors()` returns one row per population-level coefficient,
#'     group-level standard deviation and correlation, with the prior that
#'     applies to it. A parameter fixed to a constant has the distribution
#'     `"constant"` and the constant as its location.
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
#' @param parameters A regular expression the names of the returned
#'   parameters must match.
#' @param combine,source,verbose Ignored, for consistency with the generics.
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
  formulas <- find_formula.bmmfit(x)
  fixed <- lapply(formulas[!endsWith(names(formulas), "_random")], function(f) {
    setdiff(all.vars(f[[length(f)]]), names(formulas))
  })
  # every model parameter is a location parameter as well as a distributional
  # one; a component bmm models do not have, such as zi, has no predictors
  fixed <- fixed[switch(component,
    all = names(fixed),
    conditional = , location = "conditional",
    distributional = , auxiliary = setdiff(names(fixed), "conditional"),
    character(0)
  )]
  random <- if (length(fixed)) insight::find_random(x, split_nested = TRUE)
  out <- Filter(length, switch(effects, fixed = fixed, random = random, c(fixed, random)))
  if (!length(out)) {
    return(NULL)
  }
  if (flatten) unique(unlist(out, use.names = FALSE)) else out
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::find_response
find_response.bmmfit <- function(x, combine = TRUE, ...) {
  # a model fitted to counts names the columns of its response matrix, which
  # need not be its resp_vars (sdt_cdp() stores a column prefix there)
  response <- x$data[[all.vars(brms::brmsterms(x$formula)$respform)[1]]]
  colnames(response) %||% unlist(x$bmm$model$resp_vars, use.names = FALSE)
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::get_data
get_data.bmmfit <- function(x, effects = "all", component = "all",
                            source = "frame", verbose = FALSE, ...) {
  data <- x$data[setdiff(names(x$data), built_data_columns(x$bmm$model))]
  data <- do.call(cbind, lapply(names(data), function(col) {
    if (is.matrix(data[[col]])) as.data.frame(data[[col]]) else data[col]
  }))
  groups <- insight::find_random(x, split_nested = TRUE, flatten = TRUE)
  switch(effects,
    random = data[intersect(groups, names(data))],
    fixed = data[setdiff(names(data), setdiff(groups, insight::find_predictors(x, flatten = TRUE)))],
    data
  )
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::find_parameters
find_parameters.bmmfit <- function(x, effects = "all", component = "all",
                                   flatten = FALSE, parameters = NULL, ...) {
  draws <- .insight_draws(x)
  classes <- switch(effects,
    fixed = "b", random = c("sd", "cor", "r"), random_variance = c("sd", "cor"),
    grouplevel = "r", full = c("b", "sd", "cor", "r"), c("b", "sd", "cor")
  )
  pars <- union(intersect(names(x$bmm$model$parameters), draws$parameter), draws$parameter)
  pars <- switch(component, all = , conditional = , location = pars, intersect(component, pars))
  draws <- draws[draws$class %in% classes & draws$parameter %in% pars, , drop = FALSE]
  if (!is.null(parameters)) {
    draws <- draws[grepl(parameters, draws$variable, perl = TRUE), , drop = FALSE]
  }
  draws$element <- ifelse(draws$class == "b", draws$parameter, paste0(draws$parameter, "_random"))
  out <- split(draws$variable, factor(draws$element, unique(c(pars, paste0(pars, "_random")))))
  out <- Filter(length, out)
  if (flatten) unique(unlist(out, use.names = FALSE)) else out
}


#' @rdname easystats-bmmfit
#' @exportS3Method insight::get_priors
get_priors.bmmfit <- function(x, verbose = TRUE, ...) {
  draws <- .insight_draws(x)
  draws <- draws[match(find_parameters.bmmfit(x, flatten = TRUE), draws$variable), , drop = FALSE]
  provenance <- prior_provenance(x)
  prior <- vapply(seq_len(nrow(draws)), function(i) {
    .applicable_prior(draws[i, ], provenance)
  }, character(1))
  moments <- do.call(rbind, lapply(prior, .prior_moments))
  data.frame(
    Parameter = draws$variable,
    Distribution = ifelse(nzchar(prior), sub("\\(.*", "", prior), "uniform"),
    Location = moments[, 1], Scale = moments[, 2], df = moments[, 3],
    stringsAsFactors = FALSE
  )
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
# which says nothing about either, nor `/`, which brms expands.
.split_bar_terms <- function(formula) {
  terms <- stats::terms(formula)
  labels <- attr(terms, "term.labels")
  is_bar <- grepl("|", labels, fixed = TRUE)
  env <- environment(formula)
  random <- lapply(.expand_nested_groups(sub("||", "|", labels[is_bar], fixed = TRUE)), function(label) {
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

# brms reads `(1 | a/b)` as `(1 | a) + (1 | a:b)`
.expand_nested_groups <- function(labels) {
  unlist(lapply(labels, function(label) {
    effects <- sub("\\|[^|]*$", "| ", label)
    groups <- trimws(strsplit(sub("^.*\\|", "", label), "/", fixed = TRUE)[[1]])
    paste0(effects, vapply(seq_along(groups), function(i) {
      paste(groups[seq_len(i)], collapse = ":")
    }, character(1)))
  }))
}

# One row per posterior draw of a coefficient, with the brms parameter it
# belongs to and its class. brms names the draws b_<par>_<coef>,
# sd_<group>__<par>_<coef>, cor_<group>__<par>_<coef>__<par>_<coef> and
# r_<group>__<par>[<level>,<coef>], leaving out <par>_ for the main dpar. The
# parameter is the dpar or nlpar name the remainder starts with; brms forbids
# `_` in these names, so at most one matches. A population-level predictor of
# the main dpar may itself start with <par>_, so the main dpar's coefficients,
# as the fit's prior lists them, are matched first.
.insight_draws <- function(x) {
  variable <- brms::variables(x)
  class <- sub("^(b|bs|bsp|bcs)_.*", "b", sub("^(sd|cor|r)_.*", "\\1", variable))
  keep <- class %in% c("b", "sd", "cor", "r")
  variable <- variable[keep]
  class <- class[keep]
  grouped <- class != "b"
  rest <- sub("^[a-z]+_", "", variable)
  rest[grouped] <- sub("^.*?__", "", rest[grouped], perl = TRUE)
  terms <- brms::brmsterms(x$formula)
  main <- names(terms$dpars)[1]
  # a mixture names the coefficients of its first dpar, mu1, with the prefix
  candidates <- names(c(terms$dpars, terms$nlpars))
  prefix <- vapply(rest, function(r) {
    c(candidates[startsWith(r, paste0(candidates, "_")) | startsWith(r, paste0(candidates, "["))], NA)[1]
  }, character(1), USE.NAMES = FALSE)
  prior <- x$prior
  prefix[class == "b" & rest %in% prior$coef[prior$class == "b" & !nzchar(prior$dpar) & !nzchar(prior$nlpar)]] <- NA
  draws <- data.frame(
    variable, class,
    parameter = ifelse(is.na(prefix), main, prefix),
    coef = ifelse(is.na(prefix), rest, substring(rest, nchar(prefix) + 2)),
    group = ifelse(grouped, sub("^[a-z]+_(.*?)__.*", "\\1", variable, perl = TRUE), ""),
    stringsAsFactors = FALSE
  )
  drop_technical_parameters(draws, x$bmm$model)
}

# The prior brms applies to a coefficient: the one set for the coefficient,
# else the one for its class in its group, else the one for its class. The
# intercept of a dpar takes its prior from class Intercept; a correlation's
# prior belongs to no parameter.
.applicable_prior <- function(draw, provenance) {
  keys <- list(
    c("Intercept", "", ""),
    c(draw$class, draw$coef, draw$group),
    c(draw$class, "", draw$group),
    c(draw$class, "", "")
  )
  if (draw$class != "b" || draw$coef != "Intercept") keys[[1]] <- NULL
  same_parameter <- provenance$parameter %in% draw$parameter | draw$class == "cor"
  for (key in keys) {
    hit <- same_parameter & provenance$class == key[1] &
      provenance$coef == key[2] & provenance$group == key[3]
    if (any(hit)) return(provenance$prior[which(hit)[1]])
  }
  ""
}

# Location, scale and degrees of freedom of a prior, as insight's get_priors()
# reports them; NA for distributions without a location-scale form
.prior_moments <- function(prior) {
  args <- suppressWarnings(as.numeric(strsplit(sub("^[^(]*\\((.*)\\)$", "\\1", prior), ",")[[1]]))
  switch(sub("\\(.*", "", prior),
    normal = , cauchy = , logistic = c(args[1:2], NA),
    student_t = c(args[2:3], args[1]),
    constant = , lkj = c(args[1], NA, NA),
    c(NA, NA, NA)
  )
}
