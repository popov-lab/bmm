#' report and performance support for bmmfit objects
#'
#' @description
#' These S3 methods let [report::report()] describe a bmm fit as the model you
#' specified with [bmm()], and let [performance::model_performance()] return
#' only the indices that are defined for bmm models. They build on the
#' \pkg{insight} methods described in [easystats-bmmfit].
#'
#' @details
#' \itemize{
#'   \item `format_model()` names the model and its version, e.g. "Bayesian
#'     Signal Discrimination Model (SDM) by Oberauer (2023) (simple)".
#'   \item `report_model()` describes the model in one sentence, with all
#'     response variables of models such as [ezdm()] or [m3()] in one list,
#'     the formula of each model parameter as `summary()` prints it, and the
#'     grouping factors of the random effects.
#'   \item `report_priors()` gives the text of `prior_info(fit, format =
#'     "text")`: the priors in the names of the model parameters, with a note
#'     that they apply on the sampling (link) scale.
#'   \item `model_performance()` returns LOOIC and WAIC (with ELPD and the
#'     standard errors). R2, RMSE and the other indices of
#'     \pkg{performance} need the expected response or residuals, which bmm
#'     models do not define; requested ones are left out with a message.
#' }
#'
#' The parameter summaries and the paragraph on convergence come from
#' \pkg{report} unchanged.
#'
#' @param x,model A `bmmfit` object (created by [bmm()]).
#' @param table Ignored, for consistency with the generic.
#' @param metrics The indices to return: `"LOOIC"`, `"WAIC"`, or `"all"` and
#'   `"common"` for both.
#' @param verbose Should a message name the requested indices that are left
#'   out?
#' @param ... Passed on to the methods of \pkg{report} and \pkg{performance}
#'   for \pkg{brms} fits.
#'
#' @return See the generics, e.g. [report::report_model()] and
#'   [performance::model_performance()].
#'
#' @seealso [prior_info()], [easystats-bmmfit]
#'
#' @examples
#' \dontrun{
#' report::report(fit)
#' performance::model_performance(fit)
#' }
#'
#' @name easystats-report-bmmfit
NULL


#' @rdname easystats-report-bmmfit
#' @exportS3Method report::format_model
format_model.bmmfit <- function(x) {
  model <- x$bmm$model
  version <- if (is.null(model$version) || model$version %in% c(NA, "NA", "")) "" else glue(" ({model$version})")
  # model names such as "The Multinomial / Memory Measurement Model" follow "a"
  paste0("Bayesian ", sub("^The ", "", model$name), version)
}


# the method for brms fits writes one sentence per response variable, as it
# pastes the vector of find_response() into its text, and gives the formula of
# the responses on all predictors ("dev_rad ~ set_size" for an sdm with
# c ~ 0 + set_size), not the model's; the bootstrap flag it reads from
# report_table() is never set for a bmm fit, so the table is not computed
#' @rdname easystats-report-bmmfit
#' @exportS3Method report::report_model
report_model.bmmfit <- function(x, table = NULL, ...) {
  model_text <- report::format_model(x)
  algorithm <- report::format_algorithm(x)
  text_full <- if (nzchar(algorithm)) paste0(model_text, " (estimated using ", algorithm, ")") else model_text
  predict_text <- paste0(" to predict ", collapse_and(insight::find_response(x)))
  if (!suppressWarnings(insight::is_nullmodel(x))) {
    predict_text <- paste0(
      predict_text, " with ", collapse_and(insight::find_predictors(x, effects = "fixed", flatten = TRUE))
    )
  }
  formula <- paste(summarise_formula(x$bmm$user_formula, x$bmm$model), collapse = "; ")
  text_full <- paste0(text_full, predict_text, " (formula: ", formula, ")")
  model_text <- paste0(model_text, predict_text)
  if (insight::model_info(x)$is_mixed) {
    random <- report::report_random(x)
    text_full <- paste0(text_full, ". ", random)
    model_text <- paste0(model_text, ". ", summary(random))
  }
  # report_text() appends the period after the model sentence
  priors <- report::report_priors(x)
  text_full <- paste0(text_full, ". ", sub("\\.$", "", priors))
  model_text <- paste0(model_text, ". ", sub("\\.$", "", summary(priors)))
  report::as.report_model(text_full, summary = model_text, ...)
}


# the method for brms fits reads the random terms from insight::find_terms(),
# which has none for a bmm fit, and then fails on an undefined variable
#' @rdname easystats-report-bmmfit
#' @exportS3Method report::report_random
report_random.bmmfit <- function(x, ...) {
  groups <- unique(insight::find_random(x, flatten = TRUE))
  text <- paste0(
    "The model included ", collapse_and(groups), " as random effect", if (length(groups) > 1) "s"
  )
  report::as.report_random(text, summary = text, ...)
}


#' @rdname easystats-report-bmmfit
#' @exportS3Method report::report_priors
report_priors.bmmfit <- function(x, ...) {
  text <- prior_report_text(prior_info(x))
  report::as.report_priors(text, summary = sub("(?<=\\.)\\s.*$", "", text, perl = TRUE))
}


#' @rdname easystats-report-bmmfit
#' @exportS3Method performance::model_performance
model_performance.bmmfit <- function(model, metrics = c("LOOIC", "WAIC"), verbose = TRUE, ...) {
  defined <- c("LOOIC", "WAIC")
  requested <- if (all(metrics %in% c("all", "common"))) defined else toupper(metrics)
  undefined <- setdiff(requested, defined)
  if (length(undefined) && verbose) {
    message2(
      "{collapse_and(undefined)} {if (length(undefined) > 1) 'are' else 'is'} not defined \\
      for bmm models and left out."
    )
  }
  kept <- intersect(requested, defined)
  NextMethod(metrics = if (length(kept)) kept else defined)
}
