#' @title Update a bmm model
#' @description Update an existing bmm mode. This function calls
#'   [brms::update.brmsfit()], but it applies the necessary bmm postprocessing
#'   to the model object before and after the update.
#' @param object An object of class `bmmfit`
#' @param formula. A [bmmformula()]. If missing, the original formula
#'  is used. Currently you have to specify a full `bmmformula`
#' @param newdata An optional data frame containing the variables in the model
#' @param recompile Logical, indicating whether the Stan model should be recompiled. If
#'   NULL (the default), update tries to figure out internally, if recompilation
#'   is necessary. Setting it to FALSE will cause all Stan code changing
#'   arguments to be ignored.
#' @param ... Further arguments passed to [brms::update.brmsfit()]
#' @return An updated `bmmfit` object refit to the new data and/or formula
#' @details When updating a brmsfit created with the cmdstanr backend in a
#'   different R session, a recompilation will be triggered because by default,
#'   cmdstanr writes the model executable to a temporary directory. To avoid
#'   that, set option "cmdstanr_write_stan_file_dir" to a nontemporary path of
#'   your choice before creating the original bmmfit.
#'
#'   For more information and examples, see [brms::update.brmsfit()]
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # generate artificial data from the Signal Discrimination Model
#' # generate artificial data from the Signal Discrimination Model
#' dat <- data.frame(y = rsdm(2000))
#'
#' # define formula
#' ff <- bmf(c ~ 1, kappa ~ 1)
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = sdm(resp_error = "y"),
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # update the model
#' fit <- update(fit, newdata = data.frame(y = rsdm(2000, kappa = 5)))
#'
update.bmmfit <- function(object, formula., newdata = NULL, recompile = NULL, ...) {
  dots <- list(...)
  local_brms_threads(dots)
  stopif(
    isTRUE(object$version$bmm < "0.3.0"),
    "Updating bmm models works only with models fitted with version 0.3.0 or higher"
  )
  stopif(
    "data" %in% names(dots),
    "Please use argument 'newdata' to update the data."
  )
  stopif(
    "model" %in% names(dots),
    "You cannot update with a different model. Create a new fit with 'bmm()' instead."
  )

  object <- restructure(object)

  model <- object$bmm$model
  old_user_formula <- object$bmm$user_formula
  olddata <- object$data
  configure_opts <- object$bmm$configure_opts

  # revert some postprocessing changes to brmsfit from postprocess_brm
  object <- revert_postprocess_brm(model, object)

  # use the new configure_opts if they are provided
  if (any(names(dots) %in% names(configure_opts))) {
    new_opts <- names(dots)[names(dots) %in% names(configure_opts)]
    configure_opts[new_opts] <- dots[new_opts]
  }
  opts <- configure_options(configure_opts)

  # reuse or replace formula and data
  if (missing(formula.)) {
    user_formula <- old_user_formula
  } else {
    user_formula <- formula.
  }
  if (is.null(newdata)) {
    data <- check_data(model, olddata, user_formula)
    attr(data, "data_name") <- attr(olddata, "data_name")
  } else {
    data <- check_data(model, newdata, user_formula)
    attr(data, "data_name") <- substitute_name(newdata)
  }

  # standard bmm checks and transformations
  formula <- check_formula(model, data, user_formula)
  config_args <- configure_model(model, data, formula)
  prior <- configure_prior(model, data, config_args$formula, object$prior)
  prior <- combine_prior(prior, dots$prior)
  dots$prior <- NULL
  new_fit_args <- combine_args(nlist(config_args, dots, prior))

  # construct the new formula and data only if they have changed
  if (!identical(new_fit_args$formula, object$formula)) {
    formula. <- new_fit_args$formula
  }
  if (!identical(new_fit_args$data, olddata)) {
    newdata <- new_fit_args$data
  }

  # pass back to brms::update.brmsfit
  object <- NextMethod("update", object,
    formula = formula., newdata = newdata,
    prior = prior, recompile = recompile, ...
  )

  # bmm postprocessing
  postprocess_brm(model, object,
    fit_args = new_fit_args, user_formula = user_formula,
    configure_opts = configure_opts
  )
}
