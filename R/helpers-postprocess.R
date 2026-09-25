#' @title Generic S3 method for postprocessing the fitted brm model
#' @description Called by bmm() to automatically perform some type of postprocessing
#'   depending on the model type. It will call the appropriate postprocess_brm.*
#'   methods based on the list of classes defined in the .model_* functions. For
#'   models with several classes listed, it will call the functions in the order
#'   they are listed. Thus, any operations that are common to a group of models
#'   should be defined in the appropriate postprocess_brm.* function, where \*
#'   corresponds to the shared class. For example, for the sdm model, the
#'   postprocessing involves setting the link function for the c parameter to "log",
#'   because it was coded manually in the stan code, but it was specified as "identity"
#'   in the brms custom family. If your model requires no postprocessing, you can
#'   skip this method, and the default method will be used (which returns the same
#'   brmsfit object that was passed to it).
#' @param model A model list object returned from check_model()
#' @param fit the fitted brm model
#' @param ... Additional arguments passed to the method
#' @return An object of class brmsfit, with any necessary postprocessing applied
#' @export
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' fit <- readRDS("my_saved_fit.rds")
#' postprocessed_fit <- prostprocess_brm(fit)
#'
#' @keywords internal developer
postprocess_brm <- function(model, fit, ...) {
  UseMethod("postprocess_brm")
}

#' @export
postprocess_brm.bmmodel <- function(model, fit, ...) {
  dots <- list(...)
  class(fit) <- c("bmmfit", "brmsfit")
  fit$version$bmm <- utils::packageVersion("bmm")
  fit$bmm <- nlist(
    model,
    user_formula = dots$user_formula,
    configure_opts = dots$configure_opts
  )
  attr(fit$data, "data_name") <- attr(dots$fit_args$data, "data_name")

  # add bmm version to the stancode
  fit$model <- add_bmm_version_to_stancode(fit$model)
  reset_env(NextMethod("postprocess_brm"))
}

#' @export
postprocess_brm.default <- function(model, fit, ...) {
  fit
}

get_mu_pars <- function(object) {
  bterms <- brms::brmsterms(object$formula)
  dpars <- bterms$dpars
  if ("mu" %in% names(dpars)) {
    X <- get_model_matrix(dpars$mu$fe, object$data)
    return(colnames(X))
  }
  NULL
}

#' @title Generic S3 method for reverting any postprocessing of the fitted brm model
#' @description Called by update.bmmfit() to automatically revert some of the postprocessing
#'   depending on the model type. It will call the appropriate revert_postprocess_brm.*
#'   methods based on the list of classes defined in the .model_* functions. For
#'   models with several classes listed, it will call the functions in the order
#'   they are listed. For example, for the sdm model, the
#'   postprocessing involves setting the link function for the c parameter to "log",
#'   because it was coded manually in the stan code, but it was specified as "identity"
#'   in the brms custom family. However, during the update process, the link function
#'   should be set back to "identity". Only use this if you have a specific reason to
#'   revert the postprocessing (if otherwise the update method would produce incorrect
#'   results).
#' @param model A model list object returned from check_model()
#' @param fit the fitted brm model
#' @param ... Additional arguments passed to the method
#' @return An object of class brmsfit, with any necessary postprocessing applied
#' @export
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' fit <- readRDS("my_saved_fit.rds")
#' postprocessed_fit <- prostprocess_brm(fit)
#' reverted_fit <- revert_postprocess_brm(postprocessed_fit)
#'
#' @keywords internal developer
revert_postprocess_brm <- function(model, fit, ...) {
  UseMethod("revert_postprocess_brm")
}

#' @export
revert_postprocess_brm.default <- function(model, fit, ...) {
  fit
}
