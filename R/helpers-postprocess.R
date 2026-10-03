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

#' @title Expected response of a fitted bmm model
#'
#' @description `posterior_epred()` returns draws of the expected value of the
#'   response for each observation: the mean of what [brms::posterior_predict()]
#'   simulates for it. What that is depends on the model:
#'
#'   | Model | Expected response |
#'   |-------|-------------------|
#'   | `ddm()`, `cswald()` | Mean response time, averaged over both responses, under the diffusion process that `posterior_predict()` simulates from |
#'   | `ezdm()`, version `"3par"` | Mean response time of the cell |
#'   | `ezdm()`, version `"4par"` | Mean response time of the cell's upper-boundary responses |
#'   | `sdt_yn()` | Number of "old"/"signal" responses, the number of trials times their probability |
#'   | `sdt_mafc()` | Number of correct responses, the number of trials times their probability |
#'   | `m3()`, `utility()`, `sdt_rating()`, `sdt_ranking()`, `sdt_cdp()` | Expected count of each response category, from the multinomial family of \pkg{brms} |
#'   | `sdm()`, `mixture2p()`, `mixture3p()`, `imm()` | Not defined: the mean of a circular response error is not a useful quantity, so these models stop with an error |
#'
#'   With `dpar` or `nlpar`, `posterior_epred()` returns draws of that model
#'   parameter for every model, as in \pkg{brms}: a `dpar` on its native scale,
#'   an `nlpar` on the scale of its link. [native_parameters()] returns the
#'   model parameters on their native scale over a grid of predictor values.
#'
#' @param object A `bmmfit` object.
#' @param ... Further arguments passed to [brms::posterior_epred()], such as
#'   `newdata` or `ndraws`.
#' @param dpar,nlpar Name of a distributional or non-linear parameter whose
#'   draws are returned instead of the expected response.
#'
#' @return A draws by observations matrix (an array with a third dimension for
#'   the response categories of the multinomial models).
#' @seealso [brms::posterior_epred()], [native_parameters()]
#' @importFrom brms posterior_epred
#' @export
#' @examples
#' \dontrun{
#' fit <- bmm(
#'   bmf(drift ~ 1, bound ~ 1, ndt ~ 1),
#'   data = rddm(200, drift = 1.5, bound = 1.2, ndt = 0.3),
#'   model = ddm(rt = "rt", response = "response"),
#'   backend = "cmdstanr"
#' )
#' # expected response time of each observation, one row per draw
#' epred <- posterior_epred(fit)
#' }
posterior_epred.bmmfit <- function(object, ..., dpar = NULL, nlpar = NULL) {
  if (is.null(dpar) && is.null(nlpar)) {
    refuse_undefined_epred(object)
  }
  NextMethod()
}

# brms computes fitted() from a prep object, so posterior_epred.bmmfit() is
# never reached and a circular mixture would return brms's number
#' @rdname posterior_epred.bmmfit
#' @param scale As in [brms::fitted.brmsfit()]: `"response"` is the expected
#'   response, `"linear"` the linear predictor of `mu`.
#' @export
fitted.bmmfit <- function(object, ..., scale = c("response", "linear"),
                          dpar = NULL, nlpar = NULL) {
  if (match.arg(scale) == "response" && is.null(dpar) && is.null(nlpar)) {
    refuse_undefined_epred(object)
  }
  NextMethod()
}

# Stops for a fit whose expected response is not defined. The custom families
# refuse through the function stored in the family; this also covers the
# models built on a native brms family, for which brms returns a number of its
# own, wherever bmm hands such a fit to brms for an expected response.
refuse_undefined_epred <- function(object) {
  if (expected_response_defined(object$bmm$model)) {
    return(invisible())
  }
  model_name <- intersect(class(object$bmm$model), supported_models(print_call = FALSE))
  posterior_epred_undefined(model_name[1])()
}

# Whether the model's response has an expected value worth returning
expected_response_defined <- function(model) {
  UseMethod("expected_response_defined")
}

#' @export
expected_response_defined.default <- function(model) {
  TRUE
}

# brms averages the von Mises locations of the mixtures linearly, which is no
# circular mean, and for an unbiased model the circular mean is 0 anyway
#' @export
expected_response_defined.circular <- function(model) {
  FALSE
}

# the class of the circular models before bmm 1.0.1 (#216), which restructure()
# keeps on old fits
#' @export
expected_response_defined.vwm <- function(model) {
  FALSE
}

# The posterior_epred function of a custom family whose response has no
# useful expected value. brms finds a custom family's posterior_epred where the
# family stores it, and otherwise looks up posterior_epred_<family name> in the
# family environment, which reset_env() points at the global environment, so a
# family without one fails with "object not found".
posterior_epred_undefined <- function(model_name) {
  force(model_name)
  function(prep) {
    stop2("The expected response is not defined for the {model_name} model; \\
          use native_parameters() for the model parameters.")
  }
}

# A custom posterior_epred returns a draws x observations matrix. The moment
# helpers it calls flatten their arguments column by column, so their result
# comes back in that order.
.epred_matrix <- function(x, prep) {
  matrix(x, nrow = prep$ndraws, ncol = prep$nobs)
}

# A data column lined up with the draws x observations dpar matrices, which
# hold the draws of one observation in one column
.epred_data <- function(x, prep) {
  rep(x, each = prep$ndraws)
}
