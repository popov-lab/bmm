#' emmeans support for bmmfit objects
#'
#' @description
#' These S3 methods allow [emmeans::emmeans()] to work seamlessly with bmmfit
#' objects by automatically resolving whether a parameter is a distributional
#' (`dpar`) or non-linear (`nlpar`) parameter in brms. Users can always use
#' `dpar = "parameter_name"` regardless of internal classification.
#'
#' @details
#' bmm models use two types of parameters internally in brms:
#' \itemize{
#'   \item **Distributional parameters (`dpar`)**: Used in models with
#'     `custom_family(dpars = ...)` (e.g., SDM, EZDM)
#'   \item **Non-linear parameters (`nlpar`)**: Used in models with `bmf2bf()`
#'     + `nlf()` (e.g., mixture2p, mixture3p, IMM, M3)
#' }
#'
#' Users should not need to know this distinction. These methods intercept the
#' `dpar` argument and, if the parameter is actually an nlpar, silently
#' re-route it to `nlpar`.
#'
#' @param object A bmmfit object (created by [bmm()])
#' @param dpar Character string. Name of the model parameter (e.g., `"kappa"`,
#'   `"c"`, `"thetat"`). Automatically resolved to the correct brms parameter
#'   type.
#' @param nlpar Character string. Explicit non-linear parameter name. If
#'   provided, takes precedence over `dpar`.
#' @param ... Additional arguments passed to the brmsfit methods.
#' @param trms,xlev,grid Arguments passed by emmeans internally.
#'
#' @return See [emmeans::recover_data()] and [emmeans::emm_basis()] for return
#'   values.
#'
#' @seealso [emmeans::emmeans()]
#'
#' @examples
#' \dontrun{
#' # Works for any bmm model — no need to know dpar vs nlpar
#' em <- emmeans(fit, ~ condition, dpar = "kappa")
#' pairs(em)
#' confint(em)
#' }
#'
#' @name emmeans-bmmfit
NULL


#' @rdname emmeans-bmmfit
#' @exportS3Method emmeans::recover_data
recover_data.bmmfit <- function(object, ..., dpar = NULL, nlpar = NULL) {
  resolved <- .bmmfit_resolve_par(object, dpar, nlpar)
  class(object) <- class(object)[class(object) != "bmmfit"]
  emmeans::recover_data(object, ..., dpar = resolved$dpar, nlpar = resolved$nlpar)
}


#' @rdname emmeans-bmmfit
#' @exportS3Method emmeans::emm_basis
emm_basis.bmmfit <- function(object, trms, xlev, grid, ...,
                              dpar = NULL, nlpar = NULL) {
  resolved <- .bmmfit_resolve_par(object, dpar, nlpar)
  class(object) <- class(object)[class(object) != "bmmfit"]
  emmeans::emm_basis(object, trms, xlev, grid, ...,
                     dpar = resolved$dpar, nlpar = resolved$nlpar)
}


.bmmfit_resolve_par <- function(object, dpar, nlpar) {
  if (!is.null(dpar) && is.null(nlpar) && !is.null(object$bmm$model)) {
    par_info <- .get_parameter_info(object, dpar)
    if (par_info$type == "nlpar") {
      nlpar <- dpar
      dpar  <- NULL
    }
  }
  list(dpar = dpar, nlpar = nlpar)
}
