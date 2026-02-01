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
#'   (see `names(x$bmm$model$parameters)`). If `NULL`, the function will error
#'   and ask you to specify a parameter explicitly.
#' @param scale Character. Scale on which to show the parameter:
#'   \describe{
#'     \item{`"native"` (default)}{Show on natural scale using inverse link transformation.
#'       For example, `kappa` with log link shown on exp scale, `thetat` with
#'       identity link shown on probability scale (between 0 and 1).}
#'     \item{`"parameter"` or `"sampling"`}{Show on the sampling scale (as used during MCMC).
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
#'     \item `method`: Method for computing effects ("posterior_predict" or "posterior_epred")
#'   }
#'
#' @return A `brmsfit_conditional_effects` object (from brms), which can be:
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
#'   \item **Non-linear parameters (`nlpar`)**: Core model parameters like `kappa`, `c`, `a`, `thetat`
#'   \item **Distributional parameters (`dpar`)**: Derived parameters used in brms mixture distributions
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
#' @seealso [brms::conditional_effects()] for the underlying brms function
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a mixture model with set size effect on kappa
#' fit <- bmm(
#'   formula = bmf(kappa ~ set_size, thetat ~ 1),
#'   data = zhang_luck_2008,
#'   model = mixture3p(
#'     resp_error = "err_rad",
#'     nt_features = paste0("nt_", 1:7, "_err_rad"),
#'     set_size = "set_size"
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
#' ce_specific <- conditional_effects(fit, par = "kappa", effects = "set_size")
#'
#' # Combine with other brms options
#' ce_detailed <- conditional_effects(
#'   fit,
#'   par = "kappa",
#'   effects = "set_size",
#'   spaghetti = TRUE,
#'   ndraws = 100
#' )
#' }
conditional_effects.bmmfit <- function(x,
                                       par = NULL,
                                       scale = c("native", "parameter", "sampling"),
                                       ...) {
  # Validate inputs
  stopif(!inherits(x, "bmmfit"), "x must be a bmmfit object")

  scale <- match.arg(scale)

  # Require explicit parameter specification
  if (is.null(par)) {
    model_pars <- names(x$bmm$model$parameters)
    stop2(
      "Argument 'par' is required. Please specify which parameter to plot.\n",
      "Available parameters: {paste(model_pars, collapse = ', ')}"
    )
  }

  # Validate that par is a single character string
  stopif(
    !is.character(par) || length(par) != 1,
    "Argument 'par' must be a single character string"
  )

  # Get parameter information from model
  par_info <- .get_parameter_info(x, par)

  # Route to brms::conditional_effects with appropriate arguments
  # Use brms:::conditional_effects.brmsfit directly to avoid infinite recursion
  # since x has class c("bmmfit", "brmsfit") and would dispatch back here
  dots <- list(...)

  if (par_info$type == "dpar") {
    ce_result <- brms:::conditional_effects.brmsfit(x, dpar = par_info$brms_name, ...)
  } else if (par_info$type == "nlpar") {
    ce_result <- brms:::conditional_effects.brmsfit(x, nlpar = par_info$brms_name, ...)
  } else {
    stop2("Internal error: parameter type must be 'dpar' or 'nlpar'")
  }

  # Apply transformations based on parameter type and requested scale
  # - nlpars: brms returns on sampling scale, apply inverse link for native scale
  # - dpars: brms returns on native scale, apply forward link for sampling scale
  if (par_info$link != "identity") {
    if (par_info$type == "nlpar" && scale == "native") {
      # nlpar: sampling → native (apply inverse link)
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = TRUE)
    } else if (par_info$type == "dpar" && scale %in% c("parameter", "sampling")) {
      # dpar: native → sampling (apply forward link)
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = FALSE)
    }
  }

  return(ce_result)
}


#' Get parameter information from bmmfit object
#'
#' @description
#' Internal helper function that extracts parameter type (dpar/nlpar),
#' link function, and brms parameter name from a bmmfit object.
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name from model specification
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`type`}{Either "dpar" or "nlpar"}
#'     \item{`model_name`}{Original parameter name from model specification}
#'     \item{`brms_name`}{Parameter name as used in brms (may be modified)}
#'     \item{`link`}{Link function for this parameter}
#'   }
#'
#' @keywords internal
#' @noRd
.get_parameter_info <- function(bmmfit, par) {
  model <- bmmfit$bmm$model
  model_pars <- names(model$parameters)

  # First check if it's a distributional parameter in brmsterms
  # This takes precedence because dpars can share names with model parameters
  bterms <- brms::brmsterms(bmmfit$formula)
  
  if (!is.null(bterms$dpars) && par %in% names(bterms$dpars)) {
    # Found as dpar - extract link from brmsterms
    dpar_info <- bterms$dpars[[par]]
    link <- .extract_dpar_link(dpar_info)

    return(list(
      type = "dpar",
      model_name = par,
      brms_name = par,
      link = link
    ))
  }

  # Then check if par is in model parameters (these are nlpars in bmm)
  if (par %in% model_pars) {
    # This is a non-linear parameter
    link <- model$links[[par]] %||% "identity"
    return(list(
      type = "nlpar",
      model_name = par,
      brms_name = par,
      link = link
    ))
  }

  # If we get here, parameter not found
  stop2(
    "Parameter '{par}' not found in model.\n",
    "Available parameters: {paste(model_pars, collapse = ', ')}"
  )
}


#' Extract link function from distributional parameter info
#'
#' @description
#' Internal function to extract link function from brmsterms dpar structure.
#' For dpars with a family object, extracts the link from the family.
#' For non-linear dpars (btnl class), returns "identity" as they don't have links.
#'
#' @param dpar_info A dpar element from brmsterms$dpars
#'
#' @return Character string. Link function name
#'
#' @keywords internal
#' @noRd
.extract_dpar_link <- function(dpar_info) {
  # For dpars (both btl and btnl classes), extract link from family
  if ("family" %in% names(dpar_info)) {
    link <- dpar_info$family$link
    if (!is.null(link) && nzchar(link)) {
      return(link)
    }
  }

  # Default to identity if no link found
  return("identity")
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
    return(ce_object) # No transformation needed
  }

  # Apply transformation to each data frame in the list
  # conditional_effects returns a list of data frames, one per effect
  for (i in seq_along(ce_object)) {
    df <- ce_object[[i]]

    # Transform estimate and confidence intervals using link_transform()
    # brms stores these as estimate__, lower__, upper__
    if ("estimate__" %in% names(df)) {
      df$estimate__ <- link_transform(df$estimate__, link, inverse = inverse)
    }
    if ("lower__" %in% names(df)) {
      df$lower__ <- link_transform(df$lower__, link, inverse = inverse)
    }
    if ("upper__" %in% names(df)) {
      df$upper__ <- link_transform(df$upper__, link, inverse = inverse)
    }

    ce_object[[i]] <- df
  }

  return(ce_object)
}
