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
#' @method conditional_effects bmmfit
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

  # For mixture models with multinomial logit, compute softmax manually
  if (par_info$multinomial && scale == "native") {
    softmax_result <- .compute_softmax_conditional_effects(x, par, ...)
    if (!is.null(softmax_result)) {
      return(softmax_result)
    } else {
      # Fallback to warning if softmax computation not available
      warning2(
        "Parameter '{par}' uses multinomial logit transformation.\n",
        "Native scale display not available for this model configuration.\n",
        "Showing on sampling scale instead."
      )
      scale <- "sampling"
    }
  }

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
#'     \item{`multinomial`}{Logical indicating if parameter uses multinomial logit}
#'   }
#'
#' @keywords internal
#' @noRd
.get_parameter_info <- function(bmmfit, par) {
  model <- bmmfit$bmm$model
  model_pars <- names(model$parameters)
  bterms <- brms::brmsterms(bmmfit$formula)
  
  # Check if parameter exists in model
  if (!par %in% model_pars) {
    stop2(
      "Parameter '{par}' not found in model.\n",
      "Available parameters: {paste(model_pars, collapse = ', ')}"
    )
  }
  
  # Get link from bmm model (authoritative source for user-specified parameters)
  link <- model$links[[par]] %||% "identity"
  
  # Check if this is a multinomial logit parameter
  # These appear in mixture3p (thetat, thetant) and imm models
  # and cannot be meaningfully transformed in isolation
  multinomial <- .is_multinomial_param(par, model)
  
  # Determine if parameter appears as dpar or nlpar in brms
  if (!is.null(bterms$dpars) && par %in% names(bterms$dpars)) {
    type <- "dpar"
  } else if (!is.null(bterms$nlpars) && par %in% names(bterms$nlpars)) {
    type <- "nlpar"
  } else {
    # Parameter exists in model but not in brmsterms - try nlpar
    type <- "nlpar"
  }

  return(list(
    type = type,
    model_name = par,
    brms_name = par,
    link = link,
    multinomial = multinomial
  ))
}


#' Check if parameter uses multinomial logit transformation
#'
#' @description
#' Determines if a parameter uses multinomial logit (softmax) transformation.
#' These parameters cannot be meaningfully transformed individually because
#' they are interdependent.
#'
#' @param par Character string. Parameter name
#' @param model A bmmodel object
#'
#' @return Logical. TRUE if parameter uses multinomial logit
#'
#' @keywords internal
#' @noRd
.is_multinomial_param <- function(par, model) {
  model_class <- class(model)
  
  # mixture3p: thetat and thetant use multinomial logit
  if ("mixture3p" %in% model_class && par %in% c("thetat", "thetant")) {
    return(TRUE)
  }
  
  # IMM models don't use multinomial logit for their theta-like parameters
  # They use activation-based mixture weights
  
  return(FALSE)
}


#' Compute softmax transformation for multinomial parameters
#'
#' @description
#' For mixture models with multinomial logit (softmax), manually computes
#' the softmax transformation by extracting conditional effects for all
#' relevant nlpars and applying the softmax formula.
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name to return (thetat or thetant)
#' @param ... Additional arguments passed to brms::conditional_effects()
#'
#' @return A brms_conditional_effects object with softmax-transformed values
#'
#' @keywords internal
#' @noRd
.compute_softmax_conditional_effects <- function(bmmfit, par, ...) {
  model <- bmmfit$bmm$model
  model_class <- class(model)
  
  # Currently only implemented for mixture3p
  if (!"mixture3p" %in% model_class) {
    return(NULL)
  }
  
  # Extract conditional effects for all theta parameters on sampling scale
  # thetat and thetant are nlpars, thetaguess is fixed to 0
  ce_thetat <- brms:::conditional_effects.brmsfit(bmmfit, nlpar = "thetat", ...)
  ce_thetant <- brms:::conditional_effects.brmsfit(bmmfit, nlpar = "thetant", ...)
  
  # Apply softmax: prob_i = exp(theta_i) / sum(exp(theta_j))
  # For mixture3p: thetat, thetant, thetaguess (fixed at 0)
  # softmax denominator = exp(thetat) + exp(thetant) + exp(0) = exp(thetat) + exp(thetant) + 1
  
  # Process each effect in the conditional_effects list
  for (i in seq_along(ce_thetat)) {
    df_thetat <- ce_thetat[[i]]
    df_thetant <- ce_thetant[[i]]
    
    # Check that both have same structure
    if (nrow(df_thetat) != nrow(df_thetant)) {
      stop2("Mismatch in conditional_effects structure between thetat and thetant")
    }
    
    # Extract values on sampling scale
    thetat_est <- df_thetat$estimate__
    thetat_lower <- if ("lower__" %in% names(df_thetat)) df_thetat$lower__ else NULL
    thetat_upper <- if ("upper__" %in% names(df_thetat)) df_thetat$upper__ else NULL
    
    thetant_est <- df_thetant$estimate__
    thetant_lower <- if ("lower__" %in% names(df_thetant)) df_thetant$lower__ else NULL
    thetant_upper <- if ("upper__" %in% names(df_thetant)) df_thetant$upper__ else NULL
    
    # Compute softmax denominator (sum of exp for all components)
    denom_est <- exp(thetat_est) + exp(thetant_est) + 1  # 1 = exp(0) for thetaguess
    denom_lower <- if (!is.null(thetat_lower) && !is.null(thetant_lower)) {
      exp(thetat_lower) + exp(thetant_lower) + 1
    } else NULL
    denom_upper <- if (!is.null(thetat_upper) && !is.null(thetant_upper)) {
      exp(thetat_upper) + exp(thetant_upper) + 1
    } else NULL
    
    # Select which parameter to return based on 'par'
    if (par == "thetat") {
      # Probability of target response
      df_result <- df_thetat
      df_result$estimate__ <- exp(thetat_est) / denom_est
      if (!is.null(thetat_lower)) {
        df_result$lower__ <- exp(thetat_lower) / denom_lower
      }
      if (!is.null(thetat_upper)) {
        df_result$upper__ <- exp(thetat_upper) / denom_upper
      }
    } else if (par == "thetant") {
      # Probability of swapping to non-targets
      # Note: In mixture3p, thetant already includes log(1/(set_size-1))
      # So exp(thetant) gives us the TOTAL swap probability directly
      # We just need to apply softmax
      df_result <- df_thetant
      df_result$estimate__ <- exp(thetant_est) / denom_est
      if (!is.null(thetant_lower)) {
        df_result$lower__ <- exp(thetant_lower) / denom_lower
      }
      if (!is.null(thetant_upper)) {
        df_result$upper__ <- exp(thetant_upper) / denom_upper
      }
    } else {
      stop2("Unknown parameter: {par}")
    }
    
    ce_thetat[[i]] <- df_result
  }
  
  return(ce_thetat)
}


#' Map nlpar to corresponding dpar for native scale display
#'
#' @description
#' For mixture models, nlpars like thetat/thetant use multinomial logit and
#' can't be transformed individually. However, brms creates corresponding dpars
#' (theta1, theta2, etc.) that represent the mixture weights on probability scale.
#' This function maps user-facing nlpars to their corresponding dpars and
#' optionally provides a post-processing function.
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name (nlpar)
#'
#' @return List with elements:
#'   - dpar: Character string with dpar name
#'   - postprocess: Optional function(ce_result, bmmfit) for post-processing
#'   Returns NULL if no mapping available.
#'
#' @keywords internal
#' @noRd
.map_nlpar_to_dpar <- function(bmmfit, par) {
  model <- bmmfit$bmm$model
  model_class <- class(model)
  
  # mixture3p: map to theta dpars
  if ("mixture3p" %in% model_class) {
    # theta1 = target responses (from thetat)
    # theta2 = guessing responses (fixed at 0)
    # theta3+ = non-target responses (from thetant)
    if (par == "thetat") {
      return(list(dpar = "theta1", postprocess = NULL))
    } else if (par == "thetant") {
      # theta3 represents probability of swapping to ONE specific lure
      # Need to multiply by (set_size - 1) to get total swap probability
      return(list(
        dpar = "theta3",
        postprocess = .postprocess_thetant_mixture3p
      ))
    }
  }
  
  # Add more mappings here for other mixture models if needed
  
  return(NULL)
}


#' Post-process thetant for mixture3p models
#'
#' @description
#' Multiplies theta3 (probability of swapping to one specific lure) by
#' (set_size - 1) to get the total probability of swapping to any non-target.
#'
#' @param ce_result Conditional effects object from brms
#' @param bmmfit A bmmfit object
#'
#' @return Modified conditional effects object
#'
#' @keywords internal
#' @noRd
.postprocess_thetant_mixture3p <- function(ce_result, bmmfit) {
  set_size_var <- bmmfit$bmm$model$other_vars$set_size
  
  # Process each effect in the list
  for (i in seq_along(ce_result)) {
    df <- ce_result[[i]]
    
    # Determine set_size for each row
    if (set_size_var %in% names(df)) {
      # set_size varies in the data
      ss_values <- df[[set_size_var]]
    } else if (set_size_var %in% names(attr(df, "conditions"))) {
      # set_size is in conditions (constant across this effect)
      ss_values <- rep(attr(df, "conditions")[[set_size_var]], nrow(df))
    } else {
      # Try to get from cond__ or other attributes
      # Fall back to extracting from bmmfit data
      data <- bmmfit$data
      if (set_size_var %in% names(data)) {
        # Use the unique set_size values or mean if variable
        ss_unique <- unique(data[[set_size_var]])
        if (length(ss_unique) == 1) {
          ss_values <- rep(ss_unique, nrow(df))
        } else {
          warning2(
            "set_size varies in data but not specified in conditional_effects.\n",
            "Using mean set_size for thetant scaling. Consider specifying set_size explicitly."
          )
          ss_values <- rep(mean(as.numeric(as.character(data[[set_size_var]])), na.rm = TRUE), nrow(df))
        }
      } else {
        stop2(
          "Cannot determine set_size for thetant post-processing.\n",
          "Expected variable '{set_size_var}' not found in data or conditions."
        )
      }
    }
    
    # Convert set_size to numeric if it's a factor
    if (is.factor(ss_values)) {
      ss_values <- as.numeric(as.character(ss_values))
    } else if (!is.numeric(ss_values)) {
      ss_values <- as.numeric(ss_values)
    }
    
    # Multiply probability columns by (set_size - 1)
    multiplier <- ss_values - 1
    df$estimate__ <- df$estimate__ * multiplier
    if ("lower__" %in% names(df)) {
      df$lower__ <- df$lower__ * multiplier
    }
    if ("upper__" %in% names(df)) {
      df$upper__ <- df$upper__ * multiplier
    }
    if ("se__" %in% names(df)) {
      df$se__ <- df$se__ * multiplier
    }
    
    ce_result[[i]] <- df
  }
  
  return(ce_result)
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
