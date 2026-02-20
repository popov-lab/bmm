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
#'   (see `names(x$bmm$model$parameters)`). If `NULL` (the default), conditional
#'   effects are computed for all estimated (non-fixed) parameters.
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
#' @return A `brms_conditional_effects` object (from brms), which can be:
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
#' @aliases conditional_effects
#' @method conditional_effects bmmfit
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a mixture model with set size effect on kappa
#' fit <- bmm(
#'   formula = bmf(kappa ~ setsize, thetat ~ 1),
#'   data = zhang_luck_2008,
#'   model = mixture3p(
#'     resp_error = "response_error",
#'     nt_features = paste0("col_lure", 1:5),
#'     set_size = "setsize"
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
#' ce_specific <- conditional_effects(fit, par = "kappa", effects = "setsize")
#'
#' # Combine with other brms options
#' ce_detailed <- conditional_effects(
#'   fit,
#'   par = "kappa",
#'   effects = "setsize",
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

  # If par is NULL, compute conditional effects for all estimated parameters
  if (is.null(par)) {
    model <- x$bmm$model
    estimated_pars <- setdiff(names(model$parameters),
                              names(model$fixed_parameters))
    all_effects <- list()
    for (p in estimated_pars) {
      ce <- conditional_effects(x, par = p, scale = scale, ...)
      if (length(ce) > 0) {
        names(ce) <- paste0(p, ".", names(ce))
        all_effects <- c(all_effects, ce)
      }
    }
    class(all_effects) <- c("brms_conditional_effects", "list")
    return(all_effects)
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
      softmax_result <- .filter_internal_effects(softmax_result, x)
      return(softmax_result)
    } else {
      warning2(
        "Parameter '{par}' uses multinomial logit transformation.\n",
        "Native scale display not available for this model configuration.\n",
        "Showing on sampling scale instead."
      )
      scale <- "sampling"
    }
  }

  # For multinomial family models (m3), brms requires categorical = TRUE
  # even when requesting a specific nlpar, which breaks the workflow.
  # Bypass by computing conditional effects via posterior_linpred directly.
  if ("m3" %in% class(x$bmm$model)) {
    ce_result <- .compute_multinomial_conditional_effects(x, par, ...)
  } else if (par_info$type == "dpar") {
    ce_result <- .brms_conditional_effects(x, dpar = par_info$brms_name, ...)
  } else if (par_info$type == "nlpar") {
    ce_result <- .brms_conditional_effects(x, nlpar = par_info$brms_name, ...)
  } else {
    stop2("Internal error: parameter type must be 'dpar' or 'nlpar'")
  }

  # Apply transformations based on parameter type and requested scale
  # - nlpars: brms returns on sampling scale, apply inverse link for native scale
  # - dpars: brms returns on native scale, apply forward link for sampling scale
  if (par_info$link != "identity") {
    if (par_info$type == "nlpar" && scale == "native") {
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = TRUE)
    } else if (par_info$type == "dpar" && scale %in% c("parameter", "sampling")) {
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = FALSE)
    }
  }

  # Filter out internal variables before returning
  ce_result <- .filter_internal_effects(ce_result, x)

  return(ce_result)
}


#' Call brms conditional_effects without infinite recursion
#'
#' @description
#' Strips the `"bmmfit"` class so that S3 dispatch reaches
#' `brms::conditional_effects.brmsfit()` instead of recursing back to
#' `conditional_effects.bmmfit()`.
#'
#' @param x A bmmfit object
#' @param ... Arguments forwarded to [brms::conditional_effects()]
#'
#' @return A `brms_conditional_effects` object
#'
#' @keywords internal
#' @noRd
.brms_conditional_effects <- function(x, ...) {
  class(x) <- class(x)[class(x) != "bmmfit"]
  conditional_effects(x, ...)
}


#' Filter internal variables from conditional_effects results
#'
#' @description
#' Removes conditional effects plots for internal model variables
#' (like LureIdx, Idx_*, inv_ss, etc.) that are created during data
#' preprocessing but are not part of the user's formula.
#'
#' @param ce_result A brms_conditional_effects object
#' @param bmmfit A bmmfit object
#'
#' @return Filtered conditional_effects object with only user-specified predictors
#'
#' @keywords internal
#' @noRd
.filter_internal_effects <- function(ce_result, bmmfit) {
  # Get user-specified variables from the data
  # Internal variables typically match patterns like:
  # - LureIdx* (mixture models)
  # - Idx_* (m3 models)
  # - inv_ss (mixture models)
  # - Item*_Col_rad, Item*_Pos_rad (internal feature variables)
  # - nt_features and nt_distances from model (mixture3p, imm, sdm models)
  
  internal_patterns <- c(
    "^LureIdx",
    "^Idx_",
    "^inv_ss$",
    "^Item[0-9]+_",
    "^expS$"
  )
  
  # Extract nt_features and nt_distances from model if available
  model <- bmmfit$bmm$model
  if (!is.null(model$other_vars$nt_features)) {
    # Escape special regex characters and add to patterns
    nt_features <- model$other_vars$nt_features
    internal_patterns <- c(internal_patterns, paste0("^", gsub("([.|()\\^{}+$*?[]\\\\])", "\\\\\\1", nt_features), "$"))
  }
  if (!is.null(model$other_vars$nt_distances)) {
    # Escape special regex characters and add to patterns
    nt_distances <- model$other_vars$nt_distances
    internal_patterns <- c(internal_patterns, paste0("^", gsub("([.|()\\^{}+$*?[]\\\\])", "\\\\\\1", nt_distances), "$"))
  }
  
  # Get the names of effects in the conditional_effects list
  effect_names <- names(ce_result)
  
  # Filter out effects that match internal patterns
  keep_effects <- sapply(effect_names, function(name) {
    # Extract the variable name from the effect name
    # Effect names can be like "variable", "var1:var2", etc.
    vars_in_effect <- strsplit(name, ":")[[1]]
    
    # Check if any variable matches internal patterns
    any_internal <- any(sapply(vars_in_effect, function(var) {
      any(sapply(internal_patterns, function(pattern) {
        grepl(pattern, var)
      }))
    }))
    
    # Keep if not internal
    !any_internal
  })
  
  # Filter the list
  if (any(keep_effects)) {
    ce_result <- ce_result[keep_effects]
    
    # Restore class and attributes
    class(ce_result) <- c("brms_conditional_effects", "list")
  }
  
  ce_result
}


#' Extract grouping variable names from random effects in a formula
#'
#' @description
#' Parses the RHS of a formula to identify random-effects grouping variables
#' that should be excluded from conditional effects. Handles all brms grouping
#' specifications:
#' \itemize{
#'   \item Bare names: `(1 | id)`, `(1 || id)`
#'   \item Correlation IDs: `(1 |ID1| id)` — excludes both `ID1` and `id`
#'   \item `gr()`: `(1 | gr(id, by = exp))` — extracts `id`, not `exp`
#'   \item `mm()`: `(1 | mm(g1, g2))` — extracts all positional args
#'   \item Crossed: `(1 | id:group)` — extracts both `id` and `group`
#' }
#'
#' @param formula A formula object
#'
#' @return Character vector of grouping variable names to exclude
#'
#' @keywords internal
#' @noRd
.extract_re_grouping_vars <- function(formula) {
  rhs_str <- paste(deparse(formula[[length(formula)]]), collapse = " ")

  # Match text after each | that is not itself | or )
  # This captures: bare grouping vars, correlation IDs, and gr()/mm() calls
  bar_parts <- regmatches(
    rhs_str, gregexpr("(?<=\\|)[^|)]+", rhs_str, perl = TRUE)
  )[[1]]
  bar_parts <- trimws(bar_parts)
  bar_parts <- bar_parts[nchar(bar_parts) > 0]

  if (length(bar_parts) == 0) return(character(0))

  unlist(lapply(bar_parts, function(part) {
    # gr(id, ...) — first argument is the grouping variable
    if (grepl("^gr\\s*\\(", part)) {
      inner <- sub("^gr\\s*\\(\\s*", "", part)
      return(trimws(sub("[,)]+.*", "", inner)))
    }
    # mm(g1, g2, ...) — positional args (before named args) are grouping vars
    if (grepl("^mm\\s*\\(", part)) {
      inner <- sub("^mm\\s*\\(\\s*", "", part)
      args <- trimws(strsplit(inner, ",")[[1]])
      return(args[!grepl("=", args)])
    }
    # Bare variable name(s) or correlation ID — split on : only
    trimws(strsplit(part, ":")[[1]])
  }))
}


#' Build a prediction grid for conditional effects
#'
#' @description
#' Constructs a prediction grid for computing conditional effects via
#' [brms::posterior_linpred()]. For each effect variable, creates a data frame
#' where that variable varies over its range (numeric) or levels (factor) while
#' all other columns are held at reference values (mean for numeric, first level
#' for factor).
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name whose formula determines the
#'   predictor variables.
#' @param effects Character vector. Specific effect variables to include. If
#'   `NULL`, all RHS variables from the parameter's formula are used.
#' @param resolution Integer. Number of points for numeric predictors (default
#'   100).
#'
#' @return A named list of data frames, one per effect variable. Empty list if
#'   no effects are found.
#'
#' @keywords internal
#' @noRd
.ce_prediction_grid <- function(bmmfit, par, effects = NULL, resolution = 100) {
  user_formula <- bmmfit$bmm$user_formula
  par_formula <- user_formula[[par]]
  if (is.null(par_formula)) return(list())

  # Extract fixed-effect predictors only, excluding random-effects grouping
  # variables and correlation IDs from all brms grouping patterns
  f <- stats::formula(par_formula)
  re_groups <- .extract_re_grouping_vars(f)
  rhs_vars <- all.vars(f[-2])
  rhs_vars <- setdiff(rhs_vars, c("0", "1", re_groups))

  if (is.null(effects)) {
    effect_vars <- rhs_vars
  } else {
    effect_vars <- unlist(strsplit(as.character(effects), ":"))
    effect_vars <- intersect(effect_vars, rhs_vars)
  }

  if (length(effect_vars) == 0) return(list())

  orig_data <- bmmfit$data
  grids <- list()

  for (var in effect_vars) {
    col <- orig_data[[var]]
    if (is.factor(col) || is.character(col)) {
      varying <- sort(unique(col))
    } else {
      rng <- range(col, na.rm = TRUE)
      varying <- seq(rng[1], rng[2], length.out = resolution)
    }

    newdata <- data.frame(x__ = varying)
    names(newdata) <- var

    for (v in setdiff(names(orig_data), var)) {
      cv <- orig_data[[v]]
      if (is.matrix(cv)) next
      if (is.factor(cv)) {
        newdata[[v]] <- factor(levels(cv)[1], levels = levels(cv))
      } else if (is.character(cv)) {
        newdata[[v]] <- cv[1]
      } else if (is.integer(cv)) {
        newdata[[v]] <- as.integer(round(stats::median(cv, na.rm = TRUE)))
      } else if (is.numeric(cv)) {
        newdata[[v]] <- mean(cv, na.rm = TRUE)
      } else {
        newdata[[v]] <- cv[1]
      }
    }

    grids[[var]] <- newdata
  }

  grids
}


#' Summarize posterior draws into conditional-effect statistics
#'
#' @description
#' Takes a draws matrix (n_draws x n_grid_points) and computes summary
#' statistics suitable for `brms_conditional_effects` data frames.
#'
#' @param draws Matrix. Posterior draws (rows = draws, columns = grid points).
#' @param prob Numeric. Probability mass for credible intervals (default 0.95).
#' @param robust Logical. If `TRUE`, use median/MAD instead of mean/SD.
#'
#' @return A list with elements `estimate`, `lower`, `upper`, `se` — each a
#'   numeric vector of length `ncol(draws)`.
#'
#' @keywords internal
#' @noRd
.ce_summarize_draws <- function(draws, prob = 0.95, robust = FALSE) {
  probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  if (robust) {
    estimate <- apply(draws, 2, stats::median)
    se <- apply(draws, 2, stats::mad)
  } else {
    estimate <- colMeans(draws)
    se <- apply(draws, 2, stats::sd)
  }
  lower <- apply(draws, 2, stats::quantile, probs = probs[1])
  upper <- apply(draws, 2, stats::quantile, probs = probs[2])
  list(estimate = estimate, lower = lower, upper = upper, se = se)
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
  # Currently only implemented for mixture3p
  if (!"mixture3p" %in% class(bmmfit$bmm$model)) {
    return(NULL)
  }

  # Get conditional_effects for the requested par to obtain the grid
  ce_par <- .brms_conditional_effects(bmmfit, nlpar = par, ...)

  # Extract options from ... for posterior computation
  dots <- list(...)
  prob <- dots$prob %||% 0.95
  robust <- dots$robust %||% FALSE
  re_formula <- dots$re_formula %||% NA
  ndraws <- dots$ndraws

  # For each effect, compute softmax-transformed posterior draws
  result <- lapply(ce_par, function(df) {
    # Extract newdata grid (remove brms-internal *__ columns)
    internal_cols <- grep("__$", names(df), value = TRUE)
    newdata <- df[, !names(df) %in% internal_cols, drop = FALSE]

    # Get posterior draws for both theta params on sampling scale
    linpred_args <- list(
      object = bmmfit,
      newdata = newdata,
      re_formula = re_formula,
      allow_new_levels = TRUE
    )
    if (!is.null(ndraws)) linpred_args$ndraws <- ndraws

    draws_t <- do.call(
      brms::posterior_linpred,
      c(linpred_args, list(nlpar = "thetat"))
    )
    draws_nt <- do.call(
      brms::posterior_linpred,
      c(linpred_args, list(nlpar = "thetant"))
    )

    # Apply softmax jointly per draw
    # softmax denom = exp(thetat) + exp(thetant) + exp(0)
    denom <- exp(draws_t) + exp(draws_nt) + 1
    if (par == "thetat") {
      softmax_draws <- exp(draws_t) / denom
    } else {
      softmax_draws <- exp(draws_nt) / denom
    }

    summ <- .ce_summarize_draws(softmax_draws, prob = prob, robust = robust)
    df$estimate__ <- summ$estimate
    df$lower__ <- summ$lower
    df$upper__ <- summ$upper

    df
  })

  names(result) <- names(ce_par)
  class(result) <- class(ce_par)
  result
}


#' Compute conditional effects for multinomial family models
#'
#' @description
#' For models using `brms::multinomial()` family (e.g., m3), brms requires
#' `categorical = TRUE` even when requesting a specific nlpar, which conflicts
#' with nlpar-level computation. This helper bypasses that check by using
#' `brms::posterior_linpred()` directly with a manually constructed prediction
#' grid.
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name (nlpar) to compute effects for
#' @param ... Additional arguments (prob, robust, re_formula, ndraws, effects,
#'   resolution)
#'
#' @return A `brms_conditional_effects` object with one element per effect
#'
#' @keywords internal
#' @noRd
.compute_multinomial_conditional_effects <- function(bmmfit, par, ...) {
  dots <- list(...)
  prob <- dots$prob %||% 0.95
  robust <- dots$robust %||% FALSE
  re_formula <- dots$re_formula %||% NA
  ndraws <- dots$ndraws
  resolution <- dots$resolution %||% 100
  effects <- dots$effects

  # Build prediction grids using shared helper
  grids <- .ce_prediction_grid(bmmfit, par,
                               effects = effects,
                               resolution = resolution)
  if (length(grids) == 0) {
    return(structure(list(), class = c("brms_conditional_effects", "list")))
  }

  result <- list()

  for (var in names(grids)) {
    newdata <- grids[[var]]

    # Get posterior draws for the nlpar
    linpred_args <- list(
      object = bmmfit,
      newdata = newdata,
      nlpar = par,
      re_formula = re_formula,
      allow_new_levels = TRUE
    )
    if (!is.null(ndraws)) linpred_args$ndraws <- ndraws

    draws <- do.call(brms::posterior_linpred, linpred_args)

    # Summarize draws using shared helper
    summ <- .ce_summarize_draws(draws, prob = prob, robust = robust)
    newdata$estimate__ <- summ$estimate
    newdata$lower__ <- summ$lower
    newdata$upper__ <- summ$upper
    newdata$se__ <- summ$se
    newdata$effect1__ <- newdata[[var]]
    newdata$cond__ <- factor("1")

    attr(newdata, "effects") <- var
    attr(newdata, "response") <- par

    result[[var]] <- newdata
  }

  class(result) <- c("brms_conditional_effects", "list")
  result
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

  # Apply transformation to each data frame in the list using vectorized lapply
  # conditional_effects returns a list of data frames, one per effect
  result <- lapply(ce_object, function(df) {
    # Transform estimate and confidence intervals using link_transform()
    # brms stores these as estimate__, lower__, upper__
    df$estimate__ <- link_transform(df$estimate__, link, inverse = inverse)
    df$lower__ <- link_transform(df$lower__, link, inverse = inverse)
    df$upper__ <- link_transform(df$upper__, link, inverse = inverse)
    df
  })
  
  names(result) <- names(ce_object)
  class(result) <- class(ce_object)
  result
}
