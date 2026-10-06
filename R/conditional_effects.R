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
#'       logit (mixture2p) or softmax (mixture3p) link shown on the probability scale.}
#'     \item{`"sampling"`}{Show on the sampling scale (as used during MCMC).
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
#'   formula = bmf(kappa ~ 0 + setsize, thetat ~ 1),
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
                                       scale = c("native", "sampling"),
                                       ...) {
  x <- restructure(x)
  scale <- match.arg(scale)

  if (is.null(par)) {
    .ce_all_parameters(x, scale, ...)
  } else {
    stopif(
      !is.character(par) || length(par) != 1,
      "Argument 'par' must be a single character string"
    )
    .ce_single_parameter(x, par, scale, ...)
  }
}


.ce_all_parameters <- function(x, scale, ...) {
  model <- x$bmm$model
  estimated_pars <- setdiff(names(model$parameters),
                            names(model$fixed_parameters))
  all_effects <- list()
  for (p in estimated_pars) {
    ce <- .ce_single_parameter(x, par = p, scale = scale, ...)
    if (length(ce) > 0) {
      names(ce) <- paste0(p, ".", names(ce))
      all_effects <- c(all_effects, ce)
    }
  }
  class(all_effects) <- c("brms_conditional_effects")
  all_effects
}


.ce_single_parameter <- function(x, par, scale, ...) {
  par_info <- .get_parameter_info(x, par)

  if (par_info$softmax && scale == "native") {
    softmax_result <- .compute_softmax_conditional_effects(x, par, ...)
    if (!is.null(softmax_result)) {
      .filter_internal_effects(softmax_result, x)
    } else {
      warning2(
        "Parameter '{par}' uses softmax transformation.\n",
        "Native scale display not available for this model configuration.\n",
        "Showing on sampling scale instead."
      )
      .ce_compute_and_transform(x, par, par_info, "sampling", ...)
    }
  } else {
    .ce_compute_and_transform(x, par, par_info, scale, ...)
  }
}


.ce_compute_and_transform <- function(x, par, par_info, scale, ...) {
  ce_result <- if (par_info$type == "dpar") {
    .brms_conditional_effects(x, dpar = par_info$brms_name, ...)
  } else if (par_info$type == "nlpar" && .has_category_dpars(x)) {
    .ce_nlpar_category_family(x, par, par_info$brms_name, ...)
  } else if (par_info$type == "nlpar") {
    .brms_conditional_effects(x, nlpar = par_info$brms_name, ...)
  } else {
    stop2("Internal error: parameter type must be 'dpar' or 'nlpar'")
  }

  # nlpars: brms returns on sampling scale; dpars: on native scale
  if (par_info$link != "identity") {
    if (par_info$type == "nlpar" && scale == "native") {
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = TRUE)
    } else if (par_info$type == "dpar" && scale == "sampling") {
      ce_result <- .apply_link_transform(ce_result, par_info$link, inverse = FALSE)
    }
  }

  .filter_internal_effects(ce_result, x)
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
  internal_patterns <- c(
    "^LureIdx",
    "^Idx_",
    "^inv_ss$",
    "^Item[0-9]+_",
    "^expS$"
  )
  
  model <- bmmfit$bmm$model
  if (!is.null(model$other_vars$nt_features)) {
    nt_features <- model$other_vars$nt_features
    escaped <- gsub("([][(){}^$*+?.|\\\\])", "\\\\\\1", nt_features)
    internal_patterns <- c(internal_patterns, paste0("^", escaped, "$"))
  }
  if (!is.null(model$other_vars$nt_distances)) {
    nt_distances <- model$other_vars$nt_distances
    escaped <- gsub("([][(){}^$*+?.|\\\\])", "\\\\\\1", nt_distances)
    internal_patterns <- c(internal_patterns, paste0("^", escaped, "$"))
  }
  
  effect_names <- names(ce_result)
  combined_pattern <- paste(internal_patterns, collapse = "|")
  keep_effects <- !vapply(effect_names, function(name) {
    vars <- strsplit(name, ":")[[1]]
    any(grepl(combined_pattern, vars))
  }, logical(1))

  if (any(keep_effects)) {
    ce_result <- ce_result[keep_effects]
    class(ce_result) <- c("brms_conditional_effects")
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

  if (length(bar_parts) == 0) {
    character(0)
  } else {
    unlist(lapply(bar_parts, function(part) {
      if (grepl("^gr\\s*\\(", part)) {
        # gr(id, ...) — first argument is the grouping variable
        inner <- sub("^gr\\s*\\(\\s*", "", part)
        trimws(sub("[,)]+.*", "", inner))
      } else if (grepl("^mm\\s*\\(", part)) {
        # mm(g1, g2, ...) — positional args (before named args) are grouping vars
        inner <- sub("^mm\\s*\\(\\s*", "", part)
        args <- trimws(strsplit(inner, ",")[[1]])
        args[!grepl("=", args)]
      } else {
        # Bare variable name(s) or correlation ID — split on : only
        trimws(strsplit(part, ":")[[1]])
      }
    }))
  }
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
  if (!"mixture3p" %in% class(bmmfit$bmm$model)) {
    NULL
  } else {
    .compute_softmax_ce_inner(bmmfit, par, ...)
  }
}


.compute_softmax_ce_inner <- function(bmmfit, par, ...) {
  ce_par <- .brms_conditional_effects(bmmfit, nlpar = par, ...)

  dots <- list(...)
  prob <- dots$prob %||% 0.95
  robust <- dots$robust %||% FALSE
  re_formula <- dots$re_formula %||% NA
  ndraws <- dots$ndraws

  result <- lapply(ce_par, function(df) {
    internal_cols <- grep("__$", names(df), value = TRUE)
    newdata <- df[, !names(df) %in% internal_cols, drop = FALSE]

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

    # numerically stable softmax: subtract max before exponentiating
    shift <- pmax(draws_t, draws_nt, 0)
    exp_t <- exp(draws_t - shift)
    exp_nt <- exp(draws_nt - shift)
    exp_0 <- exp(-shift)
    denom <- exp_t + exp_nt + exp_0
    if (par == "thetat") {
      softmax_draws <- exp_t / denom
    } else {
      softmax_draws <- exp_nt / denom
    }

    summ <- .ce_summarize_draws(softmax_draws, prob = prob, robust = robust)
    df$estimate__ <- summ$estimate
    df$lower__ <- summ$lower
    df$upper__ <- summ$upper
    df$se__ <- summ$se

    df
  })

  names(result) <- names(ce_par)
  class(result) <- class(ce_par)
  result
}


#' Does the fit's family have one distributional parameter per category?
#'
#' @description
#' These are the families for which [brms::conditional_effects()] refuses a
#' request for a non-linear parameter and asks for `categorical = TRUE`
#' (`conv_cats_dpars()` in brms 2.23.0): **m3** and the **sdt_rating**,
#' **sdt_cdp** and **sdt_ranking** models, which all use `brms::multinomial()`.
#'
#' @param x A bmmfit object
#'
#' @return Logical scalar
#'
#' @keywords internal
#' @noRd
.has_category_dpars <- function(x) {
  x$family$family %in% c("categorical", "multinomial", "dirichlet",
                         "dirichlet2", "logistic_normal")
}


#' Conditional effects of a non-linear parameter in a family with category dpars
#'
#' @description
#' brms computes the predictions of an `nlpar` for these families and only then
#' refuses them, because it expects the categories to be plotted. A non-linear
#' parameter does not depend on the categories, so brms is asked for one of the
#' category dpars instead, at a single draw, which gives the conditions grid
#' with `effects`, `conditions`, `int_conditions`, `re_formula`, `surface` and
#' `resolution` applied. The parameter is then evaluated on that grid and
#' summarised as brms would.
#'
#' Without `effects`, brms returns the effects of every formula in the model;
#' only those built from the parameter's own predictors are kept.
#'
#' @param x A bmmfit object
#' @param par Character string. The bmm parameter name.
#' @param nlpar Character string. The brms name of the parameter.
#' @param ... Arguments of [brms::conditional_effects()]; those it does not
#'   name itself reach [brms::posterior_linpred()], as in brms.
#'
#' @return A `brms_conditional_effects` object
#'
#' @keywords internal
#' @noRd
.ce_nlpar_category_family <- function(x, par, nlpar, effects = NULL,
                                      conditions = NULL, int_conditions = NULL,
                                      re_formula = NA, spaghetti = FALSE,
                                      surface = FALSE, resolution = 100,
                                      select_points = 0, too_far = 0,
                                      prob = 0.95, probs = NULL, robust = TRUE,
                                      ndraws = NULL, draw_ids = NULL, ...) {
  grid <- .brms_conditional_effects(
    x,
    dpar = names(brms::brmsterms(x$formula)$dpars)[1],
    effects = effects, conditions = conditions,
    int_conditions = int_conditions, re_formula = re_formula,
    spaghetti = spaghetti, surface = surface, resolution = resolution,
    select_points = select_points, too_far = too_far, draw_ids = 1, ...
  )

  if (is.null(effects)) {
    own_vars <- .np_grid_vars(x, par, re_formula = NA)
    grid <- grid[vapply(grid, function(ce) all(attr(ce, "effects") %in% own_vars),
                        logical(1))]
  }

  probs <- probs %||% c((1 - prob) / 2, 1 - (1 - prob) / 2)
  out <- lapply(grid, function(ce) {
    cond_data <- ce[setdiff(names(ce), c("estimate__", "se__", "lower__", "upper__"))]
    draws <- brms::posterior_linpred(
      x, newdata = cond_data, nlpar = nlpar, re_formula = re_formula,
      allow_new_levels = TRUE, ndraws = ndraws, draw_ids = draw_ids, ...
    )
    summ <- brms::posterior_summary(draws, probs = probs, robust = robust)
    ce$estimate__ <- summ[, 1]
    ce$se__ <- summ[, 2]
    ce$lower__ <- summ[, 3]
    ce$upper__ <- summ[, 4]
    attr(ce, "response") <- nlpar
    if (!is.null(attr(ce, "spaghetti"))) {
      attr(ce, "spaghetti") <- .ce_spaghetti(cond_data, draws, attr(ce, "effects"))
    }
    ce
  })
  structure(out, class = "brms_conditional_effects")
}


#' Spaghetti lines from the draws of a conditional effect
#'
#' @description
#' Builds the `spaghetti` attribute the way `brms` does: one line per draw,
#' and per level of the second effect when there is one.
#'
#' @param cond_data The conditions grid, without the summary columns
#' @param draws Matrix of draws (rows = draws, columns = grid rows)
#' @param effects Character vector of the effect names
#'
#' @return A data frame with the grid repeated per draw, plus `estimate__` and
#'   `sample__`
#'
#' @keywords internal
#' @noRd
.ce_spaghetti <- function(cond_data, draws, effects) {
  sample <- rep(seq_len(nrow(draws)), each = ncol(draws))
  if (length(effects) == 2L) {
    sample <- paste0(sample, "_", cond_data[[effects[2]]])
  }
  cbind(
    cond_data[rep(seq_len(nrow(cond_data)), times = nrow(draws)), , drop = FALSE],
    estimate__ = as.numeric(t(draws)),
    sample__ = factor(sample)
  )
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
    ce_object
  } else {
    .apply_link_transform_inner(ce_object, link, inverse)
  }
}


.apply_link_transform_inner <- function(ce_object, link, inverse) {
  result <- lapply(ce_object, function(df) {
    df$estimate__ <- link_transform(df$estimate__, link, inverse = inverse)
    df$lower__ <- link_transform(df$lower__, link, inverse = inverse)
    df$upper__ <- link_transform(df$upper__, link, inverse = inverse)
    df
  })
  
  names(result) <- names(ce_object)
  class(result) <- class(ce_object)
  result
}
