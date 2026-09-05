############################################################################# !
# METHODS FOR MULTIVARIATE BMM FITS                                      ####
############################################################################# !

#' @export
summary.mvbmmfit <- function(object, priors = FALSE, prob = 0.95, robust = FALSE,
                             mc_se = FALSE, ..., backend = "bmm") {
  object <- restructure(object)
  backend <- match.arg(backend, c("bmm", "brms"))

  # summary.bmmfit assumes a single model, so the brms summary is obtained
  # directly from the underlying brmsfit
  brms_fit <- object
  class(brms_fit) <- "brmsfit"
  out <- summary(brms_fit,
    priors = priors, prob = prob, robust = robust,
    mc_se = mc_se, ...
  )
  if (backend == "brms") {
    return(out)
  }

  out$components <- lapply(object$bmm$components, function(comp) {
    nlist(
      model = comp$model,
      formula = comp$user_formula,
      resp_name = comp$resp_name,
      data_name = comp$data_name,
      nobs = sum(object$data[[comp$subset_var]])
    )
  })
  class(out) <- "mvbmmsummary"
  out
}

#' @export
print.mvbmmsummary <- function(x, digits = 2,
                               color = getOption("bmm.color_summary", TRUE), ...) {
  options(bmm.color_summary = color)
  cat(style("green")(paste0(
    "Multivariate bmm model with ", length(x$components), " components\n\n"
  )))
  for (i in seq_along(x$components)) {
    comp <- x$components[[i]]
    cat(style("green")(paste0("Component ", i, " [", comp$resp_name, "]\n")))
    print_model_header(comp$model, comp$formula, comp$data_name %||% "", comp$nobs)
    cat("\n")
  }
  print_summary_draws(x)
  print_summary_random(x, digits)
  if (nrow(x$fixed)) {
    for (comp in x$components) {
      rows <- x$fixed[mv_component_rows(rownames(x$fixed), comp), , drop = FALSE]
      if (nrow(rows)) {
        print_summary_fixed(rows, digits, label = paste0(" [", comp$resp_name, "]"))
      }
    }
  }
  print_summary_footer(x)
  invisible(x)
}

# rows of the brms fixed-effects summary belonging to one component: brms
# prefixes coefficients of non-linear parameters with <resp>_<par> and
# suffixes distributional parameters with <par>_<resp>; the printed rows are
# additionally filtered to the component's model parameters (plus its
# response, which labels the mu coefficients of family components)
mv_component_rows <- function(rows, comp) {
  pars <- union(names(comp$model$parameters), names(comp$formula)[!is_nl(comp$formula)])
  in_component <- grepl(paste0("^", comp$resp_name, "_"), rows) |
    grepl(paste0("_", comp$resp_name, "_"), rows)
  par_pattern <- paste0("(^|_)(", paste(pars, collapse = "|"), ")_")
  in_component & grepl(par_pattern, rows)
}

#' @export
update.mvbmmfit <- function(object, ...) {
  stop2(
    "update() is not yet supported for multivariate bmm models. Please \\
    modify the bmm_component() specifications and refit the model with bmm()."
  )
}

#' @export
conditional_effects.mvbmmfit <- function(x, ...) {
  stop2(
    "conditional_effects() is not yet supported for multivariate bmm models, \\
    because bmm's conditional effects are computed on the model parameters \\
    and this is not yet implemented across components. You can obtain \\
    response-scale effects for one component from the underlying brmsfit:
      brmsfit <- your_fit
      class(brmsfit) <- 'brmsfit'
      brms::conditional_effects(brmsfit, resp = 'response_name')"
  )
}

#' @export
parameters.mvbmmfit <- function(x, ...) {
  x <- restructure(x)
  tables <- lapply(x$bmm$components, function(comp) {
    # family components without predicted distributional parameters have
    # nothing to report
    if (length(comp$model$parameters) == 0) {
      return(NULL)
    }
    table <- parameters(comp$model, formula = comp$user_formula, ...)
    table$response <- rep(comp$resp, nrow(table))
    table
  })
  tables <- Filter(Negate(is.null), tables)
  out <- do.call(rbind, tables)
  attr(out, "model_name") <- glue("Multivariate bmm model with {length(x$bmm$components)} components")
  out
}

#' @export
pp_check.mvbmmfit <- function(object, ..., resp = NULL) {
  resps <- vapply(object$bmm$components, function(x) x$resp_name, character(1))
  stopif(
    is.null(resp),
    "For multivariate bmm models, pp_check() requires the 'resp' argument to \\
    select a component. Available responses: {collapse_comma(resps)}"
  )
  stopif(
    not_in(resp, resps),
    "Unknown response '{resp}'. Available responses: {collapse_comma(resps)}"
  )
  comp <- object$bmm$components[[match(resp, resps)]]
  stopif(
    identical(comp$model$family$family %||% "", "multinomial") ||
      inherits(comp$model, "m3"),
    "pp_check() is not yet supported for multinomial components of \\
    multivariate bmm models."
  )
  class(object) <- "brmsfit"
  pp_check(object, resp = resp, ...)
}

#' @export
reset_env.mvbmmfit <- function(object, env = globalenv(), ...) {
  object$formula <- reset_env(object$formula, env)
  for (i in seq_along(object$bmm$components)) {
    object$bmm$components[[i]]$user_formula <-
      reset_env(object$bmm$components[[i]]$user_formula, env)
  }
  object
}

#' @export
reset_env.mvbrmsformula <- function(object, env = globalenv(), ...) {
  for (i in seq_along(object$forms)) {
    object$forms[[i]] <- reset_env(object$forms[[i]], env)
  }
  object
}
