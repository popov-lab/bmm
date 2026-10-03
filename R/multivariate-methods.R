############################################################################# !
# METHODS FOR MULTIVARIATE BMM FITS                                      ####
############################################################################# !

# summary.bmmfit assumes a single model, so a multivariate fit is summarised
# as the underlying brmsfit
#' @export
summary.mvbmmfit <- function(object, ...) {
  object <- restructure(object)
  class(object) <- "brmsfit"
  summary(object, ...)
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

# Functions that read the single model in fit$bmm$model stop here: a
# multivariate fit stores its models in fit$bmm$components, and reporting on
# the missing model would describe no component or the wrong one
refuse_mvbmmfit <- function(x, fun) {
  stopif(
    inherits(x, "mvbmmfit"),
    "{fun}() is not yet supported for multivariate bmm models."
  )
}

#' @export
model_citation.mvbmmfit <- function(x, ...) {
  stop2(
    "model_citation() is not yet supported for multivariate bmm models. Call \\
    it on the model of each component instead, e.g. \\
    model_citation(fit$bmm$components[[1]]$model)."
  )
}

# The expected response of a component is refused as it is for the same model
# fitted alone; brms itself requires `resp` for a model that uses subset()
#' @export
posterior_epred.mvbmmfit <- function(object, ..., dpar = NULL, nlpar = NULL) {
  if (is.null(dpar) && is.null(nlpar)) {
    refuse_undefined_component_epred(object, list(...)$resp)
  }
  NextMethod()
}

#' @export
fitted.mvbmmfit <- function(object, ..., scale = c("response", "linear"),
                            dpar = NULL, nlpar = NULL) {
  if (match.arg(scale) == "response" && is.null(dpar) && is.null(nlpar)) {
    refuse_undefined_component_epred(object, list(...)$resp)
  }
  NextMethod()
}

refuse_undefined_component_epred <- function(object, resp) {
  for (comp in object$bmm$components) {
    if (is.null(resp) || comp$resp_name %in% resp) {
      refuse_undefined_epred(list(bmm = list(model = comp$model)))
    }
  }
}

#' @export
parameter_info.mvbmmfit <- function(x, ...) {
  x <- restructure(x)
  tables <- lapply(x$bmm$components, function(comp) {
    # family components without predicted distributional parameters have
    # nothing to report
    if (length(comp$model$parameters) == 0) {
      return(NULL)
    }
    table <- parameter_info(comp$model, formula = comp$user_formula, ...)
    table$response <- rep(comp$resp, nrow(table))
    table
  })
  tables <- Filter(Negate(is.null), tables)
  out <- do.call(rbind, tables)
  attr(out, "model_name") <- glue("Multivariate bmm model with {length(x$bmm$components)} components")
  out
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
