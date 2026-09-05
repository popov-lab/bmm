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
