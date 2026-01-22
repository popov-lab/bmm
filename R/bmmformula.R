#' @title Create formula for predicting parameters of a `bmmodel`
#' @description This function is used to specify the formulas predicting the
#' different parameters of a `bmmodel`.
#'
#' @aliases bmf
#'
#' @details # General formula structure
#'
#' The formula argument accepts formulas of the following syntax:
#'
#'   ``` r
#'   parameter ~ fixed_effects + (random_effects | grouping_variable)
#'   ```
#'
#' `bmm` formulas are built on `brms` formulas and function in nearly the same
#' way, so you can use most of the `brms` formula syntax. The main differences
#' is that in `bmm` formulas, the response variable is not specified in the
#' formula. Instead, each parameter of the model is explicitly specified as the
#' left-hand side of the formula. In `brms`, the response variable is always
#' specified as the left-hand side of the first formula, which implicitly means
#' that any predictors in the first formula are predictors of the `mu` parameter
#' of the model. In general, measurement models do not all have a `mu`
#' parameter, therefore it is more straigthforward to explicitely predict each
#' parameter of the model.
#'
#' For example, in the following `brms` formula for the drift diffusion model,
#' the first line corresponds to the drift rate parameter, but this is not
#' explicitely stated.
#'
#'           brmsformula(rt | dec(response) ~ condition + (condition | id),
#'                       bs ~ 1 + (1 | id),
#'                       ndt ~ 1 + (1 | id),
#'                       bias ~ 1 + (1 | id))
#'
#' In `bmm`, the same formula would be written as:
#'
#'           bmmformula(drift ~ condition + (condition | id),
#'                      bs ~ 1 + (1 | id),
#'                      ndt ~ 1 + (1 | id),
#'                      bias ~ 1 + (1 | id))
#'
#' and the rt and response variables would be specified in the model argument of
#' the `bmm()` function.
#'
#' Aside from that, the `bmm` formula syntax is the same as the `brms` formula
#' syntax. For more information on the `brms` formula syntax, see
#' [brms::brmsformula()].
#'
#' You can also use the `bmf()` function as a shorthand for `bmmformula()`.
#'
#' You can also set parameters to a constant value by using par = value syntax:
#'
#'           bmf(drift ~ condition + (condition | id),
#'               bs ~ 1 + (1 | id),
#'               ndt ~ 1 + (1 | id),
#'               bias = 0.5)
#'
#' in which case the bias parameter will not be estimated but it will be fixed to 0.5
#'
#' @param ... Formulas for predicting a `bmmodel` parameter. Each formula for a
#' parameter should be specified as a separate argument, separated by commas
#' @return A list of formulas for each parameters being predicted
#' @export
#' @examples
#' imm_formula <- bmmformula(
#'   c ~ 0 + set_size + (0 + set_size | id),
#'   a ~ 1,
#'   kappa ~ 0 + set_size + (0 + set_size | id)
#' )
#'
#' # or use the shorter alias 'bmf'
#' imm_formula2 <- bmf(
#'   c ~ 0 + set_size + (0 + set_size | id),
#'   a ~ 1,
#'   kappa ~ 0 + set_size + (0 + set_size | id)
#' )
#' identical(imm_formula, imm_formula2)
bmmformula <- function(...) {
  out <- list(...)
  stopif(!all(vapply(out, is_form_or_num, logical(1))), "Arguments must be formulas or numeric assignments.")
  out <- out[!vapply(out, is_null_formula, logical(1))]

  par_names <- sapply(seq_along(out), function(i) {
    if (is_formula(out[[i]])) lhs_vars(out[[i]]) else names(out[i])
  })

  stopif(any(vapply(par_names, length, integer(1)) == 0), "Formulas must have a left-hand-side variable")
  stopif(any(duplicated(par_names)), "Duplicate formula for parameter(s) {par_names[duplicates]}")
  new_bmmformula(setNames(out, par_names))
}

#' @rdname bmmformula
#' @export
bmf <- bmmformula

# an internal bmmformula constructor
new_bmmformula <- function(x = nlist()) {
  stopifnot(is_namedlist(x))
  class(x) <- "bmmformula"
  assign_nl_attr(x)
}

#' @export
`[.bmmformula` <- function(object, pars) {
  new_bmmformula(NextMethod())
}

#' @export
`[<-.bmmformula` <- function(object, pars, value) {
  if (!is.list(value)) value <- list(value)
  new_bmmformula(NextMethod())
}

#' @export
"+.bmmformula" <- function(f1, f2) {
  stopif(!is_bmmformula(f1), "The first argument must be a bmmformula.")
  f2_not_a_formula <- !(is_formula(f2) || is_bmmformula(f2)) && !is.null(f2)
  stopif(f2_not_a_formula, "The second argument must be a formula or a bmmformula.")
  stopif(length(f2) != 0 && length(lhs_vars(f2)) == 0, "Formulas must have a left-hand-side variable") 

  if (is_formula(f2) && length(f2) > 0) {
    f2 <- setNames(list(f2), lhs_vars(f2))
  }

  duplicates <- intersect(names(f1), names(f2))
  if (length(duplicates) > 0) {
    message2("Duplicate parameter(s): {collapse_comma(duplicates)}. Overwriting the initial formula.")
  }
  f1[names(f2)] <- f2
  f1
}


#' Generic S3 method for checking if the formula is valid for the specified model
#' @param model a model list object returned from check_model()
#' @param data user supplied data
#' @param formula user supplied formula
#' @return the formula object
#' @keywords internal developer
check_formula <- function(model, data, formula) {
  UseMethod("check_formula")
}

#' @export
check_formula.bmmodel <- function(model, data, formula) {
  stopif(
    !is_bmmformula(formula),
    "The provided formula is not a bmm formula. Please use the bmf() function. E.g.:
    bmmformula(kappa ~ 1, thetat ~ 1) or bmf(kappa ~ 1, thetat ~ 1)"
  )

  wpar <- unrecognized_parameters(model, formula)
  stopif(length(wpar), "Unrecognized model parameters: {collapse_comma(wpar)}")

  formula <- add_missing_parameters(model, formula)
  NextMethod("check_formula")
}

#' @export
check_formula.default <- function(model, data, formula) {
  formula
}

#' @export
check_formula.non_targets <- function(model, data, formula) {
  set_size_var <- model$other_vars$set_size
  pred_list <- rhs_vars(formula, collapse = FALSE)
  has_set_size <- vapply(pred_list, function(x) set_size_var %in% x, logical(1))
  ss_forms <- formula[has_set_size]
  intercepts <- vapply(ss_forms, has_intercept, logical(1))
  stopif(
    any(intercepts),
    "The formula for parameter(s) {names(ss_forms)[intercepts]} contains \\
    an intercept and also uses set_size as a predictor. This model requires \\
    that the intercept is supressed when set_size is used as predictor. \\
    Try using 0 + {set_size_var} instead."
  )
  NextMethod("check_formula")
}

#' @title Convert `bmmformula` objects to `brmsformula` objects
#' @description
#'  Called by [configure_model()] inside [bmm()] to convert the `bmmformula` into a
#'  `brmsformula` based on information in the model object. It will call the
#'  appropriate bmf2bf.\* methods based on the classes defined in the model_\* function.
#' @param model The model object defining one of the supported `bmmodels``
#' @param formula The `bmmformula` that should be converted to a `brmsformula`
#' @return A `brmsformula` defining the response variables and the additional parameter
#'   formulas for the specified `bmmodel`
#' @keywords internal developer
#' @examples
#' model <- mixture2p(resp_error = "error")
#'
#' formula <- bmmformula(
#'   thetat ~ 0 + set_size + (0 + set_size | id),
#'   kappa ~ 1 + (1 | id)
#' )
#'
#' brms_formula <- bmf2bf(model, formula)
#' @export
bmf2bf <- function(model, formula) {
  UseMethod("bmf2bf")
}

# If no model specific methods exist, it will construct a base brms formula from the first resp_vars
# We do it this way because most models require just one response variable
#' @export
bmf2bf.bmmodel <- function(model, formula = bmmformula()) {
  brms_formula <- NextMethod("bmf2bf") %||% brms::bf(glue("{model$resp_vars[[1]]} ~ 1"))
  components <- lapply(formula, function(x) if (is_nl(x)) brms::nlf(x) else brms::lf(x))
  Reduce(`+`, components, init = brms_formula)
}

#' @export
bmf2bf.default <- function(model, formula) {
  NULL
}

add_missing_parameters <- function(model, formula = bmmformula()) {
  formula_pars <- names(formula)
  model_pars <- names(model$parameters)
  fixed_pars <- names(model$fixed_parameters)

  estimable_formula_pars <- setdiff(formula_pars, fixed_pars)
  missing_pars <- setdiff(model_pars, estimable_formula_pars)

  for (mpar in missing_pars) {
    add <- stats::as.formula(paste(mpar, "~ 1"))
    if (mpar %in% fixed_pars) {
      attr(add, "constant") <- TRUE
    } else {
      message2("No formula for parameter {mpar} provided. Only a fixed intercept will be estimated.")
    }
    formula[mpar] <- add
  }

  # Ensure formula list matches model parameters order
  all_pars <- unique(c(model_pars, formula_pars))
  formula[all_pars]
}

unrecognized_parameters <- function(model, formula) {
  predicted_pars <- names(formula)
  possible_pars <- c(names(model$parameters), rhs_vars(formula), get_resp_vars(model))
  setdiff(predicted_pars, possible_pars)
}

get_resp_vars <- function(object, ...) {
  UseMethod("get_resp_vars")
}

#' @export
get_resp_vars.bmmodel <- function(object, ...) {
  vars <- unlist(object$resp_vars, use.names = FALSE)
  stopif(!is.character(vars) && !is.na(vars), "`resp_vars` cannot be coerced to a character vector")
  vars
}

has_intercept <- function(object) {
  UseMethod("has_intercept")
}

#' @export
has_intercept.terms <- function(object) {
  as.logical(attr(object, "intercept"))
}

#' @export
has_intercept.formula <- function(object) {
  try_terms <- try(terms(object), silent = TRUE)
  if (is_try_error(try_terms)) {
    return(FALSE)
  }
  has_intercept(try_terms)
}

rhs_vars <- function(object, ...) {
  UseMethod("rhs_vars")
}

# @param collapse Logical. Should it return a single vector with all the variables
#  or a list with the variables for each parameter?
#' @export
rhs_vars.bmmformula <- function(object, collapse = TRUE, ...) {
  vars <- lapply(object, rhs_vars)
  if (!collapse) {
    return(vars)
  }
  out <- unlist(vars, use.names = F)
  unique(out)
}

#' @export
rhs_vars.formula <- function(object, ...) {
  if (length(object) == 0) {
    warning2("The formula contains no predictors")
    return(character(0))
  }
  rhs <- object[[length(object)]]
  all.vars(rhs)
}

#' @export
rhs_vars.default <- function(object, ...) {
  character(0)
}

lhs_vars <- function(object, ...) {
  UseMethod("lhs_vars")
}

#' @export
lhs_vars.default <- function(object, ...) {
  character(0)
}

#' @export
lhs_vars.bmmformula <- function(object, ...) {
  names(object)
}

#' @export
lhs_vars.brmsformula <- function(object, ...) {
  lhs_vars(brms::brmsterms(object), ...)
}

#' @export
lhs_vars.brmsterms <- function(object, resp = FALSE, ...) {
  names(c(object$dpars, object$nlpars))
}

#' @export
lhs_vars.formula <- function(object, ...) {
  if (length(object) == 3) {
    return(as.character(object[[2]]))
  }
  character(0)
}

# adds an attribute nl to each component of the the formula indicating if the
# any of the predictors of the component are also predicted in another component
assign_nl_attr <- function(object) {
  UseMethod("assign_nl_attr")
}

#' @export
assign_nl_attr.default <- function(object) {
  stop2("Don't know how to assign nl attributes to a {class(formula)[1]} object")
}

#' @export
assign_nl_attr.bmmformula <- function(object) {
  lhs <- lhs_vars(object)
  predictors <- rhs_vars(object, collapse = FALSE)
  for (dpar in lhs) {
    attr(object[[dpar]], "nl") <- any(predictors[[dpar]] %in% lhs)
  }
  object
}

is_nl <- function(object, ...) {
  UseMethod("is_nl")
}

#' @export
is_nl.bmmformula <- function(object, ...) {
  vapply(object, is_nl, logical(1))
}

#' @export
is_nl.default <- function(object, ...) {
  isTRUE(attr(object, "nl"))
}

is_formula <- function(x) inherits(x, "formula")
is_bmmformula <- function(x) inherits(x, "bmmformula")
is_brmsformula <- function(x) inherits(x, "brmsformula")
is_form_or_num <- function(x) is_formula(x) || (is.numeric(x) && length(x) == 1)
is_null_formula <- function(x) is_formula(x) && length(x) == 0

is_constant <- function(x) {
  UseMethod("is_constant")
}

#' @export
is_constant.bmmformula <- function(x) {
  vapply(x, is_constant, logical(1))
}

#' @export
is_constant.default <- function(x) {
  isTRUE(attr(x, "constant")) || is.numeric(x)
}

#' @title Apply link functions for parameters in a formula or bmmformula
#' @description
#'   This function applies the specified link functions in the list of `links` to the
#'   `formula` or `bmmformula` that is passed to it. This function is mostly used internally for configuring
#'   `bmmodels`.
#' @param formula A `formula` or `bmmformula` that the links should be applied to
#' @param links A list of `links` that should be applied to the formula. Each element in this list
#'   should be named using the parameter labels the links should be applied for and contain
#'   a character variable specifying the link to be applied. Currently implemented links are:
#'   "log", "logit", "probit", and "identity".
#' @return The `formula` or `bmmformula` the links have been applied to
#'
#' @examples
#' # specify a bmmformula
#' form <- bmf(x ~ a + c, kappa ~ 1, a ~ 1, c ~ 1)
#' links <- list(a = "log", c = "logit")
#'
#' apply_links(form, links)
#' @keywords developer
#' @export
apply_links <- function(formula, links = nlist()) {
  if (length(links) == 0) {
    return(formula)
  }
  stopifnot(is_namedlist(links) || is.null(links), all(vapply(links, is.character, logical(1))))
  UseMethod("apply_links")
}

#' @export
apply_links.bmmformula <- function(formula, links = nlist()) {
  which_nl <- is_nl(formula)
  formula[which_nl] <- lapply(formula[which_nl], function(f) apply_links(f, links))
  formula
}

#' @export
apply_links.formula <- function(formula, links = nlist()) {
  shared_parameters <- names(links) %in% rhs_vars(formula)
  if (!any(shared_parameters)) {
    return(formula)
  }
  inv_list <- links[shared_parameters]
  inv_list[] <- lapply(names(inv_list), function(nm) inv_link(nm, inv_list[[nm]]))
  eval(methods::substituteDirect(formula, inv_list), envir = environment(formula))
}

#' The inverse of a specified link function applied to a parameter.
#'
#' @param par A character string representing the parameter to which the inverse link function will be applied.
#' @param link A character string specifying the link function. Options are "log", "logit", "probit", and "identity".
#' @return An expression representing the inverse link function applied to the parameter.
#' @examples
#' inv_link("x", "log")
#' identical(inv_link("x", "log"), quote(exp(x))) # TRUE
#' @noRd
inv_link <- function(par = character(0), link = c("log", "logit", "probit", "identity")) {
  stopifnot(is.character(par), length(par) == 1)
  par <- as.symbol(par)
  link <- match.arg(link)
  switch(link,
    log = substitute(exp(par)),
    logit = substitute(inv_logit(par)),
    probit = substitute(Phi(par)),
    identity = par
  )
}
