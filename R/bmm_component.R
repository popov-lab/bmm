############################################################################# !
# COMPONENT SPECIFICATION                                                ####
############################################################################# !

#' @title Define components of a multivariate bmm model
#' @description `bmm_component()` bundles a `bmmformula`, a model (or plain
#'   `brms` distributional family) and a dataset into a single component of a
#'   multivariate model. Components are combined with `+` into a multivariate
#'   specification that can be passed as the `formula` argument to [bmm()],
#'   which estimates all components jointly.
#'
#'   The main purpose of joint estimation is to obtain the correlation matrix
#'   of subject-level parameters across components: random effects that share
#'   an ID within `|ID|` bars (e.g. `(1 | p | id)`) across components are
#'   placed in a single correlation matrix spanning all components, which can
#'   be extracted from the fitted model with [brms::VarCorr()]. This estimates
#'   individual-difference correlations between measurement model parameters
#'   (or between model parameters and external covariates) jointly, avoiding
#'   the attenuation inherent in two-step approaches that correlate point
#'   estimates from separately fitted models.
#'
#' @details Each component carries its own dataset. The datasets are stacked
#'   internally and each component's likelihood is evaluated only on its own
#'   rows via `subset()` addition terms, so the components may have different
#'   numbers of observations. All components that should share random-effect
#'   correlations must use the same grouping variable (e.g. `id`), with
#'   identical values identifying the same subject across datasets.
#'
#'   For components using a `bmm` model, `formula`, `model` and `data` behave
#'   exactly as in a univariate [bmm()] call. For components using a plain
#'   distributional `family` (e.g. `brms::lognormal()`), the first formula
#'   must have the response variable on the left-hand side (as in `brms`), and
#'   any further formulas may predict the family's distributional parameters
#'   (e.g. `sigma`).
#'
#'   Response variable names must be unique across components, as they identify
#'   the components in the joint model (e.g. in priors and post-processing).
#'
#'   Current limitations:
#'   * the `sdm` model is not yet supported in multivariate models
#'   * families that require response addition terms (e.g. `binomial`, which
#'     needs `trials()`) cannot be used as a component family
#'   * random-effect correlations are estimated on the link scale of each
#'     parameter (e.g. log kappa, logit thetat)
#'
#' @param formula A `bmmformula` created by [bmf()]. For family components,
#'   the first formula must have the response variable on the left-hand side.
#' @param model A `bmmodel` object, e.g. [mixture2p()]. Exactly one of `model`
#'   or `family` must be provided.
#' @param family A `brms` distributional family, e.g. [brms::lognormal()] or
#'   `gaussian()`.
#' @param data A data.frame containing the data for this component.
#' @param prior An optional `brmsprior` object with priors for this component,
#'   specified as in a univariate [bmm()] call. bmm automatically assigns the
#'   priors to this component in the joint model.
#' @param e1,e2 Objects created by `bmm_component()`, or multivariate
#'   specifications resulting from adding such objects.
#'
#' @return An object of class `bmm_component`. Adding two components (or a
#'   component and an existing specification) returns an object of class
#'   `mvbmmformula` that can be passed to [bmm()].
#'
#' @seealso [bmm()], [bmf()], [brms::VarCorr()]
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simulate a visual working memory task and a speed task for the same subjects
#' dat_vwm <- data.frame(
#'   id = factor(rep(1:20, each = 50)),
#'   error = rmixture2p(1000, kappa = 8, p_mem = 0.8)
#' )
#' dat_speed <- data.frame(
#'   id = factor(rep(1:20, each = 50)),
#'   rt = rlnorm(1000, meanlog = -0.5, sdlog = 0.3)
#' )
#'
#' # each component has its own formula, model/family and data; the shared
#' # |p| ID places all tagged random effects in one correlation matrix
#' joint <- bmm_component(
#'   bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + (1 | p | id)),
#'   model = mixture2p(resp_error = "error"),
#'   data = dat_vwm
#' ) +
#'   bmm_component(
#'     bmf(rt ~ 1 + (1 | p | id)),
#'     family = brms::lognormal(),
#'     data = dat_speed
#'   )
#'
#' fit <- bmm(joint, cores = 4, backend = "cmdstanr")
#'
#' # the cross-component correlation matrix of subject-level parameters
#' brms::VarCorr(fit)
#' @export
bmm_component <- function(formula, model = NULL, family = NULL, data, prior = NULL) {
  stop_missing_args()
  data_name <- substitute_name(data)
  stopif(
    is.null(model) == is.null(family),
    "Exactly one of 'model' or 'family' must be provided to bmm_component(). \\
    Use 'model' for bmm measurement models and 'family' for plain brms \\
    distributional families."
  )
  stopif(
    !is_bmmformula(formula),
    "The 'formula' argument must be a bmm formula. Please use the bmf() \\
    function. E.g.: bmf(kappa ~ 1, thetat ~ 1)"
  )
  stopif(
    !is.null(prior) && !brms::is.brmsprior(prior),
    "The 'prior' argument must be a brmsprior object created by \\
    brms::prior() and related functions."
  )
  stopif(
    !is.null(model) && !is_bmmodel(model),
    "The 'model' argument must be a bmmodel object created by one of the \\
    model functions (see supported_models()). To use a plain brms \\
    distributional family, use the 'family' argument instead."
  )
  stopif(
    inherits(model, "sdm"),
    "The sdm model is not yet supported in multivariate bmm models, because \\
    its Stan implementation cannot be combined with other response variables."
  )
  if (!is.null(family)) {
    family <- validate_component_family(family, formula)
    model <- brms_family_model(family, formula)
  }
  structure(nlist(formula, model, data, prior, data_name), class = "bmm_component")
}

#' @rdname bmm_component
#' @export
"+.bmm_component" <- function(e1, e2) {
  if (missing(e2)) {
    return(e1)
  }
  new_mvbmmformula(c(component_list(e1), component_list(e2)))
}

#' @rdname bmm_component
#' @export
"+.mvbmmformula" <- `+.bmm_component`

new_mvbmmformula <- function(components) {
  structure(components, class = "mvbmmformula")
}

is_mvbmmformula <- function(x) {
  inherits(x, "mvbmmformula")
}

is_bmm_component <- function(x) {
  inherits(x, "bmm_component")
}

component_list <- function(x) {
  if (is_bmm_component(x)) {
    return(list(x))
  }
  if (is_mvbmmformula(x)) {
    return(unclass(x))
  }
  stop2(
    "Objects of class '{class(x)[1]}' cannot be part of a multivariate bmm \\
    model. Combine components created by bmm_component() with `+`."
  )
}

#' @export
print.bmm_component <- function(x, ...) {
  cat(style("purple1")("  Model: "))
  cat(summarise_model(x$model, newline = TRUE, wsp = 9), "\n")
  cat(style("purple1")("Formula: "))
  print(x$formula, newline = TRUE, wsp = 9)
  cat("\n")
  cat(
    style("purple1")("   Data:"), x$data_name %||% "",
    "(Number of observations:", paste0(nrow(x$data), ")\n")
  )
  if (!is.null(x$prior)) {
    cat(style("purple1")("  Prior:"), "user-specified for this component\n")
  }
  invisible(x)
}

#' @export
print.mvbmmformula <- function(x, ...) {
  cat(style("green")(paste0(
    "Multivariate bmm model with ", length(x), " components:\n\n"
  )))
  for (i in seq_along(x)) {
    cat(style("green")(paste0("Component ", i, ":\n")))
    print(x[[i]])
    cat("\n")
  }
  invisible(x)
}

############################################################################# !
# DISTRIBUTIONAL FAMILY COMPONENTS                                       ####
############################################################################# !

# family components run through the same pipeline as bmm models by wrapping
# the brms family in a minimal internal bmmodel; `parameters` contains only
# the dpars the user predicts, so add_missing_parameters() will not silently
# add free distributional parameters the user did not ask for
brms_family_model <- function(family, formula) {
  resp <- names(formula)[1]
  predicted <- intersect(names(formula), setdiff(family$dpars, "mu"))
  dpar_links <- lapply(predicted, function(dpar) family[[paste0("link_", dpar)]])
  structure(
    list(
      resp_vars = nlist(resp),
      other_vars = nlist(),
      name = glue("{family$family} distributional model"),
      family = family,
      parameters = setNames(
        as.list(glue("Distributional parameter '{predicted}' of the {family$family} family")),
        predicted
      ),
      links = c(list(mu = family$link), setNames(dpar_links, predicted)),
      fixed_parameters = list(),
      default_priors = list()
    ),
    class = c("bmmodel", "brms_family"),
    # a synthetic call so construct_model_call() prints the family instead of
    # the internal wrapper
    call = as.call(c(as.name(family$family), list(link = family$link)))
  )
}

validate_component_family <- function(family, formula) {
  if (is.function(family)) {
    family <- family()
  }
  stopif(
    !inherits(family, "family"),
    "The 'family' argument must be a distributional family, e.g. gaussian() \\
    or brms::lognormal()."
  )
  stopif(
    inherits(family, c("mixfamily", "customfamily")),
    "Mixture and custom families are not supported in bmm_component(). \\
    Models based on such families must be specified via the 'model' argument."
  )
  if (!inherits(family, "brmsfamily")) {
    family <- brms::brmsfamily(family$family, link = family$link)
  }
  stopif(
    family$family %in% c("binomial", "beta_binomial", "multinomial"),
    "The '{family$family}' family requires response addition terms (e.g. \\
    trials()), which cannot be expressed in a bmm_component() yet."
  )
  stopif(
    names(formula)[1] %in% family$dpars,
    "The first formula of a family component must have the response variable \\
    on the left-hand side, e.g. bmf(rt ~ 1, sigma ~ 1) for a lognormal model \\
    of the variable 'rt'."
  )
  family
}

#' @export
check_model.brms_family <- function(model, data = NULL, formula = NULL) {
  model
}

#' @export
check_data.brms_family <- function(model, data, formula) {
  resp_name <- model$resp_vars[[1]]
  stopif(
    not_in(resp_name, colnames(data)),
    "The response variable '{resp_name}' is not present in the data."
  )
  NextMethod("check_data")
}

#' @export
configure_model.brms_family <- function(model, data, formula) {
  resp <- model$resp_vars[[1]]
  pforms <- formula[names(formula) != resp]
  formula <- Reduce(`+`, lapply(pforms, brms::lf), init = brms::bf(formula[[resp]]))
  formula$family <- model$family
  nlist(formula, data)
}
