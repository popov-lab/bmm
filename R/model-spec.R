#' The specification of a measurement model
#'
#' @description `model_spec()` returns everything bmm states about a model in
#'   one list of plain R data: its name and version, its parameters with their
#'   links and fixed values, the default priors, the ranges initial values are
#'   drawn from, the data columns it reads, and the names of its density and
#'   generator functions. Packages that build on bmm read a model here, not
#'   from the fields of the model object.
#'
#' @param x A `bmmodel` object, such as `sdm(resp_error = "y")`, or a `bmmfit`
#'   object returned by [bmm()].
#' @param ... Unused.
#'
#' @return A list of class `bmm_model_spec` with these elements, in this order:
#'   \describe{
#'     \item{`name`}{The model's name.}
#'     \item{`package`}{The package that defines the model: `"bmm"` for the
#'       models bmm ships.}
#'     \item{`version`}{The model's version, such as `"3par"` for [ezdm()], or
#'       `NA` for a model without versions. This is not the version of bmm.}
#'     \item{`domain`, `task`}{The research domains and the tasks the model is
#'       for, as listed by [bmm_models()].}
#'     \item{`citation`}{The references for the model, one per element.}
#'     \item{`parameters`}{A data frame with one row per parameter and the
#'       columns `parameter`, `description`, `link` (the link function, see
#'       [link_transform()]), `fixed` (logical), `value` (the value of a fixed
#'       parameter on the link scale, else `NA`), `alternative_group` (the
#'       parameters of which a formula may use only one, joined by `"|"`, such
#'       as `"criterion|bias"`, else `NA`) and `softmax` (`TRUE` for mixture
#'       weights that are transformed together through a softmax, whose
#'       `link` then does not apply on its own). For [m3()] with
#'       `version = "custom"` the formula defines the parameters, so the
#'       specification of the constructor lists `b` only and the
#'       specification of a fit lists them all.}
#'     \item{`default_priors`}{A data frame with the columns `parameter`,
#'       `main`, `effects` and `sd`: the default prior on a parameter's
#'       intercept, on its other regression coefficients, and on its random
#'       effect standard deviations, as Stan strings.}
#'     \item{`init_ranges`}{A data frame with the columns `parameter`, `lower`
#'       and `upper`: the range initial values are drawn from, on the native
#'       scale (softmax weights on the sampling scale). It has no rows for
#'       models that build their initial values otherwise ([m3()], [mpt()]).}
#'     \item{`resp_vars`, `other_vars`}{Named lists of the arguments of the
#'       model's constructor, as the model was created with them: the
#'       response variables, and the other variables and settings. Values the
#'       constructor stores beyond its arguments are not listed. Expressions
#'       in them, such as the branches of [mpt()] trees, are given as
#'       strings.}
#'     \item{`data_columns`}{A data frame with the columns `argument` and
#'       `role`: what the data column named by each argument of the
#'       constructor must hold.}
#'     \item{`functions`}{A list with the names of the model's exported
#'       `density` and `generator` (random generation) functions, such as
#'       `"dsdm"` and `"rsdm"`; `NA` where a model has none.}
#'     \item{`contract_version`}{The version of this specification, `"0.9"`
#'       until bmm 1.6.0 fixes it as `"1.0"`.}
#'   }
#'   For a `bmmfit` the list ends with one more element, `bmm_version`: the
#'   version of bmm the model was fitted with. Its `name`, `domain` and `task`
#'   are those of the installed bmm (the vocabulary of [bmm_models()]); the
#'   rest is what the fit was built with.
#'
#'   Every element is a character, logical or numeric vector, a data frame, or
#'   a list of these, so the list converts to JSON as it is. `citation`,
#'   `domain`, `task` and the entries of `resp_vars` and `other_vars` are
#'   vectors, which `auto_unbox = TRUE` writes as a single value when they
#'   have one element and as an array otherwise. An `NA` in a data frame row
#'   drops the key unless `na = "null"` is given. Use
#'   `jsonlite::toJSON(spec, auto_unbox = TRUE, na = "null")`.
#'
#' @section Stability: The specification is the stable surface of a bmm model
#'   for other packages (decision record 0005, "The model contract for
#'   companion and extension packages"): "The contract is exactly what
#'   `model_spec()` returns." and "The lifecycle rules apply to these items as
#'   they apply to exported functions. Reading the model object with `$`, its
#'   attributes, and internal functions stay outside the contract." A change
#'   to an element is announced in the NEWS under "Companion packages" before
#'   the release that makes it. Fields of a `bmmodel` object that
#'   `model_spec()` does not return can change in any release.
#'
#' @family model information
#' @export
#'
#' @examples
#' spec <- model_spec(sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"))
#' spec
#' spec$default_priors
#' spec$functions
#'
#' # mixture weights transformed together through a softmax
#' model_spec(mixture3p(resp_error = "y", nt_features = "nt", set_size = 4))$parameters
model_spec <- function(x, ...) {
  UseMethod("model_spec")
}

#' @rdname model_spec
#' @export
model_spec.bmmodel <- function(x, ...) {
  structure(
    list(
      name = x$name,
      package = "bmm",
      version = if (x$version %in% c("NA", "")) NA_character_ else x$version,
      domain = x$domain,
      task = x$task,
      citation = as.character(x$citation),
      parameters = spec_parameters(x),
      default_priors = spec_default_priors(x$default_priors),
      init_ranges = spec_init_ranges(x$init_ranges, names(x$parameters)),
      resp_vars = plain_data(named_list(x$resp_vars)),
      other_vars = plain_data(named_list(spec_other_vars(x))),
      data_columns = spec_data_columns(x),
      functions = list(
        density = x$functions$density %||% NA_character_,
        generator = x$functions$generator %||% NA_character_
      ),
      contract_version = "0.9"
    ),
    class = c("bmm_model_spec", "list")
  )
}

#' @rdname model_spec
#' @export
model_spec.bmmfit <- function(x, ...) {
  model <- restructure(x)$bmm$model
  current <- current_constructor(model)
  # name, domain and task are bmm's vocabulary (#554), not a property of the
  # fit; fits made before bmm 1.5.0 carry a model object without `functions`
  if (!is.null(current)) {
    model[c("name", "domain", "task")] <- current[c("name", "domain", "task")]
  }
  model$functions <- model$functions %||% current$functions
  spec <- model_spec(model)
  spec$bmm_version <- as.character(x$version$bmm)
  spec
}

#' @export
print.bmm_model_spec <- function(x, ...) {
  cat(
    x$name, "\n",
    "Package: ", x$package, "\n",
    "Version: ", if (is.na(x$version)) "none" else x$version, "\n",
    "Parameters:\n",
    sep = ""
  )
  print(x$parameters[c("parameter", "link", "fixed", "value", "alternative_group")], row.names = FALSE)
  cat("Use `$` for the default priors, initial value ranges, data columns and functions.\n")
  invisible(x)
}

spec_parameters <- function(model) {
  pars <- names(model$parameters) %||% character(0)
  fixed <- pars %in% names(model$fixed_parameters)
  value <- rep(NA_real_, length(pars))
  value[fixed] <- vapply(model$fixed_parameters[pars[fixed]], as.numeric, numeric(1))
  groups <- vapply(model$alternative_parameters, paste, "", collapse = "|")
  data.frame(
    parameter = pars,
    description = as.character(unlist(model$parameters)),
    link = vapply(pars, function(par) model$links[[par]] %||% "identity", "", USE.NAMES = FALSE),
    fixed = fixed,
    value = value,
    alternative_group = vapply(pars, function(par) {
      groups[vapply(model$alternative_parameters, function(group) par %in% group, TRUE)][1]
    }, "", USE.NAMES = FALSE),
    softmax = vapply(pars, .is_softmax_param, TRUE, model = model, USE.NAMES = FALSE)
  )
}

spec_default_priors <- function(priors) {
  component <- function(name) {
    vapply(priors, function(prior) prior[[name]] %||% NA_character_, "", USE.NAMES = FALSE)
  }
  data.frame(
    parameter = names(priors) %||% character(0),
    main = component("main"),
    effects = component("effects"),
    sd = component("sd")
  )
}

# `mu` of the response-time models has an init range but is brms plumbing,
# not a declared parameter
spec_init_ranges <- function(ranges, pars) {
  ranges <- ranges[intersect(names(ranges), pars)]
  data.frame(
    parameter = names(ranges) %||% character(0),
    lower = vapply(ranges, `[`, 0, 1, USE.NAMES = FALSE),
    upper = vapply(ranges, `[`, 0, 2, USE.NAMES = FALSE)
  )
}

# NA for a model registered by another package
spec_model_name <- function(model) {
  intersect(rev(class(model)), model_names())[1]
}

# constructors store bookkeeping next to their arguments (mpt() its indicators,
# sdt_cdp() integer codes); only the arguments are contract items. A model of
# another package has no constructor here, so its entries are kept as stored.
spec_other_vars <- function(model) {
  name <- spec_model_name(model)
  if (is.na(name)) {
    model$other_vars
  } else {
    model$other_vars[intersect(names(model$other_vars), names(formals(get_model2(name))))]
  }
}

# a model registered by another package has no row in data_column_roles
spec_data_columns <- function(model) {
  name <- spec_model_name(model)
  roles <- if (is.na(name)) NULL else column_roles(name, model$version)
  roles <- roles[!is.na(roles)]
  data.frame(argument = names(roles) %||% character(0), role = unname(roles) %||% character(0))
}

# an empty named list is a JSON object, as a filled one is; list() would be an
# array
named_list <- function(x) {
  if (length(x) == 0) setNames(list(), character(0)) else x
}

# the arguments of a model can hold classed lists and expressions, such as the
# branches of mpt() trees; jsonlite refuses the one and deparses neither
plain_data <- function(x) {
  if (is.list(x)) return(lapply(unclass(x), plain_data))
  if (is.language(x)) return(deparse1(x))
  x
}
