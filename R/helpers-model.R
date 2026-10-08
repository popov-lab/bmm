#' @title Generic S3 method for configuring the model to be fit by brms
#' @description Called by bmm() to automatically construct the model
#'   formula, family objects and default priors for the model specified by the
#'   user. It will call the appropriate configure_model.* functions based on the
#'   list of classes defined in the .model_* functions. Currently, we have a
#'   method only for the last class listed in the .model_* functions. This is to
#'   keep model configuration as simple as possible. In the future we may add
#'   shared methods for classes of models that share the same configuration.
#' @param model A model list object returned from check_model()
#' @param data The user supplied data.frame containing the data to be checked
#' @param formula The user supplied formula
#' @return A named list containing at minimum the following elements:
#'
#'  - formula: An object of class `brmsformula`. The constructed model formula
#'  - data: the user supplied data.frame, preprocessed by check_data
#'  - family: the brms family object
#'  - prior: the brms prior object
#'  - stanvars: (optional) An object of class `stanvars` (for custom families).
#'   See [brms::custom_family()] for more details.
#'
#' @details A bare bones configure_model.* method should look like this:
#'
#'  ``` r
#'  configure_model.newmodel <- function(model, data, formula) {
#'
#'     # preprocessing - e.g. extract arguments from data check, construct new variables
#'     <preprocessing code>
#'
#'     # construct the formula
#'     formula <- bmf2bf(formula, model)
#'
#'     # construct the family
#'     family <- <code for new family>
#'
#'     # construct the default prior
#'     prior <- <code for new prior>
#'
#'     # return the list
#'     nlist(formula, data, family, prior)
#'  }
#'  ```
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' configure_model.mixture3p <- function(model, data, formula) {
#'   # retrieve arguments from the data check
#'   max_set_size <- attr(data, "max_set_size")
#'   lure_idx <- attr(data, "lure_idx_vars")
#'   nt_features <- model$other_vars$nt_features
#'   set_size_var <- model$other_vars$set_size
#'
#'   # construct initial brms formula
#'   formula <- bmf2bf(model, formula) +
#'     brms::lf(kappa2 ~ 1) +
#'     brms::lf(mu2 ~ 1) +
#'     brms::nlf(theta1 ~ thetat) +
#'     brms::nlf(kappa1 ~ kappa)
#'
#'   # additional internal terms for the mixture model formula
#'   kappa_nts <- paste0("kappa", 3:(max_set_size + 1))
#'   theta_nts <- paste0("theta", 3:(max_set_size + 1))
#'   mu_nts <- paste0("mu", 3:(max_set_size + 1))
#'
#'   for (i in 1:(max_set_size - 1)) {
#'     formula <- formula +
#'       glue_nlf("{kappa_nts[i]} ~ kappa") +
#'       glue_nlf(
#'         "{theta_nts[i]} ~ {lure_idx[i]} * (thetant + log(inv_ss)) + ",
#'         "(1 - {lure_idx[i]}) * (-100)"
#'       ) +
#'       glue_nlf("{mu_nts[i]} ~ {nt_features[i]}")
#'   }
#'
#'   # define mixture family
#'   vm_list <- lapply(1:(max_set_size + 1), function(x) brms::von_mises(link = "identity"))
#'   vm_list$order <- "none"
#'   formula$family <- brms::do_call(brms::mixture, vm_list)
#'
#'   nlist(formula, data)
#' }
#'
#' @export
#' @keywords internal developer
configure_model <- function(model, data, formula) {
  UseMethod("configure_model")
}

############################################################################# !
# CHECK_MODEL methods                                                    ####
############################################################################# !

#' Generic S3 method for checking if the model is supported and model preprocessing
#'
#' In addition for validating the model, specific methods might add information
#' to the model object based on the provided data and formula
#'
#' @param model the model argument supplied by the user
#' @param data the data argument supplied by the user
#' @param formula the formula argument supplied by the user
#'
#' @return An object of type 'bmmodel'
#' @keywords internal developer
check_model <- function(model, data = NULL, formula = NULL) {
  UseMethod("check_model")
}

#' @export
check_model.default <- function(model, data = NULL, formula = NULL) {
  if (is.function(model)) {
    fun_name <- as.character(substitute(model))
    stopif(
      fun_name %in% model_names(),
      "Did you forget to provide the required arguments to the model function?
      See ?{fun_name} for details on properly specifying the model argument"
    )
  }

  stopif(
    !is_supported_bmmodel(model),
    "You provided an object of class `{class(model)}` to the model argument.
    The model argument should be a `bmmodel` function.
    You can see the list of supported models by running `bmm_models()`

    {models_text()}"
  )
  model
}

#' @export
check_model.bmmodel <- function(model, data = NULL, formula = NULL) {
  model <- check_links(model)
  model <- replace_regex_variables(model, data)
  model <- update_model_fixed_parameters(model, formula)
  NextMethod("check_model")
}

# check if the user has provided a regular expression for any model variables and
# replace the regular expression with the actual variables
replace_regex_variables <- function(model, data) {
  regex <- isTRUE(attr(model, "regex"))
  regex_vars <- attr(model, "regex_vars")

  # check if the regex transformation has already been applied (e.g., if
  # updating a previously fit model)
  regex_applied <- isTRUE(attr(model, "regex_applied"))
  if (regex_applied || !regex || length(regex_vars) == 0) {
    return(model)
  }

  data_cols <- names(data)
  # save original user-provided variables
  user_vars <- c(model$resp_vars, model$other_vars)
  attr(model, "user_vars") <- user_vars

  for (var in regex_vars) {
    var_type <- if (var %in% names(model$other_vars)) "other_vars" else "resp_vars"
    model[[var_type]][[var]] <- get_variables(model[[var_type]][[var]], data_cols, regex)
  }

  attr(model, "regex_applied") <- regex
  model
}

# if the user has provided a constant in the bmmformula, add that info to the
# model object; if they have predicted a parameter that is constant by default,
# remove it from the model object
update_model_fixed_parameters <- function(model, formula) {
  constants <- names(formula)[is_constant(formula)]
  free <- names(formula)[!is_constant(formula)]
  # add new constants to the model object
  if (length(constants) > 0) {
    model$fixed_parameters[constants] <- strip_attributes(formula[constants],
      protect = "names",
      recursive = TRUE
    )
  }
  overwrite <- intersect(names(model$fixed_parameters), free)
  if (length(overwrite) > 0) {
    model$fixed_parameters[overwrite] <- NULL
  }
  model
}

############################################################################# !
# LINKS                                                                  ####
############################################################################# !

# The links bmm can apply, with the range each one confines its parameter to on
# the native scale. A link is the model's statement about where the parameter
# lives, which is why widening that range is worth a warning: the sampler is
# then free to propose values the likelihood is not defined for, and the
# default priors, which are written on the link scale, no longer imply the same
# range. Keep in sync with link_transform().
.link_ranges <- list(
  identity = c(-Inf, Inf),
  log = c(0, Inf),
  softplus = c(0, Inf),
  # sqrt admits exactly 0, which log does not (link_transform(0, "sqrt",
  # inverse = TRUE) is 0). The bound below is the closure of that range, so a
  # log -> sqrt swap does not warn; a single point is not a sampling hazard.
  sqrt = c(0, Inf),
  # 1/eta on an unbounded linear predictor is negative for every eta < 0
  inverse = c(-Inf, Inf),
  log1p = c(-1, Inf),
  logm1 = c(1, Inf),
  logit = c(0, 1),
  probit = c(0, 1),
  cloglog = c(0, 1),
  loglog = c(0, 1),
  softmax = c(0, 1),
  tan_half = c(-pi, pi)
)

# The parameters whose link a model passes on to the fit. Models that build
# their family from fixed links override this with character(0), so that
# setting one is refused rather than silently applied to print(), the initial
# values and the prior scale while the sampler keeps the default. NULL means
# the model has no closed set of parameters (m3), so no name can be refused.
# The methods are named for the model rather than for `bmmodel`, which comes
# first in every class vector and would shadow them.
settable_links <- function(model) {
  UseMethod("settable_links")
}

#' @exportS3Method
settable_links.default <- function(model) {
  names(model$links)
}

# The link functions a model can carry through to its likelihood. .link_ranges
# is the wider vocabulary of links bmm knows a range for, including the ones a
# model may declare as its own default but a user cannot ask for: brms writes
# no inverse-link code for loglog or softmax, so a custom family built with
# either dies in stancode() with "argument is of length zero". Models that
# apply their links themselves rather than through a brms family override this.
settable_link_functions <- function(model) {
  UseMethod("settable_link_functions")
}

#' @exportS3Method
settable_link_functions.default <- function(model) {
  setdiff(names(.link_ranges), c("loglog", "softmax"))
}

# Assign user-supplied links onto the model's defaults. A name-indexed
# assignment appends an unrecognized name instead of refusing it, so a typo
# used to advertise a parameter the model does not have while the link the user
# meant was never applied (#420). Every constructor assigns through here.
set_links <- function(out, links) {
  attr(out, "links_default") <- out$links
  attr(out, "links_checked") <- out$links
  if (length(links) == 0) {
    return(out)
  }
  links <- validate_links(links, out)
  out$links[names(links)] <- links
  attr(out, "links_checked") <- out$links
  out
}

# Links can also reach a model after construction -- `model$links <- list(...)`
# is the documented idiom for m3 -- so the pipeline re-checks whatever differs
# from the state set_links() last signed off on. Diffing against `links_default`
# instead would re-validate the user's constructor argument and warn a second
# time. Parameters in `links_fixed` are excluded because the pipeline sets those
# itself (see resolve_fixed_links). A model that never went through set_links()
# -- a custom m3, a fit from an older bmm -- carries no attribute and is left
# alone; for a custom m3 the missing-links check in check_model.m3_custom is
# what catches a garbage list.
check_links <- function(model) {
  default <- attr(model, "links_default")
  if (is.null(default)) {
    return(model)
  }
  checked <- attr(model, "links_checked") %||% default
  pars <- setdiff(names(model$links), names(model$links_fixed))
  changed <- pars[!vapply(pars, function(p) {
    identical(model$links[[p]], checked[[p]])
  }, logical(1))]
  if (length(changed) == 0) {
    return(model)
  }
  # validated against the links the model declared, so that a target the user
  # just added is not offered back as settable
  declared <- model
  declared$links <- default
  links <- validate_links(model$links[changed], declared)
  model$links[setdiff(changed, names(links))] <- NULL
  model$links[names(links)] <- links
  model
}

# One set of rules for every path by which a link reaches a model: the name
# identifies a parameter (a single typo is repaired, with a warning), the model
# can pass the link on to the fit, the link is one bmm implements, and a link
# whose range is wider than the declared one is a sampling hazard rather than
# an error.
validate_links <- function(links, model) {
  stopif(
    !is_namedlist(links) ||
      !all(vapply(links, function(l) is.character(l) && length(l) == 1, logical(1))),
    'The `links` argument must be a named list of link functions, \\
     e.g. links = list(kappa = "log")'
  )
  defaults <- attr(model, "links_default") %||% model$links
  settable <- settable_links(model)
  # a model built by use_model_template() is not in bmm_models() yet
  model_name <- c(
    intersect(class(model), model_names()),
    class(model)[2]
  )[1]
  given <- names(links)

  if (!is.null(settable)) {
    known <- unique(c(names(model$parameters), names(defaults)))
    names(links) <- vapply(given, match_link_target, character(1), known = known)
    unmatched <- given[is.na(names(links))]
    stopif(
      length(unmatched) > 0,
      "Unrecognized link target(s): {collapse_comma(unmatched)}. \\
       {model_name}() takes links for {collapse_comma(known)}"
    )
    duplicates <- given[names(links) %in% names(links)[duplicated(names(links))]]
    stopif(
      anyDuplicated(names(links)) > 0,
      "Several entries of `links` name the same parameter: \\
       {collapse_comma(duplicates)}"
    )
  }

  # a repair is reported before anything is refused, so that a user whose typo
  # resolved to a parameter they cannot set learns both halves of what happened
  repaired <- names(links) != given
  warnif(
    any(repaired),
    "Link target(s) {collapse_comma(given[repaired])} read as \\
     {collapse_comma(names(links)[repaired])}. Check the spelling of your \\
     `links` argument"
  )

  # naming the link the model already uses asks for no change, so it is a no-op
  # and neither the refusals below nor the allow-list applies to it
  asked <- names(links)[!vapply(names(links), function(p) {
    identical(links[[p]], defaults[[p]])
  }, logical(1))]

  if (!is.null(settable)) {
    settable_str <- if (length(settable) > 0) {
      glue("Links can be set for {collapse_comma(settable)}")
    } else {
      glue("No link of {model_name}() can be set")
    }
    scaling <- setdiff(asked, names(defaults))
    stopif(
      length(scaling) > 0,
      "{collapse_comma(scaling)} has no link in {model_name}(): the parameter \\
       is fixed for scaling. {settable_str}"
    )
    fixed <- setdiff(setdiff(asked, settable), scaling)
    stopif(
      length(fixed) > 0,
      "The link of {collapse_comma(fixed)} cannot be changed in {model_name}(): \\
       the model is written around {summarise_links(defaults[fixed])} -- its \\
       likelihood, the values it fixes, or both are expressed on that scale -- \\
       so another link would not give the parameter the meaning \\
       {model_name}() documents for it. {settable_str}"
    )
  }

  offered <- settable_link_functions(model)
  unsupported <- setdiff(unlist(links[asked]), offered)
  stopif(
    length(unsupported) > 0,
    "Unknown link function(s): {collapse_comma(unsupported)}. \\
     {model_name}() takes {collapse_comma(offered)}"
  )

  warn_link_range(links, defaults)
  links
}

# "kapa" is not a prefix of "kappa", so R's partial matching does not see it,
# but a single edit does. A repair must be unique: imm's a, c and s are each
# one edit apart, so a typo among them identifies no parameter.
match_link_target <- function(name, known) {
  if (name %in% known) {
    return(name)
  }
  hit <- known[startsWith(known, name)]
  if (length(hit) != 1) {
    distance <- utils::adist(name, known, ignore.case = TRUE)[1, ]
    hit <- known[distance == min(distance) & distance <= 1]
  }
  if (length(hit) == 1) hit else NA_character_
}

# A link that admits values the default one excludes (log -> identity for a
# positive parameter) is legal -- it is how a parameter is freed from a bound
# the model assumes by default -- but the likelihood is written for the default
# range, so it is flagged.
warn_link_range <- function(links, defaults) {
  pars <- intersect(names(links), names(defaults))
  wider <- vapply(pars, function(p) {
    given <- .link_ranges[[links[[p]]]]
    default <- .link_ranges[[defaults[[p]]]]
    if (is.null(given) || is.null(default)) {
      return(FALSE)
    }
    given[1] < default[1] || given[2] > default[2]
  }, logical(1))
  warnif(
    any(wider),
    "The link(s) {summarise_links(links[pars[wider]])} allow values that the \\
     model's default {summarise_links(defaults[pars[wider]])} exclude. \\
     Sampling a bounded parameter on a wider scale can push the likelihood out \\
     of its domain, so check that your priors keep \\
     {collapse_comma(pars[wider])} in range"
  )
}

# kept central rather than as a field in each model constructor so the console
# annotations stay short, uniform and reviewable in one place
response_annotations <- function(model) {
  if (inherits(model, "circular")) {
    return(list(resp_error = "radians in [-pi, pi]"))
  }
  if (inherits(model, "ddm") || inherits(model, "cswald")) {
    return(list(
      rt = "seconds",
      response = "0/1 or logical; 1 = upper boundary"
    ))
  }
  if (inherits(model, "ezdm")) {
    return(list(
      mean_rt = "seconds",
      var_rt = "seconds^2",
      n_upper = "count of upper-boundary responses"
    ))
  }
  if (inherits(model, "m3") || inherits(model, "mpt")) {
    return(list(resp_cats = "counts per response category"))
  }
  if (inherits(model, "sdt_yn")) {
    return(list(response = "count of 'old'/'signal' responses per cell"))
  }
  list()
}

#' @export
print.bmmodel <- function(x, ...) {
  cat(construct_model_call(x), "\n")
  annotations <- response_annotations(x)
  resp_str <- sapply(names(x$resp_vars), function(var) {
    annot <- annotations[[var]]
    paste0(
      var, " = ", paste(x$resp_vars[[var]], collapse = ", "),
      if (!is.null(annot)) paste0(" (", annot, ")")
    )
  })
  cat("Response:  ", paste(resp_str, collapse = "\n            "), "\n")
  cat("Parameters:", paste(names(x$parameters), collapse = ", "), "\n")
  if (length(x$links) > 0) {
    cat("Links:     ", summarise_links(x$links), "\n")
  }
  if (length(x$fixed_parameters) > 0) {
    fixed_str <- paste(
      names(x$fixed_parameters), "=", x$fixed_parameters,
      collapse = ", "
    )
    cat("Fixed:     ", fixed_str, "\n")
  }
  print_model_details(x)
  cat("Use parameter_info() for more details.\n")
  invisible(x)
}

# print.bmmodel dispatches before any model-specific print method (classes
# are ordered general to specific), so model-specific lines are added here
#' @keywords internal
print_model_details <- function(model, ...) {
  UseMethod("print_model_details")
}

#' @export
print_model_details.default <- function(model, ...) {
  invisible(NULL)
}


############################################################################# !
# HELPER FUNCTIONS                                                       ####
############################################################################# !

# maps the `domain` field of each `.model_*()` constructor to the task group
# shown by bmm_models(); a domain not listed here prints as its own group,
# so a new model never disappears from the list
model_groups <- c(
  "Visual working memory" = "Continuous reproduction",
  "Working Memory (categorical), Categorical Decision Making" = "Categorical recall and n-AFC decisions",
  "Categorical decision making, memory, and reasoning" = "Processing-tree models",
  "Perception & Recognition Memory" = "Detection, recognition and confidence judgments",
  "Recognition Memory" = "Detection, recognition and confidence judgments",
  "Decision Making / Response times" = "Choices and response times"
)

model_group <- function(domain) {
  group <- unname(model_groups[domain])
  group[is.na(group)] <- domain[is.na(group)]
  group[group == ""] <- "Other models"
  group
}

model_registry <- function(models = model_names()) {
  specs <- lapply(models, function(m) get_model(m)())
  registry <- data.frame(
    model = models,
    name = sub("\\.$", "", vapply(specs, `[[`, "", "name")),
    domain = vapply(specs, `[[`, "", "domain")
  )
  registry$group <- model_group(registry$domain)
  known <- unique(unname(model_groups))
  group_levels <- c(known, setdiff(registry$group, known))
  registry[order(match(registry$group, group_levels), registry$model), ]
}

reference_url <- "https://popov-lab.github.io/bmm/reference/"

format_model_list <- function(registry, style = "text", headers = TRUE) {
  blocks <- lapply(unique(registry$group), function(group) {
    rows <- registry[registry$group == group, ]
    header <- if (!headers) NULL else if (style == "md") glue("**{group}**") else group
    items <- if (style == "md") {
      glue("- [`{rows$model}()`]({reference_url}{rows$model}.html): {rows$name}")
    } else {
      glue("- {rows$model}(): {rows$name}")
    }
    c(header, if (headers) "", items, "")
  })
  unlist(blocks)
}

# what each data argument of a constructor holds, for the model overview in the
# Get started article. Kept central for the same reason as
# response_annotations(); the `@param` text is too long for a table cell. A
# `<model>_<version>` entry replaces the model's entry for that version, and NA
# marks an argument that names no data column. A test requires an entry for
# every argument without a default, so a new model needs its labels here
data_column_roles <- list(
  cswald = c(
    rt = "response time in seconds",
    response = "choice, 0 = lower and 1 = upper boundary"
  ),
  ddm = c(
    rt = "response time in seconds",
    response = "choice, 0 = lower and 1 = upper boundary"
  ),
  ezdm = c(
    mean_rt = "mean response time in seconds",
    var_rt = "variance of the response times in seconds\u00b2",
    n_upper = "number of upper-boundary responses",
    n_trials = "number of trials"
  ),
  ezdm_4par = c(
    mean_rt = "mean response time in seconds, one column per boundary (upper, lower)",
    var_rt = "variance of the response times in seconds\u00b2, one column per boundary (upper, lower)",
    n_upper = "number of upper-boundary responses",
    n_trials = "number of trials"
  ),
  imm = c(
    resp_error = "response error relative to the target, in radians",
    nt_features = "non-target features relative to the target, in radians, one column per non-target",
    nt_distances = "distance of each non-target to the target, one column per non-target",
    set_size = "set size (a column, or one number)"
  ),
  imm_abc = c(
    resp_error = "response error relative to the target, in radians",
    nt_features = "non-target features relative to the target, in radians, one column per non-target",
    nt_distances = NA,
    set_size = "set size (a column, or one number)"
  ),
  m3 = c(
    resp_cats = "number of responses in each response category, one column per category",
    num_options = "number of candidates in each category (columns, or one number per category)"
  ),
  m3_ss = c(
    resp_cats = "number of responses in each of 3 categories, in the order correct, other list item, not-presented lure (one column each)",
    num_options = "number of candidates in each category (columns, or one number per category)"
  ),
  m3_cs = c(
    resp_cats = "number of responses in each of 5 categories, in the order correct, distractor close in context, other list item, other distractor, not-presented lure (one column each)",
    num_options = "number of candidates in each category (columns, or one number per category)"
  ),
  mixture2p = c(
    resp_error = "response error relative to the target, in radians"
  ),
  mixture3p = c(
    resp_error = "response error relative to the target, in radians",
    nt_features = "non-target features relative to the target, in radians, one column per non-target",
    set_size = "set size (a column, or one number)"
  ),
  mpt = c(
    trees = "number of responses in each response category, one column per category, named after the branches of the trees",
    tree_id = "the tree each row belongs to (models with several trees)"
  ),
  sdm = c(
    resp_error = "response error relative to the target, in radians"
  ),
  sdt_cdp = c(
    response = "prefix of the count columns `new<k>`, `know<k>`, `remember<k>` and optionally `guess<k>`, one per confidence level (default: no prefix)",
    stimulus = "stimulus type, 0 = new and 1 = old",
    n_new = NA,
    n_old = NA
  ),
  sdt_mafc = c(
    response = "number of correct responses",
    n_trials = "number of trials",
    m = "number of alternatives (a column, or one number)"
  ),
  sdt_ranking = c(
    response = "number of trials with the target at each rank, one column per rank, from rank 1 (most likely target) to rank `m`",
    m = "number of ranked items (a column, or one number)"
  ),
  sdt_rating = c(
    response = "number of responses in each rating category, one column per category, ordered from 'definitely noise' to 'definitely signal'",
    stimulus = "stimulus type, 0 = noise/new and 1 = signal/old"
  ),
  sdt_yn = c(
    response = "number of 'old'/'signal' responses",
    stimulus = "stimulus type, 0 = noise/new and 1 = signal/old",
    n_trials = "number of trials"
  )
)

model_versions <- function(model) {
  version <- formals(get_model2(model))$version
  if (is.null(version)) NA_character_ else eval(version)
}

column_roles <- function(model, version) {
  data_column_roles[[paste0(model, "_", version)]] %||% data_column_roles[[model]]
}

format_data_columns <- function(roles) {
  roles <- roles[!is.na(roles)]
  paste0("`", names(roles), "`: ", roles, collapse = "<br>")
}

# the part of a parameter description before its first ": ", " = " or ". ",
# e.g. "Drift rate" from "Drift rate = Average rate of evidence accumulation"
parameter_label <- function(description) {
  sub("(: | = |\\. ).*$", "", description)
}

format_key_parameters <- function(spec) {
  estimated <- setdiff(names(spec$parameters), names(spec$fixed_parameters))
  fixed <- intersect(names(spec$fixed_parameters), names(spec$parameters))
  lines <- if (length(estimated) == 0) {
    "None by default: your formula defines them"
  } else {
    descriptions <- vapply(spec$parameters[estimated], as.character, "")
    paste0("`", estimated, "`: ", parameter_label(descriptions))
  }
  if (length(fixed) > 0) {
    lines <- c(lines, paste0("Fixed by default: ", paste0("`", fixed, "`", collapse = ", ")))
  }
  paste(lines, collapse = "<br>")
}

# one row per model version, with the data columns it needs and the parameters
# it estimates; used for the tables of the Get started article
model_overview <- function(group = NULL) {
  registry <- model_registry()
  if (!is.null(group)) {
    registry <- registry[registry$group %in% group, ]
  }
  rows <- lapply(seq_len(nrow(registry)), function(i) {
    model <- registry$model[i]
    versions <- model_versions(model)
    specs <- lapply(versions, function(version) {
      if (is.na(version)) get_model(model)() else get_model(model)(version = version)
    })
    label <- glue("[`{model}()`]({reference_url}{model}.html)")
    if (length(versions) > 1) {
      label <- glue("{label}, version `{versions}`")
    }
    data.frame(
      Model = paste0(label, "<br>", registry$name[i]),
      `Data columns` = vapply(versions, function(version) {
        format_data_columns(column_roles(model, version))
      }, "", USE.NAMES = FALSE),
      `Key parameters` = vapply(specs, format_key_parameters, ""),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Measurement models available in `bmm`
#'
#' @details Printed, the result lists the models grouped by the task they are
#'   meant for, one line per model with its constructor and full name. The
#'   groups are: continuous reproduction; categorical recall and n-AFC
#'   decisions; detection, recognition and confidence judgments; choices and
#'   response times; processing-tree models. Type `?modelname` (for example
#'   `?imm`) for the arguments of a model.
#' @return A character vector of model names with class `bmm_models`, which
#'   prints as the grouped list. Use it as a character vector in base R or
#'   dplyr, e.g. `"imm" %in% bmm_models()`; `as.character()` gives a plain one
#'   where a function refuses the class.
#' @export
#'
#' @examples
#' bmm_models()
#' "imm" %in% bmm_models()
bmm_models <- function() {
  structure(model_names(), class = "bmm_models")
}

#' @export
print.bmm_models <- function(x, ...) {
  # base functions such as sub() keep the class on vectors that are no longer
  # model names, so those print as what they are
  if (!all(x %in% model_names())) {
    print(unclass(x), ...)
    return(invisible(x))
  }
  cat(models_text(unclass(x)))
  invisible(x)
}

# vctrs refuses to combine an unknown class with character, so without these
# dplyr joins, binds, if_else() and assignment fail on a bmm_models() column
#' @exportS3Method vctrs::vec_ptype2 bmm_models.bmm_models
vec_ptype2.bmm_models.bmm_models <- function(x, y, ...) character()

#' @exportS3Method vctrs::vec_ptype2 bmm_models.character
vec_ptype2.bmm_models.character <- function(x, y, ...) character()

#' @exportS3Method vctrs::vec_ptype2 character.bmm_models
vec_ptype2.character.bmm_models <- function(x, y, ...) character()

#' @exportS3Method vctrs::vec_cast character.bmm_models
vec_cast.character.bmm_models <- function(x, to, ...) unclass(x)

#' @exportS3Method vctrs::vec_cast bmm_models.character
vec_cast.bmm_models.character <- function(x, to, ...) structure(x, class = "bmm_models")

#' @exportS3Method vctrs::vec_ptype2 bmm_models.factor
vec_ptype2.bmm_models.factor <- function(x, y, ...) character()

#' @exportS3Method vctrs::vec_ptype2 factor.bmm_models
vec_ptype2.factor.bmm_models <- function(x, y, ...) character()

#' @exportS3Method vctrs::vec_ptype2 bmm_models.ordered
vec_ptype2.bmm_models.ordered <- function(x, y, ...) character()

#' @exportS3Method vctrs::vec_ptype2 ordered.bmm_models
vec_ptype2.ordered.bmm_models <- function(x, y, ...) character()

# one string, not one per model, so it is not a format() method: tibble and
# print.data.frame call format() on columns and expect one string per element
models_text <- function(models = model_names()) {
  out <- paste(
    c(
      "The following models are supported:", "",
      format_model_list(model_registry(models), "text"),
      "Type `?modelname` to get information about a specific model, e.g. `?imm`", ""
    ),
    collapse = "\n"
  )
  gsub("`", " ", out)
}

# the registry behind bmm_models(), as plain names for internal lookups
model_names <- function() {
  sub("^\\.model_", "", lsp("bmm", pattern = "^\\.model_"))
}

#' Deprecated: use `bmm_models()`
#'
#' @description `supported_models()` is deprecated as of bmm 1.4.0 and will be
#'   removed in bmm 1.6.0. It shares its name with `insight::supported_models()`,
#'   which the **parameters** package re-exports, so whichever package is
#'   attached last decides what `supported_models()` returns. Replace both
#'   `supported_models()` and `supported_models(print_call = FALSE)` with
#'   [bmm_models()]; wrap it in `as.character()` if you need a plain
#'   character vector.
#'
#' @param print_call Logical. If `TRUE` (default), returns the output of
#'   [bmm_models()], which prints the models grouped by task. If `FALSE`,
#'   returns the model names as a plain character vector.
#' @return The output of [bmm_models()], or the model names as a plain
#'   character vector if `print_call = FALSE`. Before bmm 1.4.0 the default
#'   returned the printed list as one string; use
#'   `capture.output(print(bmm_models()))` for the printed list as text,
#'   one line per element.
#' @keywords internal
#' @export
supported_models <- function(print_call = TRUE) {
  warning2("`supported_models()` is deprecated as of bmm 1.4.0; use `bmm_models()`. \\
            It will be removed in 1.6.0.")
  if (print_call) bmm_models() else model_names()
}


#' @title Generate a markdown list of the measurement models available in `bmm`
#' @description Used internally to populate the README and the "Get started"
#'   article. Models are grouped as in [bmm_models()], and every model
#'   links to its reference page on the website.
#' @param group Optional character vector of group labels as printed by
#'   [bmm_models()]. Only those groups are listed and the group headers
#'   are omitted, so a document can add its own text per group.
#' @return Markdown code for printing the list of measurement models available
#'   in `bmm`
#' @export
#'
#' @examples
#' print_pretty_models_md()
#' print_pretty_models_md(group = "Continuous reproduction")
#'
#' @keywords internal
print_pretty_models_md <- function(group = NULL) {
  registry <- model_registry()
  stopif(
    !is.null(group) && length(group) == 0,
    "`group` must not be empty; omit it to print all groups."
  )
  stopif(
    !all(group %in% registry$group),
    "Unknown model group(s): {collapse_comma(setdiff(group, registry$group))}"
  )
  if (!is.null(group)) {
    registry <- registry[registry$group %in% group, ]
  }
  cat(format_model_list(registry, "md", headers = is.null(group)), sep = "\n")
}

# used to extract well formatted information from the model object to print
# in the @details section for the documentation of each model
model_docs <- function(model, components = "all") {
  UseMethod("model_docs")
}


#' @export
model_docs.bmmodel <- function(model, components = "all") {
  pars <- model$parameters
  par_info <- ""
  if (length(pars) > 0) {
    for (par in names(pars)) {
      par_info <- paste0(par_info, "   - `", par, "`: ", pars[[par]], "\n")
    }
  }

  fixed_pars <- model$fixed_parameters
  fixed_par_info <- ""
  if (length(fixed_pars) > 0) {
    for (fixed_par in names(fixed_pars)) {
      fixed_par_info <- paste0(
        fixed_par_info, "   - `", fixed_par,
        "` = ", fixed_pars[[fixed_par]], "\n"
      )
    }
  }

  links <- model$links
  links_info <- summarise_links(links)

  priors <- model$default_priors
  priors_info <- summarise_default_prior(priors)

  info_all <- list(
    domain = paste0("* **Domain:** ", model$domain, "\n\n"),
    task = paste0("* **Task:** ", model$task, "\n\n"),
    name = paste0("* **Name:** ", model$name, "\n\n"),
    citation = paste0("* **Citation:** \n\n", collapse(paste0("   - ", model$citation, "\n")), "\n"),
    version = paste0("* **Version:** ", model$version, "\n\n"),
    requirements = paste0("* **Requirements:** \n\n  ", model$requirements, "\n\n"),
    parameters = paste0("* **Parameters:** \n\n  ", par_info, "\n"),
    fixed_parameters = paste0("* **Fixed parameters:** \n\n  ", fixed_par_info, "\n"),
    links = paste0("* **Default parameter links:** \n\n     - ", links_info, "\n\n"),
    prior = paste0("* **Default priors:** \n\n", priors_info, "\n")
  )

  if (length(components) == 1 && components == "all") {
    components <- names(info_all)
  }

  if (model$version == "NA" || model$version == "") {
    components <- components[components != "version"]
  }

  # return only the specified components
  collapse(info_all[components])
}



#' @param model A string with the name of the model supplied by the user
#' @return A function of type .model_*
#' @details the returned object is a function. To get the model object, call the
#'   returned function, e.g. `get_model("mixture2p")()`
#' @noRd
get_model <- function(model) {
  get(paste0(".model_", model), mode = "function")
}

# same as get_model2, but with the new model structure for the user facing alias
get_model2 <- function(model) {
  get(model, mode = "function")
}

#' Create a file with a template for adding a new model (for developers)
#'
#' @param model_name A string with the name of the model. The file will be named
#'  `model_model_name.R` and all necessary functions will be created with
#'  the appropriate names and structure. The file will be saved in the `R/`
#'  directory
#' @param versions An optional character vector naming the model versions. If
#'  `NULL` (default), the template generates a single flat `.{model_name}_defaults`
#'  specification (like `ddm`). If supplied, it generates a
#'  `.{model_name}_version_table` with one entry per version and a versioned
#'  user-facing alias validated with `match.arg()` (like `cswald`).
#' @param testing Logical; If TRUE, the function will return the file content but
#'  will not save the file. If FALSE (default), the function will save the file
#' @param custom_family Logical; Do you plan to define a brms::custom_family()?
#'  If TRUE the function will add a section for the custom family, placeholders
#'  for the stan_vars and corresponding empty .stan files in
#'  `inst/stan_chunks/`, that you can fill For an example, see the sdm
#'  model in `/R/model_sdm.R`. If FALSE (default) the function will
#'  not add the custom family section nor stan files.
#' @param stanvar_blocks A character vector with the names of the blocks that
#'  will be added to the custom family section. See [brms::stanvar()] for more
#'  details. The default lists all the possible blocks, but it is unlikely that
#'  you will need all of them. You can specify a vector of only those that you
#'  need. The function will add a section for each block in the list
#' @param open_files Logical; If TRUE (default), the function will open the
#'  template files that were created in RStudio
#'
#' @return If `testing` is TRUE, the function will return the file content as a
#'  string. If `testing` is FALSE, the function will return NULL
#'
#' @details If you get a warning during check() about non-ASCII characters, this
#'  is often due to the citation field. You can find what the problem is by
#'  running
#'  ```r
#'  remotes::install_github("eddelbuettel/dang")
#'  dang::checkPackageAsciiCode(dir = ".")
#'  ```
#'  usually rewriting the numbers (issue, page numbers) manually fixes it
#' @keywords internal developer
#' @export
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' library(usethis)
#'
#'
#' # create a new model file with a brms::custom_family, three .stan files in
#' # inst/stan_chunks/ and open the files
#' use_model_template("abc",
#'   custom_family = TRUE,
#'   stanvar_blocks = c("functions", "likelihood", "tdata")
#' )
#'
use_model_template <- function(model_name,
                               versions = NULL,
                               custom_family = FALSE,
                               stanvar_blocks = c(
                                 "data", "tdata", "parameters",
                                 "tparameters", "model", "likelihood",
                                 "genquant", "functions"
                               ),
                               open_files = TRUE,
                               testing = FALSE) {
  file_name <- paste0("model_", model_name, ".R")

  # check if model exists
  if (model_name %in% model_names()) {
    stop2("Model {model_name} already exists")
  }
  if (file.exists(paste0("R/", file_name))) {
    stop2("File {file_name} already exists")
  }

  versioned <- !is.null(versions)

  model_header <- glue(
    "#############################################################################!
     # MODELS                                                                 ####
     #############################################################################!
     # see 'R/model_ddm.R' (flat defaults) or 'R/model_cswald.R' (versioned) for examples
     #
     # Besides this file, a new model needs entries in:
     # - `data_column_roles` in R/helpers-model.R (a test requires it)
     # - `stored_frame_cases()` in tests/testthat/test-update.R (a test requires it)
     #   plus a `revert_check_data()` method in R/update.R if `check_data()`
     #   consumes or creates columns
     # - `response_annotations()` in R/helpers-model.R, if the response columns
     #   need a unit or a coding note in the console output
     #
     # In this file, `citation` needs at least one reference, one per element, each
     # on a single line without a \"- \" bullet and ending in \".\" or a
     # https://doi.org/ URL (a test requires it)\n\n\n"
  )


  check_data_header <- glue(
    "#############################################################################!
     # CHECK_DATA S3 methods                                                  ####
     #############################################################################!
     # A check_data.* function should be defined for each class of the model.
     # If a model shares methods with other models, the shared methods should be
     # defined in helpers-data.R. Put here only the methods that are specific to
     # the model. See ?check_data for details.
     # (YOU CAN DELETE THIS SECTION IF YOU DO NOT REQUIRE ADDITIONAL DATA CHECKS)\n\n\n"
  )

  bmf2bf_header <- glue(
    "#############################################################################!
     # Convert bmmformula to brmsformla methods                               ####
     #############################################################################!
     # A bmf2bf.* function should be defined if the default method for constructing
     # the brmsformula from the bmmformula does not apply (e.g if aterms are required).
     # The shared method for all `bmmodels` is defined in bmmformula.R.
     # See ?bmf2bf for details.
     # (YOU CAN DELETE THIS SECTION IF YOUR MODEL USES A STANDARD FORMULA WITH 1 RESPONSE VARIABLE)\n\n\n"
  )

  configure_model_header <- glue(
    "#############################################################################!
     # CONFIGURE_MODEL S3 METHODS                                             ####
     #############################################################################!
     # Each model should have a corresponding configure_model.* function. See
     # ?configure_model for more information.\n\n\n"
  )

  postprocess_brm_header <- glue(
    "#############################################################################!
     # POSTPROCESS METHODS                                                    ####
     #############################################################################!
     # A postprocess_brm.* function should be defined for the model class. See
     # ?postprocess_brm for details\n\n\n"
  )


  # the parameter specification block. paste0 (not glue) because the nested
  # list() calls are full of parentheses and commas that glue mishandles.
  spec_body <- paste(
    "  parameters = list(",
    '    par1 = "Parameter 1 = description of parameter 1",',
    '    par2 = "Parameter 2 = description of parameter 2"',
    "  ),",
    "  links = list(",
    '    par1 = "identity",',
    '    par2 = "log"',
    "  ),",
    "  fixed_parameters = list(",
    "    mu = 0",
    "  ),",
    "  # pick the sd rate by meaning: 1 for sensitivity, strength, concentration",
    "  # (kappa), mixing weights and identity-linked drift; 2 for criteria,",
    "  # thresholds, boundary, ndt, start point, log-linked drift, diffusion",
    "  # constant, log ratios, log SDs and correlations; 4 for circular bias",
    "  priors = list(",
    '    par1 = list(main = "normal(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),',
    '    par2 = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)", sd = "exponential(2)")',
    "  ),",
    "  init_ranges = list(",
    "    par1 = c(-1, 1),",
    "    par2 = c(0.5, 1.5)",
    "  )",
    sep = "\n"
  )

  if (versioned) {
    indented_body <- gsub("(^|\n)", "\\1  ", spec_body)
    version_entries <- vapply(versions, function(v) {
      paste0("  ", v, " = list(\n", indented_body, "\n  )")
    }, character(1))
    defaults_block <- paste0(
      ".", model_name, "_version_table <- list(\n",
      paste(version_entries, collapse = ",\n"), "\n)\n\n\n"
    )
  } else {
    defaults_block <- paste0(
      ".", model_name, "_defaults <- list(\n", spec_body, "\n)\n\n\n"
    )
  }

  if (versioned) {
    model_object <- glue('
      .model_<<model_name>> <- function(resp_var1 = NULL, required_arg1 = NULL, required_arg2 = NULL,
                                        links = NULL, version = "<<versions[1]>>", call = NULL, ...) {
        out <- structure(
          list(
            resp_vars = nlist(resp_var1),
            other_vars = nlist(required_arg1, required_arg2),
            domain = "",
            task = "",
            name = "",
            citation = character(),
            version = version,
            requirements = "",
            parameters = .<<model_name>>_version_table[[version]][["parameters"]],
            links = .<<model_name>>_version_table[[version]][["links"]],
            fixed_parameters = .<<model_name>>_version_table[[version]][["fixed_parameters"]],
            default_priors = .<<model_name>>_version_table[[version]][["priors"]],
            init_ranges = .<<model_name>>_version_table[[version]][["init_ranges"]]
          ),
          class = c("bmmodel", "<<model_name>>", paste0("<<model_name>>_", version)),
          call = call
        )
        out <- set_links(out, links)
        out
      }

      # uncomment if configure_model() builds the links into the family or into
      # the non-linear formulas rather than reading them from the list above,
      # so that a link set by the user is refused instead of silently ignored:
      # #\' @exportS3Method
      # settable_links.<<model_name>> <- function(model) character(0)',
      .open = "<<", .close = ">>"
    )
  } else {
    model_object <- glue('
      .model_<<model_name>> <- function(resp_var1 = NULL, required_arg1 = NULL, required_arg2 = NULL,
                                        links = NULL, call = NULL, ...) {
        out <- structure(
          list(
            resp_vars = nlist(resp_var1),
            other_vars = nlist(required_arg1, required_arg2),
            domain = "",
            task = "",
            name = "",
            citation = character(),
            version = "NA",
            requirements = "",
            parameters = .<<model_name>>_defaults[["parameters"]],
            links = .<<model_name>>_defaults[["links"]],
            fixed_parameters = .<<model_name>>_defaults[["fixed_parameters"]],
            default_priors = .<<model_name>>_defaults[["priors"]],
            init_ranges = .<<model_name>>_defaults[["init_ranges"]]
          ),
          class = c("bmmodel", "<<model_name>>"),
          call = call
        )
        out <- set_links(out, links)
        out
      }

      # uncomment if configure_model() builds the links into the family or into
      # the non-linear formulas rather than reading them from the list above,
      # so that a link set by the user is refused instead of silently ignored:
      # #\' @exportS3Method
      # settable_links.<<model_name>> <- function(model) character(0)',
      .open = "<<", .close = ">>"
    )
  }

  # the dependency check is commented out; uncomment and adapt if the model
  # requires a specific backend (see e.g. 'R/model_ddm.R', which needs cmdstanr)
  dependency_check <- paste(
    "   # uncomment if your model requires a specific backend:",
    '   # stopif(!requireNamespace("cmdstanr", quietly = TRUE),',
    "   #        'The \"cmdstanr\" package is required for this model.')",
    sep = "\n"
  )

  param_docs <- c(
    "#' @param resp_var1 A description of the response variable",
    "#' @param required_arg1 A description of the required argument",
    "#' @param required_arg2 A description of the required argument",
    "#' @param links A list of links for the model parameters."
  )

  if (versioned) {
    param_docs <- c(param_docs, glue(
      "#' @param version A character string selecting the model version. \\
       One of <<collapse_comma(versions)>>.",
      .open = "<<", .close = ">>"
    ))
    version_formal <- glue(", version = c(<<collapse_comma(versions)>>)",
      .open = "<<", .close = ">>"
    )
    version_validation <- "   version <- match.arg(version)\n"
    version_pass <- "version = version, "
    alias_example_ref <- "R/model_cswald.R"
  } else {
    version_formal <- ""
    version_validation <- ""
    version_pass <- ""
    alias_example_ref <- "R/model_ddm.R"
  }

  # assemble the @param block as one contiguous chunk so an absent version
  # parameter never leaves a blank line that would split the roxygen block
  params_doc <- paste(
    c(param_docs, "#' @param ... used internally for testing, ignore it"),
    collapse = "\n"
  )

  user_facing_alias <- glue("
    # user facing alias
    # information in the title and details sections will be filled in
    # automatically based on the information in the .model_<<model_name>>()$info\n
    #\' @title `r .model_<<model_name>>()$name`
    #\' @name <<model_name>>
    #\' @details `r model_docs(.model_<<model_name>>())`
    <<params_doc>>
    #\' @return An object of class `bmmodel`
    #\' @keywords bmmodel
    #\' @export
    #\' @examples
    #\' \\dontrun{
    #\' # put a full example here (see '<<alias_example_ref>>' for an example)
    #\' }
    <<model_name>> <- function(resp_var1, required_arg1, required_arg2, links = NULL<<version_formal>>, ...) {
       call <- match.call()
       stop_missing_args()
    <<version_validation>><<dependency_check>>
       .model_<<model_name>>(resp_var1 = resp_var1, required_arg1 = required_arg1, required_arg2 = required_arg2,
                    links = links, <<version_pass>>call = call, ...)
    }\n\n\n",
    .open = "<<", .close = ">>"
  )

  check_data_method <- glue(
    "#' @export
    check_data.<<model_name>> <- function(model, data, formula) {
       # retrieve required arguments
       required_arg1 <- model$other_vars$required_arg1
       required_arg2 <- model$other_vars$required_arg2\n
       # check the data (required)\n
       # compute any necessary transformations (optional)\n
       # save some variables as attributes of the data for later use (optional)\n
       NextMethod('check_data')
    }\n\n\n",
    .open = "<<", .close = ">>"
  )

  # add bmf2bf method if necessary
  bmf2bf_method <- glue("#' @export
    bmf2bf.<<model_name>> <- function(model, formula) {
       # retrieve the variables the formula needs
       resp_var1 <- model$resp_vars$resp_var1
       required_arg1 <- model$other_vars$required_arg1\n
       # set the base brmsformula with the response and its addition terms
       brms_formula <- brms::bf(paste0(resp_var1, \" | vreal(\", required_arg1, \") ~ 1\"))\n
       # return the brms_formula to add the remaining bmmformulas to it.
       brms_formula
    }\n\n\n",
    .open = "<<", .close = ">>"
  )


  # add custom family section if custom_family is TRUE
  # PS: do not try to replace with glue - already wasted enough time, it doesn't work well
  if (custom_family) {
    family_template <- paste0(
      "   <<model_name>>_family <- brms::custom_family(\n",
      "     '<<model_name>>',\n",
      "     dpars = c(),\n",
      "     links = c(),\n",
      "     lb = c(), # lower bounds for parameters\n",
      "     ub = c(), # upper bounds for parameters\n",
      "     type = '', # real for continuous dv, int for discrete dv\n",
      "     loop = TRUE, # FALSE if the Stan likelihood is vectorized over observations\n",
      "     log_lik = log_lik_<<model_name>>,\n",
      "     posterior_predict = posterior_predict_<<model_name>>,\n",
      "     posterior_epred = posterior_epred_<<model_name>>\n",
      "   )\n   formula$family <- <<model_name>>_family\n\n"
    )

    stan_vars_template <- paste0(
      "   # prepare initial stanvars to pass to brms, model formula and priors\n",
      "   sc_path <- system.file('stan_chunks', package='bmm')\n"
    )

    for (stanvar_block in stanvar_blocks) {
      stan_vars_file <- glue("inst/stan_chunks/{model_name}_{stanvar_block}.stan")
      if (!testing) {
        file.create(stan_vars_file)
        if (open_files) {
          usethis::edit_file(stan_vars_file)
        }
      }
      # PS: do not try to replace with glue - already wasted enough time, it doesn't work well
      stan_vars_template <- paste0(
        stan_vars_template,
        "   stan_", stanvar_block, " <- read_lines2(paste0(sc_path, '/", model_name, "_", stanvar_block, ".stan'))\n"
      )
    }
    stan_vars_template <- paste0(stan_vars_template, "\n   stanvars <- ")
    i <- 1
    for (stanvar_block in stanvar_blocks) {
      if (i < length(stanvar_blocks)) {
        stan_vars_template <- paste0(stan_vars_template, "brms::stanvar(scode = stan_", stanvar_block, ", block = '", stanvar_block, "') +\n      ")
        i <- i + 1
      } else {
        stan_vars_template <- paste0(stan_vars_template, "brms::stanvar(scode = stan_", stanvar_block, ", block = '", stanvar_block, "')\n\n")
      }
    }
    out_template <- "   nlist(formula, data, stanvars)\n"

    family_functions <- paste0(
      "\n\n#############################################################################!\n",
      "# LOG_LIK, POSTERIOR_PREDICT & POSTERIOR_EPRED                           ####\n",
      "#############################################################################!\n",
      "# see posterior_epred_sdt_yn() in 'R/model_sdt_yn.R' and posterior_epred_ddm()\n",
      "# in 'R/model_ddm.R' for examples\n\n",
      "# returns one log-likelihood value per posterior draw for observation i\n",
      "log_lik_", model_name, " <- function(i, prep) {\n}\n\n",
      "# returns one simulated response per posterior draw for observation i\n",
      "posterior_predict_", model_name, " <- function(i, prep, ...) {\n}\n\n",
      "# returns a draws x observations matrix: build it with .epred_matrix() and\n",
      "# line up data columns with the draws using .epred_data()\n",
      "posterior_epred_", model_name, " <- function(prep) {\n}\n",
      "# if the model's expected response is not meaningful, use instead:\n",
      "# posterior_epred_", model_name, " <- posterior_epred_undefined(\"", model_name, "\")\n\n\n"
    )
  } else {
    family_functions <- ""
    stan_vars_template <- ""
    family_template <- "   formula$family <- NULL\n\n"
    out_template <- "   nlist(formula, data)\n"
  }


  family_comment <- ifelse(custom_family,
    "   # construct the family & add to formula object\n",
    "   # add family to formula object\n"
  )

  # PS: do not try to replace with glue - already wasted enough time, it doesn't work well
  configure_model_method <- glue::glue("#' @export\n",
    "configure_model.<<model_name>> <- function(model, data, formula) {\n",
    "   # retrieve required arguments\n",
    "   required_arg1 <- model$other_vars$required_arg1\n",
    "   required_arg2 <- model$other_vars$required_arg2\n\n",
    "   # retrieve arguments from the data check\n",
    "   my_precomputed_var <- attr(data, 'my_precomputed_var')\n\n",
    "   # construct brms formula from the bmm formula\n",
    "   formula <- bmf2bf(model, formula)\n\n",
    family_comment,
    family_template,
    stan_vars_template,
    "   # return the list\n",
    out_template,
    "}\n\n",
    .open = "<<", .close = ">>"
  )

  postprocess_brm_method <- glue(
    "#' @export
    postprocess_brm.<<model_name>> <- function(model, fit, ...) {
       # any required postprocessing (if none, delete this section)
       fit
    }\n",
    .open = "<<", .close = ">>"
  )

  file_content <- paste0(
    model_header,
    defaults_block,
    model_object,
    "\n\n",
    user_facing_alias,
    check_data_header,
    check_data_method,
    bmf2bf_header,
    bmf2bf_method,
    configure_model_header,
    configure_model_method,
    family_functions,
    postprocess_brm_header,
    postprocess_brm_method
  )

  if (!testing) {
    writeLines(file_content, paste0("R/", file_name))
    if (open_files) {
      usethis::edit_file(paste0("R/", file_name))
    }
  } else {
    cat(file_content)
  }
}


#' @title Generate Stan code for bmm models
#' @description Given the `model`, the `data` and the `formula` for the model,
#'   this function will return the combined stan code generated by `bmm` and
#'   `brms`
#'
#' @inheritParams bmm
#' @aliases stancode
#' @param object A `bmmformula` object
#' @param ... Further arguments passed to [brms::stancode()]. See the
#'   description of [brms::stancode()] for more details
#'
#' @return A character string containing the fully commented Stan code to fit a
#'   bmm model.
#'
#' @seealso [bmm_models()], [brms::stancode()]
#' @keywords extract_info
#' @examples
#' scode1 <- stancode(bmf(c ~ 1, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#' cat(scode1)
#' @importFrom brms stancode
#' @export
stancode.bmmformula <- function(object, data, model, prior = NULL, ...) {
  withr::local_options(bmm.sort_data = FALSE)
  dots <- list(...)
  local_brms_threads(dots)
  add_bmm_version_to_stancode(call_brms_extractor(
    brms::stancode,
    configure_fit(object, data, model, prior, until = "prior", frame_args = brms_frame_args(dots)),
    dots
  ))
}

# Everything brm() needs to fit a specification: the checked model, the brms
# formula with its stanvars, the prior, and the initial values, plus the model
# and user formula that postprocess_brm() stores in the fit. Shared by bmm()
# and the stancode/standata/default_prior extractors. `until` says how far to
# go, because each step costs a brms frame build: "init" for bmm(), "prior"
# for stancode() and default_prior(), "model" for standata(), whose Stan data
# depends on neither. `frame_args` are the brm() arguments besides the formula
# and data that shape the model frame (see brms_frame_args())
configure_fit <- function(formula, data = NULL, model = NULL, prior = NULL, until = "init",
                          frame_args = list()) {
  UseMethod("configure_fit")
}

#' @export
configure_fit.bmmformula <- function(formula, data = NULL, model = NULL, prior = NULL,
                                     until = "init", frame_args = list()) {
  user_formula <- formula
  # the user's expression for the data lives in the caller of the generic
  # (bmm() or an extractor). It is set before check_data(), whose methods print
  # it (order_data_query()), and again after, because some methods rebuild the
  # data frame without it
  data_name <- substitute_name(data, envir = parent.frame())
  if (!missing(data) && !is.null(data) && (is.list(data) || is.atomic(data))) {
    attr(data, "data_name") <- data_name
  }
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  attr(data, "data_name") <- data_name
  formula <- check_formula(model, data, formula)
  config_args <- configure_model(model, data, formula)
  if (until == "model") {
    return(nlist(config_args, model, user_formula))
  }
  prior <- brms::do_call(
    configure_prior, c(list(model, data, config_args$formula, prior), frame_args)
  )
  # the prior decides which parameters exist, so the inits are built from it
  if (until == "init") {
    config_args$init <- brms::do_call(
      create_initfun, c(list(model, data, config_args$formula, prior), frame_args)
    )
  }
  nlist(config_args, prior, model, user_formula)
}

#' @export
configure_fit.default <- function(formula, data = NULL, model = NULL, ...) {
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  check_formula(model, data, formula)
}

# The brms extractors take the formula as `object`
call_brms_extractor <- function(fun, cfg, dots) {
  args <- combine_args(nlist(config_args = cfg$config_args, dots, prior = cfg$prior))
  args$object <- args$formula
  args$formula <- NULL
  brms::do_call(fun, args)
}

add_bmm_version_to_stancode <- function(stancode) {
  version <- packageVersion("bmm")
  text <- paste0("and bmm ", version)
  brms_comp <- regexpr("brms.*(?=\\n)", stancode, perl = T)
  insert_loc <- brms_comp + attr(brms_comp, "match.length") - 1
  new_stancode <- paste0(
    substr(stancode, 1, insert_loc),
    " ", text,
    substr(stancode, insert_loc + 1, nchar(stancode))
  )
  class(new_stancode) <- class(stancode)
  new_stancode
}

############################################################################# !
# StanCode Helper FUNCTIONS                                              ####
############################################################################# !

# Return integer positions of all occurrences of a fixed pattern in x
# simple wrapper around gregexpr that return a 0-length vector if no matches
# instead of -1
which_positions <- function(x, pat) {
  m <- gregexpr(pat, x, fixed = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L) integer(0) else as.integer(m)
}


#' Find the matching closing brace for an opening brace position
#'
#' Given a string containing brace-like delimiters, return the index of the
#' closing brace that matches the opening brace at \code{open_pos}, accounting
#' for nested braces.
#'
#' The function operates on indices of all occurrences of \code{open_brace} and
#' \code{close_brace} in \code{x} and uses a cumulative-sum depth counter to
#' identify the first position at which nesting depth returns to zero.
#'
#' @param x A length-1 character string.
#' @param open_pos A 1-based integer index into \code{x} indicating the position
#'   of an opening brace character.
#' @param open_brace A length-1 character string giving the opening delimiter
#'   (default \code{"\{"}).
#' @param close_brace A length-1 character string giving the closing delimiter
#'   (default \code{"\}"}).
#'
#' @returns A single integer: the 1-based index of the matching closing brace in
#'   \code{x}.
#' @noRd
#' @examples
#' find_matching_brace("{}", 1L)
#'
#' x <- "{{}{}}"
#' find_matching_brace(x, 1L) # outer brace -> 6
#' find_matching_brace(x, 2L) # inner brace -> 3
#' find_matching_brace(x, 4L) # inner brace -> 5
#'
#' y <- "abc{def{ghi}jkl}mno"
#' find_matching_brace(y, 4L) # -> 16
#' find_matching_brace(y, 8L) # -> 12
find_matching_brace <- function(x, open_pos,
                                open_brace = "{", close_brace = "}") {
  stopifnot(is.character(x), length(x) == 1L)
  n <- nchar(x)
  stopifnot(length(open_pos) == 1L, open_pos >= 1L, open_pos <= n)
  stopif(
    substr(x, open_pos, open_pos) != open_brace,
    "Character at open_pos {open_pos} is not open_brace."
  )

  opens <- which_positions(x, open_brace)
  closes <- which_positions(x, close_brace)

  # combine + sort
  pos <- c(opens, closes)
  stopif(!length(pos), "No braces found.")
  delta <- c(rep.int(1L, length(opens)), rep.int(-1L, length(closes)))

  delta <- delta[order(pos)]
  pos <- pos[order(pos)]

  # locate the opening brace occurrence in the brace stream
  k0 <- match(as.integer(open_pos), pos)
  stopif(is.na(k0), "open_pos was not found among opening brace indices (unexpected).")
  stopif(delta[k0] != 1L, "open_pos does not correspond to an opening brace in the stream.")

  # cumulative sum from that opening brace
  cs <- cumsum(delta[k0:length(delta)])

  # first return to 0 gives the matching close
  j <- which(cs == 0L)[1]
  stopif(is.na(j), "No matching closing brace found (unbalanced braces?).")

  pos[k0 + j - 1L]
}

#' @title Extract code from different STAN program blocks
#'
#' @description
#'   This function extracts the code from the different program blocks of a STAN
#'   program. This can be used in combination with the `stancode` function to
#'   access information about the STAN code generated by `brms` and `bmm`.
#'
#' @param stan_code The STAN code for which the elements should be extracted
#' @param blocks A character vector specifying for which program blocks
#'   the code should be extracted. The default extracks all standard blocks:
#'   "functions", "data", "transformed data", "parameters", "transformed parameters",
#'   "model", and "generated quantities"
#'
#' @return A named list with each element containing the code of one of the STAN
#'   program blocks. If a block
#'
#' @keywords extract_info
#'
#' @examples
#' # generate simple stan code from brms
#' stan_code <- stancode(brms::bf(x ~ 1), data = data.frame(x = rnorm(100)))
#'
#' extracted_program_blocks <- extract_stan_blocks(stan_code)
#'
#' @export
extract_stan_blocks <- function(stan_code, blocks = c("functions", "data", "transformed data", "parameters", "transformed parameters", "model", "generated quantities")) {
  stopifnot(is.character(stan_code), length(stan_code) == 1L)
  blocks <- match.arg(blocks, several.ok = TRUE)
  out <- lapply(blocks, function(b) .extract_stan_block(stan_code, b))
  names(out) <- blocks
  out
}

.extract_stan_block <- function(stan_code, block,
                                include_braces = FALSE,
                                trim = TRUE) {
  stopifnot(is.character(block), length(block) == 1L)

  # Anchor to start-of-line (multiline mode), allow leading spaces/tabs.
  # Match the exact block name, then optional whitespace, then '{'.
  header_pat <- paste0("(?m)^[[:space:]]*", block, "[[:space:]]*\\{")
  m <- regexpr(header_pat, stan_code, perl = TRUE)
  if (m[1] == -1L) {
    return(NULL)
  }

  open_pos <- m[1] + attr(m, "match.length") - 1L # points at '{'
  close_pos <- find_matching_brace(stan_code, open_pos)

  out <- if (include_braces) {
    substr(stan_code, open_pos, close_pos)
  } else {
    substr(stan_code, open_pos + 1L, close_pos - 1L)
  }

  if (trim) out <- sub("^\\s+|\\s+$", "", out)
  out
}


#' @title Extract dimension from parameters in STAN parameter block
#'
#' @description
#'  This functions extracts the names, dimensions, and types from a compiled STAN
#'  parameters blocks generated by bmm or brms. This function is used to specify
#'  initial values for bmm models.
#'
#' @param parameters_block The parameters block extracted via `extract_stan_blocks`
#'
#' @return A list of all parameters, their types, and dimensions as as specified in
#'   the STAN data generated by bmm and brms. A declaration whose type this
#'   function does not model is returned with its name and `NA` for its type and
#'   dimensions, so that only that parameter is left without initial values
#'
#' @keywords extract_info
#'
#' @examples
#' # generate simple stan code from brms
#' stan_code <- stancode(brms::bf(x ~ 1 + cond + (1 + cond | ID)),
#'   data = data.frame(
#'     x = rnorm(100),
#'     ID = rep(1:50, each = 2),
#'     cond = rep(1:2, times = 50)
#'   )
#' )
#'
#' extracted_program_blocks <- extract_stan_blocks(stan_code)
#'
#' par_dims <- extract_parameter_dimensions(extracted_program_blocks$parameters) #'
#'
#' @export
extract_parameter_dimensions <- function(parameters_block) {
  lines <- unlist(strsplit(parameters_block, "\n"))
  lines <- trimws(lines)
  lines <- gsub("//.*", "", lines)
  lines <- lines[nzchar(lines)]

  res <- lapply(lines, parse_parameters_line)
  names(res) <- vapply(res, `[[`, character(1), "name")
  res
}

parse_parameters_line <- function(x) {
  # strip trailing comments and semicolon; normalize whitespace
  x <- sub("//.*$", "", x)
  x <- trimws(x)
  x <- sub(";$", "", x)
  if (!nzchar(x)) stop2("Empty or comment-only line.")

  # 1) optional leading array[...] prefix
  array_dims <- character(0)
  if (grepl("^array\\s*\\[", x, perl = TRUE)) {
    taken <- take_dims(x)
    array_dims <- taken$dims
    x <- trimws(taken$rest)
  }
  is_array <- length(array_dims) > 0

  # 2) detect base type
  base_types <- c(
    "cholesky_factor_corr", "cholesky_factor_cov",
    "corr_matrix", "cov_matrix",
    "row_vector", "unit_vector", "positive_ordered", "simplex", "ordered",
    "vector", "matrix", "real", "int"
  )
  base_type <- NULL
  for (bt in base_types) {
    if (grepl(paste0("^", bt, "\\b"), x, perl = TRUE)) {
      base_type <- bt
      break
    }
  }
  # a declaration of a type Stan has and this parser does not (sum_to_zero_vector,
  # complex, tuple) is reported without one, so that it costs the sampler's own
  # initial value for that parameter and not the whole init list
  if (is.null(base_type)) {
    return(list(
      name = sub(".*[ \\]]", "", x, perl = TRUE), type = NA_character_,
      types = NA_character_, dims = NA_character_, bounds = NULL
    ))
  }

  # consume the base type token
  x_after_bt <- sub(paste0("^", base_type), "", x, perl = TRUE)

  # 3) constraints <...>
  bounds <- parse_bounds(x_after_bt)
  if (!is.null(bounds)) {
    x_after_bt <- sub("^\\s*<[^>]*>\\s*", "", x_after_bt, perl = TRUE)
  }

  # 4) base-type dims [ ... ] (needed for most non-scalars)
  base_dims <- character(0)
  if (grepl("^\\s*\\[", x_after_bt, perl = TRUE)) {
    taken <- take_dims(x_after_bt)
    base_dims <- taken$dims
    x_after_bt <- trimws(taken$rest)
  } else {
    if (base_type %in% c(
      "vector", "row_vector", "matrix",
      "simplex", "unit_vector", "ordered", "positive_ordered",
      "corr_matrix", "cov_matrix", "cholesky_factor_corr", "cholesky_factor_cov"
    )) {
      stop2("Missing dimensions for type {base_type}.")
    }
  }

  # 5) name
  name <- trimws(x_after_bt)
  if (!nzchar(name)) stop2("Missing parameter name.")

  # 6) normalize dims by base type
  dims_by_type <- switch(base_type,
    real                 = 1,
    int                  = 1,
    vector               = base_dims[1],
    row_vector           = base_dims[1],
    simplex              = base_dims[1],
    unit_vector          = base_dims[1],
    ordered              = base_dims[1],
    positive_ordered     = base_dims[1],
    matrix               = base_dims[1:2],
    corr_matrix          = base_dims[1],
    cov_matrix           = base_dims[1],
    cholesky_factor_corr = base_dims[1],
    cholesky_factor_cov  = base_dims[1],
    character(0)
  )
  dims <- c(array_dims, dims_by_type)

  # NEW: dual typing
  types <- if (is_array) c("array", base_type) else base_type

  list(
    name    = name,
    type    = base_type, # backward-compat: base type only
    types   = types, # NEW: includes "array" when applicable
    dims    = dims,
    bounds  = bounds
  )
}

# The dimensions in the leading [...] of a declaration and the text after it.
# brms sizes the parameters of mo() and s() by an element of a data array
# (simplex[Jmo_c[1]], vector[knots_kappa_1[1]]), so the closing bracket is the
# one that balances the opening one, and only a top-level comma separates
# dimensions
take_dims <- function(x) {
  open <- regexpr("[", x, fixed = TRUE)
  close <- find_matching_brace(x, open, "[", "]")
  dims <- trimws(strsplit(substr(x, open + 1L, close - 1L), ",(?![^\\[]*\\])", perl = TRUE)[[1]])
  list(dims = dims[nzchar(dims)], rest = substring(x, close + 1L))
}

# helper: parse <...> constraints into a named list
parse_bounds <- function(s) {
  if (!grepl("<[^>]*>", s, perl = TRUE)) {
    return(NULL)
  }
  inside <- sub(".*?<([^>]*)>.*", "\\1", s, perl = TRUE)
  parts <- trimws(unlist(strsplit(inside, ",")))
  kvs <- lapply(parts, function(p) {
    if (!grepl("=", p, fixed = TRUE)) {
      return(NULL)
    }
    sp <- strsplit(p, "=", fixed = TRUE)[[1]]
    setNames(list(trimws(sp[2])), trimws(sp[1]))
  })
  # merge into a single named list
  kvs <- Filter(Negate(is.null), kvs)
  if (!length(kvs)) {
    return(list())
  }
  Reduce(function(a, b) c(a, b), kvs)
}
