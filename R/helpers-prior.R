#' @title Get Default priors for Measurement Models specified in BMM
#' @description Obtain the default priors for a Bayesian multilevel measurement
#'   model, as well as information for which parameters priors can be specified.
#'   Given the `model`, the `data` and the `formula` for the model, this
#'   function will return the default priors that would be used to estimate the
#'   model. Additionally, it will return all model parameters that have no prior
#'   specified (flat priors). This can help to get an idea about which priors
#'   need to be specified and also know which priors were used if no
#'   user-specified priors were passed to the [bmm()] function.
#'
#'   The default priors in `bmm` tend to be more informative than the default
#'   priors in `brms`, as we use domain knowledge to specify the priors.
#'
#' @inheritParams bmm
#' @aliases default_prior
#' @param object A `bmmformula` object
#' @param ... Further arguments passed to [brms::default_prior()]
#'
#' @return A data.frame with columns specifying the `prior`, the `class`, the
#'   `coef` and `group` for each of the priors specified. Separate rows contain
#'   the information on the parameters (or parameter classes) for which priors
#'   can be specified.
#'
#' @seealso [supported_models()], [brms::default_prior()]
#'
#' @keywords extract_info
#'
#' @examples
#' default_prior(bmf(c ~ 1, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#' @export
default_prior.bmmformula <- function(object, data, model, formula = object, ...) {
  withr::local_options(bmm.sort_data = FALSE)

  formula <- object
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)
  config_args <- configure_model(model, data, formula)
  prior <- configure_prior(model, data, config_args$formula, user_prior = NULL)

  dots <- list(...)
  prior_args <- combine_args(nlist(config_args, dots, prior))
  prior_args$object <- prior_args$formula
  prior_args$formula <- NULL

  brms_priors <- brms::do_call(brms::default_prior, prior_args)

  combine_prior(brms_priors, prior_args$prior)
}

#' @title Report the priors used in a fitted bmm model
#' @description For each parameter of a fitted bmm model, reports the link
#'   function, the prior actually applied on the sampling (link) scale, and
#'   where that prior came from: a bmm default, a brms default, or a
#'   user-specified prior. Parameters without a proper prior are flagged as
#'   flat.
#' @param fit A `bmmfit` object returned by [bmm()]
#' @param format Character. `"table"` (default) prints the report as a table;
#'   `"text"` prints sentences ready for a methods section.
#' @details The provenance of each prior is determined by re-deriving the
#'   default priors from the model, formula and data stored in the fit. A
#'   user-specified prior that is identical to the bmm default is therefore
#'   reported as a default. Coefficients that inherit their prior from a more
#'   general class are collapsed into the row of the prior they inherit from.
#'
#'   Flat priors are flagged because they are improper: Bayes factors via
#'   bridge sampling are undefined when any parameter has an improper prior.
#' @return A `data.frame` of class `bmm_report_priors` with columns
#'   `parameter`, `link`, `class`, `coef`, `group`, `prior` and `source`.
#' @seealso [default_prior()], [parameters()]
#' @keywords extract_info
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' fit <- bmm(
#'   bmf(c ~ 0 + set_size, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#' report_priors(fit)
#' report_priors(fit, format = "text")
#' @export
report_priors <- function(fit, format = "table") {
  stopif(!inherits(fit, "bmmfit"), "The fit argument must be a bmmfit object returned by bmm()")
  stopif(is.null(fit$prior), "The fit object contains no prior information")
  format <- match.arg(format, c("table", "text"))
  fit <- restructure(fit)
  structure(
    prior_provenance(fit),
    class = c("bmm_report_priors", "data.frame"),
    model_name = fit$bmm$model$name,
    par_labels = unlist(fit$bmm$model$parameters),
    format = format
  )
}

#' @export
print.bmm_report_priors <- function(x, ...) {
  if (identical(attr(x, "format"), "text")) {
    cat(strwrap(prior_report_text(x), width = 80), sep = "\n")
    return(invisible(x))
  }

  model_name <- attr(x, "model_name")
  if (!is.null(model_name) && nzchar(model_name)) {
    cat(style("purple1")("Model: "), model_name, "\n\n")
  }

  print_df <- as.data.frame(x)
  print_df$prior[x$source == "flat"] <- "(flat)"
  all_empty <- vapply(print_df, function(col) all(is.na(col) | !nzchar(col)), logical(1))
  print_df <- print_df[, !all_empty, drop = FALSE]
  for (col in names(print_df)) {
    print_df[[col]][is.na(print_df[[col]]) | !nzchar(print_df[[col]])] <- "--"
  }
  print.data.frame(print_df, right = FALSE, row.names = FALSE)

  if (any(x$source == "flat")) {
    flat_note <- "Flat priors are improper: Bayes factors via bridge sampling are
      undefined unless you specify proper priors for these coefficients."
    cat("\n", strwrap(gsub("\\s+", " ", flat_note), width = 80), sep = "\n")
  }
  invisible(x)
}

# fit -> tidy classified prior table; the shared extraction layer for #194
prior_provenance <- function(fit) {
  # force defaults on so reconstruction is deterministic regardless of session
  # options; fits made with bmm.default_priors = FALSE still classify correctly
  # because their flat rows can never match a non-empty default
  withr::local_options(bmm.default_priors = TRUE)
  model <- fit$bmm$model
  defaults <- suppressWarnings(suppressMessages({
    # reconstruct from the post-pipeline formula and model frame stored on the
    # fit instead of re-running the data pipeline: brms drops raw response
    # columns from the model frame for some models (e.g. m3), so check_data
    # cannot be re-run there; it is still tried because it restores helper
    # columns that model-specific configure_prior methods inspect (ss_numeric)
    data <- tryCatch(
      check_data(model, fit$data, fit$bmm$user_formula),
      error = function(e) fit$data
    )
    combine_prior(
      brms::default_prior(fit$formula, data = fit$data),
      configure_prior(model, data, fit$formula, user_prior = NULL)
    )
  }))
  out <- classify_priors(fit$prior, defaults, links = model$links)
  # constants the user set in the formula are folded into the model object by
  # check_model, so the reconstruction reports them as defaults
  user_fixed <- names(fit$bmm$user_formula)[is_constant(fit$bmm$user_formula)]
  out$source[out$parameter %in% user_fixed & out$source == "bmm default"] <- "user"
  out
}

# compare the prior table of a fit against a freshly reconstructed default
# prior table to determine the provenance of each prior; brms stamps every
# prior passed to brm(prior = ...) as "user", so the source column of the fit
# cannot distinguish bmm defaults from user-set priors
classify_priors <- function(prior, defaults, links = list()) {
  key_cols <- c("class", "dpar", "nlpar", "coef", "group", "resp")
  eff <- resolve_effective_prior(prior)
  keys <- do.call(paste, prior[key_cols])
  def_keys <- do.call(paste, defaults[key_cols])
  def_eff <- resolve_effective_prior(defaults)[match(keys, def_keys)]
  def_eff[is.na(def_eff)] <- ""
  from_bmm <- keys %in% def_keys[defaults$source == "user"]

  source <- ifelse(
    !nzchar(eff), "flat",
    ifelse(eff == def_eff & from_bmm, "bmm default",
      ifelse(eff == def_eff, "brms default", "user")
    )
  )

  parameter <- ifelse(
    nzchar(prior$nlpar), prior$nlpar,
    ifelse(nzchar(prior$dpar), prior$dpar,
      ifelse(prior$class %in% c("b", "Intercept"), "mu", NA)
    )
  )
  link <- vapply(parameter, function(p) {
    if (is.na(p)) NA_character_ else links[[p]] %||% "identity"
  }, character(1), USE.NAMES = FALSE)

  out <- data.frame(
    parameter, link,
    class = prior$class, coef = prior$coef, group = prior$group,
    prior = prior$prior, source,
    stringsAsFactors = FALSE
  )
  keep <- nzchar(prior$prior) | (!has_parent_prior(prior) & !vacuous_flat_prior(prior, eff))
  out <- out[keep, , drop = FALSE]
  row.names(out) <- NULL
  out
}

# brms semantics: an empty prior string inherits from the row with the same
# class/dpar/nlpar/resp at coef = "" (first within the same group, then at
# group = ""); a row that stays empty after resolution has a flat prior
resolve_effective_prior <- function(prior) {
  base <- paste(prior$class, prior$dpar, prior$nlpar, prior$resp)
  eff <- prior$prior
  parent_prior <- function(i, group) {
    hit <- which(base == base[i] & prior$group == group & !nzchar(prior$coef) & nzchar(prior$prior))
    if (length(hit)) prior$prior[hit[1]] else ""
  }
  for (i in which(!nzchar(eff))) {
    eff[i] <- parent_prior(i, prior$group[i])
    if (!nzchar(eff[i])) eff[i] <- parent_prior(i, "")
  }
  eff
}

# rows that inherit from a parent row shown in the report are redundant
has_parent_prior <- function(prior) {
  base <- paste(prior$class, prior$dpar, prior$nlpar, prior$resp)
  vapply(seq_len(nrow(prior)), function(i) {
    parents <- base == base[i] & !nzchar(prior$coef) &
      (prior$group == prior$group[i] | !nzchar(prior$group))
    parents[i] <- FALSE
    any(parents)
  }, logical(1))
}

# a flat class- or group-level row is vacuous when every row that could
# inherit from it resolved to its own prior (brms emits a class-level b row
# even when the only coefficient below it has an explicit prior)
vacuous_flat_prior <- function(prior, eff) {
  base <- paste(prior$class, prior$dpar, prior$nlpar, prior$resp)
  vapply(seq_len(nrow(prior)), function(i) {
    if (nzchar(eff[i]) || nzchar(prior$coef[i])) {
      return(FALSE)
    }
    children <- base == base[i] &
      (prior$group == prior$group[i] | !nzchar(prior$group[i]))
    children[i] <- FALSE
    any(children) && all(nzchar(eff[children]))
  }, logical(1))
}

prior_report_text <- function(x) {
  par_labels <- attr(x, "par_labels")
  in_par <- !is.na(x$parameter)
  sentences <- unlist(lapply(unique(x$parameter[in_par]), function(par) {
    parameter_prior_sentences(x[in_par & x$parameter == par, ], par_labels)
  }))
  standalone <- vapply(which(!in_par), function(i) {
    scope <- prior_scope_phrase(x[i, ])
    if (!nzchar(x$prior[i])) {
      glue("No prior was specified for {scope}; an improper flat prior was used.")
    } else {
      glue("{upfirst(scope)} received a {x$prior[i]} prior ({x$source[i]}).")
    }
  }, character(1))
  paste(
    c(
      "Priors were specified on the sampling (link) scale of each parameter.",
      sentences, standalone
    ),
    collapse = " "
  )
}

parameter_prior_sentences <- function(rows, par_labels) {
  subject <- prior_subject_phrase(rows$parameter[1], par_labels)
  is_const <- grepl("^constant\\(", rows$prior)
  sentences <- character(0)
  if (any(is_const)) {
    value <- sub("^constant\\((.*)\\)$", "\\1", rows$prior[is_const])
    sentences <- glue("{upfirst(subject)} was fixed to {value}.")
  }
  rest <- rows[!is_const, , drop = FALSE]
  if (nrow(rest) > 0) {
    clauses <- vapply(seq_len(nrow(rest)), function(i) {
      scope <- prior_scope_phrase(rest[i, ])
      if (!nzchar(rest$prior[i])) {
        glue("{scope} received no prior (improper flat prior)")
      } else {
        glue("{scope} received a {rest$prior[i]} prior ({rest$source[i]})")
      }
    }, character(1))
    sentences <- c(
      sentences,
      glue("For {subject} ({rest$link[1]} link), {collapse_and(clauses)}.")
    )
  }
  sentences
}

# use the model's parameter description only when it names the parameter in
# the form "<label> parameter ..."; descriptions vary too much across models
# to be quoted whole
prior_subject_phrase <- function(par, par_labels) {
  desc <- if (par %in% names(par_labels)) par_labels[[par]] else ""
  label <- trimws(sub("\\s*[Pp]arameter.*$", "", desc))
  if (identical(label, trimws(desc)) || !nzchar(label) ||
    nchar(label) > 40 || tolower(label) == tolower(par)) {
    return(glue("parameter {par}"))
  }
  glue("the {tolower(label)} parameter {par}")
}

collapse_and <- function(x) {
  if (length(x) < 3) {
    return(paste(x, collapse = " and "))
  }
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}

prior_scope_phrase <- function(row) {
  grouping <- if (nzchar(row$group)) glue(" (grouping: {row$group})") else ""
  if (row$class == "b" && row$coef == "Intercept") {
    "the intercept"
  } else if (row$class == "b" && nzchar(row$coef)) {
    glue("the coefficient {row$coef}")
  } else if (row$class == "b") {
    "all population-level coefficients"
  } else if (row$class == "Intercept") {
    "the intercept"
  } else if (row$class == "sd") {
    glue("the group-level standard deviations{grouping}")
  } else if (row$class %in% c("cor", "L")) {
    glue("the correlations among group-level effects{grouping}")
  } else {
    glue("the {row$class} parameters")
  }
}

upfirst <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substring(x, 2))
}

#' Generic S3 method for configuring the default prior for a bmmodel
#'
#' Called by bmm() to automatically construct the priors for a given
#' model, data and formula, and combine it with the prior given by the user. The
#' first method executed is configure_prior.bmmodel, which will build the prior
#' based on information from the model object such as fixed_parameters,
#' default_priors, etc. Thus it is important to define these values in the model
#' object. The function will also recognize if the user has specified that some
#' parameters should be fixed to a constant and put the appropriate constant
#' priors. Any additional priors that a developer wants to specify, which are
#' not based on information in the model object, can be defined in the
#' configure_prior.* method for the model. See configure_prior.imm_full for an
#' example.
#' @param model A `bmmodel` object
#' @param data A data.frame containing the data used in the model
#' @param formula A `brmsformula` object returned from configure_model()
#' @param user_prior A `brmsprior` object given by the user as an argument to
#'  bmm()
#' @param ... Additional arguments passed to the method
#'
#' @return A `brmsprior` object containing the default priors for the model
#'
#' @export
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' configure_prior.mixture3p <- function(model, data, formula, user_prior, ...) {
#'   # if there is set_size 1 in the data, set constant prior over thetant for set_size1
#'   prior <- brms::empty_prior()
#'   set_size_var <- model$other_vars$set_size
#'   prior_cond <- any(data$ss_numeric == 1) && !is.numeric(data[[set_size_var]])
#'
#'   thetant_preds <- rhs_vars(formula$pforms$thetant)
#'   if (prior_cond && set_size_var %in% thetant_preds) {
#'     prior <- prior + brms::prior_("constant(-100)",
#'       class = "b",
#'       coef = paste0(set_size_var, 1),
#'       nlpar = "thetant"
#'     )
#'   }
#'   # check if there is a random effect on theetant that include set_size as predictor
#'   bterms <- brms::brmsterms(formula$pforms$thetant)
#'   re_terms <- bterms$dpars$mu$re
#'   if (!is.null(re_terms)) {
#'     for (i in 1:nrow(re_terms)) {
#'       group <- re_terms$group[[i]]
#'       form <- re_terms$form[[i]]
#'       thetant_preds <- rhs_vars(form)
#'
#'       if (prior_cond && set_size_var %in% thetant_preds) {
#'         prior <- prior + brms::prior_("constant(1e-8)",
#'           class = "sd",
#'           coef = paste0(set_size_var, 1),
#'           group = group,
#'           nlpar = "thetant"
#'         )
#'       }
#'     }
#'   }
#'
#'   prior
#' }
#'
#' @keywords internal developer
configure_prior <- function(model, data, formula, user_prior, ...) {
  UseMethod("configure_prior")
}

#' @export
configure_prior.default <- function(model, data, formula, user_prior, ...) {
  NULL
}

#' @export
configure_prior.bmmodel <- function(model, data, formula, user_prior = NULL, ...) {
  prior <- fixed_pars_priors(model, formula)
  default_prior <- set_default_prior(model, data, formula)
  prior <- combine_prior(default_prior, prior)
  prior <- combine_prior(prior, user_prior)
  additional_prior <- NextMethod("configure_prior")
  combine_prior(prior, additional_prior)
}

#' @title Construct constant priors for fixed model parameters
#' @param model a `bmmodel` object
#' @param formula a `brmsformula` object
#' @param additional_pars a list of name=value pairs to fix additional
#'   parameters where the name is the parameter name and the value is the fixed
#'   value
#'
#' @return an object of class brmsprior of the form prior("constant(value)",
#'   class="Intercept", dpar=parameter_name) for all fixed parameters in the
#'   model
#' @noRd
fixed_pars_priors <- function(model, formula, additional_pars = list()) {
  fix_pars <- model$fixed_parameters
  if (length(fix_pars) == 0) {
    return(brms::empty_prior())
  }

  # construct parameter names and prior values
  par_list <- c(model$fixed_parameters, additional_pars)
  pars <- names(par_list)
  values <- unlist(par_list)
  priors <- glue("constant({values})")

  # determine type of parameters
  bterms <- brms::brmsterms(formula)
  dpars <- names(bterms$dpars)
  nlpars <- names(bterms$nlpars)

  # internal consistency check: a fixed parameter that configure_model never
  # wires into the formula would collapse to a malformed b_Intercept prior
  missing_pars <- pars[!pars %in% c(dpars, nlpars)]
  stopif(
    length(missing_pars) > 0,
    "Fixed parameter(s) {collapse_comma(missing_pars)} are not part of the model \\
    formula (neither a distributional nor a non-linear parameter). This is a \\
    model-definition error: configure_model() must wire every fixed parameter \\
    into the formula."
  )

  # flexibly set the variables for set_prior
  classes <- ifelse(pars %in% dpars, "Intercept", "b")
  coefs <- ifelse(pars %in% dpars, "", "Intercept")
  dpars <- ifelse(pars %in% dpars, pars, "")
  nlpars <- ifelse(pars %in% nlpars, pars, "")
  brms::set_prior(priors, class = classes, coef = coefs, dpar = dpars, nlpar = nlpars)
}

#' Set default priors for a bmmodel
#'
#' Specify default priors flexibly regardless of the formula the user has supplied.
#' The function will automatically recognize when intercepts are present or suppressed,
#' and will set the default priors accordingly. You can specify priors of intercepts/main
#' levels of a predictor with suppressed intercept, and priors on the effects of the
#' predictors relative to the intercept.
#'
#' @param model A `bmmodel` object
#' @param formula A `brmsformula` object
#' @param data A data.frame containing the data used in the model
#'
#' @noRd
#' @keywords internal developer
set_default_prior <- function(model, data, formula) {
  if (isFALSE(getOption("bmm.default_priors", TRUE))) {
    return(NULL)
  }

  default_priors <- validate_default_priors(model, formula)
  bterms <- brms::brmsterms(formula)
  pars <- intersect(lhs_vars(bterms), names(default_priors))

  priors <- lapply(pars, function(par) {
    construct_default_priors_list(par, bterms, default_priors, data)
  })
  priors <- unnest_list(priors)
  Reduce(combine_prior, priors, init = brms::empty_prior())
}

#' Removes default priors for parameters that are predicted by a non-linear formula.
#' @param model A `bmmodel` object
#' @param formula A `brmsformula` object
#' @return A list of valid default priors
#' @noRd
validate_default_priors <- function(model, formula) {
  default_priors <- model$default_priors
  stopif(
    !is.list(default_priors) || !all(sapply(default_priors, is.list)),
    "The default_priors should be a list of lists"
  )
  dpars_predicted_by_nlpars <- names(which(sapply(formula$pforms, is_nl)))
  default_priors[dpars_predicted_by_nlpars] <- NULL
  warnif(
    any(dpars_predicted_by_nlpars %in% names(model$parameter)),
    "Non-linear transformations of model parameters detected in the formula.
    Consider specifying priors for better estimation; otherwise, flat priors will be used."
  )
  default_priors
}

#' Determine priors for a single model parameter based on formula structure
#'
#' This function evaluates the model terms and determines the appropriate priors based on
#' whether an intercept is present, whether predictors are involved in interactions, and
#' whether priors should be assigned to levels of categorical predictors.
#'
#' @param par The model parameter name
#' @param bterms The parsed `brms` model terms
#' @param default_priors The prior descriptions for the parameter
#' @param data The data frame used in the model
#' @return A list of prior objects
#' @noRd
construct_default_priors_list <- function(par, bterms, default_priors, data) {
  bterms$allpars <- c(bterms$dpars, bterms$nlpars)
  terms <- stats::terms(bterms$allpars[[par]]$fe)
  prior_desc <- default_priors[[par]]
  has_effects_prior <- !is.null(prior_desc$effects)
  fixed_effects_count <- sum(attr(terms, "order") == 1)
  interactions_count <- sum(attr(terms, "order") > 1)
  interaction_only <- fixed_effects_count == 0 && interactions_count > 0

  priors <- list()

  # priors on fixed effects
  if (has_effects_prior && fixed_effects_count > 0) {
    fixed_effects_prior <- .build_prior(prior_desc$effects, "b", par = par, bterms = bterms)
    priors <- c(priors, list(fixed_effects_prior))
  }

  # priors on intercept; unfortunately too convoluted to use the create_prior function
  if (has_intercept(terms)) {
    intercept_prior <- .build_prior(prior_desc$main, "Intercept", par = par, bterms = bterms)
    priors <- c(priors, list(intercept_prior))
    return(priors)
  }

  # priors when intercept is supressed and all levels are explicit
  if ((fixed_effects_count == 1 && interactions_count == 0) || interaction_only) {
    levels_only_prior <- .build_prior(prior_desc[[1]], "b", par = par, bterms = bterms)
    priors <- c(priors, list(levels_only_prior))
    return(priors)
  }

  # edge case: with multiple predictors and no intercept, set the main prior on
  # the levels of the first predictor
  first_predictor_coefs <- paste0(rhs_vars(terms)[1], levels(data[[rhs_vars(terms)[1]]]))
  for (coef in first_predictor_coefs) {
    first_predictor_prior <- .build_prior(prior_desc[[1]], "b", par, coef = coef, bterms = bterms)
    priors <- c(priors, list(first_predictor_prior))
  }
  priors
}

# Helper function to create a prior object conditional on parameter type
.build_prior <- function(prior_desc, class, par, bterms, ...) {
  args <- c(list(prior = prior_desc, class = class), list(...))
  if (par %in% names(bterms$nlpars)) {
    if (class == "Intercept") {
      args$class <- "b" # Intercept priors in brms are stored in `b` for nlpars
      args$coef <- "Intercept"
    }
    args$nlpar <- par
  } else {
    args$dpar <- par
  }
  do.call(brms::prior_, args)
}

# internal function to combine two priors (e.g. the default prior with the user
# given prior) parts present in prior2 will overwrite the corresponding parts in
# prior1
combine_prior <- function(prior1, prior2) {
  if (is.null(prior2)) {
    return(prior1)
  }

  cols <- c("class", "dpar", "nlpar", "coef", "group", "resp")
  prior1_types <- do.call(paste, prior1[, cols])
  prior2_types <- do.call(paste, prior2[, cols])
  is_duplicate <- prior1_types %in% prior2_types
  prior <- prior1[!is_duplicate, ] + prior2
  row.names(prior) <- 1:nrow(prior)
  prior
}

summarise_default_prior <- function(prior_list) {
  pars <- names(prior_list)
  prior_info <- ""
  for (par in pars) {
    prior_info <- paste0(prior_info, "   - `", par, "`:\n")
    types <- names(prior_list[[par]])
    for (type in types) {
      prior <- prior_list[[par]][[type]]
      prior_info <- paste0(prior_info, "      - `", type, "`: ", prior, "\n")
    }
  }
  prior_info
}

constrain_set_size1_fixef <- function(formula, nlpars, set_size_var, prior_value) {
  nl_pforms <- formula$pforms[nlpars]
  has_setsize <- vapply(nl_pforms, function(x) set_size_var %in% rhs_vars(x), logical(1))
  
  if (!any(has_setsize)) {
    return(NULL)
  }
    
  brms::prior_(prior_value,
    class = "b",
    coef = paste0(set_size_var, 1),
    nlpar = nlpars[has_setsize]
  )
}

constrain_set_size1_ranef <- function(formula, nlpars, set_size_var, prior_value) {
  re_terms_list <- lapply(nlpars, function(nlp) {
    bterms <- brms::brmsterms(formula$pforms[[nlp]])
    re <- bterms$dpars$mu$re
    if (NROW(re) > 0) {
      data.frame(nlpar = nlp, group = re$group, form = I(re$form), stringsAsFactors = FALSE)
    }
  })
  
  re_terms_df <- do.call(rbind, re_terms_list)
  if (is.null(re_terms_df)) {
    return(NULL)
  }
  
  has_setsize <- vapply(re_terms_df$form, function(x) set_size_var %in% rhs_vars(x), logical(1))
  
  if (!any(has_setsize)) {
    return(NULL)
  }

  brms::prior_(prior_value,
    class = "sd",
    coef = paste0(set_size_var, 1),
    group = re_terms_df$group[has_setsize],
    nlpar = re_terms_df$nlpar[has_setsize]
  )
}
