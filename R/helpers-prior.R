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

  # Skip fixed parameters the user has overridden in the formula
  # (e.g., sdratio ~ 1 overrides fixed_parameters = list(sdratio = 0))
  # Also skip fixed parameters that configure_model stripped from the formula
  # (e.g., sdratio removed for sdt_ranking_ev family).
  bterms <- brms::brmsterms(formula)
  active_pars <- c(names(bterms$dpars), names(bterms$nlpars))
  fix_pars <- fix_pars[names(fix_pars) %in% active_pars]
  if (length(fix_pars) == 0) return(brms::empty_prior())

  overridden <- intersect(names(fix_pars), active_pars)
  overridden <- Filter(function(p) {
    !isTRUE(attr(formula$pforms[[p]], "constant"))
  }, overridden)
  fix_pars <- fix_pars[!names(fix_pars) %in% overridden]
  if (length(fix_pars) == 0) {
    return(brms::empty_prior())
  }

  # construct parameter names and prior values
  par_list <- c(fix_pars, additional_pars)
  pars <- names(par_list)
  values <- unlist(par_list)
  priors <- glue("constant({values})")

  # determine type of parameters (reuse bterms from above)
  dpars <- names(bterms$dpars)
  nlpars <- names(bterms$nlpars)

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
