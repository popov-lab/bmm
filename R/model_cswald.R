#############################################################################!
# MODELS                                                                 ####
#############################################################################!

.cswald_version_table <- list(
  simple = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary seperation",
      ndt = "non-decision time",
      s = "diffusion constant"
    ),
    links = list(
      drift = "log",
      bound = "log",
      ndt = "log",
      s = "log"
    ),
    fixed_parameters = list(
      s = 0
    ),
    priors = list(
      drift = list(main = "normal(0,1)", effects = "normal(0,0.3)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)"),
      s = list(main = "normal(0,0.3)", effects = "normal(0,0.2)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(1,2),
      bound = c(1.5,2),
      ndt = c(0.025, 0.05),
      s = c(0.95,1.05)
    )
  ),
  crisk = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary seperation",
      ndt = "non-decision time",
      zr = "relative starting point",
      s = "diffusion constant"
    ),
    links = list(
      drift = "identity",
      bound = "log",
      ndt = "log",
      zr = "logit",
      s = "log"
    ),
    fixed_parameters = list(
      zr = 0,
      s = 0
    ),
    priors = list(
      drift = list(main = "normal(1,1)", effects = "normal(0,0.3)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)"),
      zr = list(main = "normal(0,0.3)", effects = "normal(0,0.2)"),
      s = list(main = "normal(0,0.5)", effects = "normal(0,0.2)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(1,2),
      bound = c(1.5,2),
      ndt = c(0.025, 0.05),
      zr = c(0.45, 0.55),
      s = c(0.95,1.05)
    )
  )
)


.model_cswald <- function(
    rt = NULL,
    response = NULL,
    links = NULL,
    version = "simple",
    call = NULL,
    ...
) {
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = nlist(),
      domain = "Processing Speed, Decision Making",
      task = "Choice Reaction Time tasks (with few errors)",
      name = "Censored-Shifted Wald Model",
      citation = "Miller, R., Scherbaum, S., Heck, D. W., Goschke, T., & Enge, S. (2017).
        On the Relation Between the (Censored) Shifted Wald and the Wiener Distribution as Measurement Models
        for Choice Response Times. Applied Psychological Measurement, 42(2), 116-135. https://doi.org/10.1177/0146621617710465",
      version = version,
      requirements = "",
      parameters = .cswald_version_table[[version]][["parameters"]],
      links = .cswald_version_table[[version]][["links"]],
      fixed_parameters = .cswald_version_table[[version]][["fixed_parameters"]],
      default_priors = .cswald_version_table[[version]][["priors"]],
      init_ranges = .cswald_version_table[[version]][["init_ranges"]],
      void_mu = TRUE
    ),
    class = c("bmmodel", "cswald", paste0("cswald_", version)),
    call = call
  )

  if (!is.null(version)) {
    class(out) <- c(class(out))
  }

  out$links[names(links)] <- links
  out
}

# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_cswald()$info

#' @title `r .model_cswald()$name`
#' @name cswald
#' @details `r model_info(.model_cswald())`
#' @param rt A description of the response variable
#' @param response A description of the response variable
#' @param links A list of links for the parameters.
#' @param version A character string which version of the cswald model to use.
#'   The standard cswald model is referred to as "simple" version, and the competing risk
#'   version is referred to as "crisk". For more details on these model versions, please look
#'   into Miller et al. (2017) given as the reference in the `cswald` model object.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @export
#' @examples
#' \dontrun{
#' # put a full example here (see 'R/model_mixture3p.R' for an example)
#' }
cswald <- function(rt, response, links = NULL, version = "simple", ...) {
  call <- match.call()
  stop_missing_args()
  .model_cswald(
    rt = rt,
    response = response,
    links = links,
    version = version,
    call = call,
    ...
  )
}

#############################################################################!
# CHECK_DATA S3 methods                                                  ####
#############################################################################!
# A check_data.* function should be defined for each class of the model.
# If a model shares methods with other models, the shared methods should be
# defined in helpers-data.R. Put here only the methods that are specific to
# the model. See ?check_data for details.
# (YOU CAN DELETE THIS SECTION IF YOU DO NOT REQUIRE ADDITIONAL DATA CHECKS)

#' @export
check_data.cswald <- function(model, data, formula) {
  # retrieve required arguments
  rt_var <- model$resp_vars$rt
  response_var <- model$resp_vars$response

  # stop due to missing information
  stopif(
    not_in(rt_var, colnames(data)),
    "The response variable '{rt_var}' is not present in the data."
  )

  stopif(
    not_in(response_var, colnames(data)),
    "The response variable '{response_var}' is not present in the data."
  )

  # checks for rt_var
  if (typeof(data[, rt_var]) %in% c("double", "numerical")) {
    stopif(
      any(data[, rt_var] < 0),
      glue("Some reaction times are lower than zero, please check your data.")
    )

    warnif(any(data[,rt_var] > 10),
           glue::glue("Your data contains reaction times larger then 10.\n",
                      "Either you have passed reaction times in milli-seconds, then please recode them to seconds and rerun the model.\n",
                      "Or you have very long RTs in your data in which case you might want to consider outlier deletion."))

    warnif(any(data[,rt_var] < .100),
           glue::glue("Your data contains reaction times smaller the 0.100 seconds.\n",
                      "It is likely that the model will not be able to sample with the current settings of the inital values.\n",
                      "Either pass your own initial value function or consider filtering reaction times below 0.100 seconds"))
  } else {
    stop(glue(
      "The rt variable: ",
      rt_var,
      " needs to be of type double or numerical."
    ))
  }

  # checks for response_var
  if (typeof(data[, response_var]) %in% c("integer", "double", "numerical")) {
    stopif(any(!data[, response_var] %in% c(0, 1)), glue("Only values of zero"))
  } else if (typeof(data[, response_var]) == "logical") {
    warning(glue::glue(
      "The response variable you provided is boolean, it will be internally transformed ",
      "to an integer variable with values 0 for FALSE and 1 for TRUE."
    ))
    data[, response_var] <- ifelse(data[, response_var], 1, 0)
  } else if (typeof(data[, response_var]) == "character") {
    data[,response_var] <- tolower(data[,response_var])
    stopif(any(!data[,response_var] %in% c("upper","lower")),
           glue::glue("You have passed a character variable as response variable containing invalid responses.\n",
                      "Please pass only upper or lower responses as response variables either coded numerically (0 = \"lower\" and 1 = \"upper\")\n",
                      "or characters that match \"upper\" and \"lower\"."))
    warning(glue::glue(
      "The response variable you provided is a character variable, it will be internally transformed ",
      "to an integer variable with values 0 for \"lower\" and 1 for \"upper\"."
    ))
    data[,response_var] <- ifelse(data[,response_var] == "upper",1,0)
  }else {
    stop(glue(
      "The response variable: ",
      response_var,
      " needs to be of type integer, numerical, or logical."
    ))
  }

  stopif(any(!data[,response_var] %in% c(0,1)),
         glue::glue("Invalid values in the response variable passed to ddm.\n",
                    "Please pass either numeric or character variables that only contain 0 and 1 or  \"upper\" and \"lower\"."))

  data
}

#############################################################################!
# Convert bmmformula to brmsformla methods                               ####
#############################################################################!
# A bmf2bf.* function should be defined if the default method for constructing
# the brmsformula from the bmmformula does not apply (e.g if aterms are required).
# The shared method for all `bmmodels` is defined in bmmformula.R.
# See ?bmf2bf for details.
# (YOU CAN DELETE THIS SECTION IF YOUR MODEL USES A STANDARD FORMULA WITH 1 RESPONSE VARIABLE)

#' @export
bmf2bf.cswald <- function(model, formula) {
  # retrieve required response arguments
  rt_var <- model$resp_vars$rt
  response_var <- model$resp_vars$response

  # set the base brmsformula based
  brms_formula <- brms::bf(glue(rt_var, " | dec(", response_var, ") ~ 1"))

  # return the brms_formula to add the remaining bmmformulas to it.
  brms_formula
}

#############################################################################!
# CONFIGURE_MODEL S3 METHODS                                             ####
#############################################################################!
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.cswald_simple <- function(model, data, formula) {
  # retrieve required arguments
  rt_var <- model$other_vars$rt
  response_var <- model$other_vars$response
  trials <- data[, response_var]
  links <- model$links

  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)

  # construct the family & add to formula object
  cswald_family <- function (link_drift, link_bound, link_ndt, link_s) {
    brms::custom_family(
      'cswald',
      dpars = c("mu", "drift", "bound", "ndt", "s"),
      links = c("identity", link_drift, link_bound, link_ndt, link_s),
      ub = c(NA, NA, NA, NA, NA), # upper bounds for parameters
      lb = c(NA, 0, 0, 0, 0), # lower bounds for parameters
      type = 'real', # real for continous dv, int for discrete dv
      vars = 'dec[n]',
      loop = TRUE, # is the likelihood vectorized
    )
  }

  formula$family <- cswald_family(
    link_drift = links$drift,
    link_bound = links$bound,
    link_ndt = links$ndt,
    link_s = links$s
  )

  # prepare initial stanvars to pass to brms, model formula and priors
  sc_path <- system.file('stan_chunks', package = 'bmm')
  stan_helpers <- read_lines2(paste0(sc_path, "/cswald_helper_functions.stan"))
  stan_functions <- read_lines2(paste0(sc_path, '/cswald_simple_functions.stan'))

  stanvars <- brms::stanvar(scode = stan_helpers, block = 'functions') +
    brms::stanvar(scode = stan_functions, block = 'functions')

  # return the list
  nlist(formula, data, stanvars)
}

posterior_predict_cswald_simple <- function(i, prep, ...) {

}

log_lik_cswald_simple <- function(i, prep) {

}

#' @export
configure_model.cswald_crisk <- function(model, data, formula) {
  env_cswald_crisk <- new.env(parent = baseenv())

  # retrieve required arguments
  rt_var <- model$other_vars$rt
  response_var <- model$other_vars$response
  trials <- data[, response_var]
  links <- model$links

  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)

  # construct the family & add to formula object
  cswald_crisk_family <- function(link_drift, link_bound, link_ndt, link_zr, link_s) {
    brms::custom_family(
      'cswald_crisk',
      dpars = c("mu", "drift", "bound", "ndt", "zr", "s"),
      links = c("identity", link_drift, link_bound, link_ndt, link_zr, link_s),
      ub = c(NA, NA, NA, NA, 1, NA), # upper bounds for parameters
      lb = c(NA, 0, 0, 0, 0, 0), # lower bounds for parameters
      type = 'real', # real for continous dv, int for discrete dv
      vars = 'dec[n]',
      loop = TRUE, # is the likelihood vectorized
    )
  }
  formula$family <- cswald_crisk_family(
    link_drift = links$drift,
    link_bound = links$bound,
    link_ndt = links$ndt,
    link_zr = links$zr,
    link_s = links$s
  )

  # prepare initial stanvars to pass to brms, model formula and priors
  sc_path <- system.file('stan_chunks', package = 'bmm')
  stan_helpers <- read_lines2(paste0(sc_path, "/cswald_helper_functions.stan"))
  stan_functions <- read_lines2(paste0(sc_path, '/cswald_crisk_functions.stan'))

  stanvars <- brms::stanvar(scode = stan_helpers, block = 'functions') +
    brms::stanvar(scode = stan_functions, block = 'functions')

  # return the list
  nlist(formula, data, stanvars)
}

log_lik_cswald_crisk <- function(i, prep) {

}

posterior_predict_cswald_crisk <- function(i, prep, ...) {

}
