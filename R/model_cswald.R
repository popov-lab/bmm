#############################################################################!
# MODELS                                                                 ####
#############################################################################!

.cswald_version_table <- list(
  simple = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary (distance from starting point to correct boundary)",
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
      bound = "boundary separation (total distance between boundaries)",
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
      drift = list(main = "normal(0,1)", effects = "normal(0,0.5)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)"),
      zr = list(main = "normal(0,0.3)", effects = "normal(0,0.2)"),
      s = list(main = "normal(0,0.5)", effects = "normal(0,0.2)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(-0.5, 0.5),
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
      requirements = glue(
        "- Reaction times should be passed in seconds","\n",
        "- The response variable should be passed numerically: 0 = lower response, 1 = upper response"
      ),
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
#' @param rt The name of the variable in the dataset containing the response
#'   times. Response times should be coded in seconds (not milliseconds).
#' @param response The name of the variable in the dataset containing the
#'   response/decision. Responses should be coded as 0 (lower boundary) or
#'   1 (upper boundary). Alternatively, character values "lower" and "upper"
#'   or logical values (FALSE/TRUE) are accepted and will be converted
#'   automatically.
#' @param links A named list of link functions for the model parameters.
#'   Available parameters depend on the version: "simple" has `drift`, `bound`,
#'   `ndt`, and `s`; "crisk" additionally has `zr`. Default links are "log" for
#'   most parameters and "logit" for `zr`.
#' @param version A character string specifying which version of the cswald
#'   model to use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): The standard censored shifted Wald model,
#'       which treats error responses as censored correct responses. Best suited
#'       for tasks with few errors (<20%). **Note:** The `bound` parameter in
#'       the simple version represents the distance from the starting point to
#'       the correct boundary, which is half the total boundary separation
#'       in the diffusion model (assuming an unbiased starting point). To
#'       convert to the full boundary separation (as in DDM or crisk), multiply
#'       by 2.
#'     \item `"crisk"`: The competing risks version, which models both response
#'       types as arising from racing accumulators toward opposite boundaries.
#'       Better suited for tasks with substantial error rates. The `bound`
#'       parameter represents the total boundary separation, consistent with
#'       the diffusion model parameterization.
#'   }
#'   For more details, see Miller et al. (2017).
#' @param ... Additional arguments passed internally (for testing purposes).
#' @return An object of class `bmmodel`
#' @export
#' @keywords bmmodel
#' @seealso [dcswald()] and [rcswald()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # generate simulated data from the diffusion model
#' dat <- rcswald(n = 500, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5, s = 1)
#'
#' # specify the model
#' model <- cswald(rt = "rt", response = "response", version = "simple")
#'
#' # specify the formula
#' formula <- bmf(
#'   drift ~ 1,
#'   bound ~ 1,
#'   ndt ~ 1
#' )
#'
#' # fit the model
#' fit <- bmm(
#'   formula = formula,
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
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

  # check that required variables are present
 stopif(
    not_in(rt_var, colnames(data)),
    "The RT variable '{rt_var}' is not present in the data."
  )

  stopif(
    not_in(response_var, colnames(data)),
    "The response variable '{response_var}' is not present in the data."
  )

  # check for NA values
  n_na_rt <- sum(is.na(data[, rt_var]))
  n_na_resp <- sum(is.na(data[, response_var]))

  stopif(
    n_na_rt > 0,
    "The RT variable '{rt_var}' contains {n_na_rt} NA values. \\
    Please remove or impute missing values before fitting the model."
  )

  stopif(
    n_na_resp > 0,
    "The response variable '{response_var}' contains {n_na_resp} NA values. \\
    Please remove or impute missing values before fitting the model."
  )

  # checks for rt_var
  if (typeof(data[, rt_var]) %in% c("double", "integer")) {
    stopif(
      any(data[, rt_var] < 0),
      "Some reaction times are lower than zero, please check your data."
    )

    warnif(
      any(data[, rt_var] > 10),
      "Your data contains reaction times larger than 10 seconds.\n
      Either you have passed reaction times in milliseconds, then please \\
      recode them to seconds and rerun the model.\n
      Or you have very long RTs in your data in which case you might want \\
      to consider outlier filtering."
    )

    warnif(
      any(data[, rt_var] < 0.100),
      "Your data contains reaction times smaller than 0.100 seconds.\n
      It is likely that the model will not be able to sample with the \\
      current settings of the initial values.\n
      Either pass your own initial value function or consider filtering \\
      reaction times below 0.100 seconds."
    )
  } else {
    stop(glue(
      "The RT variable '{rt_var}' needs to be of type double or integer."
    ))
  }

  # checks for response_var
  resp_type <- typeof(data[, response_var])

  # handle factor variables by converting to character first
  if (is.factor(data[, response_var])) {
    data[, response_var] <- as.character(data[, response_var])
    resp_type <- "character"
  }

  if (resp_type %in% c("integer", "double")) {
    stopif(
      any(!data[, response_var] %in% c(0, 1)),
      "The response variable '{response_var}' contains values other than 0 and 1.\n
      Please pass responses coded as 0 (lower boundary) and 1 (upper boundary)."
    )
  } else if (resp_type == "logical") {
    warning(glue::glue(
      "The response variable is boolean and will be internally transformed ",
      "to an integer variable with values 0 for FALSE and 1 for TRUE."
    ))
    data[, response_var] <- as.integer(data[, response_var])
  } else if (resp_type == "character") {
    data[, response_var] <- tolower(data[, response_var])
    stopif(
      any(!data[, response_var] %in% c("upper", "lower")),
      "The response variable '{response_var}' contains invalid character values.\n
      Please pass only 'upper' or 'lower' as response values, or use \\
      numeric coding (0 = lower, 1 = upper)."
    )
    warning(glue::glue(
      "The response variable is a character variable and will be internally ",
      "transformed to an integer variable with 0 for 'lower' and 1 for 'upper'."
    ))
    data[, response_var] <- ifelse(data[, response_var] == "upper", 1L, 0L)
  } else {
    stop(glue(
      "The response variable '{response_var}' is of type '{resp_type}'.\n
      Please provide responses as integer (0/1), logical, character \\
      ('upper'/'lower'), or factor."
    ))
  }

  # final validation of response values
  stopif(
    any(!data[, response_var] %in% c(0, 1)),
    "Invalid values in the response variable '{response_var}'.\n
    After processing, responses must be coded as 0 (lower) or 1 (upper)."
  )

  # warn about high error rates for simple version (assumes few errors)
  if ("cswald_simple" %in% class(model)) {
    error_rate <- mean(data[, response_var] == 0)
    warnif(
      error_rate > 0.20,
      "Your data has an error rate of {round(error_rate * 100, 1)}%.\n
      The simple censored shifted Wald model assumes few errors. \\
      Consider using version = 'crisk' (competing risks) for data with \\
      substantial error rates."
    )
  }

  # warn about small sample sizes
  n_obs <- nrow(data)
  warnif(
    n_obs < 50,
    "Your data contains only {n_obs} observations.\n
    Parameter estimation may be unreliable with very small sample sizes."
  )

  NextMethod("check_data")
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
  response_var <- model$other_vars$response
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
      log_lik = log_lik_cswald_simple,
      posterior_predict = posterior_predict_cswald_simple
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

  # extract parameters from posterior draws

  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # Simple version estimates distance to upper bound only (single-boundary Wald).
  # To generate data from the full DDM (via rcswald), convert to total boundary
  # separation: a = 2 * bound, with unbiased starting point zr = 0.5
  out <- rcswald(
    n = length(drift),
    drift = drift,
    bound = bound * 2,
    ndt = ndt,
    zr = 0.5,
    s = s
  )

  # handle negative_rt option (consistent with brms::wiener)
  dots <- list(...)
  if (!is.null(dots$negative_rt) && dots$negative_rt) {
    # encode lower bound responses as negative RTs
    out$rt * ifelse(out$response == 1, 1, -1)
  } else {
    out$rt
  }
}

log_lik_cswald_simple <- function(i, prep) {
  # extract parameters from posterior draws
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # get observed data for observation i
  rt <- prep$data$Y[i]
  response <- prep$data$dec[i]

  # number of posterior draws
  ndraws <- length(drift)

  # replicate rt and response to match number of posterior draws
  rt <- rep(rt, ndraws)
  response <- rep(response, ndraws)

  # compute log-likelihood for each posterior draw
  dcswald(rt, response, drift, bound, ndt, s = s, version = "simple", log = TRUE)
}

#' @export
configure_model.cswald_crisk <- function(model, data, formula) {
  # retrieve required arguments
  response_var <- model$other_vars$response
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
      lb = c(NA, NA, 0, 0, 0, 0), # lower bounds (drift can be negative)
      type = 'real', # real for continous dv, int for discrete dv
      vars = 'dec[n]',
      loop = TRUE, # is the likelihood vectorized
      log_lik = log_lik_cswald_crisk,
      posterior_predict = posterior_predict_cswald_crisk
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
  # extract parameters from posterior draws
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # get observed data for observation i
  rt <- prep$data$Y[i]
  response <- prep$data$dec[i]

  # number of posterior draws
  ndraws <- length(drift)

  # replicate rt and response to match number of posterior draws
  rt <- rep(rt, ndraws)
  response <- rep(response, ndraws)

  # compute log-likelihood for each posterior draw
  dcswald(rt, response, drift, bound, ndt, zr = zr, s = s, version = "crisk", log = TRUE)
}

posterior_predict_cswald_crisk <- function(i, prep, ...) {
  # extract parameters from posterior draws
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # Crisk version estimates full boundary separation and relative starting point
  out <- rcswald(
    n = length(drift),
    drift = drift,
    bound = bound,
    ndt = ndt,
    zr = zr,
    s = s
  )

  # handle negative_rt option (consistent with brms::wiener)
  dots <- list(...)
  if (!is.null(dots$negative_rt) && dots$negative_rt) {
    # encode lower bound responses as negative RTs
    out$rt * ifelse(out$response == 1, 1, -1)
  } else {
    out$rt
  }
}
