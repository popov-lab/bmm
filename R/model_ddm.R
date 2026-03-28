#############################################################################!
# MODELS                                                                 ####
#############################################################################!

.ddm_defaults <- list(
  parameters = list(
    drift = "Drift rate = Average rate of evidence accumulation of the decision processes",
    bound = "Boundary separation = Distance between the decision boundaries that need to be reached",
    ndt   = "Non-decision time = Additional time required beyond the evidence accumulation process",
    zr    = "Relative starting point = Starting point between the decision thresholds relative to the upper bound."
  ),
  links = list(
    drift = "identity", bound = "log", ndt = "log", zr = "logit"
  ),
  fixed_parameters = list(
    zr = 0,
    mu = 0
  ),
  priors = list(
    drift = list(main = "cauchy(0,1)", effects = "normal(0,0.5)"),
    bound = list(main = "normal(0,0.5)", effects = "normal(0,0.5)"),
    ndt   = list(main = "normal(-1.5,0.5)", effects = "normal(0,0.3)"),
    zr    = list(main = "normal(0,0.5)", effects = "normal(0,0.3)")
  ),
  init_ranges = list(
    mu = c(-0.1,0.1),
    drift  = c(-1,1),
    bound  = c(1,2),
    ndt    = c(0.01,0.05),
    zr     = c(0.45,0.55)
  )
)

.model_ddm <- function(rt = NULL, response = NULL, links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = nlist(),
      domain = "Decision Making / Processing Speed",
      task = "Two-Alternative Force Choice RT",
      name = "Diffusion Decision Model",
      version = "NA",
      citation = glue(
        "Ratcliff, R. (1978). A theory of memory retrieval. Psychological Review, 85(2), 59-108. https://doi.org/10/fjwm2f;"
      ),
      requirements = glue(
        "- The response time should be in seconds and \\
          represent the time between onset of the target stimulus until the response execution
          - The response should be coded numerically: \\
          0 = lower response, 1 = upper response"
      ),
      parameters = .ddm_defaults[["parameters"]],
      links = .ddm_defaults[["links"]],
      fixed_parameters = .ddm_defaults[["fixed_parameters"]],
      default_priors = .ddm_defaults[["priors"]],
      init_ranges = .ddm_defaults[["init_ranges"]],
      void_mu = TRUE
    ),
    class = c("bmmodel", "ddm"),
    call = call
  )
  out$links[names(links)] <- links
  out
}
# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_ddm()$info

#' @title `r .model_ddm()$name`
#' @name ddm
#' @details `r model_info(.model_ddm())`
#' @param rt Name of the reaction time variable coding reaction time in seconds in the data.
#' @param response Name of the response variable coding the response numerically (0 = lower response / incorrect, 1 = upper response / correct)
#' @param links A list of links for the parameters.
#' @section Default behavior:
#' By default, `zr` is fixed at 0. If you want to estimate `zr`, add a formula
#' for `zr` in your `bmf()` call.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # Simulate data for one subject using rddm
#' set.seed(123)
#' n_trials <- 500
#' 
#' # Simulate DDM data with fixed parameters
#' sim_data <- rddm(
#'   n = n_trials,
#'   drift = 1.5,    # drift rate
#'   bound = 1.2,    # boundary separation
#'   ndt = 0.3,      # non-decision time
#'   zr = 0.5        # relative starting point
#' )
#' 
#' # Prepare data frame
#' dat <- data.frame(
#'   rt = sim_data$rt,
#'   response = sim_data$response
#' )
#' 
#' # Define formula (intercept-only model)
#' ff <- bmmformula(
#'   drift ~ 1,
#'   bound ~ 1,
#'   ndt ~ 1
#' )
#' 
#' # Specify the DDM model
#' model <- ddm(rt = "rt", response = "response")
#' 
#' # Fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   iter = 1000,
#'   backend = "cmdstanr"
#' )
#' 
#' # Check parameter recovery
#' summary(fit)
ddm <- function(rt, response, links = NULL, ...) {
  stopif(
    !requireNamespace("cmdstanr", quietly = TRUE),
    'The "cmdstanr" package is required for this functionality'
  )

  call <- match.call()
  stop_missing_args()
  .model_ddm(rt = rt, response = response,
             links = links, call = call, ...)
}

#############################################################################!
# CHECK_DATA S3 methods                                                  ####
#############################################################################!

#' @export
check_data.ddm <- function(model, data, formula) {
  rt_var <- model$resp_vars$rt
  response_var <- model$resp_vars$response

  stopif(
    not_in(rt_var, colnames(data)),
    "The RT variable '{rt_var}' is not present in the data."
  )

  stopif(
    not_in(response_var, colnames(data)),
    "The response variable '{response_var}' is not present in the data."
  )

  if (any(is.na(data[, rt_var]))) {
    data <- data[!is.na(data[, rt_var]), ]
    warning2("Some values in {rt_var} were NA. These were removed from the analysis.")
  }

  if (any(is.na(data[, response_var]))) {
    data <- data[!is.na(data[, response_var]), ]
    warning2("Some values in {response_var} were NA. These were removed from the analysis.")
  }

  if (typeof(data[, rt_var]) %in% c("double", "integer")) {
    stopif(
      any(data[, rt_var] < 0, na.rm = TRUE),
      "Some reaction times are lower than zero, please check your data."
    )

    warnif(
      any(data[, rt_var] > 10, na.rm = TRUE),
      "Your data contains reaction times larger than 10 seconds.\n
      Either you have passed reaction times in milliseconds, then please \\
      recode them to seconds and rerun the model.\n
      Or you have very long RTs in your data in which case you might want \\
      to consider outlier filtering."
    )

    warnif(
      any(data[, rt_var] < 0.100, na.rm = TRUE),
      "Your data contains reaction times smaller than 0.100 seconds.\n
      It is likely that the model will not be able to sample with the \\
      current settings of the initial values.\n
      Either pass your own initial value function or consider filtering \\
      reaction times below 0.100 seconds."
    )
  } else {
    stop2("The RT variable '{rt_var}' needs to be of type double or integer.")
  }

  resp_type <- typeof(data[, response_var])

  if (resp_type %in% c("integer", "double", "numerical")) {
    stopif(
      any(!data[, response_var] %in% c(0, 1), na.rm = TRUE),
      "The response variable '{response_var}' should only contain values of 0 and 1."
    )
  } else if (resp_type == "logical") {
    warning2(
      "The response variable is boolean and will be internally transformed ",
      "to an integer variable with values 0 for FALSE and 1 for TRUE."
    )
    data[, response_var] <- ifelse(data[, response_var], 1, 0)
  } else if (resp_type == "character") {
    data[, response_var] <- tolower(data[, response_var])
    stopif(
      any(!data[, response_var] %in% c("upper", "lower")),
      "The response variable '{response_var}' contains invalid character values.\n
      Please pass only 'upper' or 'lower' as response values, or use \\
      numeric coding (0 = lower, 1 = upper)."
    )
    warning2(
      "The response variable is a character variable and will be internally ",
      "transformed to an integer variable with 0 for 'lower' and 1 for 'upper'."
    )
    data[, response_var] <- ifelse(data[, response_var] == "upper", 1, 0)
  } else {
    stop2(
      "The response variable '{response_var}' is of type '{resp_type}'.\n
      Please provide responses as integer (0/1), logical, or character \\
      ('upper'/'lower')."
    )
  }

  stopif(
    any(!data[, response_var] %in% c(0, 1)),
    "Invalid values in the response variable '{response_var}'.\n
    After processing, responses must be coded as 0 (lower) or 1 (upper)."
  )

  NextMethod("check_data")
}

#############################################################################!
# Convert bmmformula to brmsformla methods                               ####
#############################################################################!

#' @export
bmf2bf.ddm <- function(model, formula) {
  rt <- model$resp_vars$rt
  response <- model$resp_vars$response
  brms::bf(paste0(rt, " | dec(", response, ") ~ 1"))
}

#############################################################################!
# CONFIGURE_MODEL S3 METHODS                                             ####
#############################################################################!

#' @export
configure_model.ddm <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  ddm_family <- function(link_drift, link_bound, link_ndt, link_zr) {
    brms::custom_family(
      "ddm",
      dpars = c("mu", "drift", "bound", "ndt", "zr"),
      links = c("identity", link_drift, link_bound, link_ndt, link_zr),
      lb = c(NA, NA, 0.1, 0, 0),
      ub = c(NA, NA, NA, NA, 1),
      type = "real",
      vars = "dec[n]",
      loop = TRUE,
      log_lik = log_lik_ddm,
      posterior_predict = posterior_predict_ddm
    )
  }

  formula$family <- ddm_family(
    link_drift = model$links$drift,
    link_bound = model$links$bound,
    link_ndt = model$links$ndt,
    link_zr = model$links$zr
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_functions <- read_lines2(paste0(sc_path, "/ddm_functions.stan"))
  stanvars <- brms::stanvar(scode = stan_functions, block = "functions")

  nlist(formula, data, stanvars)
}

log_lik_ddm <- function(i, prep) {
  dddm(
    rt = prep$data$Y[i],
    response = prep$data[["dec"]][i],
    drift = brms::get_dpar(prep, "drift", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    zr = brms::get_dpar(prep, "zr", i = i),
    log = TRUE
  )
}

posterior_predict_ddm <- function(i, prep, ...) {
  out <- rddm(
    n = length(brms::get_dpar(prep, "drift", i = i)),
    drift = brms::get_dpar(prep, "drift", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    zr = brms::get_dpar(prep, "zr", i = i)
  )

  dots <- list(...)
  if (!is.null(dots$negative_rt) && dots$negative_rt) {
    out[["rt"]] * ifelse(out[["response"]] == 1, 1, -1)
  } else {
    out[["rt"]]
  }
}
