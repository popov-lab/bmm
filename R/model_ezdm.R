#############################################################################!
# MODELS                                                                 ####
#############################################################################!

.ezdm_version_table <- list(
  "3par" = list(
    parameters = list(
      drift = "Drift rate = Average rate of evidence accumulation of the decision processes",
      bound = "Boundary separation = Distance between the decision boundaries that need to be reached",
      ndt = "Non-decision time = Additional time required beyond the evidence accumulation process",
      s = "The diffusion constant, that is the standard deviation of the Gaussian noise during sampling"
    ),
    links = list(
      drift = "identity", bound = "log", ndt = "log", s = "log"
    ),
    fixed_parameters = list(
      s = 0
    ),
    priors = list(
      drift = list(main = "cauchy(0,1)", effects = "normal(0,0.5)"),
      bound = list(main = "normal(0,0.5)", effects = "normal(0,0.5)"),
      ndt = list(main = "normal(-1.5,0.5)", effects = "normal(0,0.3)"),
      s = list(main = "normal(0,1)", effects = "normal(0,0.3)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(-1,1),
      bound = c(1,2),
      ndt = c(0.25, 0.5),
      s = c(0.99,1.01)
    )
  ),
  "4par" = list(
    parameters = list(
      drift = "Drift rate = Average rate of evidence accumulation of the decision processes",
      bound = "Boundary separation = Distance between the decision boundaries that need to be reached",
      ndt = "Non-decision time = Additional time required beyond the evidence accumulation process",
      zr = "Relative starting point = Starting point between the decision thresholds relative to the upper bound.",
      s = "The diffusion constant, that is the standard deviation of the Gaussian noise during sampling"
    ),
    links = list(
      drift = "identity", bound = "log", ndt = "log", zr = "logit", s = "log"
    ),
    fixed_parameters = list(
      s = 0
    ),
    priors = list(
      drift = list(main = "cauchy(0,1)", effects = "normal(0,0.5)"),
      bound = list(main = "normal(0,0.5)", effects = "normal(0,0.5)"),
      ndt = list(main = "normal(-1.5,0.5)", effects = "normal(0,0.3)"),
      zr = list(main = "normal(0,0.5)", effects = "normal(0,0.3)"),
      s = list(main = "normal(0,1)", effects = "normal(0,0.3)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(-1,1),
      bound = c(1,2),
      ndt = c(0.25, 0.5),
      zr = c(0.45, 0.55),
      s = c(0.99,1.01)
    )
  )
)


.model_ezdm <- function(mean_rt = NULL, var_rt = NULL, n_upper = NULL, n_trials = NULL, version = "3par", links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(mean_rt, var_rt, n_upper),
      other_vars = nlist(n_trials),
      domain = "Processing Speed, Decision Making",
      task = "Choice Reaction Time tasks",
      name = "EZ-Diffusion Model",
      citation = glue(
        "Wagenmakers, E.-J., Van Der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. Psychonomic Bulletin & Review, 14(1), 3-22. https://doi.org/10/fk447c","\n",
        "- Ch\u00e1vez De la Pe\u00f1a, A. F., & Vandekerckhove, J. (2025). An EZ Bayesian hierarchical drift diffusion model for response time and accuracy. Psychonomic Bulletin & Review. https://doi.org/10.3758/s13423-025-02729-y"
      ),
      version = version,
      requirements = glue(
        "Provide aggregated statistics for each subject and condition that model parameters should vary over:","\n\n",
        "  - Mean reaction times (mean_rt) in seconds","\n",
        "  - Variance of reaction times (var_rt) in seconds","\n",
        "  - Number of responses to the upper decision threshold (n_upper)", "\n",
        "  - Total number of trials used to calculate aggregated statistics (n_trials)"
      ),
      parameters = .ezdm_version_table[[version]][["parameters"]],
      links = .ezdm_version_table[[version]][["links"]],
      fixed_parameters = .ezdm_version_table[[version]][["fixed_parameters"]],
      default_priors = .ezdm_version_table[[version]][["priors"]],
      init_ranges = .ezdm_version_table[[version]][["init_ranges"]],
      void_mu = TRUE
    ),
    class = c("bmmodel", "ezdm"),
    call = call
  )
  if(!is.null(version)) class(out) <- c(class(out), paste0("ezdm_",version))
  out$links[names(links)] <- links
  out
}
# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_ezdm()$info

#' @title `r .model_ezdm()$name`
#' @name ezdm
#' @details `r model_info(.model_ezdm(version = "4par"))`
#' 
#' **Note:** The drift parameter uses an identity link (not log) to allow
#' estimation of negative drift rates for below-chance performance. A soft
#' absolute value approximation ensures smooth gradients for MCMC sampling.
#' 
#' @param mean_rt The variable or variables (for 4par version) coding the mean reaction time in seconds in the data.
#' @param var_rt The variable or variables (for 4par version) coding the variance of the reaction time in seconds in the data
#' @param n_upper The variable coding the number of responses that hit the upper response threshold (typically the number of correct responses) in the data.
#' @param n_trials The variable coding the number of trials that was used to calculate the aggregated statistics.
#' @param links A list of links for the parameters.
#' @param version A character label for the version of the model. There is a three-parameter version
#'   (version = "3par") of the `ezdm` that fixes the relative starting point `zr` to 0.5, and a
#'   four parameter version (version = "4par"), that allows to freely estimate the starting point.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # Minimal parameter recovery example with 3-parameter EZDM
#' 
#' # Simulate data from known parameters
#' set.seed(123)
#' sim_data <- rezdm(
#'   n = 10,
#'   n_trials = 100,
#'   drift = 2,
#'   bound = 1.5,
#'   ndt = 0.3,
#'   version = "3par"
#' )
#' 
#' # Add subject ID
#' sim_data$id <- 1:10
#' 
#' # Specify model
#' model <- ezdm(
#'   mean_rt = "mean_rt",
#'   var_rt = "var_rt", 
#'   n_upper = "n_upper",
#'   n_trials = "n_trials",
#'   version = "3par"
#' )
#' 
#' # Specify formula with random effects
#' formula <- bmf(
#'   drift ~ 1 + (1 | id),
#'   bound ~ 1 + (1 | id),
#'   ndt ~ 1
#' )
#' 
#' # Fit model (using cmdstanr backend)
#' fit <- bmm(
#'   formula = formula,
#'   data = sim_data,
#'   model = model,
#'   backend = "cmdstanr",
#'   cores = 4,
#'   chains = 4,
#'   iter = 2000,
#'   warmup = 1000
#' )
#' 
#' # Check parameter recovery
#' summary(fit)
#' 
#' # Extract population-level effects
#' # True values: drift = 2, bound = 1.5, ndt = 0.3 (on log scale for drift/bound)
#' exp(fixef(fit))
#' }
ezdm <- function(mean_rt, var_rt, n_upper, n_trials, links = NULL, version = "3par", ...) {
  call <- match.call()
  stop_missing_args()
  .model_ezdm(mean_rt = mean_rt, var_rt = var_rt, n_upper = n_upper, n_trials = n_trials,
              links = links, version = version, call = call, ...)
}

#############################################################################!
# CHECK_DATA S3 methods                                                  ####
#############################################################################!

#' @export
check_data.ezdm <- function(model, data, formula) {
  # retrieve required variable names
  mean_rt <- model$resp_vars$mean_rt
  var_rt <- model$resp_vars$var_rt
  n_upper <- model$resp_vars$n_upper
  n_trials <- model$other_vars$n_trials

  # check that required model arguments were provided
  stopif(
    is.null(mean_rt) || is.null(var_rt) || is.null(n_upper) || is.null(n_trials),
    "ezdm model requires mean_rt, var_rt, n_upper, and n_trials arguments"
  )

  # validate length of mean_rt and var_rt dependent on version
  if(model$version == "3par") {
    stopif(
      length(mean_rt) != 1,
      "For ezdm version '3par', mean_rt must be a single variable name."
    )
    stopif(
      length(var_rt) != 1,
      "For ezdm version '3par', var_rt must be a single variable name."
    )
  } else if(model$version == "4par") {
    stopif(
      length(mean_rt) != 2,
      "For ezdm version '4par', mean_rt must be a vector of two variable names: c(mean_rt_upper, mean_rt_lower)."
    )
    stopif(
      length(var_rt) != 2,
      "For ezdm version '4par', var_rt must be a vector of two variable names: c(var_rt_upper, var_rt_lower)."
    )
  } else {
    stop2("Unknown ezdm version: {model$version}. Supported versions are '3par' and '4par'.")
  }

  # check that all required variables exist in data
  required_vars <- c(mean_rt, var_rt, n_upper, n_trials)
  missing_vars <- setdiff(required_vars, colnames(data))
  stopif(
    length(missing_vars) > 0,
    "The following required variables are missing from the data: {paste(missing_vars, collapse = ', ')}"
  )

  # check that mean RT values are plausible (warn if likely in milliseconds)
  # typical RTs in seconds are 0.2-3s; values > 10 suggest milliseconds
  mean_rt_values <- unlist(data[mean_rt])
  warnif(
    any(mean_rt_values > 10, na.rm = TRUE),
    "Some mean RT values are greater than 10. If your reaction times are in
    milliseconds, please convert them to seconds before fitting the model.
    The model assumes reaction times are measured in seconds."
  )

  # check that mean RT values are positive
  stopif(
    any(mean_rt_values <= 0, na.rm = TRUE),
    "Mean RT values must be positive. Found non-positive values in the data."
  )

  # check that variance values are positive
  var_rt_values <- unlist(data[var_rt])
  stopif(
    any(var_rt_values <= 0, na.rm = TRUE),
    "Variance of RT must be positive. Found non-positive values in the data."
  )

  # check that n_trials is a positive integer
  n_trials_values <- data[[n_trials]]
  stopif(
    any(n_trials_values <= 2, na.rm = TRUE),
    "Number of trials (n_trials) must be larger than two."
  )
  warnif(
    any(n_trials_values != round(n_trials_values), na.rm = TRUE),
    "Number of trials (n_trials) should be whole numbers. Found non-integer values."
  )

  # check that n_upper is a non-negative integer
  n_upper_values <- data[[n_upper]]
  stopif(
    any(n_upper_values < 0, na.rm = TRUE),
    "Number of upper boundary responses (n_upper) needs to be positive."
  )
  warnif(
    any(n_upper_values != round(n_upper_values), na.rm = TRUE),
    "Number of upper boundary responses (n_upper) should be whole numbers."
  )

  # check that n_upper <= n_trials (proportion correct between 0 and 1)
  stopif(
    any(n_upper_values > n_trials_values, na.rm = TRUE),
    "Number of upper boundary responses (n_upper) cannot exceed total trials (n_trials)."
  )

  NextMethod('check_data')
}

#############################################################################!
# Convert bmmformula to brmsformla methods                               ####
#############################################################################!

#' @export
bmf2bf.ezdm_3par <- function(model, formula) {
  # retrieve required response arguments
  mean_rt <- model$resp_vars$mean_rt
  var_rt <- model$resp_vars$var_rt
  n_upper <- model$resp_vars$n_upper
  n_trials <- model$other_vars$n_trials

  # set the base brmsformula based
  brms_formula <- brms::bf(paste0(mean_rt, " | vreal(", var_rt, ") + vint(", n_upper, ") + trials(", n_trials, ") ~ 1"))

  # return the brms_formula to add the remaining bmmformulas to it.
  brms_formula
}

#' @export
bmf2bf.ezdm_4par <- function(model, formula) {
  # retrieve required response arguments
  mean_rt <- model$resp_vars$mean_rt
  var_rt <- model$resp_vars$var_rt
  n_upper <- model$resp_vars$n_upper
  n_trials <- model$other_vars$n_trials

  # set the base brmsformula based
  brms_formula <- brms::bf(glue::glue("{mean_rt[1]} | vreal({mean_rt[2]}, {var_rt[1]}, {var_rt[2]}) + vint({n_upper}, {n_trials}) ~ 1"))

  # return the brms_formula to add the remaining bmmformulas to it.
  brms_formula
}

#############################################################################!
# CONFIGURE_MODEL S3 METHODS                                             ####
#############################################################################!

#' @export
configure_model.ezdm_3par <- function(model, data, formula) {
  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)
  links <- model$links

  # construct the family & add to formula object
  ezdm_3par_family <- function(link_drift, link_bound, link_ndt, link_s) {
    brms::custom_family(
      'ezdm_3par',
      dpars = c("mu","drift","bound","ndt","s"),
      links = c("identity",link_drift,link_bound,link_ndt,link_s),
      lb = c(NA,0,0,0,0), # lower bounds for parameters
      ub = c(NA,NA,NA,NA,NA), # upper bounds for parameters
      type = 'real', # real for continous dv, int for discrete dv
      log_lik = log_lik_ezdm_3par,
      posterior_predict = posterior_predict_ezdm_3par,
      loop = TRUE, # is the likelihood vectorized
      vars = c('vreal1[n]','vint1[n]','trials[n]')
    )
  }
  formula$family <- ezdm_3par_family(link_drift = links$drift, link_bound = links$bound, link_ndt = links$ndt, link_s = links$s)

  # prepare initial stanvars to pass to brms, model formula and priors
  sc_path <- system.file('stan_chunks', package='bmm')
  stan_functions <- read_lines2(paste0(sc_path, '/ezdm_3par_functions.stan'))

  stanvars <- brms::stanvar(scode = stan_functions, block = 'functions')

  # return the list
  nlist(formula, data, stanvars)
}

# log_lik for 3par ezdm
log_lik_ezdm_3par <- function(i, prep) {

  # extract posterior samples for distributional parameters
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # extract observed data for observation i
  mean_rt <- prep$data$Y[i]
  var_rt <- prep$data$vreal1[i]
  n_upper <- prep$data$vint1[i]
  n_trials <- prep$data$trials[i]

  # compute log-likelihood using dezdm (vectorized over posterior samples)
  dezdm(
    mean_rt = mean_rt,
    var_rt = var_rt,
    n_upper = n_upper,
    n_trials = n_trials,
    drift = drift,
    bound = bound,
    ndt = ndt,
    s = s,
    version = "3par",
    log = TRUE
  )
}

# posterior_predict for pp_check
# Returns predictions for EZDM 3par dependent variables
# By default returns mean_rt (the primary response Y)
# Use dv argument to select other variables: "var_rt", "n_upper"
# Usage: posterior_predict(fit, dv = "var_rt")
posterior_predict_ezdm_3par <- function(i, prep, ...) {
  dots <- list(...)
  dv <- if(is.null(dots$dv)) "mean_rt" else dots$dv

  # validate dv argument
  valid_dvs <- c("mean_rt", "var_rt", "n_upper")
  if (!dv %in% valid_dvs) {
    stop2("dv must be one of: {paste(valid_dvs, collapse = ', ')}")
  }

  # extract posterior samples for distributional parameters
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # extract n_trials for this observation
  n_trials <- prep$data$trials[i]

  # number of posterior draws
  n_draws <- length(drift)

  # generate predictions for each posterior draw
  pred <- sapply(seq_len(n_draws), function(j) {
    sim <- rezdm(
      n = 1,
      n_trials = n_trials,
      drift = drift[j],
      bound = bound[j],
      ndt = ndt[j],
      s = s[j],
      version = "3par"
    )
    sim[[dv]]
  })

  pred
}

#' @export
configure_model.ezdm_4par <- function(model, data, formula) {
  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)
  links <- model$links

  # construct the family & add to formula object
  ezdm_4par_family <- function(link_drift, link_bound, link_ndt, link_zr, link_s) {
    brms::custom_family(
      'ezdm_4par',
      dpars = c("mu","drift","bound","ndt","zr","s"),
      links = c("identity",link_drift,link_bound,link_ndt,link_zr,link_s),
      lb = c(NA,0,0,0,0,0), # lower bounds for parameters
      ub = c(NA,NA,NA,NA,1,NA), # upper bounds for parameters
      type = 'real', # real for continous dv, int for discrete dv
      log_lik = log_lik_ezdm_4par,
      posterior_predict = posterior_predict_ezdm_4par,
      loop = TRUE, # is the likelihood vectorized
      vars = c('vreal1[n]','vreal2[n]','vreal3[n]','vint1[n]','vint2[n]')
    )
  }

  formula$family <- ezdm_4par_family(link_drift = links$drift, link_bound = links$bound, link_ndt = links$ndt, link_zr = links$zr, link_s = links$s)

  # prepare initial stanvars to pass to brms, model formula and priors
  sc_path <- system.file('stan_chunks', package='bmm')
  stan_functions <- read_lines2(paste0(sc_path, '/ezdm_4par_functions.stan'))

  stanvars <- brms::stanvar(scode = stan_functions, block = 'functions')

  # return the list
  nlist(formula, data, stanvars)
}

# log_lik for 4par ezdm
log_lik_ezdm_4par <- function(i, prep) {

  # extract posterior samples for distributional parameters
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # extract observed data for observation i
  # based on bmf2bf.ezdm_4par formula:
  # Y = mean_rt_upper, vreal1 = mean_rt_lower
  # vreal2 = var_rt_upper, vreal3 = var_rt_lower
  # vint1 = n_upper, vint2 = n_trials
  mean_rt_upper <- prep$data$Y[i]
  mean_rt_lower <- prep$data$vreal1[i]
  var_rt_upper <- prep$data$vreal2[i]
  var_rt_lower <- prep$data$vreal3[i]
  n_upper <- prep$data$vint1[i]
  n_trials <- prep$data$vint2[i]

  # compute log-likelihood using dezdm (vectorized over posterior samples)
  dezdm(
    mean_rt = c(mean_rt_upper, mean_rt_lower),
    var_rt = c(var_rt_upper, var_rt_lower),
    n_upper = n_upper,
    n_trials = n_trials,
    drift = drift,
    bound = bound,
    ndt = ndt,
    zr = zr,
    s = s,
    version = "4par",
    log = TRUE
  )
}

# posterior_predict for pp_check
# Returns predictions for EZDM 4par dependent variables
# By default returns mean_rt_upper (the primary response Y)
# Use dv argument to select other variables:
#   "mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower", "n_upper"
# Usage: posterior_predict(fit, dv = "var_rt_upper")
posterior_predict_ezdm_4par <- function(i, prep, ...) {
  dots <- list(...)
  dv <- if(is.null(dots$dv)) "mean_rt_upper" else dots$dv

  # validate dv argument
  valid_dvs <- c("mean_rt_upper", "mean_rt_lower",
                 "var_rt_upper", "var_rt_lower", "n_upper")
  if (!dv %in% valid_dvs) {
    stop2("dv must be one of: {paste(valid_dvs, collapse = ', ')}")
  }

  # extract posterior samples for distributional parameters
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # extract n_trials for this observation
  n_trials <- prep$data$vint2[i]

  # number of posterior draws
  n_draws <- length(drift)

  # generate predictions for each posterior draw
  pred <- sapply(seq_len(n_draws), function(j) {
    sim <- rezdm(
      n = 1,
      n_trials = n_trials,
      drift = drift[j],
      bound = bound[j],
      ndt = ndt[j],
      zr = zr[j],
      s = s[j],
      version = "4par"
    )
    sim[[dv]]
  })

  pred
}
