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
      drift = list(main = "normal(0,1)", effects = "normal(0,0.5)"),
      bound = list(main = "normal(0,0.5)", effects = "normal(0,0.5)"),
      ndt = list(main = "normal(-1.5,0.5)", effects = "normal(0,0.3)"),
      s = list(main = "normal(0,1)", effects = "normal(0,0.3)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(0.5,1.5),
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
      zr = "Relative startin point = Starting point between the decision thresholds relative to the upper bound.",
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
      drift = c(-0.5,0.5),
      bound = c(1,2),
      ndt = c(0.25, 0.5),
      zr = c(0.45, 0.55),
      s = c(0.99,1.01)
    )
  )
)


.model_ezdm <- function(mean_rt = NULL, var_rt = NULL, n_upper = NULL, n_trials = NULL, response = NULL, version = "4par", links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(mean_rt, var_rt, n_upper, response),
      other_vars = nlist(n_trials),
      domain = "Processing Speed, Decision Making",
      task = "Choice Reaction Time tasks",
      name = "EZ-Diffusion Model",
      citation = glue(
        "Wagenmakers, E.-J., Van Der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. Psychonomic Bulletin & Review, 14(1), 3–22. https://doi.org/10/fk447c","\n",
        "- Chávez De la Peña, A. F., & Vandekerckhove, J. (2025). An EZ Bayesian hierarchical drift diffusion model for response time and accuracy. Psychonomic Bulletin & Review. https://doi.org/10.3758/s13423-025-02729-y"
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
#' @name ezdm,
#' @details `r model_info(.model_ezdm())`
#' @param mean_rt The variable coding the mean reaction time in seconds in the data.
#' @param var_rt The variable coding the variance of the reaction time in seconds in the data
#' @param n_upper The variable coding the number of responses that hit the upper response threshold (typically the number of correct responses) in the data.
#' @param n_trials The variable coding the number of trials that was used to calculated the aggregated statistics.
#' @param links A list of links for the parameters.
#' @param version A character label for the version of the model. There is a three-parameter version
#'   (version = "3par")of the `ezdm` that fixes the relativ starting point `zr` to 0.5, and a
#'   four parameter verison (version = "4par"), that allows to freely estimate the starting point.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @export
#' @examples
#' \dontrun{
#' # put a full example here (see 'R/model_mixture3p.R' for an example)
#' }
ezdm <- function(mean_rt, var_rt, n_upper, n_trials, response = NULL, links = NULL, version = "3par", ...) {
  call <- match.call()
  stop_missing_args()
  .model_ezdm(mean_rt = mean_rt, var_rt = var_rt, n_upper = n_upper, n_trials = n_trials, response = response,
              links = links, version = version, call = call, ...)
}

#############################################################################!
# CHECK_DATA S3 methods                                                  ####
#############################################################################!

#' @export
check_data.ezdm <- function(model, data, formula) {
  # retrieve required arguments

  # check the data (required)

  # compute any necessary transformations (optional)

  # save some variables as attributes of the data for later use (optional)

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
  response <- model$resp_vars$response

  # set the base brmsformula based
  brms_formula <- brms::bf(paste0(mean_rt, " | vreal(", var_rt, ") + vint(", n_upper, ") + trials(", n_trials, ") + dec(",response,") ~ 1"))

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

log_lik_ezdm_3par <- function(i, prep) {

  out
}

posterior_predict_ezdm_3par <- function(i, prep, ...) {
  dots <- list(...)

  # out
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
      vars = c('vreal1[n]','vint1[n]','trials[n]','dec[n]')
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



log_lik_ezdm_4par <- function(i, prep) {

  out
}

posterior_predict_ezdm_4par <- function(i, prep, ...) {
  dots <- list(...)

  # out
}
