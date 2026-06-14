############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_binary <- function(response = NULL, stimulus = NULL,
                              n_trials = NULL, dist = "normal",
                              links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)

  parameters <- list(
    dprime = "Sensitivity: distance between signal and noise distributions",
    criterion = "Response bias: location of decision boundary",
    sdratio = paste0(
      "Log SD ratio: log ratio of signal to noise standard deviations ",
      "(exp(sdratio) is the natural SD ratio, 0 = equal variance)"
    )
  )
  default_priors <- list(
    dprime = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
    criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)"),
    sdratio = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)")
  )
  param_links <- list(dprime = "identity", criterion = "identity",
                      sdratio = "identity")

  requirements <- glue(
    "Provide pre-aggregated data with the following columns:", "\n\n",
    "  - Response counts (n_old): number of 'old'/'signal' responses", "\n",
    "  - Stimulus type (stimulus): 0 = noise, 1 = signal", "\n",
    "  - Number of trials (n_trials): total trials per cell"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(stimulus, n_trials, dist, dist_int),
      domain = "Perception & Recognition Memory",
      task = "Signal/Noise or Old/New Recognition",
      name = "Signal Detection Theory (Binary)",
      citation = glue(
        "Green, D. M., & Swets, J. A. (1966). Signal detection theory ",
        "and psychophysics. Wiley."
      ),
      version = "binary",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = list(sdratio = 0),
      default_priors = default_priors,
      init_ranges = list(
        dprime    = c(0.5, 1.5),
        criterion = c(-0.5, 0.5),
        sdratio   = c(-0.3, 0.3)
      ),
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdt", "sdt_binary"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Binary Signal Detection Theory Model
#' @name sdt_binary
#' @details `r model_info(.model_sdt_binary())`
#' @param response A single string naming the column with "old"/"signal"
#'   response counts.
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (noise/new) and 1 (signal/old).
#' @param n_trials The name of the variable containing the total number of
#'   trials per cell.
#' @param dist Character. The noise distribution to use:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT
#'     \item "logistic": Logistic SDT
#'     \item "gumbel_min": Extreme-value (Gumbel minimum) SDT
#'     \item "gumbel_max": Gumbel maximum SDT
#'   }
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Green, D. M., & Swets, J. A. (1966). \emph{Signal detection theory and
#'   psychophysics}. Wiley.
#'
#' DeCarlo, L. T. (1998). Signal detection theory and generalized linear
#'   models. \emph{Psychological Methods}, \emph{3}(2), 186--205.
#'   \doi{10.1037/1082-989X.3.2.186}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' dat <- rsdt_binary(n_per_cell = 100, n_subjects = 20,
#'                    dprime = 1.5, criterion = 0.2)
#'
#' model <- sdt_binary(
#'   response = "n_old",
#'   stimulus = "stimulus",
#'   n_trials = "n_trials"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Unequal-variance binary SDT
#' fit_uv <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, sdratio ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' }
sdt_binary <- function(response, stimulus, n_trials,
                       dist = c("normal", "logistic", "gumbel_min", "gumbel_max"),
                       links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)

  .model_sdt_binary(response = response, stimulus = stimulus,
                    n_trials = n_trials, dist = dist,
                    links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_binary <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stim_vals <- unique(data[[stim_var]])
  stopif(!all(stim_vals %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (noise) and 1 (signal)")

  resp_var <- model$resp_vars$response
  n_trials_var <- model$other_vars$n_trials

  .validate_sdt_counts(data, resp_var, n_trials_var)

  warnif(any(data[[n_trials_var]] != round(data[[n_trials_var]]), na.rm = TRUE),
         "Variable '{n_trials_var}' should contain integer counts")

  data$dist_type <- model$other_vars$dist_int

  NextMethod("check_data")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt_binary <- function(model, formula) {
  resp_var <- model$resp_vars$response
  stim_var <- model$other_vars$stimulus
  n_trials_var <- model$other_vars$n_trials

  brms::bf(paste0(
    resp_var, " | vint(", stim_var, ", dist_type) + trials(",
    n_trials_var, ") ~ 0"
  ))
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_binary <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  formula$family <- brms::custom_family(
    "sdt_binary",
    dpars = c("mu", "dprime", "criterion", "sdratio"),
    links = c("identity", model$links$dprime, model$links$criterion,
              model$links$sdratio),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_sdt_binary,
    posterior_predict = posterior_predict_sdt_binary,
    vars = c("vint1[n]", "vint2[n]", "trials[n]")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdt_binary_funs.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

log_lik_sdt_binary <- function(i, prep) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  sdratio <- exp(brms::get_dpar(prep, "sdratio", i = i))
  stimulus <- prep$data$vint1[i]
  dist <- .sdt_dist_names[prep$data$vint2[i]]
  n_trials <- prep$data$trials[i]
  y <- prep$data$Y[i]

  eta <- .sdt_eta(dprime, criterion, stimulus, sdratio)
  p <- .sdt_cdf(eta, dist)
  stats::dbinom(y, n_trials, p, log = TRUE)
}

posterior_predict_sdt_binary <- function(i, prep, ...) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  sdratio <- exp(brms::get_dpar(prep, "sdratio", i = i))
  stimulus <- prep$data$vint1[i]
  dist <- .sdt_dist_names[prep$data$vint2[i]]
  n_trials <- prep$data$trials[i]

  eta <- .sdt_eta(dprime, criterion, stimulus, sdratio)
  p <- .sdt_cdf(eta, dist)
  stats::rbinom(length(p), n_trials, p)
}
