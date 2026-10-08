############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.ezcdm_parameters <- list(
  driftrate = "Drift rate = Length of the drift vector, the average speed of evidence accumulation",
  driftangle = "Drift angle = Direction of the drift vector in radians, the average response bias relative to the target",
  bound = "Boundary = Radius of the circular decision boundary that the accumulated evidence has to reach",
  ndt = "Non-decision time = Additional time required beyond the evidence accumulation process"
)

.ezcdm_priors <- list(
  driftrate = list(main = "normal(0.5, 1)", effects = "normal(0, 0.5)"),
  driftangle = list(main = "normal(0, 1)", effects = "normal(0, 0.5)"),
  bound = list(main = "normal(0.5, 1)", effects = "normal(0, 0.5)"),
  ndt = list(main = "normal(-1.5, 0.5)", effects = "normal(0, 0.3)")
)

.ezcdm_links <- list(driftrate = "log", driftangle = "identity", bound = "log", ndt = "log")

.ezcdm_init_ranges <- list(
  mu = c(0, 1),
  driftrate = c(1, 2),
  driftangle = c(-0.1, 0.1),
  bound = c(1, 2),
  ndt = c(0.25, 0.5)
)

.model_ezcdm <- function(mean_angle = NULL, var_angle = NULL, mean_rt = NULL, var_rt = NULL,
                         n_trials = NULL, version = "3par", links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(mean_angle, var_angle, mean_rt, var_rt),
      other_vars = nlist(n_trials),
      domain = "Decision Making / Response times",
      task = "Continuous reproduction tasks with response times",
      name = "EZ Circular Diffusion Model",
      citation = glue(
        "Qarehdaghi, H., & Amani Rad, J. (2024). EZ-CDM: Fast, simple, robust, and accurate estimation of circular diffusion model parameters. Psychonomic Bulletin & Review, 31(5), 2058-2091. https://doi.org/10.3758/s13423-024-02483-7", "\n",
        "- Smith, P. L. (2016). Diffusion theory of decision making in continuous report. Psychological Review, 123(4), 425-451. https://doi.org/10.1037/rev0000023"
      ),
      version = version,
      requirements = glue(
        "Provide aggregated statistics for each subject and condition that model parameters should vary over:", "\n\n",
        "  - Circular mean of the response angles (mean_angle) in radians", "\n",
        "  - Circular variance of the response angles (var_angle), one minus the mean resultant length", "\n",
        "  - Mean reaction time (mean_rt) in seconds", "\n",
        "  - Variance of the reaction times (var_rt) in seconds^2", "\n",
        "  - Total number of trials used to calculate the aggregated statistics (n_trials)", "\n\n",
        "ezcdm_summary_stats() computes all of them from trial-level data."
      ),
      parameters = .ezcdm_parameters,
      links = .ezcdm_links,
      fixed_parameters = if (version == "3par") list(mu = 0, driftangle = 0) else list(mu = 0),
      default_priors = .ezcdm_priors,
      init_ranges = .ezcdm_init_ranges
    ),
    class = c("bmmodel", "ezcdm", paste0("ezcdm_", version)),
    call = call
  )
  out$links[names(links)] <- links
  out
}

#' @title `r .model_ezcdm()$name`
#' @name ezcdm
#' @details `r model_docs(.model_ezcdm(version = "4par"))`
#'
#'   The response angle of the circular diffusion model follows a von Mises
#'   distribution with mean `driftangle` and concentration
#'   `bound * driftrate`, independent of the decision time. The likelihood
#'   of the angle summaries is therefore exact for any number of trials, while
#'   the variance of the reaction times enters through a Gamma approximation
#'   matched to the first two moments of the sample variance of the
#'   right-skewed decision times, and the mean reaction time through a normal
#'   approximation conditional on the variance that keeps the exact covariance
#'   of the two statistics (see [ezcdm_dist]).
#'   Angles should be coded relative to the target, so that a drift angle of 0
#'   means responses centred on the target; [ezcdm_summary_stats()] does this
#'   with its `target` argument.
#'
#' @param mean_angle The name of the variable coding the circular mean of the
#'   response angles (in radians) in the data.
#' @param var_angle The name of the variable coding the circular variance of the
#'   response angles (one minus the mean resultant length, in \[0, 1\]) in the
#'   data.
#' @param mean_rt The name of the variable coding the mean reaction time in
#'   seconds in the data.
#' @param var_rt The name of the variable coding the variance of the reaction
#'   times in seconds^2 in the data.
#' @param n_trials The name of the variable coding the number of trials that was
#'   used to calculate the aggregated statistics.
#' @param links A list of links for the parameters. For the positive parameters
#'   (`driftrate`, `bound`, `ndt`), "softplus" is available as an
#'   alternative to the default "log" link that grows linearly for large values
#'   and avoids the numerical blow-up of `exp()`. The drift angle `driftangle`
#'   has an identity link on the radian scale. Coding the responses as
#'   deviations from the target keeps it near 0, far from the wrap-around at
#'   \eqn{\pm\pi}; for data centred near \eqn{\pm\pi},
#'   `links = list(driftangle = "tan_half")` keeps the angle in
#'   \eqn{(-\pi, \pi]}.
#' @param version A character label for the version of the model. The
#'   three-parameter version (`version = "3par"`, default) fixes the drift angle
#'   `driftangle` to 0, i.e. responses centred on the target; the four-parameter
#'   version (`version = "4par"`) estimates it. Note that providing a formula
#'   for `driftangle` (e.g. `driftangle ~ 1`) with the 3par version frees the
#'   drift angle, which makes it equivalent to the 4par version.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # simulate summary statistics for 20 subjects from known parameters
#' set.seed(123)
#' sim_data <- rezcdm(
#'   n = 20, n_trials = 100,
#'   driftrate = 2, bound = 1.5, ndt = 0.3
#' )
#' sim_data$id <- 1:20
#'
#' model <- ezcdm(
#'   mean_angle = "mean_angle", var_angle = "var_angle",
#'   mean_rt = "mean_rt", var_rt = "var_rt", n_trials = "n_trials"
#' )
#'
#' formula <- bmf(
#'   driftrate ~ 1 + (1 | id),
#'   bound ~ 1 + (1 | id),
#'   ndt ~ 1
#' )
#'
#' fit <- bmm(
#'   formula = formula,
#'   data = sim_data,
#'   model = model,
#'   backend = "cmdstanr",
#'   cores = 4,
#'   chains = 4
#' )
#'
#' summary(fit)
#' }
ezcdm <- function(mean_angle, var_angle, mean_rt, var_rt, n_trials, links = NULL,
                  version = c("3par", "4par"), ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  .model_ezcdm(
    mean_angle = mean_angle, var_angle = var_angle, mean_rt = mean_rt,
    var_rt = var_rt, n_trials = n_trials, links = links, version = version,
    call = call, ...
  )
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !

#' @export
check_data.ezcdm <- function(model, data, formula) {
  resp_vars <- model$resp_vars
  n_trials <- model$other_vars$n_trials
  vars <- c(resp_vars, model$other_vars)

  stopif(
    any(lengths(vars) != 1),
    "Each of these arguments must be a single variable name: \\
    {collapse_comma(names(vars)[lengths(vars) != 1])}"
  )

  missing_vars <- setdiff(unlist(vars), colnames(data))
  stopif(
    length(missing_vars),
    "The following required variables are missing from the data: {collapse_comma(missing_vars)}"
  )

  warnif(
    any(abs(data[[resp_vars$mean_angle]]) > 2 * pi, na.rm = TRUE),
    "It appears your mean_angle variable is in degrees.
    The model requires the circular mean to be in radians.
    The model will continue to run, but the results may be compromised."
  )

  var_angle_values <- data[[resp_vars$var_angle]]
  stopif(
    any(var_angle_values < 0 | var_angle_values > 1, na.rm = TRUE),
    "Circular variance (var_angle) must be between 0 and 1."
  )

  check_rt_summary_vars(data, resp_vars$mean_rt, resp_vars$var_rt)
  check_n_trials_var(data, n_trials)

  NextMethod("check_data")
}

############################################################################# !
# Convert bmmformula to brmsformla methods                               ####
############################################################################# !

#' @export
bmf2bf.ezcdm <- function(model, formula) {
  resp <- model$resp_vars
  brms::bf(glue(
    "{resp$mean_angle} | vreal({resp$var_angle}, {resp$mean_rt}, {resp$var_rt}) + trials({model$other_vars$n_trials}) ~ 1"
  ))
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.ezcdm <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)
  links <- model$links

  formula$family <- brms::custom_family(
    "ezcdm",
    dpars = c("mu", "driftrate", "driftangle", "bound", "ndt"),
    links = c("identity", links$driftrate, links$driftangle, links$bound, links$ndt),
    lb = c(NA, 0, NA, 0, 0),
    ub = rep(NA, 5),
    type = "real",
    log_lik = log_lik_ezcdm,
    posterior_predict = posterior_predict_ezcdm,
    loop = TRUE,
    vars = c("vreal1[n]", "vreal2[n]", "vreal3[n]", "trials[n]")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stanvars <- brms::stanvar(
    scode = read_lines2(paste0(sc_path, "/ezcdm_functions.stan")),
    block = "functions"
  )

  nlist(formula, data, stanvars)
}

log_lik_ezcdm <- function(i, prep) {
  dezcdm(
    mean_angle = prep$data$Y[i],
    var_angle = prep$data$vreal1[i],
    mean_rt = prep$data$vreal2[i],
    var_rt = prep$data$vreal3[i],
    n_trials = prep$data$trials[i],
    driftrate = brms::get_dpar(prep, "driftrate", i = i),
    driftangle = brms::get_dpar(prep, "driftangle", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    log = TRUE
  )
}

# brms forwards posterior_predict() dots to prepare_predictions() only, never
# to this function; use pp_check(fit, resp_var = ...) for the other summaries
posterior_predict_ezcdm <- function(i, prep, ...) {
  rezcdm(
    n = prep$ndraws,
    n_trials = prep$data$trials[i],
    driftrate = brms::get_dpar(prep, "driftrate", i = i),
    driftangle = brms::get_dpar(prep, "driftangle", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i)
  )$mean_angle
}

############################################################################# !
# PP_CHECK OBSERVABLES                                                    ####
############################################################################# !

#' @export
pp_observables.ezcdm <- function(model) {
  list(
    observed = c(mean_angle = "Y", var_angle = "vreal1", mean_rt = "vreal2",
                 var_rt = "vreal3", n_trials = "trials"),
    checks = list(
      mean_angle = .pp_observable(function(d) d$mean_angle,
                                  label = "Circular mean of the response angles (radians)"),
      var_angle = .pp_observable(function(d) d$var_angle,
                                 label = "Circular variance of the response angles"),
      mean_rt = .pp_observable(function(d) d$mean_rt,
                               label = "Mean response time"),
      var_rt = .pp_observable(function(d) d$var_rt,
                              label = "Response time variance")
    )
  )
}

#' @export
pp_simulate.ezcdm <- function(model, prep) {
  .pp_simulate_joint(prep, rezcdm, c("driftrate", "driftangle", "bound", "ndt"),
                     n_trials = rep(prep$data$trials, each = prep$ndraws))
}
