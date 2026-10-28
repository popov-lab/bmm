############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.ezcdm_parameters <- list(
  drift = "Drift rate = Length of the drift vector, the average speed of evidence accumulation",
  driftangle = "Drift angle = Direction of the drift vector in radians, the average response bias relative to the target",
  bound = "Boundary = Radius of the circular decision boundary that the accumulated evidence has to reach",
  ndt = "Non-decision time = Additional time required beyond the evidence accumulation process"
)

# sd rates by meaning (#342): 2 for the log-linked drift rate, boundary and
# non-decision time, as in ezdm and cswald, and 2 for the drift angle, a
# circular bias like sdm's mu. mu has rate 4 on the tan_half scale, where
# angle = 2 * atan(eta) is about twice eta near 0, so the same between-person
# spread in radians is rate 2 (median 0.35 rad)
.ezcdm_priors <- list(
  drift = list(main = "normal(0.5, 1)", effects = "normal(0, 0.5)", sd = "exponential(2)"),
  driftangle = list(main = "normal(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(2)"),
  bound = list(main = "normal(0.5, 1)", effects = "normal(0, 0.5)", sd = "exponential(2)"),
  ndt = list(main = "normal(-1.5, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)")
)

.ezcdm_links <- list(drift = "log", driftangle = "identity", bound = "log", ndt = "log")

.ezcdm_init_ranges <- list(
  mu = c(0, 1),
  drift = c(1, 2),
  driftangle = c(-0.1, 0.1),
  bound = c(1, 2),
  ndt = c(0.25, 0.5)
)

.model_ezcdm <- function(mean_angle = NULL, var_angle = NULL, mean_rt = NULL, var_rt = NULL,
                         n_trials = NULL, version = "3par", links = NULL, call = NULL, ...) {
  set_links(structure(
    list(
      resp_vars = nlist(mean_angle, var_angle, mean_rt, var_rt),
      other_vars = nlist(n_trials),
      domain = "Decision Making / Response times",
      task = "Continuous reproduction tasks with response times",
      name = "EZ Circular Diffusion Model",
      citation = c(
        "Qarehdaghi, H., & Amani Rad, J. (2024). EZ-CDM: Fast, simple, robust, and accurate estimation of circular diffusion model parameters. Psychonomic Bulletin & Review, 31(5), 2058-2091. https://doi.org/10.3758/s13423-024-02483-7",
        "Smith, P. L. (2016). Diffusion theory of decision making in continuous report. Psychological Review, 123(4), 425-451. https://doi.org/10.1037/rev0000023"
      ),
      version = version,
      requirements = glue(
        "Provide aggregated statistics for each subject and condition that model parameters should vary over:", "\n\n",
        "  - Circular mean of the response angles (mean_angle) in radians", "\n",
        "  - Circular variance of the response angles (var_angle), one minus the mean resultant length", "\n",
        "  - Mean reaction time (mean_rt) in seconds", "\n",
        "  - Variance of the reaction times (var_rt) in seconds^2", "\n",
        "  - Total number of trials used to calculate the aggregated statistics (n_trials), at least 3"
      ),
      parameters = .ezcdm_parameters,
      links = .ezcdm_links,
      fixed_parameters = if (version == "3par") list(mu = 0, driftangle = 0) else list(mu = 0),
      default_priors = .ezcdm_priors,
      init_ranges = .ezcdm_init_ranges
    ),
    class = c("bmmodel", "ezcdm", paste0("ezcdm_", version)),
    call = call
  ), links)
}

#' @title `r .model_ezcdm()$name`
#' @name ezcdm
#' @details `r model_docs(.model_ezcdm())`
#'
#'   The default version "3par" fixes the drift angle `driftangle` to 0;
#'   version "4par" estimates it.
#'
#'   The within-trial noise of the two-dimensional diffusion is fixed to
#'   \eqn{\sigma = 1}, which sets the scale of the model. Unlike [ezdm()],
#'   `ezcdm()` therefore has no diffusion constant `s`. With time in seconds,
#'   `bound` is in evidence units, that is in units of the noise, and `drift`
#'   in evidence units per second. `bound` is the radius of the circular
#'   decision boundary, not a separation between two boundaries as in
#'   [ezdm()], so the `bound` estimates of the two models are not comparable.
#'
#'   The response angle of the circular diffusion model follows a von Mises
#'   distribution with mean `driftangle` and concentration
#'   `bound * drift`, independent of the decision time. The likelihood
#'   of the angle summaries is therefore exact for any number of trials, while
#'   the variance of the reaction times enters through a Gamma approximation
#'   matched to the first two moments of the sample variance of the
#'   right-skewed decision times, and the mean reaction time through a normal
#'   approximation conditional on the variance that keeps the exact covariance
#'   of the two statistics (see [ezcdm_dist]).
#'   Angles should be coded relative to the target, so that a drift angle of 0
#'   means responses centered on the target; [ezcdm_summary_stats()] computes
#'   all five summaries from trial-level data and does this with its `target`
#'   argument.
#'
#' @param mean_angle The name of the variable coding the circular mean of the
#'   response angles (in radians) in the data. Circular means coded on
#'   \[0, \eqn{2\pi}) are wrapped to \[\eqn{-\pi}, \eqn{\pi}), the range of the
#'   circular means that `posterior_predict()` and `pp_check()` simulate.
#' @param var_angle The name of the variable coding the circular variance of the
#'   response angles (one minus the mean resultant length, in \[0, 1\]) in the
#'   data.
#' @param mean_rt The name of the variable coding the mean reaction time in
#'   seconds in the data.
#' @param var_rt The name of the variable coding the variance of the reaction
#'   times in seconds^2 in the data.
#' @param n_trials The name of the variable coding the number of trials that was
#'   used to calculate the aggregated statistics: whole numbers, at least 3.
#' @param links A named list of links for the parameters, e.g.
#'   `links = list(bound = "softplus")`. For the positive parameters
#'   (`drift`, `bound`, `ndt`), "softplus" is available as an
#'   alternative to the default "log" link that grows linearly for large values
#'   and avoids the numerical blow-up of `exp()`. The drift angle `driftangle`
#'   has an identity link on the radian scale. Coding the responses as
#'   deviations from the target keeps it near 0, far from the wrap-around at
#'   \eqn{\pm\pi}; for data centered near \eqn{\pm\pi},
#'   `links = list(driftangle = "tan_half")` keeps the angle in
#'   \eqn{(-\pi, \pi]}. A name that is not a parameter of the model is an
#'   error (a close misspelling is read as the parameter, with a warning), and
#'   a link that allows values the default link excludes (e.g. "identity" for
#'   a positive parameter) is a warning.
#' @param version A character label for the version of the model. The
#'   three-parameter version (`version = "3par"`, default) fixes the drift angle
#'   `driftangle` to 0, i.e. responses centered on the target; the
#'   four-parameter version (`version = "4par"`) estimates it. A formula for
#'   `driftangle` (e.g. `driftangle ~ 1`) with version "3par" switches the
#'   model to version "4par", with a message, so that the fit, its `summary()`
#'   and `update()` carry the version that was estimated.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # simulate summary statistics for 20 subjects whose drift rate and boundary
#' # vary around known population values (drift = 2, bound = 1.5)
#' set.seed(123)
#' n_subjects <- 20
#' sim_data <- rezcdm(
#'   n = n_subjects, n_trials = 100,
#'   drift = exp(rnorm(n_subjects, log(2), 0.2)),
#'   bound = exp(rnorm(n_subjects, log(1.5), 0.2)),
#'   ndt = 0.3
#' )
#' sim_data$id <- seq_len(n_subjects)
#'
#' model <- ezcdm(
#'   mean_angle = "mean_angle", var_angle = "var_angle",
#'   mean_rt = "mean_rt", var_rt = "var_rt", n_trials = "n_trials"
#' )
#'
#' formula <- bmf(
#'   drift ~ 1 + (1 | id),
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
#'
#' # drift, bound and ndt have log links
#' exp(brms::fixef(fit)[c("drift_Intercept", "bound_Intercept", "ndt_Intercept"), "Estimate"])
#' }
ezcdm <- function(mean_angle, var_angle, mean_rt, var_rt, n_trials, links = NULL,
                  version = c("3par", "4par"), ...) {
  call <- match.call()
  .check_constructor_dots(list(...))
  stop_missing_args()
  version <- match.arg(version)
  .model_ezcdm(
    mean_angle = mean_angle, var_angle = var_angle, mean_rt = mean_rt,
    var_rt = var_rt, n_trials = n_trials, links = links, version = version,
    call = call, ...
  )
}

############################################################################# !
# CHECK_MODEL S3 methods                                                 ####
############################################################################# !

# The version only sets whether driftangle is fixed by default, so a formula
# that predicts it asks for the 4par model. The version, class and call then
# say so, which summary() prints and update() keeps; update() calls this
# after update_model_fixed_parameters() has freed the drift angle, as bmm()
# does through check_model.bmmodel().
#' @export
check_model.ezcdm_3par <- function(model, data = NULL, formula = NULL) {
  if (!"driftangle" %in% names(model$fixed_parameters)) {
    message2("The formula predicts driftangle, which ezcdm version '3par' \\
             fixes to 0, so the model is fitted as version '4par'.")
    model$version <- "4par"
    class(model)[class(model) == "ezcdm_3par"] <- "ezcdm_4par"
    if (!is.null(attr(model, "call"))) {
      attr(model, "call")$version <- "4par"
    }
  }
  NextMethod("check_model")
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !

#' @export
check_data.ezcdm <- function(model, data, formula) {
  resp_vars <- model$resp_vars
  vars <- c(resp_vars, model$other_vars)

  # an argument check, kept here as in check_data.ezdm(), where the number of
  # names depends on the version
  stopif(
    any(lengths(vars) != 1),
    "Each of these arguments must be a single variable name: \\
    {collapse_comma(names(vars)[lengths(vars) != 1])}"
  )

  columns <- unlist(vars)
  missing_vars <- setdiff(columns, colnames(data))
  stopif(
    length(missing_vars),
    "The following required variables are missing from the data: {collapse_comma(missing_vars)}"
  )

  # comparisons on a character column compare strings and on a factor warn,
  # so the range checks below would pass or fail for the wrong reason
  non_numeric <- !vapply(data[columns], is.numeric, logical(1))
  stopif(
    any(non_numeric),
    "These variables must be numeric: {collapse_comma(columns[non_numeric])}. \\
    A factor holding numbers converts with as.numeric(as.character(x))."
  )

  mean_angle <- resp_vars$mean_angle
  warnif(
    any(abs(data[[mean_angle]]) > 2 * pi, na.rm = TRUE),
    "It appears your mean_angle variable ('{mean_angle}') is in degrees.
    The model requires the circular mean to be in radians.
    The model will continue to run, but the results may be compromised."
  )
  # posterior_predict() and pp_check() simulate circular means in (-pi, pi], so
  # means coded on [0, 2 * pi) would be compared with simulations that never
  # reach them. The likelihood reads cos(mean_angle - driftangle) only, so the
  # fit is unchanged.
  outside <- which(abs(data[[mean_angle]]) > pi)
  data[[mean_angle]][outside] <- wrap(data[[mean_angle]][outside])

  var_angle <- resp_vars$var_angle
  stopif(
    any(data[[var_angle]] < 0 | data[[var_angle]] > 1, na.rm = TRUE),
    "Circular variance ('{var_angle}') must be between 0 and 1."
  )

  check_rt_summary_vars(data, resp_vars$mean_rt, resp_vars$var_rt)
  check_n_trials_var(data, model$other_vars$n_trials)

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
    dpars = c("mu", "drift", "driftangle", "bound", "ndt"),
    links = c("identity", links$drift, links$driftangle, links$bound, links$ndt),
    lb = c(NA, 0, NA, 0, 0),
    ub = rep(NA, 5),
    type = "real",
    log_lik = log_lik_ezcdm,
    posterior_predict = posterior_predict_ezcdm,
    posterior_epred = posterior_epred_undefined("ezcdm"),
    loop = TRUE,
    vars = c("vreal1[n]", "vreal2[n]", "vreal3[n]", "trials[n]")
  )

  stanvars <- brms::stanvar(
    scode = read_lines2(file.path(system.file("stan_chunks", package = "bmm"), "ezcdm_functions.stan")),
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
    drift = brms::get_dpar(prep, "drift", i = i),
    driftangle = brms::get_dpar(prep, "driftangle", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    log = TRUE
  )
}

# The brms response is the circular mean of the angles, whose expected value
# is no more useful than that of the circular models' response errors.
# posterior_epred() and fitted() on a bmmfit refuse through this method, which
# dispatches on the model class, before brms prepares any draws; the family's
# posterior_epred refuses where brms reaches the family itself.
#' @export
expected_response_defined.ezcdm <- function(model) {
  FALSE
}

# brms forwards posterior_predict() dots to prepare_predictions() only, never
# to this function; use pp_check(fit, resp_var = ...) for the other summaries
posterior_predict_ezcdm <- function(i, prep, ...) {
  rezcdm(
    n = prep$ndraws,
    n_trials = prep$data$trials[i],
    drift = brms::get_dpar(prep, "drift", i = i),
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
  .pp_simulate_joint(prep, rezcdm, c("drift", "driftangle", "bound", "ndt"),
                     n_trials = rep(prep$data$trials, each = prep$ndraws))
}
