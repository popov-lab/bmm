############################################################################# !
# MODELS                                                                 ####
############################################################################# !

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
    fixed_parameters = list(s = 0, mu = 0),
    priors = list(
      drift = list(main = "cauchy(0,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
      bound = list(main = "normal(0,0.5)", effects = "normal(0,0.5)", sd = "exponential(2)"),
      ndt = list(main = "normal(-1.5,0.5)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      s = list(main = "normal(0,1)", effects = "normal(0,0.3)", sd = "exponential(2)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(-1,1),
      bound = c(1,2),
      ndt = c(0.25, 0.5),
      s = c(0.99, 1.01)
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
    fixed_parameters = list(s = 0, mu = 0),
    priors = list(
      drift = list(main = "cauchy(0,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
      bound = list(main = "normal(0,0.5)", effects = "normal(0,0.5)", sd = "exponential(2)"),
      ndt = list(main = "normal(-1.5,0.5)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      zr = list(main = "normal(0,0.5)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      s = list(main = "normal(0,1)", effects = "normal(0,0.3)", sd = "exponential(2)")
    ),
    init_ranges = list(
      mu = c(0,1),
      drift = c(-1,1),
      bound = c(1,2),
      ndt = c(0.25, 0.5),
      zr = c(0.45, 0.55),
      s = c(0.99, 1.01)
    )
  )
)


.model_ezdm <- function(mean_rt = NULL, var_rt = NULL, n_upper = NULL, n_trials = NULL, version = "3par", links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(mean_rt, var_rt, n_upper),
      other_vars = nlist(n_trials),
      domain = "Decision Making / Response times",
      task = "Choice Reaction Time tasks",
      name = "EZ-Diffusion Model",
      citation = c(
        glue(
          "Wagenmakers, E.-J., Van Der Maas, H. L. J., & Grasman, R. P. P. P. \\
          (2007). An EZ-diffusion model for response time and accuracy. \\
          Psychonomic Bulletin & Review, 14(1), 3-22. \\
          https://doi.org/10.3758/BF03194023"
        ),
        glue(
          "Ch\u00e1vez De la Pe\u00f1a, A. F., & Vandekerckhove, J. (2025). An \\
          EZ Bayesian hierarchical drift diffusion model for response time and \\
          accuracy. Psychonomic Bulletin & Review, 32(6), 3067-3087. \\
          https://doi.org/10.3758/s13423-025-02729-y"
        )
      ),
      version = version,
      requirements = glue(
        "Provide aggregated statistics for each subject and condition that model parameters should vary over:", "\n\n",
        "  - Mean reaction times (mean_rt) in seconds", "\n",
        "  - Variance of reaction times (var_rt) in seconds", "\n",
        "  - Number of responses to the upper decision threshold (n_upper)", "\n",
        "  - Total number of trials used to calculate aggregated statistics (n_trials)"
      ),
      parameters = .ezdm_version_table[[version]][["parameters"]],
      links = .ezdm_version_table[[version]][["links"]],
      fixed_parameters = .ezdm_version_table[[version]][["fixed_parameters"]],
      default_priors = .ezdm_version_table[[version]][["priors"]],
      init_ranges = .ezdm_version_table[[version]][["init_ranges"]]
    ),
    class = c("bmmodel", "ezdm", paste0("ezdm_", version)),
    call = call
  )
  out <- set_links(out, links)
  out
}
# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_ezdm()$info

#' @title `r .model_ezdm()$name`
#' @name ezdm
#' @details `r model_docs(.model_ezdm(version = "4par"))`
#'
#'   In version "4par", a boundary reached fewer than twice, or without RT
#'   summaries (`NA`), enters the model through the response counts only.
#'   `bmm()` adds two columns to the data, `rt_used_upper` and
#'   `rt_used_lower`, which it sets to 0 for such a boundary and to 1
#'   otherwise, and replaces the summaries of such a boundary with a
#'   placeholder the likelihood never reads. A 0 already in these columns is
#'   kept, so that `update()` does not read the placeholders in `fit$data` as
#'   data.
#'
#'   `newdata` passed to `log_lik()`, `predict()` or `posterior_predict()`
#'   needs these columns too. Rows of `fit$data` have them. For raw data, set
#'   both to 1: a boundary reached fewer than twice or with `NA` summaries is
#'   then still left out of the likelihood. brms functions that compare the
#'   response with predictions, `predictive_error()` and
#'   `residuals(method = "posterior_predict")`, read the placeholders as data
#'   in cells whose `rt_used_upper` is 0, so leave those cells out first, for
#'   example with `newdata = subset(fit$data, rt_used_upper == 1)`.
#' @param mean_rt The names of the variable or variables (for 4par version) coding the mean reaction time in seconds in the data.
#' @param var_rt The names of the variable or variables (for 4par version) coding the variance of the reaction time in seconds in the data
#' @param n_upper The name of the variable coding the number of responses that hit the upper response threshold (typically the number of correct responses) in the data.
#' @param n_trials The name of the variable coding the number of trials that was used to calculated the aggregated statistics.
#' @param links A named list of links for the parameters, e.g.
#'   `links = list(bound = "softplus")`. For positive parameters
#'   (e.g. `bound`, `ndt`), "softplus" is available as an alternative to the
#'   default "log" link that grows linearly for large values and avoids the
#'   numerical blow-up of `exp()`. A name that is not a parameter of the model
#'   is an error, and a link that allows values the default link excludes
#'   (e.g. "identity" for a positive parameter) is a warning.
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
#' exp(brms::fixef(fit))
#' }
ezdm <- function(mean_rt, var_rt, n_upper, n_trials, links = NULL,
                 version = c("3par", "4par"), ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  .model_ezdm(
    mean_rt = mean_rt, var_rt = var_rt, n_upper = n_upper, n_trials = n_trials,
    links = links, version = version, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !

#' @export
check_data.ezdm <- function(model, data, formula) {
  # retrieve required variable names
  mean_rt <- model$resp_vars$mean_rt
  var_rt <- model$resp_vars$var_rt
  n_upper <- model$resp_vars$n_upper
  n_trials <- model$other_vars$n_trials


  # validate length of mean_rt and var_rt dependent on version
  if (model$version == "3par") {
    stopif(length(mean_rt) != 1, "mean_rt must be a single variable name.")
    stopif(length(var_rt) != 1, "var_rt must be a single variable name.")
  } else if (model$version == "4par") {
    stopif(length(mean_rt) != 2, "mean_rt must be a vector of two variable names: c(mean_rt_upper, mean_rt_lower).")
    stopif(length(var_rt) != 2, "var_rt must be a vector of two variable names: c(var_rt_upper, var_rt_lower).")
  } else {
    stop2("Unknown ezdm version: {model$version}. Supported versions are '3par' and '4par'.")
  }

  # check that all required variables exist in data
  required_vars <- c(mean_rt, var_rt, n_upper, n_trials)
  missing_vars <- setdiff(required_vars, colnames(data))
  stopif(
    length(missing_vars),
    "The following required variables are missing from the data: {collapse_comma(missing_vars)}"
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

  # the summaries of a 4par boundary the likelihood does not use can be
  # anything, or missing, so only the used ones are validated
  earlier <- intersect(.EZDM_RT_USED, names(data))
  stopif(
    model$version == "4par" && !all(unlist(data[earlier]) %in% c(0, 1, NA)),
    "The columns {collapse_comma(earlier)} are reserved for the ezdm 4par \\
    indicators of which boundaries have usable RT summaries and may hold \\
    only 0, 1 or NA. Rename them in your data."
  )
  rt_used <- if (model$version == "4par") {
    .ezdm_rt_used(data, mean_rt, var_rt, n_upper_values, n_trials_values)
  } else {
    TRUE
  }

  # check that mean RT values are plausible (warn if likely in milliseconds)
  # typical RTs in seconds are 0.2-3s; values > 10 suggest milliseconds
  mean_rt_values <- as.matrix(data[mean_rt])[rt_used %in% TRUE]
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
  var_rt_values <- as.matrix(data[var_rt])[rt_used %in% TRUE]
  stopif(
    any(var_rt_values <= 0, na.rm = TRUE),
    "Variance of RT must be positive. Found non-positive values in the data."
  )

  if (model$version == "4par") {
    # brms drops every row holding an NA, response counts included (#430), so
    # an unused boundary gets a placeholder instead. -1 rather than a plausible
    # value: the likelihood reads the indicator, never the placeholder, and a
    # negative variance would make it fail loudly if it ever did.
    unused <- !is.na(rt_used) & !rt_used
    for (k in 1:2) {
      data[[mean_rt[k]]][unused[, k]] <- -1
      data[[var_rt[k]]][unused[, k]] <- -1
    }
    data[.EZDM_RT_USED] <- list(as.integer(rt_used[, 1]), as.integer(rt_used[, 2]))
    n_sparse <- sum(rowSums(unused) > 0)
    # out of the cells brms keeps: a row with a missing count is dropped
    warnif(
      n_sparse > 0,
      "{n_sparse} of {sum(!is.na(rt_used[, 1]))} cells have fewer than two \\
      responses, or no RT summaries, at one or both boundaries. Such a \\
      boundary enters the model through the response counts only."
    )
  }

  NextMethod("check_data")
}

# The 4par boundaries whose RT summaries enter the likelihood: reached at least
# twice, as a sample variance needs two responses, and summarised. A row with a
# missing count stays NA, and brms drops it. update() re-checks fit$data, whose
# unused boundaries hold placeholders by then, so an indicator from an earlier
# pass keeps them unused. A missing indicator is not decided yet: brms dropped
# such rows before fitting, so they are raw cells appended to fit$data.
.ezdm_rt_used <- function(data, mean_rt, var_rt, n_upper, n_trials) {
  used <- cbind(n_upper, n_trials - n_upper) >= 2 &
    !is.na(as.matrix(data[mean_rt])) & !is.na(as.matrix(data[var_rt]))
  earlier <- .EZDM_RT_USED %in% names(data)
  previous <- as.matrix(data[.EZDM_RT_USED[earlier]])
  used[, earlier] <- used[, earlier] & (is.na(previous) | previous == 1)
  used[is.na(n_upper) | is.na(n_trials), ] <- NA
  unname(used)
}

.EZDM_RT_USED <- c("rt_used_upper", "rt_used_lower")

############################################################################# !
# Convert bmmformula to brmsformla methods                               ####
############################################################################# !

#' @export
bmf2bf.ezdm_3par <- function(model, formula) {
  # retrieve required response arguments
  mean_rt <- model$resp_vars$mean_rt
  var_rt <- model$resp_vars$var_rt
  n_upper <- model$resp_vars$n_upper
  n_trials <- model$other_vars$n_trials

  brms::bf(paste0(mean_rt, " | vreal(", var_rt, ") + vint(", n_upper, ") + trials(", n_trials, ") ~ 1"))
}

#' @export
bmf2bf.ezdm_4par <- function(model, formula) {
  # retrieve required response arguments
  mean_rt <- model$resp_vars$mean_rt
  var_rt <- model$resp_vars$var_rt
  n_upper <- model$resp_vars$n_upper
  n_trials <- model$other_vars$n_trials

  brms::bf(glue::glue(
    "{mean_rt[1]} | vreal({mean_rt[2]}, {var_rt[1]}, {var_rt[2]}) + ",
    "vint({n_upper}, {n_trials}, {.EZDM_RT_USED[1]}, {.EZDM_RT_USED[2]}) ~ 1"
  ))
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

# Stan functions of one ezdm version. The order matters: each chunk defines what
# the next one calls.
.ezdm_stan_functions <- function(version) {
  chunks <- c(
    "ezdm_series.stan", "ezdm_cumulants.stan",
    paste0("ezdm_", version, "_functions.stan")
  )
  sc_path <- system.file("stan_chunks", package = "bmm")
  paste(vapply(file.path(sc_path, chunks), read_lines2, character(1)), collapse = "\n")
}

#' @export
configure_model.ezdm_3par <- function(model, data, formula) {
  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)
  links <- model$links

  # construct the family & add to formula object
  formula$family <- brms::custom_family(
    "ezdm_3par",
    dpars = c("mu", "drift", "bound", "ndt", "s"),
    links = c("identity", links$drift, links$bound, links$ndt, links$s),
    lb = c(NA, NA, 0, 0, 0),
    ub = c(NA, NA, NA, NA, NA),
    type = "real",
    log_lik = log_lik_ezdm_3par,
    posterior_predict = posterior_predict_ezdm_3par,
    loop = TRUE,
    vars = c("vreal1[n]", "vint1[n]", "trials[n]")
  )

  # prepare initial stanvars to pass to brms, model formula and priors
  stanvars <- brms::stanvar(scode = .ezdm_stan_functions("3par"), block = "functions")

  # return the list
  nlist(formula, data, stanvars)
}


log_lik_ezdm_3par <- function(i, prep) {
  dezdm(
    mean_rt = prep$data$Y[i],
    var_rt = prep$data$vreal1[i],
    n_upper = prep$data$vint1[i],
    n_trials = prep$data$trials[i],
    drift = brms::get_dpar(prep, "drift", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    s = brms::get_dpar(prep, "s", i = i),
    version = "3par",
    log = TRUE
  )
}

# brms forwards posterior_predict() dots to prepare_predictions() only, never
# to this function, so the other observables cannot be selected here; use
# pp_check(fit, resp_var = ...) instead (#401)
posterior_predict_ezdm_3par <- function(i, prep, ...) {
  rezdm(
    n = length(brms::get_dpar(prep, "drift", i = i)),
    n_trials = prep$data$trials[i],
    drift = brms::get_dpar(prep, "drift", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    s = brms::get_dpar(prep, "s", i = i),
    version = "3par"
  )[["mean_rt"]]
}

#' @export
configure_model.ezdm_4par <- function(model, data, formula) {
  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)
  links <- model$links

  # construct the family & add to formula object
  formula$family <- brms::custom_family(
    "ezdm_4par",
    dpars = c("mu", "drift", "bound", "ndt", "zr", "s"),
    links = c("identity", links$drift, links$bound, links$ndt, links$zr, links$s),
    lb = c(NA, NA, 0, 0, 0, 0), # lower bounds for parameters
    ub = c(NA, NA, NA, NA, 1, NA), # upper bounds for parameters
    type = "real", # real for continous dv, int for discrete dv
    log_lik = log_lik_ezdm_4par,
    posterior_predict = posterior_predict_ezdm_4par,
    loop = TRUE, # is the likelihood vectorized
    vars = c("vreal1[n]", "vreal2[n]", "vreal3[n]", "vint1[n]", "vint2[n]",
             "vint3[n]", "vint4[n]")
  )

  # prepare initial stanvars to pass to brms, model formula and priors
  stanvars <- brms::stanvar(scode = .ezdm_stan_functions("4par"), block = "functions")

  # return the list
  nlist(formula, data, stanvars)
}


log_lik_ezdm_4par <- function(i, prep) {
  # based on bmf2bf.ezdm_4par formula:
  # Y = mean_rt_upper, vreal1 = mean_rt_lower
  # vreal2 = var_rt_upper, vreal3 = var_rt_lower
  # vint3, vint4 = whether the upper, lower summaries are used; where they are
  # not, they hold placeholders, and NA tells dezdm() to leave them out
  used <- c(prep$data$vint3[i], prep$data$vint4[i]) == 1
  dezdm(
    mean_rt = ifelse(used, c(prep$data$Y[i], prep$data$vreal1[i]), NA),
    var_rt = ifelse(used, c(prep$data$vreal2[i], prep$data$vreal3[i]), NA),
    n_upper = prep$data$vint1[i],
    n_trials = prep$data$vint2[i],
    drift = brms::get_dpar(prep, "drift", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    zr = brms::get_dpar(prep, "zr", i = i),
    s = brms::get_dpar(prep, "s", i = i),
    version = "4par",
    log = TRUE
  )
}

# no dv argument: see posterior_predict_ezdm_3par()
posterior_predict_ezdm_4par <- function(i, prep, ...) {
  rezdm(
    n = length(brms::get_dpar(prep, "drift", i = i)),
    n_trials = prep$data$vint2[i],
    drift = brms::get_dpar(prep, "drift", i = i),
    bound = brms::get_dpar(prep, "bound", i = i),
    ndt = brms::get_dpar(prep, "ndt", i = i),
    zr = brms::get_dpar(prep, "zr", i = i),
    s = brms::get_dpar(prep, "s", i = i),
    version = "4par"
  )[["mean_rt_upper"]]
}

############################################################################# !
# PP_CHECK OBSERVABLES                                                    ####
############################################################################# !

# Each ezdm observation is one design cell, and the checks compare the
# distribution of a statistic across cells, as the ddm checks do across trials;
# type = "intervals" gives the per-cell view instead. mean_pc rather than raw
# n_upper: n_trials varies across cells, so counts are not comparable between
# observations while proportions are. A proportion is binned rather than
# smoothed because it is bounded and takes n_trials + 1 values only.
.pp_ezdm_accuracy <- function() {
  .pp_observable(function(d) d$n_upper / d$n_trials,
                 label = "Proportion of upper responses", type = "bars_binned")
}

#' @export
pp_observables.ezdm_3par <- function(model) {
  list(
    observed = c(mean_rt = "Y", var_rt = "vreal1", n_upper = "vint1",
                 n_trials = "trials"),
    checks = list(
      mean_rt = .pp_observable(function(d) d$mean_rt,
                               label = "Mean response time"),
      var_rt = .pp_observable(function(d) d$var_rt,
                              label = "Response time variance"),
      mean_pc = .pp_ezdm_accuracy()
    )
  )
}

#' @export
pp_simulate.ezdm_3par <- function(model, prep) {
  .pp_simulate_joint(prep, .rezdm_3par, c("drift", "bound", "ndt", "s"),
                     n_trials = rep(prep$data$trials, each = prep$ndraws))
}

# The summaries of a boundary check_data() found unused are placeholders
# (#430), so the RT checks turn them into NA, which the reduction drops.
# Multiplying by 1 or NA leaves the other values exact and yrep a matrix. The
# simulated summaries are masked by the observed indicator too, but only in
# cells whose observed value is dropped anyway.
.pp_ezdm_boundary <- function(x, used) {
  x * ifelse(used == 0, NA, 1)
}

# a fit saved before the indicators existed has no placeholders either, so its
# indicators default to 1
#' @export
pp_observables.ezdm_4par <- function(model) {
  list(
    observed = c(mean_rt_upper = "Y", mean_rt_lower = "vreal1",
                 var_rt_upper = "vreal2", var_rt_lower = "vreal3",
                 n_upper = "vint1", n_trials = "vint2",
                 rt_used_upper = "vint3", rt_used_lower = "vint4"),
    defaults = c(rt_used_upper = 1L, rt_used_lower = 1L),
    y_placeholders = function(data) {
      any(data[[.EZDM_RT_USED[1]]] == 0, na.rm = TRUE)
    },
    checks = list(
      mean_rt_upper = .pp_observable(
        function(d) .pp_ezdm_boundary(d$mean_rt_upper, d$rt_used_upper),
        label = "Mean RT (upper responses)"
      ),
      mean_rt_lower = .pp_observable(
        function(d) .pp_ezdm_boundary(d$mean_rt_lower, d$rt_used_lower),
        label = "Mean RT (lower responses)"
      ),
      var_rt_upper = .pp_observable(
        function(d) .pp_ezdm_boundary(d$var_rt_upper, d$rt_used_upper),
        label = "RT variance (upper responses)"
      ),
      var_rt_lower = .pp_observable(
        function(d) .pp_ezdm_boundary(d$var_rt_lower, d$rt_used_lower),
        label = "RT variance (lower responses)"
      ),
      mean_pc = .pp_ezdm_accuracy()
    )
  )
}

#' @export
pp_simulate.ezdm_4par <- function(model, prep) {
  .pp_simulate_joint(prep, .rezdm_4par, c("drift", "bound", "ndt", "zr", "s"),
                     n_trials = rep(prep$data$vint2, each = prep$ndraws))
}
