############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# Meanlog defaults for the simple correct/error parameters. Per-category
# parameters discovered in the custom version default to the error-accumulator
# values. Note: these are locations of the finishing-time distribution (higher
# = slower), not drift rates — hence the names differ from lba/rdm's
# driftc/drifte on purpose.
.lnr_meanlog_spec <- list(
  desc = "meanlog",
  link = "identity",
  priors = list(
    correct = list(main = "normal(-1, 0.5)", effects = "normal(0, 0.3)",
                   sd = "exponential(2)"),
    error = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)",
                 sd = "exponential(2)")
  ),
  inits = list(correct = c(-1.5, -0.5), error = c(-0.5, 0.5))
)

# The ndt/s block shared by both versions, declared once.
.lnr_shared <- list(
  parameters = list(
    ndt = "non-decision time",
    s = "sdlog (shared across accumulators)"
  ),
  links = list(ndt = "log", s = "log"),
  priors = list(
    ndt = list(main = "normal(-1.5, 0.5)", effects = "normal(0, 0.3)",
               sd = "exponential(4)"),
    s = list(main = "normal(0, 0.5)", effects = "normal(0, 0.2)",
             sd = "exponential(4)")
  ),
  inits = list(mu = c(-0.5, 0.5), ndt = c(0.025, 0.05), s = c(0.8, 1.2))
)

# Compose the spec for one version. Meanlog parameters precede the shared
# block: downstream code recovers accumulator names via
# setdiff(names(parameters), c("ndt", "s")), which relies on order.
.lnr_model_spec <- function(version) {
  fixed_parameters <- list(mu = 0)

  if (version == "custom") {
    return(nlist(
      parameters = .lnr_shared$parameters,
      links = .lnr_shared$links,
      fixed_parameters,
      priors = .lnr_shared$priors,
      init_ranges = .lnr_shared$inits
    ))
  }

  nlist(
    parameters = c(
      list(
        correct = paste(.lnr_meanlog_spec$desc, "for correct accumulator"),
        error = paste(.lnr_meanlog_spec$desc, "for error accumulators")
      ),
      .lnr_shared$parameters
    ),
    links = c(
      list(correct = .lnr_meanlog_spec$link, error = .lnr_meanlog_spec$link),
      .lnr_shared$links
    ),
    fixed_parameters,
    priors = c(.lnr_meanlog_spec$priors, .lnr_shared$priors),
    init_ranges = c(.lnr_shared$inits, .lnr_meanlog_spec$inits)
  )
}

# the parameters every version has, which a custom version's response
# categories cannot be named after
.lnr_shared_pars <- c("ndt", "s")

# the identifiers .lnr_stan_code() declares, which a category would shadow
# (measured: each made a category and the program handed to stanc)
.lnr_stan_names <- c("rt", "response", "t", "m", "n", "win", "lp", "reps", "i", "j")


.model_lnr <- function(
    rt = NULL,
    response = NULL,
    n_choices = NULL,
    accumulators = NULL,
    links = NULL,
    version = "simple",
    call = NULL,
    ...) {
  vt <- .lnr_model_spec(version)
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = nlist(n_choices, accumulators),
      domain = "Decision Making / Response times",
      task = "Choice Reaction Time tasks (multi-alternative)",
      name = "Log-Normal Race Model",
      citation = glue(
        "Rouder, J. N., Province, J. M., Morey, R. D., Gomez, P., & Heathcote, A. \\
        (2015). The Lognormal Race: A Cognitive-Process Model of Choice and \\
        Latency with Desirable Psychometric Properties. Psychometrika, 80(2), \\
        491-513. https://doi.org/10.1007/s11336-013-9396-3"
      ),
      version = version,
      requirements = glue(
        "- Reaction times should be passed in seconds", "\n",
        "- For version 'simple': response variable should be integer-coded ",
        "(1 = correct, 2:K = errors)", "\n",
        "- For version 'custom': response variable should contain character ",
        "labels matching formula parameter names"
      ),
      parameters = vt[["parameters"]],
      links = vt[["links"]],
      fixed_parameters = vt[["fixed_parameters"]],
      default_priors = vt[["priors"]],
      init_ranges = vt[["init_ranges"]]
    ),
    class = c("bmmodel", "lnr", paste0("lnr_", version)),
    call = call
  )

  set_links(out, links)
}

# the accumulator parameters of the custom version are the response categories
# of the user's formula, so at construction there is no set of names to check a
# link target against (check_model.lnr_custom fills in the default link for a
# category the user left alone). The simple version's parameters are fixed.
#' @exportS3Method
settable_links.lnr <- function(model) {
  if (model$version == "custom") NULL else names(model$links)
}

#' @title `r .model_lnr()$name`
#' @name lnr
#' @details `r model_docs(.model_lnr())`
#' @param rt The name of the variable in the dataset containing the response
#'   times. Response times should be coded in seconds (not milliseconds).
#' @param response The name of the variable in the dataset containing the
#'   response/choice. For the `"simple"` version, responses should be
#'   integer-coded: 1 = correct response, 2 through K = error responses.
#'   Factor and character-digit responses are accepted and converted
#'   automatically. For the `"custom"` version, responses should be character
#'   or factor labels matching the accumulator names in the formula. Category
#'   names must not use reserved internal parameter names such as `"mu"`,
#'   `"ndt"`, or `"s"`, must not contain an underscore and must not end in a
#'   number.
#' @param n_choices An integer specifying the total number of response
#'   alternatives (K >= 2). Required for `version = "simple"`. Not used for
#'   `version = "custom"` (inferred from the formula).
#' @param accumulators For `version = "custom"` only. A named vector
#'   specifying the number of racing accumulators per response category.
#'   Can be a named integer vector of positive counts for constant numbers of
#'   accumulators (e.g., `c(correct = 1, other = 3, npl = 5)`) or a named
#'   character vector of column names whose values are positive integers for
#'   trial-varying counts (e.g., `c(correct = "n_corr", other = "n_other",
#'   npl = "n_npl")`). If omitted, defaults to 1 accumulator per category.
#' @param version A character string specifying which version of the LNR model
#'   to use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): Two meanlog parameters — `correct` for the
#'       correct accumulator (response = 1) and `error` for all error
#'       accumulators, plus one sdlog `s` shared by all of them. This covers
#'       the common case where interest is in the speed of correct vs. error
#'       processing.
#'     \item `"custom"`: Per-category meanlog parameters. Response categories
#'       are defined by the formula LHS names (e.g., `correct ~ 1, other ~ 1,
#'       npl ~ 1`). The response column must contain character labels matching
#'       these names. Category names must not be `"mu"`, `"ndt"`, `"s"`,
#'       Stan reserved words, names containing an underscore, or names ending
#'       in a number, because `brms` and `Stan` use them as identifiers of
#'       their own. Supports per-category `accumulators`.
#'   }
#' @param links A named list of link functions for the model parameters.
#'   For `"simple"`: parameters are `correct`, `error`, `ndt`, and `s`.
#'   Default links are "identity" for `correct` and `error`, "log" for
#'   `ndt` and `s`. For `"custom"`: category parameters default to "identity".
#' @param ... Additional arguments passed internally (for testing purposes).
#' @return An object of class `bmmodel`
#' @section Default behavior:
#' All four parameters are estimated. `s` is the sdlog shared by the
#' accumulators; it is not a scale constraint, because finishing times are
#' observed in seconds and the likelihood is not invariant to rescaling them.
#' Fixing it costs accuracy: on 2000 trials simulated with sdlog = 0.5, holding
#' it at 1 pushed `ndt` onto the fastest observed response time and left the
#' 95% credible intervals of `correct`, `error` and `ndt` all clear of their
#' true values. To fix it anyway — for instance to reproduce a published fit
#' that did — write `s = 0` in `bmf()`. The constant is read on the parameter's
#' own link scale, so `s = 0` means sdlog = 1.
#' @section Default priors:
#' `s` has `normal(0, 0.5)` on the log link, a median sdlog of 1 with a 95%
#' interval of roughly 0.37 to 2.7. That covers the values simulated data are
#' fitted at without contributing measurably to the posterior at realistic
#' trial counts, and is tighter than the `normal(0, 1)` that `EMC2` puts on the
#' same quantity.
#'
#' `ndt` has `normal(-1.5, 0.5)` on the log link, the prior the **ddm** model
#' uses, rather than the tighter one the two-boundary response-time models in
#' `bmm` were first tuned with: multi-alternative choice tasks carry longer
#' non-decision times than those models were calibrated on, and a prior whose
#' upper bound sits below the true value drags the meanlogs with it.
#'
#' The group-level standard deviations get one rate per kind of parameter
#' rather than one per model, shared with the other racing models, so that the
#' same quantity is given the same prior wherever it appears. A meanlog gets
#' `exponential(2)` (median 0.35) and the two log-link parameters `ndt` and `s`
#' get `exponential(4)` (median 0.17). Those medians come from hierarchical
#' fits of simulated multi-subject data, which put the between-subject standard
#' deviation of a meanlog at 0.26 to 0.43 and of `ndt` and `s` at 0.17 to 0.19.
#' @section Likelihood and `loo()`:
#' For `version = "simple"` the likelihood is the probability of the *category*
#' that was observed, not of one named accumulator: an error trial contributes
#' the density of "some error accumulator won", which for K alternatives adds
#' log(K - 1) per error trial relative to a per-accumulator likelihood. The
#' term does not depend on any parameter, so posteriors are unaffected, but
#' `loo()` and `waic()` values are shifted by that constant and are not
#' comparable with a `version = "custom"` fit, which names each accumulator, or
#' with another package that does.
#' @note Both versions describe the same response type (a categorical winner in
#'   a choice-RT race), so they live in one constructor rather than separate
#'   model functions: `"simple"` is an accuracy-coded convenience layer (correct
#'   vs. error) over the general per-accumulator case handled by `"custom"`.
#' @export
#' @keywords bmmodel
#' @seealso [dlnr()] and [rlnr()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simple version with 2 alternatives
#' dat <- rlnr(n = 500, m = c(-1, 0), s = c(1, 1), ndt = 0.2)
#' model <- lnr(rt = "rt", response = "response", n_choices = 2)
#' formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1, s ~ 1)
#' fit <- bmm(formula, dat, model, cores = 4, backend = "cmdstanr")
#'
#' # hold the sdlog at 1 instead of estimating it (s = 0 on the log link)
#' fit_fixed <- bmm(bmf(correct ~ 1, error ~ 1, ndt ~ 1, s = 0), dat, model,
#'                  cores = 4, backend = "cmdstanr")
#'
#' # custom version with named categories
#' model2 <- lnr(rt = "rt", response = "resp", version = "custom",
#'               accumulators = c(correct = 1, similar = 3, other = 5))
#' formula2 <- bmf(correct ~ 1, similar ~ 1, other ~ 1, ndt ~ 1)
lnr <- function(rt, response, n_choices = NULL,
                version = c("simple", "custom"),
                accumulators = NULL, links = NULL, ...) {
  call <- match.call()
  dots <- list(...)
  if ("n_alternatives" %in% names(dots)) {
    n_choices <- dots$n_alternatives
    warning2("The argument 'n_alternatives' is deprecated. Please use 'n_choices' instead.")
  }
  if ("num_alternatives" %in% names(dots)) {
    accumulators <- dots$num_alternatives
    warning2("The argument 'num_alternatives' is deprecated. Please use 'accumulators' instead.")
  }
  stop_missing_args()
  version <- match.arg(version)
  if (version == "simple") {
    stopif(
      !is.null(accumulators),
      "accumulators is only supported for version 'custom'."
    )
    stopif(
      is.null(n_choices) || !is.numeric(n_choices) ||
        n_choices < 2 || n_choices != round(n_choices),
      "n_choices must be an integer >= 2 for version 'simple'."
    )
    n_choices <- as.integer(n_choices)
  } else {
    stopif(
      !is.null(n_choices),
      "n_choices is only supported for version 'simple'. Use accumulators for version 'custom'."
    )
    n_choices <- NULL
  }
  .model_lnr(
    rt = rt,
    response = response,
    n_choices = n_choices,
    accumulators = accumulators,
    links = links,
    version = version,
    call = call
  )
}

############################################################################# !
# CHECK_MODEL S3 methods                                                 ####
############################################################################# !

#' @export
check_model.lnr_custom <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    cat_pars <- race_category_names(formula, .lnr_shared_pars, .lnr_stan_names)

    for (p in cat_pars) {
      model$parameters[[p]] <- paste0(.lnr_meanlog_spec$desc, " for '", p, "' accumulator")
      if (is.null(model$links[[p]])) model$links[[p]] <- .lnr_meanlog_spec$link
      if (is.null(model$default_priors[[p]])) {
        model$default_priors[[p]] <- .lnr_meanlog_spec$priors$error
      }
      if (is.null(model$init_ranges[[p]])) {
        model$init_ranges[[p]] <- .lnr_meanlog_spec$inits$error
      }
    }

    model$other_vars$resp_cats <- cat_pars
  }

  NextMethod("check_model")
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !

#' @export
check_data.lnr <- function(model, data, formula) {
  race_check_data(model, data)
  NextMethod("check_data")
}

#' @export
check_data.lnr_simple <- function(model, data, formula) {
  response_var <- model$resp_vars$response
  n_alt <- model$other_vars$n_choices
  data <- race_code_simple_response(model, data, "lnr")

  chosen <- unique(data[, response_var])
  # set by revert_check_data() on a fit's stored frame, which codes every error 2
  if (isTRUE(attr(data, "lnr_errors_pooled")) && 2L %in% chosen) {
    chosen <- seq_len(n_alt)
  }
  attr(data, "lnr_errors_pooled") <- NULL
  never_chosen <- setdiff(seq_len(n_alt), chosen)
  warnif(
    length(never_chosen) > 0,
    "Response option(s) {collapse_comma(never_chosen)} never occur in \\
    '{response_var}'. The model still estimates a meanlog for them, informed \\
    only by the trials they lost."
  )

  NextMethod("check_data")
}

#' @export
check_data.lnr_custom <- function(model, data, formula) {
  data <- race_code_custom_response(model, data, "lnr", .lnr_shared_pars)
  NextMethod("check_data")
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.lnr_simple <- function(model, formula) {
  race_bf(model, "lnr", 2)
}

#' @export
bmf2bf.lnr_custom <- function(model, formula) {
  race_bf(model, "lnr", length(model$other_vars$resp_cats))
}

############################################################################# !
# Stan code generation                                                   ####
############################################################################# !

.lnr_stan_code <- function(family_name, cat_names) {
  n_cats <- length(cat_names)
  cat_args <- paste(paste0("vector ", cat_names), collapse = ", ")
  n_args <- paste(paste0("array[] int n", seq_len(n_cats)), collapse = ", ")
  m_assignments <- paste(
    paste0("m[", seq_len(n_cats), "] = ", cat_names, "[i];"),
    collapse = "\n    "
  )
  n_assignments <- paste(
    paste0("n[", seq_len(n_cats), "] = n", seq_len(n_cats), "[i];"),
    collapse = "\n    "
  )

  glue(
    "real {family_name}_lpdf(vector rt, vector mu, {cat_args}, ",
    "vector ndt, vector s, array[] int response, {n_args}) {{\n",
    "  int N = num_elements(rt);\n",
    "  real log_lik = 0;\n",
    "  for (i in 1:N) {{\n",
    "    real t = rt[i] - ndt[i];\n",
    "    array[{n_cats}] real m;\n",
    "    array[{n_cats}] int n;\n",
    "    int win;\n",
    "    real lp;\n",
    "    if (t <= 0) return negative_infinity();\n",
    "    {m_assignments}\n",
    "    {n_assignments}\n",
    "    win = response[i];\n",
    "    lp = log(n[win]) + lognormal_lpdf(t | m[win], s[i]);\n",
    "    for (j in 1:{n_cats}) {{\n",
    "      int reps = (j == win) ? n[j] - 1 : n[j];\n",
    "      if (reps > 0) {{\n",
    "        lp += reps * lognormal_lccdf(t | m[j], s[i]);\n",
    "      }}\n",
    "    }}\n",
    "    log_lik += lp;\n",
    "  }}\n",
    "  return log_lik;\n",
    "}}"
  )
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.lnr_simple <- function(model, data, formula) {
  links <- model$links
  cat_names <- c("correct", "error")
  formula <- bmf2bf(model, formula)

  dpars <- c("mu", cat_names, "ndt", "s")
  link_vec <- c(
    "identity", links$correct, links$error, links$ndt, links$s
  )

  formula$family <- brms::custom_family(
    "lnr_simple",
    dpars = dpars,
    links = link_vec,
    ub = rep(NA, length(dpars)),
    lb = c(NA, NA, NA, 0, 0),
    type = "real",
    vars = race_family_vars(length(cat_names)),
    loop = FALSE,
    log_lik = log_lik_lnr_simple,
    posterior_predict = posterior_predict_lnr_simple,
    posterior_epred = posterior_epred_lnr_simple
  )

  stan_code <- .lnr_stan_code("lnr_simple", cat_names)
  stanvars <- brms::stanvar(scode = stan_code, block = "functions")

  nlist(formula, data, stanvars)
}

#' @export
configure_model.lnr_custom <- function(model, data, formula) {
  links <- model$links
  cat_names <- model$other_vars$resp_cats
  n_cats <- length(cat_names)
  formula <- bmf2bf(model, formula)

  dpars <- c("mu", cat_names, "ndt", "s")
  link_vec <- c(
    "identity",
    vapply(cat_names, function(p) links[[p]], character(1)),
    links$ndt, links$s
  )
  lb_vec <- c(NA, rep(NA, n_cats), 0, 0)

  formula$family <- brms::custom_family(
    "lnr_custom",
    dpars = dpars,
    links = link_vec,
    ub = rep(NA, length(dpars)),
    lb = lb_vec,
    type = "real",
    vars = race_family_vars(n_cats),
    loop = FALSE,
    log_lik = log_lik_lnr_custom,
    posterior_predict = posterior_predict_lnr_custom,
    posterior_epred = posterior_epred_lnr_custom
  )

  stan_code <- .lnr_stan_code("lnr_custom", cat_names)
  stanvars <- brms::stanvar(scode = stan_code, block = "functions")

  nlist(formula, data, stanvars)
}

############################################################################# !
# Post-processing functions (shared helpers)                             ####
############################################################################# !

.lnr_log_lik <- function(i, prep, cat_names, n_cats) {
  rt <- prep$data$Y[i]
  response <- prep$data$vint1[i]
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  t <- rt - ndt
  t[t <= 0] <- NA

  n_cat <- race_counts(prep, i, n_cats)

  m_win <- brms::get_dpar(prep, cat_names[response], i = i)
  log_lik <- log(n_cat[response]) +
    stats::dlnorm(t, meanlog = m_win, sdlog = s, log = TRUE)

  if (n_cat[response] > 1) {
    log_lik <- log_lik + (n_cat[response] - 1) *
      stats::plnorm(t, meanlog = m_win, sdlog = s,
                    lower.tail = FALSE, log.p = TRUE)
  }

  for (j in seq_len(n_cats)) {
    if (j == response || n_cat[j] == 0) next
    m_j <- brms::get_dpar(prep, cat_names[j], i = i)
    log_lik <- log_lik + n_cat[j] *
      stats::plnorm(t, meanlog = m_j, sdlog = s,
                    lower.tail = FALSE, log.p = TRUE)
  }

  log_lik[is.na(log_lik)] <- -Inf
  log_lik
}

# The per-draw parameters of one observation as the row-per-draw matrices
# lnr_race() and the survivor integrator take. get_dpar() returns a scalar for a
# dpar brms stores fixed (s under `s = 0`), so every vector is grown to ndraws,
# and the meanlogs go through matrix() because vapply() returns a bare vector
# when there is a single draw.
.lnr_draw_pars <- function(i, prep, cat_names, n_cats) {
  n_draws <- prep$ndraws
  list(
    ndt = rep_len(brms::get_dpar(prep, "ndt", i = i), n_draws),
    s = rep_len(brms::get_dpar(prep, "s", i = i), n_draws),
    m = matrix(vapply(cat_names, function(p) {
      rep_len(brms::get_dpar(prep, p, i = i), n_draws)
    }, numeric(n_draws)), nrow = n_draws),
    counts = race_counts(prep, i, n_cats)
  )
}

.lnr_posterior_predict <- function(i, prep, cat_names, n_cats, ...) {
  d <- .lnr_draw_pars(i, prep, cat_names, n_cats)
  race <- lnr_race(
    m = d$m,
    s = matrix(d$s, prep$ndraws, n_cats),
    counts = matrix(d$counts, prep$ndraws, n_cats, byrow = TRUE)
  )
  race$rt + d$ndt
}

# E[RT] = ndt + int_0^inf prod_j S_j(t)^n_j dt. A Monte-Carlo estimate of this
# integral moved by several percent between two calls on the same draws, which
# reached the user as noise on conditional_effects(); the grid is deterministic.
.lnr_posterior_epred <- function(prep, cat_names, n_cats) {
  epred <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs)
  for (i in seq_len(prep$nobs)) {
    d <- .lnr_draw_pars(i, prep, cat_names, n_cats)
    racing <- which(d$counts > 0)

    log_surv <- function(t) {
      out <- matrix(0, prep$ndraws, length(t))
      for (j in racing) {
        out <- out + d$counts[j] * lnorm_log_surv(t, d$m[, j], d$s)
      }
      out
    }

    # the race is over once the slowest single accumulator is, and it has
    # barely started before the fastest one could finish, so the two extreme
    # meanlogs bracket the product survivor for every draw
    m_lo <- matrixStats::rowMins(d$m[, racing, drop = FALSE])
    m_hi <- matrixStats::rowMaxs(d$m[, racing, drop = FALSE])
    t_hi <- max(stats::qlnorm(1e-12, m_hi, d$s, lower.tail = FALSE))
    t_lo <- max(min(stats::qlnorm(1e-12, m_lo, d$s)), t_hi * 1e-12)

    epred[, i] <- d$ndt + race_expected_time(log_surv, t_lo, t_hi)
  }
  epred
}

log_lik_lnr_simple <- function(i, prep) {
  .lnr_log_lik(i, prep, cat_names = c("correct", "error"), n_cats = 2)
}

posterior_predict_lnr_simple <- function(i, prep, ...) {
  .lnr_posterior_predict(i, prep, cat_names = c("correct", "error"),
                         n_cats = 2, ...)
}

log_lik_lnr_custom <- function(i, prep) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "ndt", "s")
  )
  .lnr_log_lik(i, prep, cat_names = cat_names, n_cats = length(cat_names))
}

posterior_predict_lnr_custom <- function(i, prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "ndt", "s")
  )
  .lnr_posterior_predict(i, prep, cat_names = cat_names,
                         n_cats = length(cat_names), ...)
}

posterior_epred_lnr_simple <- function(prep) {
  .lnr_posterior_epred(prep, cat_names = c("correct", "error"), n_cats = 2)
}

posterior_epred_lnr_custom <- function(prep) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "ndt", "s")
  )
  .lnr_posterior_epred(prep, cat_names = cat_names, n_cats = length(cat_names))
}

############################################################################# !
# PP_CHECK OBSERVABLES                                                    ####
############################################################################# !

#' @export
pp_observables.lnr <- function(model) {
  race_pp_observables(model)
}

# One method for both versions: the accumulator dpars carry their own order, so
# nothing here depends on whether they came from n_choices or from the formula
#' @export
pp_simulate.lnr <- function(model, prep) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "ndt", "s"))
  n_cats <- length(cat_names)
  n_row <- prep$ndraws * prep$nobs

  race <- lnr_race(
    m = matrix(vapply(cat_names, .pp_dpar_vector, numeric(n_row), prep = prep),
               nrow = n_row),
    s = matrix(.pp_dpar_vector(prep, "s"), n_row, n_cats),
    counts = matrix(vapply(seq_len(n_cats), function(j) {
      rep(prep$data[[paste0("vint", j + 1)]], each = prep$ndraws)
    }, integer(n_row)), nrow = n_row)
  )

  list(
    rt = matrix(race$rt + .pp_dpar_vector(prep, "ndt"), nrow = prep$ndraws),
    response = matrix(race$response, nrow = prep$ndraws)
  )
}

#' @exportS3Method
revert_check_data.lnr_simple <- function(model, data) {
  # .lnr_cat records only whether a trial was correct, so a rebuilt response
  # codes every error as option 2, and check_data() learns from the attribute
  # that the other error options are pooled rather than never chosen
  pooled <- not_in(model$resp_vars$response, colnames(data))
  data <- race_revert_check_data(model, data, "lnr")
  if (pooled) attr(data, "lnr_errors_pooled") <- TRUE
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.lnr_custom <- function(model, data) {
  data <- race_revert_check_data(model, data, "lnr")
  NextMethod("revert_check_data")
}

#' @exportS3Method
response_annotations.lnr <- function(model) {
  list(rt = "seconds")
}

#' @exportS3Method
model_column_roles.lnr <- function(model) {
  if (model$version == "custom") {
    return(c(
      rt = "response time in seconds",
      response = "chosen category, labelled by the accumulator names in the formula",
      accumulators = "number of accumulators per category (columns, or one number per category)"
    ))
  }
  c(
    rt = "response time in seconds",
    response = "chosen option, 1 = correct and 2 to `n_choices` = errors",
    n_choices = NA
  )
}
