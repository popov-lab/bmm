############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# Drift-rate defaults for the simple driftc/drifte parameters. Per-category
# parameters discovered in the custom version default to the error-accumulator
# (drifte) values.
.rdm_drift_spec <- list(
  desc = "drift rate",
  link = "log",
  priors = list(
    driftc = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)"),
    drifte = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)")
  ),
  inits = list(driftc = c(2, 4), drifte = c(1, 2.5))
)

# The gap/ndt/s/sp block shared by both versions, declared once.
.rdm_shared <- list(
  parameters = list(
    gap = "threshold gap (b = gap + sp)",
    ndt = "non-decision time",
    s = "diffusion constant",
    sp = "maximum starting point (uniform on 0 to sp)"
  ),
  links = list(gap = "log", ndt = "log", s = "log", sp = "log"),
  priors = list(
    gap = list(main = "normal(0, 0.3)", effects = "normal(0, 0.3)"),
    ndt = list(main = "normal(-2, 0.3)", effects = "normal(0, 0.3)"),
    s = list(main = "normal(0, 0.3)", effects = "normal(0, 0.2)"),
    sp = list(main = "normal(-1, 0.5)", effects = "normal(0, 0.3)")
  ),
  inits = list(
    mu = c(-0.5, 0.5), gap = c(0.8, 1.2), ndt = c(0.01, 0.05),
    s = c(0.8, 1.2), sp = c(0.2, 0.5)
  )
)

# sp = -100 on the log link (A = exp(-100), zero to double precision) is the
# constructor's "no start-point variability" default and the only value that
# takes the plain-Wald fast path in Stan. Any other constant a user writes in
# bmf() -- sp = log(0.3) -- is honoured as the start-point range it names, on
# both the Stan side and in every posterior method, which read A = sp as it is.
.rdm_sp_off <- -100

.rdm_start_var <- function(model) {
  !identical(model$fixed_parameters$sp, .rdm_sp_off)
}

# Compose the spec for one version. Drift parameters precede the shared block:
# downstream code recovers accumulator names via
# setdiff(names(parameters), c("gap", "ndt", "s", "sp")), which relies on order.
.rdm_model_spec <- function(version) {
  fixed_parameters <- list(mu = 0, s = 0, sp = .rdm_sp_off)

  if (version == "custom") {
    return(nlist(
      parameters = .rdm_shared$parameters,
      links = .rdm_shared$links,
      fixed_parameters,
      priors = .rdm_shared$priors,
      init_ranges = .rdm_shared$inits
    ))
  }

  nlist(
    parameters = c(
      list(
        driftc = paste(.rdm_drift_spec$desc, "for correct accumulator"),
        drifte = paste(.rdm_drift_spec$desc, "for error accumulators")
      ),
      .rdm_shared$parameters
    ),
    links = c(
      list(driftc = .rdm_drift_spec$link, drifte = .rdm_drift_spec$link),
      .rdm_shared$links
    ),
    fixed_parameters,
    priors = c(.rdm_drift_spec$priors, .rdm_shared$priors),
    init_ranges = c(.rdm_shared$inits, .rdm_drift_spec$inits)
  )
}

# Names the generated program already gives something else. The model's own
# dpars are there, and so is `Intercept`: brms declares `real Intercept` for the
# response's intercept, which collides with the `vector[N] Intercept` a category
# of that name would declare ("Identifier "Intercept" is already in use", stanc
# 2.40). Unlike .stan_reserved these are bmm's and brms's own names, not Stan's.
.rdm_reserved_dpars <- c("mu", "gap", "ndt", "s", "sp", "intercept")

.model_rdm <- function(
    rt = NULL,
    response = NULL,
    n_choices = NULL,
    accumulators = NULL,
    links = NULL,
    version = "simple",
    call = NULL,
    ...) {
  vt <- .rdm_model_spec(version)
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = nlist(n_choices, accumulators),
      domain = "Decision Making / Response times",
      task = "Choice Reaction Time tasks (multi-alternative)",
      name = "Racing Diffusion Model",
      citation = "Tillman, G., Van Zandt, T., & Logan, G. D. (2020).
        Sequential sampling models without random between-trial variability:
        the racing diffusion model of speeded decision making. Psychonomic
        Bulletin & Review, 27, 911-936. https://doi.org/10.3758/s13423-020-01719-6",
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
    class = c("bmmodel", "rdm", paste0("rdm_", version)),
    call = call
  )

  set_links(out, links)
}

# the accumulator parameters of the custom version are the response categories
# of the user's formula, so at construction there is no set of names to check a
# link target against (check_model.rdm_custom fills in the default link for a
# category the user left alone). The simple version's parameters are fixed.
#' @exportS3Method
settable_links.rdm <- function(model) {
  if (model$version == "custom") NULL else names(model$links)
}

# every parameter is positive and enters the Stan likelihood on the natural
# scale, so a wider link would reach print(), the initial values and the prior
# scale but not the sampler
#' @exportS3Method
settable_link_functions.rdm <- function(model) {
  "log"
}

#' @title `r .model_rdm()$name`
#' @name rdm
#' @details `r model_info(.model_rdm())`
#' @param rt The name of the variable in the dataset containing the response
#'   times. Response times should be coded in seconds (not milliseconds).
#' @param response The name of the variable in the dataset containing the
#'   response/choice. For the `"simple"` version, responses should be
#'   integer-coded: 1 = correct response, 2 through K = error responses.
#'   Factor and character-digit responses are accepted and converted
#'   automatically. For the `"custom"` version, responses should be character
#'   or factor labels matching the accumulator names in the formula.
#' @param n_choices An integer specifying the total number of response
#'   alternatives (K >= 2). Required for `version = "simple"`. Not used for
#'   `version = "custom"` (inferred from the formula).
#' @param accumulators For `version = "custom"` only. A named vector
#'   specifying the number of racing accumulators per response category.
#'   Can be a named integer vector for constant counts (e.g.,
#'   `c(correct = 1, other = 3, npl = 5)`) or a named character vector
#'   of column names for trial-varying counts (e.g.,
#'   `c(correct = "n_corr", other = "n_other", npl = "n_npl")`). If omitted,
#'   defaults to 1 accumulator per category.
#' @param version A character string specifying which version of the RDM to
#'   use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): Two drift parameters — `driftc` for the
#'       correct accumulator (response = 1) and `drifte` for all error
#'       accumulators. The diffusion constant `s` is shared and fixed by
#'       default (s = 1). The start-point range `sp` is fixed to zero by
#'       default; add `sp ~ 1` to the formula to estimate it.
#'     \item `"custom"`: Per-category drift parameters. Response categories
#'       are defined by the formula LHS names (e.g., `correct ~ 1, other ~ 1,
#'       npl ~ 1`). The response column must contain character labels matching
#'       these names. Supports per-category `accumulators`.
#'   }
#' @param links A named list of link functions for the model parameters.
#'   For `"simple"`: parameters are `driftc`, `drifte`, `gap`, `ndt`, `s`,
#'   and `sp`. All positive-valued parameters, including `ndt`, use a "log"
#'   link and only support that link.
#' @param ... Additional arguments passed internally (for testing purposes).
#' @return An object of class `bmmodel`
#' @details
#' ## Start-point variability (`sp`)
#'
#' Each accumulator starts at a point drawn uniformly from `[0, sp]` and
#' finishes when it has travelled `gap + sp` minus that start; `sp` is the
#' range of the start point and `gap` the distance from its top to the
#' threshold, so the threshold is `b = gap + sp`. By default `sp` is fixed to
#' `-100` on the log link, i.e. to zero: every accumulator starts at zero and
#' its finishing time is a plain Wald with threshold `gap`. Write `sp ~ 1` (or
#' a predictor) in `bmf()` to estimate it, or fix it to another range with a
#' constant on the log link, e.g. `bmf(..., sp = log(0.3))`, which is honoured
#' as a start point uniform on `[0, 0.3]` by the sampler and by `log_lik()`,
#' `posterior_predict()`, `posterior_epred()` and `pp_check()` alike.
#'
#' ## Likelihood convention
#'
#' `log_lik()`, and therefore `loo()` and `waic()`, use the likelihood of the
#' *category* that was observed, not of one named accumulator: the winner's
#' density carries `log(n_j)` for a category raced by `n_j` identical
#' accumulators. For `version = "simple"` the error category holds all K - 1
#' error accumulators, so every error trial sits log(K - 1) above a
#' per-response likelihood. That term does not depend on any parameter, so
#' posteriors and predictions are unaffected, but `loo()` and `waic()` values
#' are shifted by it and are not comparable with a `version = "custom"` fit,
#' which names each accumulator, or with another package that does.
#'
#' @note Both versions describe the same response type (a categorical winner in
#'   a choice-RT race), so they live in one constructor rather than separate
#'   model functions: `"simple"` is an accuracy-coded convenience layer (correct
#'   vs. error) over the general per-accumulator case handled by `"custom"`.
#' @export
#' @keywords bmmodel
#' @seealso [drdm()] and [rrdm()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simple version with 2 alternatives
#' dat <- rrdm(n = 500, drift = c(3, 1.5), gap = 1, sp = 0, ndt = 0.2)
#' model <- rdm(rt = "rt", response = "response", n_choices = 2)
#' formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
#' fit <- bmm(formula, dat, model, cores = 4, backend = "cmdstanr")
#'
#' # with starting point variability
#' formula2 <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, sp ~ 1)
#' fit2 <- bmm(formula2, dat, model, cores = 4, backend = "cmdstanr")
rdm <- function(rt, response, n_choices = NULL,
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
      is.null(n_choices) || !is.numeric(n_choices) ||
        n_choices < 2 || n_choices != round(n_choices),
      "n_choices must be an integer >= 2 for version 'simple'."
    )
    n_choices <- as.integer(n_choices)
  } else {
    n_choices <- NULL
  }
  .model_rdm(
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
check_model.rdm_custom <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    reserved_pars <- c("gap", "ndt", "s", "sp")
    formula_pars <- names(formula)
    cat_pars <- setdiff(formula_pars, reserved_pars)

    stopif(
      length(cat_pars) == 0,
      "Custom version requires at least one accumulator parameter in the formula."
    )

    bad_internal_names <- cat_pars[tolower(cat_pars) %in% .rdm_reserved_dpars]
    stopif(
      length(bad_internal_names) > 0,
      "Category names cannot use reserved internal parameter names: \\
      {collapse_comma(bad_internal_names)}."
    )

    bad_names <- intersect(cat_pars, .stan_reserved)
    stopif(
      length(bad_names) > 0,
      "Category names cannot be Stan reserved words: {collapse_comma(bad_names)}. \\
      Please rename the affected response categories."
    )

    bad_dpar_names <- cat_pars[grepl("[0-9]$", cat_pars)]
    stopif(
      length(bad_dpar_names) > 0,
      "Category names cannot end in a number because brms uses them as \\
      distributional parameters: {collapse_comma(bad_dpar_names)}. \\
      Please rename the affected response categories."
    )

    bad_underscore_names <- cat_pars[grepl("_", cat_pars, fixed = TRUE)]
    stopif(
      length(bad_underscore_names) > 0,
      "Category names cannot contain underscores because brms rejects them as \\
      distributional parameters: {collapse_comma(bad_underscore_names)}. \\
      Please rename the affected response categories."
    )

    for (p in cat_pars) {
      model$parameters[[p]] <- paste0(.rdm_drift_spec$desc, " for '", p, "' accumulator")
      if (is.null(model$links[[p]])) model$links[[p]] <- .rdm_drift_spec$link
      if (is.null(model$default_priors[[p]])) {
        model$default_priors[[p]] <- .rdm_drift_spec$priors$drifte
      }
      if (is.null(model$init_ranges[[p]])) {
        model$init_ranges[[p]] <- .rdm_drift_spec$inits$drifte
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
check_data.rdm <- function(model, data, formula) {
  rt_var <- model$resp_vars$rt
  response_var <- model$resp_vars$response

  stopif(
    not_in(rt_var, colnames(data)),
    "The RT variable '{rt_var}' is not present in the data."
  )

  n_na_rt <- sum(is.na(data[, rt_var]))
  stopif(
    n_na_rt > 0,
    "The RT variable '{rt_var}' contains {n_na_rt} NA values. \\
    Please remove or impute missing values before fitting the model."
  )

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
    stop2("The RT variable '{rt_var}' needs to be of type double or integer.")
  }

  stopif(
    not_in(response_var, colnames(data)),
    "The response variable '{response_var}' is not present in the data."
  )

  n_na_resp <- sum(is.na(data[, response_var]))
  stopif(
    n_na_resp > 0,
    "The response variable '{response_var}' contains {n_na_resp} NA values. \\
    Please remove or impute missing values before fitting the model."
  )

  NextMethod("check_data")
}

#' @export
check_data.rdm_simple <- function(model, data, formula) {
  response_var <- model$resp_vars$response
  n_alt <- model$other_vars$n_choices

  if (is.factor(data[, response_var]) || is.character(data[, response_var])) {
    labels <- as.character(data[, response_var])
    coded <- suppressWarnings(as.integer(labels))
    # as.integer() turns a label that is not a number into NA, and every check
    # below then compares against NA rather than refusing the label
    stopif(
      anyNA(coded),
      "The response variable '{response_var}' must be integer-coded 1:{n_alt} \\
      (1 = correct, 2:{n_alt} = errors) for version 'simple', but contains the \\
      non-numeric label(s) {collapse_comma(unique(labels[is.na(coded)]))}. Use \\
      version 'custom' for named response categories."
    )
    data[, response_var] <- coded
  }

  stopif(
    !is.numeric(data[, response_var]),
    "The response variable '{response_var}' must be numeric (integer-coded 1:{n_alt})."
  )

  resp_vals <- data[, response_var]
  stopif(
    any(resp_vals < 1) || any(resp_vals > n_alt) ||
      any(resp_vals != round(resp_vals)),
    "The response variable '{response_var}' must contain integers in 1:{n_alt}."
  )
  data[, response_var] <- as.integer(data[, response_var])

  data$.rdm_cat <- ifelse(data[, response_var] == 1L, 1L, 2L)
  data$.rdm_n1 <- 1L
  data$.rdm_n2 <- n_alt - 1L

  NextMethod("check_data")
}

#' @export
check_data.rdm_custom <- function(model, data, formula) {
  response_var <- model$resp_vars$response
  cat_names <- model$other_vars$resp_cats
  num_alt <- model$other_vars$accumulators
  n_cats <- length(cat_names)

  if (is.factor(data[, response_var])) {
    data[, response_var] <- as.character(data[, response_var])
  }

  stopif(
    !is.character(data[, response_var]),
    "For version 'custom', the response variable '{response_var}' must \\
    contain character labels matching the formula category names."
  )

  data_levels <- unique(data[, response_var])
  # a level named after one of the model's own parameters cannot have reached the
  # formula as a category, so it would otherwise be reported as unspecified
  bad_levels <- data_levels[tolower(data_levels) %in% .rdm_reserved_dpars]
  stopif(
    length(bad_levels) > 0,
    "Response levels cannot use reserved internal parameter names: \\
    {collapse_comma(bad_levels)}."
  )
  missing_in_formula <- setdiff(data_levels, cat_names)
  missing_in_data <- setdiff(cat_names, data_levels)
  stopif(
    length(missing_in_formula) > 0,
    "Response levels {collapse_comma(missing_in_formula)} in the data are \\
    not specified in the formula."
  )
  warnif(
    length(missing_in_data) > 0,
    "Formula categories {collapse_comma(missing_in_data)} are not present \\
    in the data."
  )

  data$.rdm_cat <- setNames(seq_along(cat_names), cat_names)[data[, response_var]]

  if (is.null(num_alt)) {
    for (i in seq_along(cat_names)) {
      data[[paste0(".rdm_n", i)]] <- 1L
    }
  } else if (is.numeric(num_alt)) {
    stopif(
      !all(cat_names %in% names(num_alt)),
      "accumulators must have names matching formula categories: \\
      {collapse_comma(cat_names)}"
    )
    counts <- num_alt[cat_names]
    invalid <- counts[!is.finite(counts) | counts < 1 | counts != round(counts)]
    stopif(
      length(invalid) > 0,
      "accumulators must contain a positive integer for each formula category. \\
      Invalid value(s): \\
      {collapse_comma(glue('{names(invalid)} = {invalid}'))}"
    )
    for (i in seq_along(cat_names)) {
      data[[paste0(".rdm_n", i)]] <- as.integer(counts[[i]])
    }
  } else if (is.character(num_alt)) {
    stopif(
      !all(cat_names %in% names(num_alt)),
      "accumulators must have names matching formula categories: \\
      {collapse_comma(cat_names)}"
    )
    missing_cols <- setdiff(num_alt, colnames(data))
    stopif(
      length(missing_cols) > 0,
      "accumulators columns {collapse_comma(missing_cols)} not found \\
      in the data."
    )
    # a trial-varying count may be zero for a category that sits that trial out
    # -- the likelihood and the simulator both skip it -- but never fractional
    # or negative, which Stan would turn into a NaN target
    for (i in seq_along(cat_names)) {
      col_name <- num_alt[cat_names[i]]
      col_vals <- data[, col_name]
      stopif(
        !is.numeric(col_vals),
        "accumulators column '{col_name}' must be numeric."
      )
      stopif(
        anyNA(col_vals) || any(!is.finite(col_vals)),
        "accumulators column '{col_name}' contains NA or non-finite values."
      )
      stopif(
        any(col_vals < 0 | col_vals != round(col_vals)),
        "accumulators column '{col_name}' must contain integers >= 0."
      )
      data[[paste0(".rdm_n", i)]] <- as.integer(col_vals)
    }
  } else {
    stop2(
      "accumulators must be NULL, a named numeric vector of positive integers, \\
      or a named character vector of column names."
    )
  }

  # the Stan likelihood opens with log(n[response]); a winner with no
  # accumulator makes that -Inf and the chain dies with no usable message
  per_trial <- as.matrix(data[paste0(".rdm_n", seq_along(cat_names))])
  won <- per_trial[cbind(seq_len(nrow(data)), data$.rdm_cat)]
  empty_winners <- unique(cat_names[data$.rdm_cat[won < 1L]])
  stopif(
    length(empty_winners) > 0,
    "Category {collapse_comma(empty_winners)} wins on trials where accumulators \\
    gives it no accumulator. A category that can be chosen needs at least one \\
    accumulator on every trial where it wins."
  )

  model$other_vars$n_choices <- n_cats
  NextMethod("check_data")
}

############################################################################# !
# CHECK_FORMULA S3 methods                                               ####
############################################################################# !

# Scaling drift, gap, sp and s by a common factor leaves the likelihood exactly
# unchanged, so an intercept for s slides along that ray together with the
# drifts and the thresholds and none of them is identified. Contrasts of s are
# identified, because the reference cell pins the scale.
#' @export
check_formula.rdm <- function(model, data, formula) {
  s_form <- formula[["s"]]
  warnif(
    is_formula(s_form) && !is_constant(s_form) && has_intercept(s_form),
    "The formula for 's' has an intercept, so 's' is estimated on the same \\
    scale as the drift rates, 'gap' and 'sp': multiplying all of them by a \\
    common factor leaves the likelihood unchanged, so only their ratios are \\
    identified. Either fix 's' (the default, or to another value with \\
    s = log(0.9) in bmf()) or suppress the intercept to estimate contrasts of \\
    's' only, as in s ~ 0 + condition."
  )
  NextMethod("check_formula")
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.rdm_simple <- function(model, formula) {
  rt_var <- model$resp_vars$rt
  brms::bf(glue("{rt_var} | vint(.rdm_cat, .rdm_n1, .rdm_n2) ~ 1"))
}

#' @export
bmf2bf.rdm_custom <- function(model, formula) {
  rt_var <- model$resp_vars$rt
  cat_names <- model$other_vars$resp_cats
  n_cols <- paste0(".rdm_n", seq_along(cat_names))
  vint_args <- paste(c(".rdm_cat", n_cols), collapse = ", ")
  brms::bf(glue("{rt_var} | vint({vint_args}) ~ 1"))
}

############################################################################# !
# Stan code generation                                                   ####
############################################################################# !

# The family's log-likelihood is unrolled over the response categories rather
# than collecting the per-row drifts and counts into arrays for a per-trial
# helper: a vector or array temporary in a per-row Stan function is allocated
# on the autodiff stack for every observation on every leapfrog step, and the
# unrolled form measured 22% (sp fixed) and 8% (sp free) less time per
# gradient with bitwise-identical log density and gradients
# (local/logs/benchmark_rdm_2026-09-25.md). The winner contributes log(n_win)
# and n_win - 1 survival copies, each loser n_j copies; the strict reps > 0
# guard keeps zero-count or underflowed survivals out of the sum (0 * -inf
# would poison the likelihood with NaN). Only the constructor's default sp
# (see .rdm_start_var) takes the plain-Wald path.
.rdm_stan_code <- function(family_name, cat_names, start_var) {
  n_cats <- length(cat_names)
  cat_args <- paste(paste0("vector ", cat_names), collapse = ", ")
  n_args <- paste(paste0("array[] int n", seq_len(n_cats)), collapse = ", ")
  pdf <- if (start_var) {
    "rdm_log_pdf(t, {cat}[i], gap[i], sp[i], s[i])"
  } else {
    "swald_lpdf(rt[i] | {cat}[i], gap[i] + sp[i], ndt[i], s[i])"
  }
  surv <- if (start_var) {
    "rdm_log_surv(t, {cat}[i], gap[i], sp[i], s[i])"
  } else {
    "swald_lccdf(rt[i] | {cat}[i], gap[i] + sp[i], ndt[i], s[i])"
  }
  per_cat <- function(template, j) glue(template, cat = cat_names[j], j = j)
  reps <- vapply(seq_len(n_cats), function(j) {
    per_cat("    int reps{j} = (response[i] == {j}) ? n{j}[i] - 1 : n{j}[i];", j)
  }, character(1))
  winner <- vapply(seq_len(n_cats), function(j) {
    head <- if (j == 1) "    if" else if (j < n_cats) "    else if" else "    else"
    cond <- if (j < n_cats) glue(" (response[i] == {j})") else ""
    per_cat(paste0(head, cond, " lp = log(n{j}[i]) + ", pdf, ";"), j)
  }, character(1))
  survivals <- vapply(seq_len(n_cats), function(j) {
    per_cat(paste0("    if (reps{j} > 0) lp += reps{j} * ", surv, ";"), j)
  }, character(1))

  paste(c(
    glue("real {family_name}_lpdf(vector rt, vector mu, {cat_args}, vector gap, ",
         "vector ndt, vector s, vector sp, array[] int response, {n_args}) {{"),
    "  int N = rows(rt);",
    "  real log_lik = 0;",
    "  for (i in 1:N) {",
    "    real t = rt[i] - ndt[i];",
    "    real lp;",
    "    if (t <= 0) return negative_infinity();",
    reps, winner, survivals,
    "    log_lik += lp;",
    "  }",
    "  return log_lik;",
    "}"
  ), collapse = "\n")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

# The vint() columns hold the winning category and the per-category accumulator
# counts, one value per observation, in the order the Stan signature expects.
# reduce_sum slices Y inside partial_log_lik but passes a custom family's `vars`
# through whole, so a bare "vint1" would pair each slice's response times with
# the top of the data. start/end only exist inside partial_log_lik, so the
# columns are sliced only where brms really threads (see brms_slices_likelihood)
.rdm_family_vars <- function(n_cats) {
  vars <- paste0("vint", seq_len(n_cats + 1))
  if (brms_slices_likelihood()) paste0(vars, "[start:end]") else vars
}

#' @export
configure_model.rdm_simple <- function(model, data, formula) {
  cat_names <- c("driftc", "drifte")
  formula <- bmf2bf(model, formula)

  formula$family <- brms::custom_family(
    "rdm_simple",
    dpars = c("mu", cat_names, "gap", "ndt", "s", "sp"),
    links = c("identity", model$links$driftc, model$links$drifte,
              model$links$gap, model$links$ndt, model$links$s,
              model$links$sp),
    ub = rep(NA, 7),
    lb = c(NA, 0, 0, 0, 0, 0, 0),
    type = "real",
    vars = .rdm_family_vars(length(cat_names)),
    loop = FALSE,
    log_lik = log_lik_rdm_simple,
    posterior_predict = posterior_predict_rdm_simple,
    posterior_epred = posterior_epred_rdm_simple
  )

  stanvars <- brms::stanvar(
    scode = read_lines2(paste0(
      system.file("stan_chunks", package = "bmm"),
      "/cswald_helper_functions.stan"
    )),
    block = "functions"
  ) + brms::stanvar(
    scode = read_lines2(paste0(
      system.file("stan_chunks", package = "bmm"),
      "/rdm_functions.stan"
    )),
    block = "functions"
  ) + brms::stanvar(
    scode = .rdm_stan_code("rdm_simple", cat_names, .rdm_start_var(model)),
    block = "functions"
  )

  nlist(formula, data, stanvars)
}

#' @export
configure_model.rdm_custom <- function(model, data, formula) {
  cat_names <- model$other_vars$resp_cats
  n_cats <- length(cat_names)
  formula <- bmf2bf(model, formula)

  n_dpars <- n_cats + 5
  formula$family <- brms::custom_family(
    "rdm_custom",
    dpars = c("mu", cat_names, "gap", "ndt", "s", "sp"),
    links = c(
      "identity",
      vapply(cat_names, function(p) model$links[[p]], character(1)),
      model$links$gap, model$links$ndt, model$links$s, model$links$sp
    ),
    ub = rep(NA, n_dpars),
    lb = c(NA, rep(0, n_cats), 0, 0, 0, 0),
    type = "real",
    vars = .rdm_family_vars(n_cats),
    loop = FALSE,
    log_lik = log_lik_rdm_custom,
    posterior_predict = posterior_predict_rdm_custom,
    posterior_epred = posterior_epred_rdm_custom
  )

  stanvars <- brms::stanvar(
    scode = read_lines2(paste0(
      system.file("stan_chunks", package = "bmm"),
      "/cswald_helper_functions.stan"
    )),
    block = "functions"
  ) + brms::stanvar(
    scode = read_lines2(paste0(
      system.file("stan_chunks", package = "bmm"),
      "/rdm_functions.stan"
    )),
    block = "functions"
  ) + brms::stanvar(
    scode = .rdm_stan_code("rdm_custom", cat_names, .rdm_start_var(model)),
    block = "functions"
  )

  nlist(formula, data, stanvars)
}

############################################################################# !
# Post-processing functions (shared helpers)                             ####
############################################################################# !

# The per-draw parameters of one observation as the row-per-draw matrices
# .rdm_race_lpdf() and .rdm_race() take. get_dpar() returns a scalar for a dpar
# brms stores fixed, so every vector is grown to ndraws, and the drifts go
# through matrix() because vapply() returns a bare vector for a single draw.
.rdm_draw_pars <- function(i, prep, cat_names, n_cats) {
  n_draws <- prep$ndraws
  list(
    rt = prep$data$Y[i],
    response = prep$data$vint1[i],
    gap = rep_len(brms::get_dpar(prep, "gap", i = i), n_draws),
    ndt = rep_len(brms::get_dpar(prep, "ndt", i = i), n_draws),
    s = rep_len(brms::get_dpar(prep, "s", i = i), n_draws),
    sp = rep_len(brms::get_dpar(prep, "sp", i = i), n_draws),
    drift = matrix(vapply(cat_names, function(p) {
      rep_len(brms::get_dpar(prep, p, i = i), n_draws)
    }, numeric(n_draws)), nrow = n_draws),
    counts = matrix(vapply(
      seq_len(n_cats),
      function(j) prep$data[[paste0("vint", j + 1)]][i],
      integer(1)
    ), n_draws, n_cats, byrow = TRUE)
  )
}

.rdm_log_lik <- function(i, prep, cat_names, n_cats) {
  d <- .rdm_draw_pars(i, prep, cat_names, n_cats)
  .rdm_race_lpdf(
    t = d$rt - d$ndt, response = d$response, drift = d$drift, counts = d$counts,
    gap = d$gap, A = d$sp, s = d$s
  )
}

.rdm_posterior_predict <- function(i, prep, cat_names, n_cats, ...) {
  d <- .rdm_draw_pars(i, prep, cat_names, n_cats)
  race <- .rdm_race(drift = d$drift, gap = d$gap, A = d$sp, s = d$s,
                    counts = d$counts)
  race$rt + d$ndt
}

# E[RT] = ndt + int_0^inf prod_j S_j(t)^n_j dt. A Monte-Carlo estimate of this
# integral moved by several percent between two calls on the same draws, which
# reached the user as noise on conditional_effects(); the grid is deterministic.
.rdm_posterior_epred <- function(prep, cat_names, n_cats, ...) {
  epred <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs)
  for (i in seq_len(prep$nobs)) {
    d <- .rdm_draw_pars(i, prep, cat_names, n_cats)
    # the counts belong to the observation, not the draw, so every row is the
    # same and the first one names the categories that race at all
    racing <- which(d$counts[1, ] > 0)

    log_surv <- function(t) {
      out <- matrix(0, prep$ndraws, length(t))
      for (j in racing) {
        out <- out + d$counts[, j] *
          wald_log_surv(t, d$drift[, j], d$gap, d$sp, d$s)
      }
      out
    }

    # the race is over once the slowest single accumulator is, so the plain
    # Wald bound S(t) <= Phi((b - v t) / (s sqrt t)) on the smallest drift and
    # the full distance b = gap + sp sets t_hi where it falls below 1e-12;
    # t_lo is where the fastest accumulator's CDF, at the shortest distance
    # gap, is still that small, so the survivor is 1 below it
    v_lo <- matrixStats::rowMins(d$drift[, racing, drop = FALSE])
    v_hi <- matrixStats::rowMaxs(d$drift[, racing, drop = FALSE])
    root <- 7.1 * d$s
    t_hi <- max(((root + sqrt(root^2 + 4 * v_lo * (d$gap + d$sp))) /
                   (2 * v_lo))^2)
    t_lo <- max(
      min(((sqrt(root^2 + 4 * v_hi * d$gap) - root) / (2 * v_hi))^2),
      t_hi * 1e-12
    )

    epred[, i] <- d$ndt + race_expected_time(log_surv, t_lo, t_hi)
  }
  epred
}

log_lik_rdm_simple <- function(i, prep) {
  .rdm_log_lik(i, prep, cat_names = c("driftc", "drifte"), n_cats = 2)
}

posterior_predict_rdm_simple <- function(i, prep, ...) {
  .rdm_posterior_predict(i, prep, cat_names = c("driftc", "drifte"),
                         n_cats = 2, ...)
}

posterior_epred_rdm_simple <- function(prep, ...) {
  .rdm_posterior_epred(prep, cat_names = c("driftc", "drifte"),
                       n_cats = 2, ...)
}

log_lik_rdm_custom <- function(i, prep) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "gap", "ndt", "s", "sp")
  )
  .rdm_log_lik(i, prep, cat_names = cat_names, n_cats = length(cat_names))
}

posterior_predict_rdm_custom <- function(i, prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "gap", "ndt", "s", "sp")
  )
  .rdm_posterior_predict(i, prep, cat_names = cat_names,
                         n_cats = length(cat_names), ...)
}

posterior_epred_rdm_custom <- function(prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "gap", "ndt", "s", "sp")
  )
  .rdm_posterior_epred(prep, cat_names = cat_names,
                       n_cats = length(cat_names), ...)
}

############################################################################# !
# PP_CHECK OBSERVABLES                                                    ####
############################################################################# !

# brms::posterior_predict() returns one matrix, so the winning category rides
# along as a second observable here instead. vint1 is the category code that
# check_data() wrote, in the order of the family's accumulator dpars.
#' @export
pp_observables.rdm <- function(model) {
  cats <- model$other_vars$resp_cats %||% c("correct", "error")
  coding <- collapse_comma(glue("{seq_along(cats)} = {cats}"))
  list(
    observed = c(rt = "Y", response = "vint1"),
    checks = list(
      rt = .pp_observable(function(d) d$rt, label = "Response time"),
      response = .pp_observable(
        function(d) d$response,
        label = glue("Response category ({coding})"),
        type = "bars"
      )
    )
  )
}

# One method for both versions: the accumulator dpars carry their own order, so
# nothing here depends on whether they came from n_choices or from the formula.
# rt and response come out of ONE race, so that a fast trial is a fast trial of
# the accumulator that actually won it.
#' @export
pp_simulate.rdm <- function(model, prep) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "gap", "ndt", "s", "sp"))
  n_cats <- length(cat_names)
  n_row <- prep$ndraws * prep$nobs

  race <- .rdm_race(
    drift = matrix(vapply(cat_names, .pp_dpar_vector, numeric(n_row),
                          prep = prep), nrow = n_row),
    gap = .pp_dpar_vector(prep, "gap"),
    A = .pp_dpar_vector(prep, "sp"),
    s = .pp_dpar_vector(prep, "s"),
    counts = matrix(vapply(seq_len(n_cats), function(j) {
      rep(prep$data[[paste0("vint", j + 1)]], each = prep$ndraws)
    }, integer(n_row)), nrow = n_row)
  )

  list(
    rt = matrix(race$rt + .pp_dpar_vector(prep, "ndt"), nrow = prep$ndraws),
    response = matrix(race$response, nrow = prep$ndraws)
  )
}
