# Per-distribution spec: everything that differs between the four drift-rate
# distributions. `drift` describes the driftc/drifte accumulator parameters; the
# gap/sp/ndt/s block lives in .lba_shared (identical across all distributions).
.lba_dist_specs <- list(
  normal = list(
    s_desc = "drift rate SD (fixed to 1 by default)",
    drift = list(
      desc = "mean drift rate",
      link = "identity",
      priors = list(
        driftc = list(main = "normal(3, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
        drifte = list(main = "normal(1, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)")
      ),
      inits = list(driftc = c(2, 4), drifte = c(0.5, 2))
    )
  ),
  gamma = list(
    s_desc = "rate parameter (fixed to 1 by default)",
    drift = list(
      desc = "drift rate (gamma shape)",
      link = "log",
      priors = list(
        driftc = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"),
        drifte = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)")
      ),
      inits = list(driftc = c(1.5, 3), drifte = c(0.8, 2))
    )
  ),
  frechet = list(
    s_desc = "scale parameter (fixed to 1 by default)",
    drift = list(
      desc = "drift rate (Frechet shape)",
      link = "log",
      priors = list(
        driftc = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"),
        drifte = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)")
      ),
      inits = list(driftc = c(1.5, 3), drifte = c(0.8, 2))
    )
  ),
  lognormal = list(
    s_desc = "sdlog parameter (fixed to 1 by default)",
    drift = list(
      desc = "drift rate (meanlog)",
      link = "identity",
      priors = list(
        driftc = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"),
        drifte = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)")
      ),
      inits = list(driftc = c(0.2, 0.8), drifte = c(-0.3, 0.3))
    )
  )
)

# The gap/sp/ndt/s block, shared by all distributions and both versions. `s` has
# a shared link but a per-distribution description (added from .lba_dist_specs).
.lba_shared <- list(
  parameters = list(
    gap = "threshold gap (b = gap + sp)",
    sp = "maximum starting point (uniform on 0 to sp)",
    ndt = "non-decision time"
  ),
  links = list(gap = "log", sp = "log", ndt = "log", s = "log"),
  priors = list(
    gap = list(main = "normal(-0.5, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"),
    sp = list(main = "normal(-1, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"),
    ndt = list(main = "normal(-1.5, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"),
    s = list(main = "normal(0, 0.3)", effects = "normal(0, 0.2)", sd = "exponential(2)")
  ),
  inits = list(
    mu = c(-0.5, 0.5), gap = c(0.3, 0.8), sp = c(0.2, 0.5),
    ndt = c(0.025, 0.05), s = c(0.8, 1.2)
  )
)

# Compose the spec for one (distribution, version). Drift parameters precede the
# shared block: downstream code recovers accumulator names via
# setdiff(names(parameters), c("gap", "sp", "ndt", "s")), which relies on order.
.lba_model_spec <- function(distribution, version) {
  ds <- .lba_dist_specs[[distribution]]
  shared_params <- c(.lba_shared$parameters, list(s = ds$s_desc))
  fixed_parameters <- list(mu = 0, s = 0)

  if (version == "custom") {
    return(nlist(
      parameters = shared_params,
      links = .lba_shared$links,
      fixed_parameters,
      priors = .lba_shared$priors,
      init_ranges = .lba_shared$inits
    ))
  }

  nlist(
    parameters = c(
      list(
        driftc = paste(ds$drift$desc, "for correct accumulator"),
        drifte = paste(ds$drift$desc, "for error accumulators")
      ),
      shared_params
    ),
    links = c(list(driftc = ds$drift$link, drifte = ds$drift$link), .lba_shared$links),
    fixed_parameters,
    priors = c(ds$drift$priors, .lba_shared$priors),
    init_ranges = c(.lba_shared$inits, ds$drift$inits)
  )
}

# Stan keywords and the names the generated model already uses (mu, Y, N, the
# shared parameters); one copy per racing model until a shared helper lands
.stan_reserved <- c(
  "int", "real", "vector", "matrix", "array", "if", "else", "for", "while",
  "return", "void", "data", "model", "target", "print", "reject", "log",
  "exp", "lower", "upper", "in", "functions", "generated", "transformed",
  "parameters", "break", "continue", "true", "false", "profile", "offset",
  "multiplier", "mu", "y", "n"
)

.model_lba <- function(
    rt = NULL,
    response = NULL,
    n_choices = NULL,
    accumulators = NULL,
    links = NULL,
    version = "simple",
    distribution = "normal",
    call = NULL,
    ...) {
  vt <- .lba_model_spec(distribution, version)
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = nlist(n_choices, accumulators),
      domain = "Decision Making / Response times",
      task = "Choice Reaction Time tasks (multi-alternative)",
      name = "Linear Ballistic Accumulator",
      citation = "Brown, S. D., & Heathcote, A. (2008). The simplest complete
        model of choice response time: Linear ballistic accumulation. Cognitive
        Psychology, 57(3), 153-178. https://doi.org/10.1016/j.cogpsych.2007.12.002",
      version = version,
      distribution = distribution,
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
    class = c("bmmodel", "lba", paste0("lba_", version)),
    call = call
  )

  set_links(out, links)
}

# The custom version learns its accumulator names from the formula, so at
# construction there is no vocabulary to validate a links argument against;
# check_model.lba_custom does it once the categories are known
#' @exportS3Method
settable_links.lba_custom <- function(model) {
  NULL
}

#' @title `r .model_lba()$name`
#' @name lba
#' @details `r model_info(.model_lba())`
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
#' @param version A character string specifying which version of the LBA model
#'   to use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): Two drift parameters, one for the correct
#'       accumulator (response = 1) and one for all error accumulators. The
#'       scale parameter `s` is shared and fixed to 1 by default (`s = 0` on
#'       the log scale). This covers the common case where interest is in the
#'       speed of correct vs. error processing.
#'     \item `"custom"`: Per-category drift parameters. Response categories
#'       are defined by the formula LHS names (e.g., `cat1 ~ 1, cat2 ~ 1`).
#'       The response column must contain character labels matching these names.
#'       Supports per-category `accumulators`.
#'   }
#' @param distribution A character string specifying the trial-to-trial drift
#'   rate distribution. All distributions use `driftc`/`drifte` as drift
#'   parameters (for the simple version), but their interpretation and link
#'   function differ:
#'   \itemize{
#'     \item `"normal"` (default): Drift rates drawn from a normal distribution.
#'       `driftc`/`drifte` = mean drift rate (identity link), `s` = SD (fixed=1).
#'     \item `"gamma"`: Drift rates drawn from a gamma distribution.
#'       `driftc`/`drifte` = shape parameter (log link), `s` = rate (fixed=1).
#'     \item `"frechet"`: Drift rates drawn from a Frechet distribution.
#'       `driftc`/`drifte` = shape parameter (log link), `s` = scale (fixed=1).
#'     \item `"lognormal"`: Drift rates drawn from a lognormal distribution.
#'       `driftc`/`drifte` = meanlog parameter, which may take any real value
#'       (identity link), `s` = sdlog (fixed=1).
#'   }
#' @param accumulators For `version = "custom"` only. A named vector
#'   specifying the number of racing accumulators per response category.
#'   Can be a named integer vector for constant counts or a named character
#'   vector of column names for trial-varying counts. If omitted, defaults to
#'   1 accumulator per category.
#' @param links A named list of link functions for the model parameters.
#' @param ... Additional arguments passed internally (for testing purposes).
#' @section Default priors: The `ndt` intercept prior is `normal(-1.5, 0.5)`
#'   on the log scale (median 0.22 s), the same as `ddm()`'s. An earlier
#'   `normal(-2, 0.3)` (median 0.135 s, 2% of its mass above 0.25 s) pulled
#'   `ndt` down and drift rates/thresholds up on a 400-trial recovery.
#'
#'   Random-effects SDs get `exponential(1)` for the identity-linked normal
#'   mean drift and `exponential(2)` for the log-linked shape drifts (gamma,
#'   Frechet), the lognormal meanlog, and for `gap`, `sp`, `ndt` and `s`,
#'   following the package-wide rule of rate 1 for identity-linked drift
#'   parameters and rate 2 for log-linked, boundary, non-decision-time and
#'   start-point parameters. These rates are anchored on the hierarchical
#'   simulation design used for the LBA recovery (between-subject SDs
#'   0.30-0.60 for `driftc`, 0.20-0.50 for `drifte`, and 0.10-0.25 on the log
#'   scale for `gap`, `sp`, `ndt`): `exponential(1)` puts 36% of its mass
#'   below 0.45, `exponential(2)` puts 29% below 0.17.
#' @section Numerical notes: For `distribution = "gamma"`, the accumulator's
#'   shape is `driftc`/`drifte` + 1, and `s * b / t` grows without bound as
#'   `ndt` approaches the fastest response time. Stan Math's gamma CDF has an
#'   approximate shape partial above shape ~5 and an undefined one once the
#'   CDF saturates (stan-dev/math #3408), making the gamma drift the slowest
#'   and least robust of the four distributions near the fastest response
#'   times; the default gamma prior keeps shape away from this region.
#'
#'   For `distribution = "frechet"`, the accumulator's moment integral has no
#'   closed form and is evaluated by a 16-point Gauss-Legendre rule, accurate
#'   to 1e-6 nats for `sp / gap <= 6`, 6e-5 nats at 10, and 8.7e-3 nats at 40;
#'   the default priors put `sp / gap` near 0.6.
#' @section Identifiability of `s`: The evidence scale of the LBA is
#'   arbitrary: multiplying every drift rate, `gap`, `sp` and `s` by the same
#'   constant leaves the likelihood unchanged. `s` is therefore fixed to 1 by
#'   convention (`0` on the log scale) to identify the rest of the
#'   parameters. `check_model()` warns if a formula frees `s`; free it only
#'   with another parameter fixed in its place.
#' @section log_lik and posterior predictions: `log_lik()`/`loo()` use the
#'   category-level likelihood: for `version = "simple"`, an error trial with
#'   `K` response categories contributes `log(K - 1)` more than a
#'   per-response likelihood would, i.e. the likelihood of "some error"
#'   rather than of the specific error observed.
#'
#'   `posterior_epred()` returns the expected response time (`ndt` plus the
#'   integral of the race survivor); it is deterministic, not
#'   simulation-based.
#'
#'   `pp_check(fit, resp_var = "response")` checks the response category
#'   (1 = correct, 2 = error for `version = "simple"`; category order for
#'   `version = "custom"`). `pp_check_vars(fit)` lists the available checks.
#' @return An object of class `bmmodel`
#' @note Both versions describe the same response type (a categorical winner in
#'   a choice-RT race), so they live in one constructor rather than separate
#'   model functions: `"simple"` is an accuracy-coded convenience layer (correct
#'   vs. error) over the general per-accumulator case handled by `"custom"`.
#' @export
#' @keywords bmmodel
#' @seealso [dlba()] and [rlba()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simple version with 2 alternatives
#' dat <- rlba(n = 500, drift = c(3, 1.5), gap = 0.5, sp = 0.5, ndt = 0.2)
#' model <- lba(rt = "rt", response = "response", n_choices = 2)
#' formula <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, sp ~ 1, ndt ~ 1)
#' fit <- bmm(formula, dat, model, cores = 4, backend = "cmdstanr")
lba <- function(rt, response, n_choices = NULL,
                version = c("simple", "custom"),
                distribution = c("normal", "gamma", "frechet", "lognormal"),
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
  distribution <- match.arg(distribution)
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
  .model_lba(
    rt = rt,
    response = response,
    n_choices = n_choices,
    accumulators = accumulators,
    links = links,
    version = version,
    distribution = distribution,
    call = call
  )
}

# The evidence scale of the LBA is arbitrary: multiplying every drift rate,
# gap, sp and s by one constant leaves the likelihood unchanged, so s = 1 is
# the convention that identifies the rest. Freeing s in the formula removes
# that constraint silently; the user may constrain another parameter instead,
# so this is a warning rather than an error.
#' @export
check_model.lba <- function(model, data = NULL, formula = NULL) {
  warnif(
    !is.null(formula) && "s" %in% names(formula) && !is_constant(formula)[["s"]],
    "The formula estimates `s`, the scale of the evidence, which the LBA \\
    likelihood cannot identify: scaling every drift rate, `gap`, `sp` and `s` \\
    by one constant leaves it unchanged. `s` is fixed to 1 by convention; free \\
    it only with another parameter fixed in its place. See ?lba."
  )
  NextMethod("check_model")
}

#' @export
check_model.lba_custom <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    cat_pars <- setdiff(names(formula), c("gap", "sp", "ndt", "s"))

    stopif(
      length(cat_pars) == 0,
      "Custom version requires at least one accumulator parameter in the formula."
    )

    bad_names <- intersect(tolower(cat_pars), .stan_reserved)
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

    known <- c(cat_pars, "gap", "sp", "ndt", "s")
    unknown <- setdiff(names(model$links), known)
    stopif(
      length(unknown) > 0,
      "Unrecognized link target(s): {collapse_comma(unknown)}. \\
      lba() takes links for {collapse_comma(known)}"
    )

    drift <- .lba_dist_specs[[model$distribution]]$drift
    for (p in cat_pars) {
      model$parameters[[p]] <- paste0(drift$desc, " for '", p, "' accumulator")
      if (is.null(model$links[[p]])) model$links[[p]] <- drift$link
      if (is.null(model$default_priors[[p]])) {
        model$default_priors[[p]] <- drift$priors$driftc
      }
      if (is.null(model$init_ranges[[p]])) {
        model$init_ranges[[p]] <- drift$inits$driftc
      }
    }

    model$other_vars$resp_cats <- cat_pars
  }

  NextMethod("check_model")
}

#' @export
check_data.lba <- function(model, data, formula) {
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
check_data.lba_simple <- function(model, data, formula) {
  response_var <- model$resp_vars$response
  n_alt <- model$other_vars$n_choices

  if (is.factor(data[, response_var]) || is.character(data[, response_var])) {
    labels <- as.character(data[, response_var])
    stopif(
      !all(grepl("^[0-9]+$", labels)),
      "The response variable '{response_var}' must be integer-coded for \\
      version 'simple' (1 = correct, 2:{n_alt} = errors); it contains the \\
      labels {collapse_comma(unique(labels[!grepl('^[0-9]+$', labels)]))}. \\
      Recode the responses or use version = 'custom' with labelled categories."
    )
    data[, response_var] <- as.integer(labels)
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

  data$.lba_cat <- ifelse(data[, response_var] == 1L, 1L, 2L)
  data$.lba_n1 <- 1L
  data$.lba_n2 <- n_alt - 1L

  NextMethod("check_data")
}

#' @export
check_data.lba_custom <- function(model, data, formula) {
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

  data$.lba_cat <- setNames(seq_along(cat_names), cat_names)[data[, response_var]]

  if (is.null(num_alt)) {
    for (i in seq_along(cat_names)) {
      data[[paste0(".lba_n", i)]] <- 1L
    }
  } else if (is.numeric(num_alt)) {
    stopif(
      !all(cat_names %in% names(num_alt)),
      "accumulators must have names matching formula categories: \\
      {collapse_comma(cat_names)}"
    )
    for (i in seq_along(cat_names)) {
      data[[paste0(".lba_n", i)]] <- as.integer(num_alt[cat_names[i]])
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
    for (i in seq_along(cat_names)) {
      data[[paste0(".lba_n", i)]] <- as.integer(data[, num_alt[cat_names[i]]])
    }
  }

  n_win <- vapply(seq_len(nrow(data)), function(i) {
    data[[paste0(".lba_n", data$.lba_cat[i])]][i]
  }, integer(1))
  stopif(
    any(n_win < 1),
    "{sum(n_win < 1)} trial(s) were answered with a response category that has \\
    no accumulator on that trial (accumulators = 0); the likelihood of such a \\
    trial is zero. Check the `accumulators` columns against the responses."
  )

  model$other_vars$n_choices <- n_cats
  NextMethod("check_data")
}

#' @export
bmf2bf.lba_simple <- function(model, formula) {
  rt_var <- model$resp_vars$rt
  brms::bf(glue("{rt_var} | vint(.lba_cat, .lba_n1, .lba_n2) ~ 1"))
}

#' @export
bmf2bf.lba_custom <- function(model, formula) {
  rt_var <- model$resp_vars$rt
  cat_names <- model$other_vars$resp_cats
  n_cols <- paste0(".lba_n", seq_along(cat_names))
  vint_args <- paste(c(".lba_cat", n_cols), collapse = ", ")
  brms::bf(glue("{rt_var} | vint({vint_args}) ~ 1"))
}

.lba_stan_code <- function(family_name, cat_names, distribution) {
  n_cats <- length(cat_names)
  cat_args <- paste(paste0("vector ", cat_names), collapse = ", ")
  n_args <- paste(paste0("array[] int n", seq_len(n_cats)), collapse = ", ")
  drift_assignments <- paste(
    paste0("drift[", seq_len(n_cats), "] = ", cat_names, "[i];"),
    collapse = "\n    "
  )
  n_assignments <- paste(
    paste0("n[", seq_len(n_cats), "] = n", seq_len(n_cats), "[i];"),
    collapse = "\n    "
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  shared_helpers <- read_lines2(file.path(sc_path, "lba_shared_functions.stan"))
  helpers <- read_lines2(
    file.path(sc_path, paste0("lba_", distribution, "_functions.stan"))
  )

  prefix <- paste0("lba_", distribution)

  main_fn <- glue(
    "real {family_name}_lpdf(vector rt, vector mu, {cat_args}, ",
    "vector gap, vector sp, vector ndt, vector s, array[] int response, ",
    "{n_args}) {{\n",
    "  int N = num_elements(rt);\n",
    "  real log_lik = 0;\n",
    "  for (i in 1:N) {{\n",
    "    real b = gap[i] + sp[i];\n",
    "    real A = sp[i];\n",
    "    real t = rt[i] - ndt[i];\n",
    "    array[{n_cats}] real drift;\n",
    "    array[{n_cats}] int n;\n",
    "    int win;\n",
    "    real lp;\n",
    "    if (t <= 0) return negative_infinity();\n",
    "    {drift_assignments}\n",
    "    {n_assignments}\n",
    "    win = response[i];\n",
    "    lp = log(n[win]) + {prefix}_single_lpdf(t | drift[win], b, A, s[i]);\n",
    "    for (j in 1:{n_cats}) {{\n",
    "      int reps = (j == win) ? n[j] - 1 : n[j];\n",
    "      if (reps > 0) {{\n",
    "        lp += reps * {prefix}_single_lccdf(t | drift[j], b, A, s[i]);\n",
    "      }}\n",
    "    }}\n",
    "    log_lik += lp;\n",
    "  }}\n",
    "  return log_lik;\n",
    "}}"
  )

  paste0(shared_helpers, "\n", helpers, "\n", main_fn)
}

# The vint() columns hold the winning category and the per-category
# accumulator counts. brms slices Y per thread but pastes a custom family's
# `vars` in whole, so under real threading the columns are sliced here;
# start/end only exist inside partial_log_lik, which is why the slicing follows
# brms_slices_likelihood() rather than the threads argument (threading(force =
# TRUE) compiles threaded but keeps the serial likelihood)
.lba_family_vars <- function(n_cats) {
  vars <- paste0("vint", seq_len(n_cats + 1))
  if (brms_slices_likelihood()) paste0(vars, "[start:end]") else vars
}

#' @export
configure_model.lba_simple <- function(model, data, formula) {
  links <- model$links
  dist <- model$distribution
  cat_names <- c("driftc", "drifte")
  formula <- bmf2bf(model, formula)

  dpars <- c("mu", cat_names, "gap", "sp", "ndt", "s")
  link_vec <- c(
    "identity",
    vapply(cat_names, function(p) links[[p]], character(1)),
    links$gap, links$sp, links$ndt,
    if (is.null(links$s)) "log" else links$s
  )

  family_name <- paste0("lba_", dist, "_simple")
  formula$family <- brms::custom_family(
    family_name,
    dpars = dpars,
    links = link_vec,
    ub = rep(NA, length(dpars)),
    lb = rep(NA, length(dpars)),
    type = "real",
    vars = .lba_family_vars(length(cat_names)),
    loop = FALSE,
    log_lik = log_lik_lba_simple,
    posterior_predict = posterior_predict_lba_simple,
    posterior_epred = posterior_epred_lba_simple
  )

  stan_code <- .lba_stan_code(family_name, cat_names, dist)
  stanvars <- brms::stanvar(scode = stan_code, block = "functions")

  nlist(formula, data, stanvars)
}

#' @export
configure_model.lba_custom <- function(model, data, formula) {
  links <- model$links
  dist <- model$distribution
  cat_names <- model$other_vars$resp_cats
  n_cats <- length(cat_names)
  formula <- bmf2bf(model, formula)

  dpars <- c("mu", cat_names, "gap", "sp", "ndt", "s")
  link_vec <- c(
    "identity",
    vapply(cat_names, function(p) links[[p]], character(1)),
    links$gap, links$sp, links$ndt,
    if (is.null(links$s)) "log" else links$s
  )

  family_name <- paste0("lba_", dist, "_custom")

  formula$family <- brms::custom_family(
    family_name,
    dpars = dpars,
    links = link_vec,
    ub = rep(NA, length(dpars)),
    lb = rep(NA, length(dpars)),
    type = "real",
    vars = .lba_family_vars(n_cats),
    loop = FALSE,
    log_lik = log_lik_lba_custom,
    posterior_predict = posterior_predict_lba_custom,
    posterior_epred = posterior_epred_lba_custom
  )

  stan_code <- .lba_stan_code(family_name, cat_names, dist)
  stanvars <- brms::stanvar(scode = stan_code, block = "functions")

  nlist(formula, data, stanvars)
}

.lba_log_lik <- function(i, prep, cat_names, n_cats) {
  rt <- prep$data$Y[i]
  response <- prep$data$vint1[i]
  gap <- brms::get_dpar(prep, "gap", i = i)
  sp <- brms::get_dpar(prep, "sp", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  dist <- .lba_dist_from_family(prep$family$name)

  t <- rt - ndt
  t[t <= 0] <- NA
  b <- gap + sp
  A <- sp

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )

  drift_win <- brms::get_dpar(prep, cat_names[response], i = i)
  log_lik <- log(n_cat[response]) + .lba_lpdf_single(t, drift_win, b, A, s, dist)
  if (n_cat[response] > 1) {
    log_lik <- log_lik + (n_cat[response] - 1) * .lba_lsurv_single(t, drift_win, b, A, s, dist)
  }
  for (j in seq_len(n_cats)) {
    if (j == response || n_cat[j] == 0) next
    drift_j <- brms::get_dpar(prep, cat_names[j], i = i)
    log_lik <- log_lik + n_cat[j] * .lba_lsurv_single(t, drift_j, b, A, s, dist)
  }

  log_lik[is.na(log_lik)] <- -Inf
  log_lik
}

.lba_posterior_predict <- function(i, prep, cat_names, n_cats, ...) {
  gap <- brms::get_dpar(prep, "gap", i = i)
  sp <- brms::get_dpar(prep, "sp", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  dist <- .lba_dist_from_family(prep$family$name)
  n_draws <- length(ndt)
  drift <- lapply(cat_names, function(p) brms::get_dpar(prep, p, i = i))
  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )

  b <- gap + sp
  min_ft <- rep(Inf, n_draws)
  for (j in seq_len(n_cats)) {
    if (n_cat[j] == 0) next
    for (k in seq_len(n_cat[j])) {
      start <- stats::runif(n_draws, 0, sp)
      dv <- .rlba_drift(dist, drift[[j]], s)
      min_ft <- pmin(min_ft, (b - start) / dv)
    }
  }
  min_ft + ndt
}

# E[RT] = ndt + int_0^inf S_race(t) dt, with S_race = prod_j S_j(t)^n_j the
# probability that no accumulator has finished by decision time t, on a log
# grid extended until the survivor's tail no longer contributes (the K = 2
# normal race decays like 1/t^2: a finite mean with an infinite variance, so a
# Monte-Carlo mean of the draws never settles)
.lba_posterior_epred <- function(prep, cat_names, n_cats, ...) {
  n_obs <- prep$nobs
  n_draws <- prep$ndraws
  dist <- .lba_dist_from_family(prep$family$name)

  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (i in seq_len(n_obs)) {
    gap <- brms::get_dpar(prep, "gap", i = i)
    sp <- brms::get_dpar(prep, "sp", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    s <- brms::get_dpar(prep, "s", i = i)
    drift <- lapply(cat_names, function(p) brms::get_dpar(prep, p, i = i))
    n_cat <- vapply(
      seq_len(n_cats),
      function(j) prep$data[[paste0("vint", j + 1)]][i],
      integer(1)
    )
    b <- gap + sp
    log_survivor <- function(t) {
      t_all <- rep(t, each = n_draws)
      Reduce(`+`, lapply(which(n_cat > 0), function(j) {
        n_cat[j] * .lba_lsurv_single(t_all, drift[[j]], b, sp, s, dist)
      }))
    }
    epred[, i] <- rep_len(ndt, n_draws) + race_mean_decision_time(log_survivor, n_draws)
  }
  epred
}

############################################################################# !
# PP_CHECK OBSERVABLES                                                    ####
############################################################################# !

#' @export
pp_observables.lba <- function(model) {
  cats <- if (model$version == "simple") c("correct", "error") else model$other_vars$resp_cats
  list(
    observed = c(rt = "Y", response = "vint1"),
    checks = list(
      rt = .pp_observable(function(d) d$rt, label = "Response time"),
      response = .pp_observable(
        function(d) d$response,
        label = paste0("Response category (", paste0(seq_along(cats), " = ", cats, collapse = ", "), ")"),
        type = "bars"
      )
    )
  )
}

# One race per draw x observation: every accumulator of every category draws
# its start point and drift, the fastest finishing time and its category are
# the observables. Category counts vary by observation (the custom version's
# accumulators columns), so the inner loop runs over the largest count and
# masks the observations where that accumulator exists.
.lba_pp_simulate <- function(model, prep) {
  cat_names <- setdiff(names(model$parameters), c("gap", "sp", "ndt", "s"))
  dist <- model$distribution
  n <- prep$ndraws * prep$nobs
  sp <- .pp_dpar_vector(prep, "sp")
  s <- .pp_dpar_vector(prep, "s")
  ndt <- .pp_dpar_vector(prep, "ndt")
  b <- .pp_dpar_vector(prep, "gap") + sp

  min_ft <- rep(Inf, n)
  winner <- rep(NA_integer_, n)
  for (j in seq_along(cat_names)) {
    n_cat <- as.vector(.pp_expand_data(prep$data[[paste0("vint", j + 1)]], prep$ndraws))
    drift <- .pp_dpar_vector(prep, cat_names[j])
    for (k in seq_len(max(n_cat))) {
      active <- which(n_cat >= k)
      ft <- (b[active] - stats::runif(length(active), 0, sp[active])) /
        .rlba_drift(dist, drift[active], s[active])
      faster <- active[ft < min_ft[active]]
      min_ft[faster] <- ft[ft < min_ft[active]]
      winner[faster] <- j
    }
  }
  list(
    rt = matrix(min_ft + ndt, nrow = prep$ndraws),
    response = matrix(winner, nrow = prep$ndraws)
  )
}

#' @export
pp_simulate.lba_simple <- function(model, prep) {
  .lba_pp_simulate(model, prep)
}

#' @export
pp_simulate.lba_custom <- function(model, prep) {
  .lba_pp_simulate(model, prep)
}

.lba_dist_from_family <- function(family_name) {
  sub("^lba_(.+)_(simple|custom)$", "\\1", family_name)
}

log_lik_lba_simple <- function(i, prep) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "gap", "sp", "ndt", "s"))
  .lba_log_lik(i, prep, cat_names = cat_names, n_cats = length(cat_names))
}

posterior_predict_lba_simple <- function(i, prep, ...) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "gap", "sp", "ndt", "s"))
  .lba_posterior_predict(i, prep, cat_names = cat_names,
                         n_cats = length(cat_names), ...)
}

log_lik_lba_custom <- function(i, prep) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "gap", "sp", "ndt", "s")
  )
  .lba_log_lik(i, prep, cat_names = cat_names, n_cats = length(cat_names))
}

posterior_predict_lba_custom <- function(i, prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "gap", "sp", "ndt", "s")
  )
  .lba_posterior_predict(i, prep, cat_names = cat_names,
                         n_cats = length(cat_names), ...)
}

posterior_epred_lba_simple <- function(prep, ...) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "gap", "sp", "ndt", "s"))
  .lba_posterior_epred(prep, cat_names = cat_names,
                       n_cats = length(cat_names), ...)
}

posterior_epred_lba_custom <- function(prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "gap", "sp", "ndt", "s")
  )
  .lba_posterior_epred(prep, cat_names = cat_names,
                       n_cats = length(cat_names), ...)
}
