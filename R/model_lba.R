.lba_dist_specs <- list(
  normal = list(
    s_desc = "drift rate SD (fixed to 1 by default)",
    drift_desc = "mean drift rate",
    drift_link = "identity",
    drift_priors = list(
      driftc = list(main = "normal(3, 1)", effects = "normal(0, 0.5)"),
      drifte = list(main = "normal(1, 1)", effects = "normal(0, 0.5)")
    ),
    drift_inits = list(driftc = c(2, 4), drifte = c(0.5, 2))
  ),
  gamma = list(
    s_desc = "rate parameter (fixed to 1 by default)",
    drift_desc = "drift rate (gamma shape)",
    drift_link = "log",
    drift_priors = list(
      driftc = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)"),
      drifte = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)")
    ),
    drift_inits = list(driftc = c(1.5, 3), drifte = c(0.8, 2))
  ),
  frechet = list(
    s_desc = "scale parameter (fixed to 1 by default)",
    drift_desc = "drift rate (Frechet shape)",
    drift_link = "log",
    drift_priors = list(
      driftc = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)"),
      drifte = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)")
    ),
    drift_inits = list(driftc = c(1.5, 3), drifte = c(0.8, 2))
  ),
  lognormal = list(
    s_desc = "sdlog parameter (fixed to 1 by default)",
    drift_desc = "drift rate (meanlog)",
    drift_link = "identity",
    drift_priors = list(
      driftc = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)"),
      drifte = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)")
    ),
    drift_inits = list(driftc = c(0.2, 0.8), drifte = c(-0.3, 0.3))
  )
)

.lba_shared_priors <- list(
  B = list(main = "normal(-0.5, 0.5)", effects = "normal(0, 0.3)"),
  A = list(main = "normal(-1, 0.5)", effects = "normal(0, 0.3)"),
  ndt = list(main = "normal(-2, 0.3)", effects = "normal(0, 0.3)"),
  s = list(main = "normal(0, 0.3)", effects = "normal(0, 0.2)")
)

.lba_shared_inits <- list(
  mu = c(-0.5, 0.5), B = c(0.3, 0.8), A = c(0.2, 0.5),
  ndt = c(0.025, 0.05), s = c(0.8, 1.2)
)

.build_lba_version <- function(dist) {
  spec <- .lba_dist_specs[[dist]]
  shared_params <- list(
    B = "threshold gap (b = A + B)",
    A = "maximum starting point (uniform on 0 to A)",
    ndt = "non-decision time",
    s = spec$s_desc
  )
  shared_links <- list(B = "log", A = "log", ndt = "log", s = "log")

  simple_params <- c(
    list(
      driftc = paste(spec$drift_desc, "for correct accumulator"),
      drifte = paste(spec$drift_desc, "for error accumulators")
    ),
    shared_params
  )
  simple_links <- c(
    list(driftc = spec$drift_link, drifte = spec$drift_link),
    shared_links
  )

  list(
    simple = list(
      parameters = simple_params,
      links = simple_links,
      fixed_parameters = list(mu = 0, s = 0),
      priors = c(spec$drift_priors, .lba_shared_priors),
      init_ranges = c(.lba_shared_inits, spec$drift_inits)
    ),
    custom = list(
      parameters = shared_params,
      links = shared_links,
      fixed_parameters = list(mu = 0, s = 0),
      priors = .lba_shared_priors,
      init_ranges = .lba_shared_inits
    )
  )
}

.lba_version_table <- Map(.build_lba_version, names(.lba_dist_specs))

.stan_reserved <- c(
  "int", "real", "vector", "matrix", "array", "if", "else", "for", "while",
  "return", "void", "data", "model", "target", "print", "reject", "log",
  "exp", "lower", "upper", "in", "functions", "generated", "transformed",
  "parameters"
)

.lba_cat_defaults <- lapply(.lba_dist_specs, function(spec) {
  list(
    link = spec$drift_link,
    prior = spec$drift_priors$driftc$main,
    init = spec$drift_inits$driftc,
    desc = spec$drift_desc
  )
})


.model_lba <- function(
    rt = NULL,
    response = NULL,
    n_alternatives = NULL,
    num_alternatives = NULL,
    links = NULL,
    version = "simple",
    distribution = "normal",
    call = NULL,
    ...) {
  vt <- .lba_version_table[[distribution]][[version]]
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = nlist(n_alternatives, num_alternatives),
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
      init_ranges = vt[["init_ranges"]],
      void_mu = TRUE
    ),
    class = c("bmmodel", "lba", paste0("lba_", version)),
    call = call
  )

  out$links[names(links)] <- links
  out
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
#' @param n_alternatives An integer specifying the total number of response
#'   alternatives (K >= 2). Required for `version = "simple"`. Not used for
#'   `version = "custom"` (inferred from the formula).
#' @param version A character string specifying which version of the LBA model
#'   to use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): Two drift parameters — one for the correct
#'       accumulator (response = 1) and one for all error accumulators. The
#'       scale parameter `s` is shared and fixed by default (s = 1). This
#'       covers the common case where interest is in the speed of correct vs.
#'       error processing.
#'     \item `"custom"`: Per-category drift parameters. Response categories
#'       are defined by the formula LHS names (e.g., `cat1 ~ 1, cat2 ~ 1`).
#'       The response column must contain character labels matching these names.
#'       Supports per-category `num_alternatives`.
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
#'       `driftc`/`drifte` = meanlog parameter (identity link), `s` = sdlog
#'       (fixed=1).
#'   }
#' @param num_alternatives For `version = "custom"` only. A named vector
#'   specifying the number of racing accumulators per response category.
#'   Can be a named integer vector for constant counts or a named character
#'   vector of column names for trial-varying counts. If omitted, defaults to
#'   1 accumulator per category.
#' @param links A named list of link functions for the model parameters.
#' @param ... Additional arguments passed internally (for testing purposes).
#' @return An object of class `bmmodel`
#' @export
#' @keywords bmmodel
#' @seealso [dlba()] and [rlba()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simple version with 2 alternatives
#' dat <- rlba(n = 500, drift = c(3, 1.5), B = 0.5, A = 0.5, ndt = 0.2)
#' model <- lba(rt = "rt", response = "response", n_alternatives = 2)
#' formula <- bmf(driftc ~ 1, drifte ~ 1, B ~ 1, A ~ 1, ndt ~ 1)
#' fit <- bmm(formula, dat, model, cores = 4, backend = "cmdstanr")
lba <- function(rt, response, n_alternatives = NULL,
                version = c("simple", "custom"),
                distribution = c("normal", "gamma", "frechet", "lognormal"),
                num_alternatives = NULL, links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  distribution <- match.arg(distribution)
  if (version == "simple") {
    stopif(
      is.null(n_alternatives) || !is.numeric(n_alternatives) ||
        n_alternatives < 2 || n_alternatives != round(n_alternatives),
      "n_alternatives must be an integer >= 2 for version 'simple'."
    )
    n_alternatives <- as.integer(n_alternatives)
  } else {
    n_alternatives <- NULL
  }
  .model_lba(
    rt = rt,
    response = response,
    n_alternatives = n_alternatives,
    num_alternatives = num_alternatives,
    links = links,
    version = version,
    distribution = distribution,
    call = call,
    ...
  )
}

#' @export
check_model.lba_custom <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    cat_pars <- setdiff(names(formula), c("B", "A", "ndt", "s"))

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

    defaults <- .lba_cat_defaults[[model$distribution]]
    for (p in cat_pars) {
      model$parameters[[p]] <- paste0(defaults$desc, " for '", p, "' accumulator")
      if (is.null(model$links[[p]])) model$links[[p]] <- defaults$link
      if (is.null(model$default_priors[[p]])) {
        model$default_priors[[p]] <- list(
          main = defaults$prior, effects = "normal(0, 0.3)"
        )
      }
      if (is.null(model$init_ranges[[p]])) {
        model$init_ranges[[p]] <- defaults$init
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
  n_alt <- model$other_vars$n_alternatives

  if (is.factor(data[, response_var])) {
    data[, response_var] <- as.integer(as.character(data[, response_var]))
  } else if (is.character(data[, response_var])) {
    data[, response_var] <- as.integer(data[, response_var])
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
  num_alt <- model$other_vars$num_alternatives
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
      "num_alternatives must have names matching formula categories: \\
      {collapse_comma(cat_names)}"
    )
    for (i in seq_along(cat_names)) {
      data[[paste0(".lba_n", i)]] <- as.integer(num_alt[cat_names[i]])
    }
  } else if (is.character(num_alt)) {
    stopif(
      !all(cat_names %in% names(num_alt)),
      "num_alternatives must have names matching formula categories: \\
      {collapse_comma(cat_names)}"
    )
    missing_cols <- setdiff(num_alt, colnames(data))
    stopif(
      length(missing_cols) > 0,
      "num_alternatives columns {collapse_comma(missing_cols)} not found \\
      in the data."
    )
    for (i in seq_along(cat_names)) {
      data[[paste0(".lba_n", i)]] <- as.integer(data[, num_alt[cat_names[i]]])
    }
  }

  model$other_vars$n_alternatives <- n_cats
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
  cat_args <- paste(paste0("real ", cat_names), collapse = ", ")
  n_args <- paste(paste0("int n", seq_len(n_cats)), collapse = ", ")
  drift_array <- paste0(
    "array[", n_cats, "] real drift = {",
    paste(cat_names, collapse = ", "), "};"
  )
  n_array <- paste0(
    "array[", n_cats, "] int n = {",
    paste(paste0("n", seq_len(n_cats)), collapse = ", "), "};"
  )

  stan_body <- switch(distribution,
    normal = .lba_stan_body_normal(n_cats),
    gamma = .lba_stan_body_gamma(n_cats),
    frechet = .lba_stan_body_frechet(n_cats),
    lognormal = .lba_stan_body_lognormal(n_cats)
  )

  glue(
    "real {family_name}_lpdf(real rt, real mu, {cat_args}, ",
    "real B, real A, real ndt, real s, int response, {n_args}) {{\n",
    "  real b = A + B;\n",
    "  real t = rt - ndt;\n",
    "  if (t <= 0) return negative_infinity();\n",
    "  {drift_array}\n",
    "  {n_array}\n",
    "{stan_body}\n",
    "}}"
  )
}

.lba_stan_body_normal <- function(n_cats) {
  "
  // Brown & Heathcote (2008) normal drift LBA
  real log_lik;
  real hi = b / t;
  real lo = (b - A) / t;
  if (A < 1e-6) {
    log_lik = log(n[response]) + log(b) - 2 * log(t)
              + normal_lpdf(hi | drift[response], s);
    for (j in 1:K) {
      real log_ccdf = normal_lccdf(hi | drift[j], s);
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log_ccdf;
      } else {
        log_lik += n[j] * log_ccdf;
      }
    }
  } else {
    real z_lo = (lo - drift[response]) / s;
    real z_hi = (hi - drift[response]) / s;
    real M_win = -drift[response] * Phi(z_lo) + s * exp(std_normal_lpdf(z_lo | ))
                 + drift[response] * Phi(z_hi) - s * exp(std_normal_lpdf(z_hi | ));
    real pdf_val = M_win / A;
    if (pdf_val <= 0) return negative_infinity();
    log_lik = log(n[response]) + log(pdf_val);
    for (j in 1:K) {
      real zl = (lo - drift[j]) / s;
      real zh = (hi - drift[j]) / s;
      real M_j = -drift[j] * Phi(zl) + s * exp(std_normal_lpdf(zl | ))
                 + drift[j] * Phi(zh) - s * exp(std_normal_lpdf(zh | ));
      real cdf_j = 1 + (1.0 / A) * (
        t * M_j - b * normal_cdf(hi | drift[j], s)
        + (b - A) * normal_cdf(lo | drift[j], s));
      real surv = 1 - cdf_j;
      if (surv < 1e-10) surv = 1e-10;
      if (surv > 1) surv = 1.0;
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log(surv);
      } else {
        log_lik += n[j] * log(surv);
      }
    }
  }
  return log_lik;"
}

.lba_stan_body_gamma <- function(n_cats) {
  "
  real log_lik;
  if (A < 1e-6) {
    log_lik = log(n[response]) + log(b) - 2 * log(t)
              + gamma_lpdf(b / t | drift[response], s);
    for (j in 1:K) {
      real log_ccdf = gamma_lccdf(b / t | drift[j], s);
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log_ccdf;
      } else {
        log_lik += n[j] * log_ccdf;
      }
    }
  } else {
    real hi = b / t;
    real lo = (b - A) / t;
    real pdf_val = (drift[response] / (s * A)) * (
      gamma_cdf(hi | drift[response] + 1, s) -
      gamma_cdf(lo | drift[response] + 1, s));
    if (pdf_val <= 0) return negative_infinity();
    log_lik = log(n[response]) + log(pdf_val);
    for (j in 1:K) {
      real M_j = (drift[j] / s) * (
        gamma_cdf(hi | drift[j] + 1, s) -
        gamma_cdf(lo | drift[j] + 1, s));
      real cdf_j = 1 + (1.0 / A) * (
        t * M_j - b * gamma_cdf(hi | drift[j], s) +
        (b - A) * gamma_cdf(lo | drift[j], s));
      real surv = 1 - cdf_j;
      if (surv < 1e-10) surv = 1e-10;
      if (surv > 1) surv = 1.0;
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log(surv);
      } else {
        log_lik += n[j] * log(surv);
      }
    }
  }
  return log_lik;"
}

.lba_stan_body_frechet <- function(n_cats) {
  "
  real log_lik;
  if (A < 1e-6) {
    real x = b / t;
    real z = x / s;
    log_lik = log(n[response]) + log(b) - 2 * log(t)
              + log(drift[response]) - log(s)
              - (1 + drift[response]) * log(z) - pow(z, -drift[response]);
    for (j in 1:K) {
      real xj = b / t;
      real zj = xj / s;
      real log_ccdf = -exp(-pow(zj, -drift[j]));
      log_ccdf = log1m_exp(log_ccdf);
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log_ccdf;
      } else {
        log_lik += n[j] * log_ccdf;
      }
    }
  } else {
    real hi = b / t;
    real lo = (b - A) / t;
    real M_win = 0;
    int n_quad = 20;
    real h = (hi - lo) / n_quad;
    for (q in 0:n_quad) {
      real u = lo + q * h;
      real z = u / s;
      real f_u = (drift[response] / s) * pow(z, -(1 + drift[response])) * exp(-pow(z, -drift[response]));
      real w = (q == 0 || q == n_quad) ? 0.5 : 1.0;
      M_win += w * u * f_u;
    }
    M_win *= h;
    real pdf_val = M_win / A;
    if (pdf_val <= 0) return negative_infinity();
    log_lik = log(n[response]) + log(pdf_val);
    for (j in 1:K) {
      real M_j = 0;
      for (q in 0:n_quad) {
        real u = lo + q * h;
        real z = u / s;
        real f_u = (drift[j] / s) * pow(z, -(1 + drift[j])) * exp(-pow(z, -drift[j]));
        real w = (q == 0 || q == n_quad) ? 0.5 : 1.0;
        M_j += w * u * f_u;
      }
      M_j *= h;
      real cdf_hi = exp(-pow(hi / s, -drift[j]));
      real cdf_lo = exp(-pow(lo / s, -drift[j]));
      real cdf_j = 1 + (1.0 / A) * (
        t * M_j - b * cdf_hi + (b - A) * cdf_lo);
      real surv = 1 - cdf_j;
      if (surv < 1e-10) surv = 1e-10;
      if (surv > 1) surv = 1.0;
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log(surv);
      } else {
        log_lik += n[j] * log(surv);
      }
    }
  }
  return log_lik;"
}

.lba_stan_body_lognormal <- function(n_cats) {
  "
  real log_lik;
  if (A < 1e-6) {
    log_lik = log(n[response]) + log(b) - 2 * log(t)
              + lognormal_lpdf(b / t | drift[response], s);
    for (j in 1:K) {
      real log_ccdf = lognormal_lccdf(b / t | drift[j], s);
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log_ccdf;
      } else {
        log_lik += n[j] * log_ccdf;
      }
    }
  } else {
    real hi = b / t;
    real lo = (b - A) / t;
    real mu_s2 = exp(drift[response] + square(s) / 2);
    real pdf_val = (mu_s2 / A) * (
      Phi((log(hi) - drift[response] - square(s)) / s) -
      Phi((log(lo) - drift[response] - square(s)) / s));
    if (pdf_val <= 0) return negative_infinity();
    log_lik = log(n[response]) + log(pdf_val);
    for (j in 1:K) {
      real mu_s2_j = exp(drift[j] + square(s) / 2);
      real M_j = mu_s2_j * (
        Phi((log(hi) - drift[j] - square(s)) / s) -
        Phi((log(lo) - drift[j] - square(s)) / s));
      real cdf_j = 1 + (1.0 / A) * (
        t * M_j - b * lognormal_cdf(hi | drift[j], s) +
        (b - A) * lognormal_cdf(lo | drift[j], s));
      real surv = 1 - cdf_j;
      if (surv < 1e-10) surv = 1e-10;
      if (surv > 1) surv = 1.0;
      if (j == response) {
        if (n[j] > 1) log_lik += (n[j] - 1) * log(surv);
      } else {
        log_lik += n[j] * log(surv);
      }
    }
  }
  return log_lik;"
}

#' @export
configure_model.lba_simple <- function(model, data, formula) {
  links <- model$links
  dist <- model$distribution
  cat_names <- setdiff(
    names(.lba_version_table[[dist]][["simple"]]$parameters),
    c("B", "A", "ndt", "s")
  )
  formula <- bmf2bf(model, formula)

  dpars <- c("mu", cat_names, "B", "A", "ndt", "s")
  link_vec <- c(
    "identity",
    vapply(cat_names, function(p) links[[p]], character(1)),
    links$B, links$A, links$ndt,
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
    vars = c("vint1[n]", "vint2[n]", "vint3[n]"),
    loop = TRUE,
    log_lik = log_lik_lba_simple,
    posterior_predict = posterior_predict_lba_simple,
    posterior_epred = posterior_epred_lba_simple
  )

  stan_code <- .lba_stan_code(family_name, cat_names, dist)
  stan_code <- gsub("\\bK\\b", length(cat_names), stan_code)
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

  dpars <- c("mu", cat_names, "B", "A", "ndt", "s")
  link_vec <- c(
    "identity",
    vapply(cat_names, function(p) links[[p]], character(1)),
    links$B, links$A, links$ndt,
    if (is.null(links$s)) "log" else links$s
  )

  vars_vec <- paste0("vint", seq_len(n_cats + 1), "[n]")
  family_name <- paste0("lba_", dist, "_custom")

  formula$family <- brms::custom_family(
    family_name,
    dpars = dpars,
    links = link_vec,
    ub = rep(NA, length(dpars)),
    lb = rep(NA, length(dpars)),
    type = "real",
    vars = vars_vec,
    loop = TRUE,
    log_lik = log_lik_lba_custom,
    posterior_predict = posterior_predict_lba_custom,
    posterior_epred = posterior_epred_lba_custom
  )

  stan_code <- .lba_stan_code(family_name, cat_names, dist)
  stan_code <- gsub("\\bK\\b", n_cats, stan_code)
  stanvars <- brms::stanvar(scode = stan_code, block = "functions")

  nlist(formula, data, stanvars)
}

.lba_log_lik <- function(i, prep, cat_names, n_cats) {
  rt <- prep$data$Y[i]
  response <- prep$data$vint1[i]
  B <- brms::get_dpar(prep, "B", i = i)
  A <- brms::get_dpar(prep, "A", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  # Determine distribution from family name
  dist <- .lba_dist_from_family(prep$family$name)

  t <- rt - ndt
  t[t <= 0] <- NA
  b <- A + B

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )

  drift_win <- brms::get_dpar(prep, cat_names[response], i = i)
  log_lik <- log(n_cat[response]) +
    .dlba_single(t, drift_win, b, A, s, dist, log = TRUE)

  if (n_cat[response] > 1) {
    surv_win <- 1 - .plba_single(t, drift_win, b, A, s, dist)
    surv_win[surv_win <= 0] <- 1e-300
    log_lik <- log_lik + (n_cat[response] - 1) * log(surv_win)
  }

  for (j in seq_len(n_cats)) {
    if (j == response) next
    drift_j <- brms::get_dpar(prep, cat_names[j], i = i)
    surv_j <- 1 - .plba_single(t, drift_j, b, A, s, dist)
    surv_j[surv_j <= 0] <- 1e-300
    log_lik <- log_lik + n_cat[j] * log(surv_j)
  }

  log_lik[is.na(log_lik)] <- -Inf
  log_lik
}

.lba_posterior_predict <- function(i, prep, cat_names, n_cats, ...) {
  B <- brms::get_dpar(prep, "B", i = i)
  A <- brms::get_dpar(prep, "A", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  dist <- .lba_dist_from_family(prep$family$name)
  n_draws <- length(ndt)

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )
  total_acc <- sum(n_cat)

  rt <- numeric(n_draws)
  for (d in seq_len(n_draws)) {
    b <- A[d] + B[d]
    ft <- numeric(total_acc)
    idx <- 1
    for (j in seq_len(n_cats)) {
      drift_j <- brms::get_dpar(prep, cat_names[j], i = i)[d]
      for (k in seq_len(n_cat[j])) {
        sp <- stats::runif(1, 0, A[d])
        dv <- switch(dist,
          normal = stats::rnorm(1, drift_j, s[d]),
          gamma = stats::rgamma(1, shape = drift_j, rate = s[d]),
          frechet = .rfrechet(1, shape = drift_j, scale = s[d]),
          lognormal = stats::rlnorm(1, meanlog = drift_j, sdlog = s[d])
        )
        ft[idx] <- if (dv > 0) (b - sp) / dv else Inf
        idx <- idx + 1
      }
    }
    rt[d] <- min(ft) + ndt[d]
  }
  rt
}

.lba_posterior_epred <- function(prep, cat_names, n_cats, ...) {
  n_obs <- prep$nobs
  n_draws <- prep$ndraws
  n_sim <- 100
  dist <- .lba_dist_from_family(prep$family$name)

  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (i in seq_len(n_obs)) {
    B <- brms::get_dpar(prep, "B", i = i)
    A <- brms::get_dpar(prep, "A", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    s <- brms::get_dpar(prep, "s", i = i)

    n_cat <- vapply(
      seq_len(n_cats),
      function(j) prep$data[[paste0("vint", j + 1)]][i],
      integer(1)
    )
    total_acc <- sum(n_cat)

    for (d in seq_len(n_draws)) {
      b <- A[d] + B[d]
      min_rts <- numeric(n_sim)
      for (sim in seq_len(n_sim)) {
        ft <- numeric(total_acc)
        idx <- 1
        for (j in seq_len(n_cats)) {
          drift_j <- brms::get_dpar(prep, cat_names[j], i = i)[d]
          for (k in seq_len(n_cat[j])) {
            sp <- stats::runif(1, 0, A[d])
            dv <- switch(dist,
              normal = stats::rnorm(1, drift_j, s[d]),
              gamma = stats::rgamma(1, shape = drift_j, rate = s[d]),
              frechet = .rfrechet(1, shape = drift_j, scale = s[d]),
              lognormal = stats::rlnorm(1, meanlog = drift_j, sdlog = s[d])
            )
            ft[idx] <- if (dv > 0) (b - sp) / dv else Inf
            idx <- idx + 1
          }
        }
        min_rts[sim] <- min(ft)
      }
      epred[d, i] <- mean(min_rts) + ndt[d]
    }
  }
  epred
}

.lba_dist_from_family <- function(family_name) {
  sub("^lba_(.+)_(simple|custom)$", "\\1", family_name)
}

log_lik_lba_simple <- function(i, prep) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "B", "A", "ndt", "s"))
  .lba_log_lik(i, prep, cat_names = cat_names, n_cats = length(cat_names))
}

posterior_predict_lba_simple <- function(i, prep, ...) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "B", "A", "ndt", "s"))
  .lba_posterior_predict(i, prep, cat_names = cat_names,
                         n_cats = length(cat_names), ...)
}

log_lik_lba_custom <- function(i, prep) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "B", "A", "ndt", "s")
  )
  .lba_log_lik(i, prep, cat_names = cat_names, n_cats = length(cat_names))
}

posterior_predict_lba_custom <- function(i, prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "B", "A", "ndt", "s")
  )
  .lba_posterior_predict(i, prep, cat_names = cat_names,
                         n_cats = length(cat_names), ...)
}

posterior_epred_lba_simple <- function(prep, ...) {
  cat_names <- setdiff(prep$family$dpars, c("mu", "B", "A", "ndt", "s"))
  .lba_posterior_epred(prep, cat_names = cat_names,
                       n_cats = length(cat_names), ...)
}

posterior_epred_lba_custom <- function(prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "B", "A", "ndt", "s")
  )
  .lba_posterior_epred(prep, cat_names = cat_names,
                       n_cats = length(cat_names), ...)
}
