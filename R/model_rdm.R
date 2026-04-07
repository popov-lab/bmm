############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.rdm_version_table <- list(
  simple = list(
    parameters = list(
      driftc = "drift rate for correct accumulator",
      drifte = "drift rate for error accumulators",
      bound = "decision threshold",
      ndt = "non-decision time",
      s = "diffusion constant",
      sp = "starting point proportion (A/bound, on logit scale)"
    ),
    links = list(
      driftc = "log",
      drifte = "log",
      bound = "log",
      ndt = "log",
      s = "log",
      sp = "logit"
    ),
    fixed_parameters = list(
      mu = 0,
      s = 0,
      sp = -100
    ),
    priors = list(
      driftc = list(main = "normal(1, 0.5)", effects = "normal(0, 0.3)"),
      drifte = list(main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)"),
      bound = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)"),
      ndt = list(main = "normal(-2, 0.3)", effects = "normal(0, 0.3)"),
      s = list(main = "normal(0, 0.3)", effects = "normal(0, 0.2)"),
      sp = list(main = "normal(-1, 1)", effects = "normal(0, 0.5)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      driftc = c(2, 4),
      drifte = c(1, 2.5),
      bound = c(0.8, 1.5),
      ndt = c(0.025, 0.05),
      s = c(0.8, 1.2),
      sp = c(0.15, 0.45)
    )
  ),
  custom = list(
    parameters = list(
      bound = "decision threshold",
      ndt = "non-decision time",
      s = "diffusion constant",
      sp = "starting point proportion (A/bound, on logit scale)"
    ),
    links = list(
      bound = "log",
      ndt = "log",
      s = "log",
      sp = "logit"
    ),
    fixed_parameters = list(
      mu = 0,
      s = 0,
      sp = -100
    ),
    priors = list(
      bound = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)"),
      ndt = list(main = "normal(-2, 0.3)", effects = "normal(0, 0.3)"),
      s = list(main = "normal(0, 0.3)", effects = "normal(0, 0.2)"),
      sp = list(main = "normal(-1, 1)", effects = "normal(0, 0.5)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      bound = c(0.8, 1.5),
      ndt = c(0.025, 0.05),
      s = c(0.8, 1.2),
      sp = c(0.15, 0.45)
    )
  )
)

.rdm_stan_reserved <- c(
  "int", "real", "vector", "matrix", "array", "if", "else", "for", "while",
  "return", "void", "data", "model", "target", "print", "reject", "log",
  "exp", "lower", "upper", "in", "functions", "generated", "transformed",
  "parameters"
)


.model_rdm <- function(
    rt = NULL,
    response = NULL,
    n_alternatives = NULL,
    num_alternatives = NULL,
    links = NULL,
    version = "simple",
    call = NULL,
    ...) {
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = nlist(n_alternatives, num_alternatives),
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
      parameters = .rdm_version_table[[version]][["parameters"]],
      links = .rdm_version_table[[version]][["links"]],
      fixed_parameters = .rdm_version_table[[version]][["fixed_parameters"]],
      default_priors = .rdm_version_table[[version]][["priors"]],
      init_ranges = .rdm_version_table[[version]][["init_ranges"]],
      void_mu = TRUE
    ),
    class = c("bmmodel", "rdm", paste0("rdm_", version)),
    call = call
  )

  out$links[names(links)] <- links
  out
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
#' @param n_alternatives An integer specifying the total number of response
#'   alternatives (K >= 2). Required for `version = "simple"`. Not used for
#'   `version = "custom"` (inferred from the formula).
#' @param num_alternatives For `version = "custom"` only. A named vector
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
#'       default (s = 1). Starting point proportion `sp` is fixed to ~0 by
#'       default; add `sp ~ 1` to the formula to estimate it.
#'     \item `"custom"`: Per-category drift parameters. Response categories
#'       are defined by the formula LHS names (e.g., `correct ~ 1, other ~ 1,
#'       npl ~ 1`). The response column must contain character labels matching
#'       these names. Supports per-category `num_alternatives`.
#'   }
#' @param links A named list of link functions for the model parameters.
#'   For `"simple"`: parameters are `driftc`, `drifte`, `bound`, `ndt`, `s`,
#'   and `sp`. Default links are "log" for all except `sp` which uses "logit".
#' @param ... Additional arguments passed internally (for testing purposes).
#' @return An object of class `bmmodel`
#' @export
#' @keywords bmmodel
#' @seealso [drdm()] and [rrdm()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simple version with 2 alternatives
#' dat <- rrdm(n = 500, drift = c(3, 1.5), bound = 1, ndt = 0.2)
#' model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
#' formula <- bmf(driftc ~ 1, drifte ~ 1, bound ~ 1, ndt ~ 1)
#' fit <- bmm(formula, dat, model, cores = 4, backend = "cmdstanr")
#'
#' # with starting point variability
#' formula2 <- bmf(driftc ~ 1, drifte ~ 1, bound ~ 1, ndt ~ 1, sp ~ 1)
#' fit2 <- bmm(formula2, dat, model, cores = 4, backend = "cmdstanr")
rdm <- function(rt, response, n_alternatives = NULL,
                version = c("simple", "custom"),
                num_alternatives = NULL, links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
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
  .model_rdm(
    rt = rt,
    response = response,
    n_alternatives = n_alternatives,
    num_alternatives = num_alternatives,
    links = links,
    version = version,
    call = call,
    ...
  )
}

############################################################################# !
# CHECK_MODEL S3 methods                                                 ####
############################################################################# !

#' @export
check_model.rdm_custom <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    reserved_pars <- c("bound", "ndt", "s", "sp")
    formula_pars <- names(formula)
    cat_pars <- setdiff(formula_pars, reserved_pars)

    stopif(
      length(cat_pars) == 0,
      "Custom version requires at least one accumulator parameter in the formula."
    )

    bad_names <- intersect(tolower(cat_pars), .rdm_stan_reserved)
    stopif(
      length(bad_names) > 0,
      "Category names cannot be Stan reserved words: {collapse_comma(bad_names)}. \\
      Please rename the affected response categories."
    )

    for (p in cat_pars) {
      model$parameters[[p]] <- paste0("drift rate for '", p, "' accumulator")
      if (is.null(model$links[[p]])) model$links[[p]] <- "log"
      if (is.null(model$default_priors[[p]])) {
        model$default_priors[[p]] <- list(
          main = "normal(0.5, 0.5)", effects = "normal(0, 0.3)"
        )
      }
      if (is.null(model$init_ranges[[p]])) {
        model$init_ranges[[p]] <- c(1, 2.5)
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

  data$.rdm_cat <- ifelse(data[, response_var] == 1L, 1L, 2L)
  data$.rdm_n1 <- 1L
  data$.rdm_n2 <- n_alt - 1L

  NextMethod("check_data")
}

#' @export
check_data.rdm_custom <- function(model, data, formula) {
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

  data$.rdm_cat <- setNames(seq_along(cat_names), cat_names)[data[, response_var]]

  if (is.null(num_alt)) {
    for (i in seq_along(cat_names)) {
      data[[paste0(".rdm_n", i)]] <- 1L
    }
  } else if (is.numeric(num_alt)) {
    stopif(
      !all(cat_names %in% names(num_alt)),
      "num_alternatives must have names matching formula categories: \\
      {collapse_comma(cat_names)}"
    )
    for (i in seq_along(cat_names)) {
      data[[paste0(".rdm_n", i)]] <- as.integer(num_alt[cat_names[i]])
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
      data[[paste0(".rdm_n", i)]] <- as.integer(data[, num_alt[cat_names[i]]])
    }
  }

  model$other_vars$n_alternatives <- n_cats
  NextMethod("check_data")
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

.rdm_stan_code <- function(family_name, cat_names, has_sp) {
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

  if (!has_sp) {
    glue(
      "real {family_name}_lpdf(real rt, real mu, {cat_args}, ",
      "real bound, real ndt, real s, real sp, int response, {n_args}) {{\n",
      "  real t = rt - ndt;\n",
      "  if (t <= 0) return negative_infinity();\n",
      "  {drift_array}\n",
      "  {n_array}\n",
      "  real log_lik = log(n[response]) + ",
      "swald_lpdf(rt | drift[response], bound, ndt, s);\n",
      "  for (j in 1:{n_cats}) {{\n",
      "    if (j == response) {{\n",
      "      if (n[j] > 1)\n",
      "        log_lik += (n[j] - 1) * ",
      "swald_lccdf(rt | drift[j], bound, ndt, s);\n",
      "    }} else {{\n",
      "      log_lik += n[j] * ",
      "swald_lccdf(rt | drift[j], bound, ndt, s);\n",
      "    }}\n",
      "  }}\n",
      "  return log_lik;\n",
      "}}"
    )
  } else {
    glue(
      "real {family_name}_lpdf(real rt, real mu, {cat_args}, ",
      "real bound, real ndt, real s, real sp, int response, {n_args}) {{\n",
      "  real b = bound;\n",
      "  real A = sp * bound;\n",
      "  real t = rt - ndt;\n",
      "  if (t <= 0) return negative_infinity();\n",
      "  {drift_array}\n",
      "  {n_array}\n",
      "  real s_sqrt_t = s * sqrt(t);\n",
      "  // Winner PDF (Tillman et al. 2020, Eq. 5)\n",
      "  real v_win = drift[response];\n",
      "  real alpha_w = (b - A - t * v_win) / s_sqrt_t;\n",
      "  real beta_w = (b - t * v_win) / s_sqrt_t;\n",
      "  real pdf_val = (1.0 / A) * (\n",
      "    -v_win * Phi(alpha_w) + (s / sqrt(t)) * exp(std_normal_lpdf(alpha_w)) +\n",
      "     v_win * Phi(beta_w)  - (s / sqrt(t)) * exp(std_normal_lpdf(beta_w))\n",
      "  );\n",
      "  if (pdf_val <= 0) return negative_infinity();\n",
      "  real log_lik = log(n[response]) + log(pdf_val);\n",
      "  // Loser survivor functions (1 - CDF from Appendix A)\n",
      "  for (j in 1:{n_cats}) {{\n",
      "    real v_j = drift[j];\n",
      "    real alpha1 = (v_j * t - b) / s_sqrt_t;\n",
      "    real alpha2 = (v_j * t - (b - A)) / s_sqrt_t;\n",
      "    real beta1 = -(v_j * t + b) / s_sqrt_t;\n",
      "    real beta2 = -(v_j * t + (b - A)) / s_sqrt_t;\n",
      "    real cdf_val = (1.0 / (2 * v_j * A)) * (Phi(alpha2) - Phi(alpha1))\n",
      "      + (s * sqrt(t) / A) * (\n",
      "          alpha2 * Phi(alpha2) - alpha1 * Phi(alpha1)\n",
      "          + exp(std_normal_lpdf(alpha2)) - exp(std_normal_lpdf(alpha1))\n",
      "        )\n",
      "      - (1.0 / (2 * v_j * A)) * (\n",
      "          exp(2 * v_j * (b - A) / square(s)) * Phi(beta2)\n",
      "          - exp(2 * v_j * b / square(s)) * Phi(beta1)\n",
      "        );\n",
      "    if (cdf_val < 0) cdf_val = 0;\n",
      "    if (cdf_val > 1) cdf_val = 1;\n",
      "    real log_surv = log1m(cdf_val);\n",
      "    if (j == response) {{\n",
      "      if (n[j] > 1) log_lik += (n[j] - 1) * log_surv;\n",
      "    }} else {{\n",
      "      log_lik += n[j] * log_surv;\n",
      "    }}\n",
      "  }}\n",
      "  return log_lik;\n",
      "}}"
    )
  }
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.rdm_simple <- function(model, data, formula) {
  cat_names <- c("driftc", "drifte")
  has_sp <- !("sp" %in% names(model$fixed_parameters))
  formula <- bmf2bf(model, formula)

  formula$family <- brms::custom_family(
    "rdm_simple",
    dpars = c("mu", cat_names, "bound", "ndt", "s", "sp"),
    links = c("identity", model$links$driftc, model$links$drifte,
              model$links$bound, model$links$ndt, model$links$s,
              model$links$sp),
    ub = rep(NA, 7),
    lb = rep(NA, 7),
    type = "real",
    vars = c("vint1[n]", "vint2[n]", "vint3[n]"),
    loop = TRUE,
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
    scode = .rdm_stan_code("rdm_simple", cat_names, has_sp),
    block = "functions"
  )

  nlist(formula, data, stanvars)
}

#' @export
configure_model.rdm_custom <- function(model, data, formula) {
  cat_names <- model$other_vars$resp_cats
  n_cats <- length(cat_names)
  has_sp <- !("sp" %in% names(model$fixed_parameters))
  formula <- bmf2bf(model, formula)

  n_dpars <- n_cats + 5
  formula$family <- brms::custom_family(
    "rdm_custom",
    dpars = c("mu", cat_names, "bound", "ndt", "s", "sp"),
    links = c(
      "identity",
      vapply(cat_names, function(p) model$links[[p]], character(1)),
      model$links$bound, model$links$ndt, model$links$s, model$links$sp
    ),
    ub = rep(NA, n_dpars),
    lb = rep(NA, n_dpars),
    type = "real",
    vars = paste0("vint", seq_len(n_cats + 1), "[n]"),
    loop = TRUE,
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
    scode = .rdm_stan_code("rdm_custom", cat_names, has_sp),
    block = "functions"
  )

  nlist(formula, data, stanvars)
}

############################################################################# !
# Post-processing functions (shared helpers)                             ####
############################################################################# !

.rdm_log_lik <- function(i, prep, cat_names, n_cats) {
  rt <- prep$data$Y[i]
  response <- prep$data$vint1[i]
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sp_val <- brms::get_dpar(prep, "sp", i = i)

  t <- rt - ndt
  t[t <= 0] <- NA
  b <- bound
  A <- sp_val * bound

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )

  drift_win <- brms::get_dpar(prep, cat_names[response], i = i)

  use_full <- any(A >= 1e-10)

  if (!use_full) {
    log_lik <- log(n_cat[response]) +
      .dwald(t, drift = drift_win, bound = b, s = s, log = TRUE)

    if (n_cat[response] > 1) {
      log_lik <- log_lik + (n_cat[response] - 1) *
        .pwald(t, drift = drift_win, bound = b, s = s,
               lower.tail = FALSE, log.p = TRUE)
    }

    for (j in seq_len(n_cats)) {
      if (j == response) next
      drift_j <- brms::get_dpar(prep, cat_names[j], i = i)
      log_lik <- log_lik + n_cat[j] *
        .pwald(t, drift = drift_j, bound = b, s = s,
               lower.tail = FALSE, log.p = TRUE)
    }
  } else {
    log_lik <- log(n_cat[response]) +
      .dwald_full(t, drift = drift_win, bound = b, A = A, s = s, log = TRUE)

    if (n_cat[response] > 1) {
      log_lik <- log_lik + (n_cat[response] - 1) *
        .pwald_full(t, drift = drift_win, bound = b, A = A, s = s,
                    lower.tail = FALSE, log.p = TRUE)
    }

    for (j in seq_len(n_cats)) {
      if (j == response) next
      drift_j <- brms::get_dpar(prep, cat_names[j], i = i)
      log_lik <- log_lik + n_cat[j] *
        .pwald_full(t, drift = drift_j, bound = b, A = A, s = s,
                    lower.tail = FALSE, log.p = TRUE)
    }
  }

  log_lik[is.na(log_lik)] <- -Inf
  log_lik
}

.rdm_posterior_predict <- function(i, prep, cat_names, n_cats, ...) {
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sp_val <- brms::get_dpar(prep, "sp", i = i)
  n_draws <- length(ndt)

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )
  total_acc <- sum(n_cat)

  rt <- numeric(n_draws)
  for (d in seq_len(n_draws)) {
    b <- bound[d]
    A <- sp_val[d] * bound[d]
    drift_vec <- unlist(lapply(seq_len(n_cats), function(j) {
      rep(brms::get_dpar(prep, cat_names[j], i = i)[d], n_cat[j])
    }))

    if (A < 1e-10) {
      ft <- .rwald_ig(total_acc, drift = drift_vec, bound = b, s = s[d])
    } else {
      start <- stats::runif(total_acc, min = 0, max = A)
      ft <- .rwald_ig(total_acc, drift = drift_vec, bound = b - start,
                      s = s[d])
    }
    rt[d] <- min(ft) + ndt[d]
  }
  rt
}

.rdm_posterior_epred <- function(prep, cat_names, n_cats, ...) {
  n_obs <- prep$nobs
  n_draws <- prep$ndraws
  n_sim <- 100

  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (i in seq_len(n_obs)) {
    bound <- brms::get_dpar(prep, "bound", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    s <- brms::get_dpar(prep, "s", i = i)
    sp_val <- brms::get_dpar(prep, "sp", i = i)

    n_cat <- vapply(
      seq_len(n_cats),
      function(j) prep$data[[paste0("vint", j + 1)]][i],
      integer(1)
    )
    total_acc <- sum(n_cat)

    for (d in seq_len(n_draws)) {
      b <- bound[d]
      A <- sp_val[d] * bound[d]
      drift_vec <- unlist(lapply(seq_len(n_cats), function(j) {
        rep(brms::get_dpar(prep, cat_names[j], i = i)[d], n_cat[j])
      }))

      if (A < 1e-10) {
        ft <- matrix(
          .rwald_ig(total_acc * n_sim, drift = drift_vec, bound = b, s = s[d]),
          nrow = n_sim, ncol = total_acc, byrow = TRUE
        )
      } else {
        start <- matrix(stats::runif(total_acc * n_sim, min = 0, max = A),
                        nrow = n_sim, ncol = total_acc)
        ft <- matrix(NA_real_, nrow = n_sim, ncol = total_acc)
        for (k in seq_len(total_acc)) {
          ft[, k] <- .rwald_ig(n_sim, drift = drift_vec[k],
                               bound = b - start[, k], s = s[d])
        }
      }
      epred[d, i] <- mean(apply(ft, 1, min)) + ndt[d]
    }
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
    prep$family$dpars, c("mu", "bound", "ndt", "s", "sp")
  )
  .rdm_log_lik(i, prep, cat_names = cat_names, n_cats = length(cat_names))
}

posterior_predict_rdm_custom <- function(i, prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "bound", "ndt", "s", "sp")
  )
  .rdm_posterior_predict(i, prep, cat_names = cat_names,
                         n_cats = length(cat_names), ...)
}

posterior_epred_rdm_custom <- function(prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "bound", "ndt", "s", "sp")
  )
  .rdm_posterior_epred(prep, cat_names = cat_names,
                       n_cats = length(cat_names), ...)
}
