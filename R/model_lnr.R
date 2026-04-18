############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.lnr_version_table <- list(
  simple = list(
    parameters = list(
      correct = "meanlog for correct accumulator",
      error = "meanlog for error accumulators",
      ndt = "non-decision time",
      s = "sdlog (shared across accumulators)"
    ),
    links = list(
      correct = "identity",
      error = "identity",
      ndt = "log",
      s = "log"
    ),
    fixed_parameters = list(
      mu = 0,
      s = 0
    ),
    priors = list(
      correct = list(main = "normal(-1, 0.5)", effects = "normal(0, 0.3)"),
      error = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)"),
      ndt = list(main = "normal(-2, 0.3)", effects = "normal(0, 0.3)"),
      s = list(main = "normal(0, 0.3)", effects = "normal(0, 0.2)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      correct = c(-1.5, -0.5),
      error = c(-0.5, 0.5),
      ndt = c(0.025, 0.05),
      s = c(0.8, 1.2)
    )
  ),
  custom = list(
    parameters = list(
      ndt = "non-decision time",
      s = "sdlog (shared across accumulators)"
    ),
    links = list(
      ndt = "log",
      s = "log"
    ),
    fixed_parameters = list(
      mu = 0,
      s = 0
    ),
    priors = list(
      ndt = list(main = "normal(-2, 0.3)", effects = "normal(0, 0.3)"),
      s = list(main = "normal(0, 0.3)", effects = "normal(0, 0.2)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      ndt = c(0.025, 0.05),
      s = c(0.8, 1.2)
    )
  )
)

# Stan reserved words that cannot be used as category names
.stan_reserved <- c(
  "int", "real", "vector", "matrix", "array", "if", "else", "for", "while",
  "return", "void", "data", "model", "target", "print", "reject", "log",
  "exp", "lower", "upper", "in", "functions", "generated", "transformed",
  "parameters"
)

.lnr_reserved_dpars <- c("mu", "ndt", "s")


.model_lnr <- function(
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
      name = "Log-Normal Race Model",
      citation = "Rouder, J. N., Province, J. M., Morey, R. D., Gomez, P., & Heathcote, A. (2015).
        The Lognormal Race: A Cognitive-Process Model of Choice and Latency with Desirable
        Psychometric Properties. Psychometrika, 80(2), 491-513. https://doi.org/10.1007/s11336-013-9396-3",
      version = version,
      requirements = glue(
        "- Reaction times should be passed in seconds", "\n",
        "- For version 'simple': response variable should be integer-coded ",
        "(1 = correct, 2:K = errors)", "\n",
        "- For version 'custom': response variable should contain character ",
        "labels matching formula parameter names"
      ),
      parameters = .lnr_version_table[[version]][["parameters"]],
      links = .lnr_version_table[[version]][["links"]],
      fixed_parameters = .lnr_version_table[[version]][["fixed_parameters"]],
      default_priors = .lnr_version_table[[version]][["priors"]],
      init_ranges = .lnr_version_table[[version]][["init_ranges"]],
      void_mu = TRUE
    ),
    class = c("bmmodel", "lnr", paste0("lnr_", version)),
    call = call
  )

  out$links[names(links)] <- links
  out
}

#' @title `r .model_lnr()$name`
#' @name lnr
#' @details `r model_info(.model_lnr())`
#' @param rt The name of the variable in the dataset containing the response
#'   times. Response times should be coded in seconds (not milliseconds).
#' @param response The name of the variable in the dataset containing the
#'   response/choice. For the `"simple"` version, responses should be
#'   integer-coded: 1 = correct response, 2 through K = error responses.
#'   Factor and character-digit responses are accepted and converted
#'   automatically. For the `"custom"` version, responses should be character
#'   or factor labels matching the accumulator names in the formula. Category
#'   names must not use reserved internal parameter names such as `"mu"`,
#'   `"ndt"`, or `"s"`.
#' @param n_alternatives An integer specifying the total number of response
#'   alternatives (K >= 2). Required for `version = "simple"`. Not used for
#'   `version = "custom"` (inferred from the formula).
#' @param num_alternatives For `version = "custom"` only. A named vector
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
#'       accumulators. The sdlog parameter `s` is shared and fixed by default
#'       (s = 1). This covers the common case where interest is in the speed
#'       of correct vs. error processing.
#'     \item `"custom"`: Per-category meanlog parameters. Response categories
#'       are defined by the formula LHS names (e.g., `correct ~ 1, other ~ 1,
#'       npl ~ 1`). The response column must contain character labels matching
#'       these names. Category names must not be `"mu"`, `"ndt"`, `"s"`,
#'       Stan reserved words, or names ending in a number. Supports
#'       per-category `num_alternatives`.
#'   }
#' @param links A named list of link functions for the model parameters.
#'   For `"simple"`: parameters are `correct`, `error`, `ndt`, and `s`.
#'   Default links are "identity" for `correct` and `error`, "log" for
#'   `ndt` and `s`. For `"custom"`: category parameters default to "identity".
#' @param ... Additional arguments passed internally (for testing purposes).
#' @return An object of class `bmmodel`
#' @export
#' @keywords bmmodel
#' @seealso [dlnr()] and [rlnr()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simple version with 2 alternatives
#' dat <- rlnr(n = 500, m = c(-1, 0), s = c(1, 1), ndt = 0.2)
#' model <- lnr(rt = "rt", response = "response", n_alternatives = 2)
#' formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)
#' fit <- bmm(formula, dat, model, cores = 4, backend = "cmdstanr")
#'
#' # custom version with named categories
#' model2 <- lnr(rt = "rt", response = "resp", version = "custom",
#'               num_alternatives = c(target = 1, similar = 3, other = 5))
#' formula2 <- bmf(target ~ 1, similar ~ 1, other ~ 1, ndt ~ 1)
lnr <- function(rt, response, n_alternatives = NULL,
                version = c("simple", "custom"),
                num_alternatives = NULL, links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  if (version == "simple") {
    stopif(
      !is.null(num_alternatives),
      "num_alternatives is only supported for version 'custom'."
    )
    stopif(
      is.null(n_alternatives) || !is.numeric(n_alternatives) ||
        n_alternatives < 2 || n_alternatives != round(n_alternatives),
      "n_alternatives must be an integer >= 2 for version 'simple'."
    )
    n_alternatives <- as.integer(n_alternatives)
  } else {
    stopif(
      !is.null(n_alternatives),
      "n_alternatives is only supported for version 'simple'. Use num_alternatives for version 'custom'."
    )
    n_alternatives <- NULL
  }
  .model_lnr(
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
check_model.lnr_custom <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    reserved_pars <- c("ndt", "s")
    formula_pars <- names(formula)
    cat_pars <- setdiff(formula_pars, reserved_pars)

    stopif(
      length(cat_pars) == 0,
      "Custom version requires at least one accumulator parameter in the formula."
    )

    bad_internal_names <- cat_pars[tolower(cat_pars) %in% .lnr_reserved_dpars]
    stopif(
      length(bad_internal_names) > 0,
      "Category names cannot use reserved internal parameter names: {collapse_comma(bad_internal_names)}."
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

    for (p in cat_pars) {
      model$parameters[[p]] <- paste0("meanlog for '", p, "' accumulator")
      if (is.null(model$links[[p]])) model$links[[p]] <- "identity"
      if (is.null(model$default_priors[[p]])) {
        model$default_priors[[p]] <- list(
          main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
        )
      }
      if (is.null(model$init_ranges[[p]])) {
        model$init_ranges[[p]] <- c(-0.5, 0.5)
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
check_data.lnr_simple <- function(model, data, formula) {
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

  data$.lnr_cat <- ifelse(data[, response_var] == 1L, 1L, 2L)
  data$.lnr_n1 <- 1L
  data$.lnr_n2 <- n_alt - 1L

  NextMethod("check_data")
}

#' @export
check_data.lnr_custom <- function(model, data, formula) {
  response_var <- model$resp_vars$response
  cat_names <- model$other_vars$resp_cats
  num_alt <- model$other_vars$num_alternatives

  if (is.factor(data[, response_var])) {
    data[, response_var] <- as.character(data[, response_var])
  }

  stopif(
    !is.character(data[, response_var]),
    "For version 'custom', the response variable '{response_var}' must \\
    contain character labels matching the formula category names."
  )

  data_levels <- unique(data[, response_var])
  bad_levels <- data_levels[tolower(data_levels) %in% .lnr_reserved_dpars]
  stopif(
    length(bad_levels) > 0,
    "Response levels cannot use reserved internal parameter names: {collapse_comma(bad_levels)}."
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

  data$.lnr_cat <- setNames(seq_along(cat_names), cat_names)[data[, response_var]]

  if (is.null(num_alt)) {
    for (i in seq_along(cat_names)) {
      data[[paste0(".lnr_n", i)]] <- 1L
    }
  } else if (is.numeric(num_alt)) {
    stopif(
      is.null(names(num_alt)) || any(names(num_alt) == "") ||
        anyDuplicated(names(num_alt)) > 0,
      "num_alternatives must be a uniquely named vector with one entry for each formula category: \\
      {collapse_comma(cat_names)}"
    )
    missing_cats <- setdiff(cat_names, names(num_alt))
    extra_cats <- setdiff(names(num_alt), cat_names)
    stopif(
      length(missing_cats) > 0 || length(extra_cats) > 0,
      "num_alternatives must have exactly the formula categories: \\
      {collapse_comma(cat_names)}"
    )
    invalid_num_alt <- num_alt[
      !is.finite(num_alt) | num_alt < 1 | num_alt != round(num_alt)
    ]
    stopif(
      length(invalid_num_alt) > 0,
      "num_alternatives must contain positive integers for each formula category. Invalid value(s): \\
      {collapse_comma(glue::glue('{names(invalid_num_alt)} = {invalid_num_alt}'))}"
    )
    for (i in seq_along(cat_names)) {
      data[[paste0(".lnr_n", i)]] <- as.integer(num_alt[cat_names[i]])
    }
  } else if (is.character(num_alt)) {
    stopif(
      is.null(names(num_alt)) || any(names(num_alt) == "") ||
        anyDuplicated(names(num_alt)) > 0,
      "num_alternatives must be a uniquely named vector with one entry for each formula category: \\
      {collapse_comma(cat_names)}"
    )
    missing_cats <- setdiff(cat_names, names(num_alt))
    extra_cats <- setdiff(names(num_alt), cat_names)
    stopif(
      length(missing_cats) > 0 || length(extra_cats) > 0,
      "num_alternatives must have exactly the formula categories: \\
      {collapse_comma(cat_names)}"
    )
    missing_cols <- setdiff(num_alt, colnames(data))
    stopif(
      length(missing_cols) > 0,
      "num_alternatives columns {collapse_comma(missing_cols)} not found \\
      in the data."
    )
    for (i in seq_along(cat_names)) {
      col_name <- num_alt[cat_names[i]]
      col_vals <- data[, col_name]
      stopif(
        !is.numeric(col_vals),
        "num_alternatives column '{col_name}' must be numeric."
      )
      stopif(
        anyNA(col_vals),
        "num_alternatives column '{col_name}' contains NA values."
      )
      stopif(
        any(!is.finite(col_vals)),
        "num_alternatives column '{col_name}' contains non-finite values."
      )
      stopif(
        any(col_vals < 1 | col_vals != round(col_vals)),
        "num_alternatives column '{col_name}' must contain integers >= 1."
      )
      data[[paste0(".lnr_n", i)]] <- as.integer(col_vals)
    }
  } else {
    stop2(
      "num_alternatives must be NULL, a named numeric vector of positive integers, \\
      or a named character vector of column names."
    )
  }

  NextMethod("check_data")
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.lnr_simple <- function(model, formula) {
  rt_var <- model$resp_vars$rt
  brms::bf(glue("{rt_var} | vint(.lnr_cat, .lnr_n1, .lnr_n2) ~ 1"))
}

#' @export
bmf2bf.lnr_custom <- function(model, formula) {
  rt_var <- model$resp_vars$rt
  cat_names <- model$other_vars$resp_cats
  n_cols <- paste0(".lnr_n", seq_along(cat_names))
  vint_args <- paste(c(".lnr_cat", n_cols), collapse = ", ")
  brms::bf(glue("{rt_var} | vint({vint_args}) ~ 1"))
}

############################################################################# !
# Stan code generation                                                   ####
############################################################################# !

.lnr_stan_code <- function(family_name, cat_names) {
  n_cats <- length(cat_names)
  cat_args <- paste(paste0("real ", cat_names), collapse = ", ")
  n_args <- paste(paste0("int n", seq_len(n_cats)), collapse = ", ")
  m_array <- paste0(
    "array[", n_cats, "] real m = {",
    paste(cat_names, collapse = ", "), "};"
  )
  n_array <- paste0(
    "array[", n_cats, "] int n = {",
    paste(paste0("n", seq_len(n_cats)), collapse = ", "), "};"
  )

  glue(
    "real {family_name}_lpdf(real rt, real mu, {cat_args}, ",
    "real ndt, real s, int response, {n_args}) {{\n",
    "  real t = rt - ndt;\n",
    "  if (t <= 0) return negative_infinity();\n",
    "  {m_array}\n",
    "  {n_array}\n",
    "  real log_lik = log(n[response]) + ",
    "lognormal_lpdf(t | m[response], s);\n",
    "  for (j in 1:{n_cats}) {{\n",
    "    if (j == response) {{\n",
    "      if (n[j] > 1)\n",
    "        log_lik += (n[j] - 1) * lognormal_lccdf(t | m[j], s);\n",
    "    }} else {{\n",
    "      log_lik += n[j] * lognormal_lccdf(t | m[j], s);\n",
    "    }}\n",
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
    vars = c("vint1[n]", "vint2[n]", "vint3[n]"),
    loop = TRUE,
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

  # brms vint() order must match Stan function signature
  vars_vec <- paste0("vint", seq_len(n_cats + 1), "[n]")

  formula$family <- brms::custom_family(
    "lnr_custom",
    dpars = dpars,
    links = link_vec,
    ub = rep(NA, length(dpars)),
    lb = lb_vec,
    type = "real",
    vars = vars_vec,
    loop = TRUE,
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

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )

  m_win <- brms::get_dpar(prep, cat_names[response], i = i)
  log_lik <- log(n_cat[response]) +
    stats::dlnorm(t, meanlog = m_win, sdlog = s, log = TRUE)

  if (n_cat[response] > 1) {
    log_lik <- log_lik + (n_cat[response] - 1) *
      stats::plnorm(t, meanlog = m_win, sdlog = s,
                    lower.tail = FALSE, log.p = TRUE)
  }

  for (j in seq_len(n_cats)) {
    if (j == response) next
    m_j <- brms::get_dpar(prep, cat_names[j], i = i)
    log_lik <- log_lik + n_cat[j] *
      stats::plnorm(t, meanlog = m_j, sdlog = s,
                    lower.tail = FALSE, log.p = TRUE)
  }

  log_lik[is.na(log_lik)] <- -Inf
  log_lik
}

.lnr_posterior_predict <- function(i, prep, cat_names, n_cats, ...) {
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  n_draws <- length(ndt)

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )
  total_acc <- sum(n_cat)

  rt <- numeric(n_draws)
  for (d in seq_len(n_draws)) {
    m_vec <- unlist(lapply(seq_len(n_cats), function(j) {
      rep(brms::get_dpar(prep, cat_names[j], i = i)[d], n_cat[j])
    }))
    s_vec <- rep(s[d], total_acc)
    ft <- stats::rlnorm(total_acc, meanlog = m_vec, sdlog = s_vec)
    rt[d] <- min(ft) + ndt[d]
  }
  rt
}

.lnr_posterior_epred <- function(prep, cat_names, n_cats, ...) {
  n_obs <- prep$nobs
  n_draws <- prep$ndraws
  n_sim <- 100

  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (i in seq_len(n_obs)) {
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    s <- brms::get_dpar(prep, "s", i = i)

    n_cat <- vapply(
      seq_len(n_cats),
      function(j) prep$data[[paste0("vint", j + 1)]][i],
      integer(1)
    )
    total_acc <- sum(n_cat)

    for (d in seq_len(n_draws)) {
      m_vec <- unlist(lapply(seq_len(n_cats), function(j) {
        rep(brms::get_dpar(prep, cat_names[j], i = i)[d], n_cat[j])
      }))
      s_vec <- rep(s[d], total_acc)
      ft <- matrix(stats::rlnorm(total_acc * n_sim, meanlog = m_vec,
                                  sdlog = s_vec),
                   nrow = n_sim, ncol = total_acc, byrow = TRUE)
      epred[d, i] <- mean(apply(ft, 1, min)) + ndt[d]
    }
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

posterior_epred_lnr_simple <- function(prep, ...) {
  .lnr_posterior_epred(prep, cat_names = c("correct", "error"),
                       n_cats = 2, ...)
}

posterior_epred_lnr_custom <- function(prep, ...) {
  cat_names <- setdiff(
    prep$family$dpars, c("mu", "ndt", "s")
  )
  .lnr_posterior_epred(prep, cat_names = cat_names,
                       n_cats = length(cat_names), ...)
}
