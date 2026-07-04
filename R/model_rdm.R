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

# Compose the spec for one version. Drift parameters precede the shared block:
# downstream code recovers accumulator names via
# setdiff(names(parameters), c("gap", "ndt", "s", "sp")), which relies on order.
# sp = -100 on the log scale fixes the starting point to ~0 by default; adding
# sp to the formula frees it and switches the likelihood to the full-Wald forms.
.rdm_model_spec <- function(version) {
  fixed_parameters <- list(mu = 0, s = 0, sp = -100)

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

.rdm_stan_reserved <- c(
  "int", "real", "vector", "matrix", "array", "if", "else", "for", "while",
  "return", "void", "data", "model", "target", "print", "reject", "log",
  "exp", "lower", "upper", "in", "functions", "generated", "transformed",
  "parameters"
)

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
      init_ranges = vt[["init_ranges"]],
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
#'       default (s = 1). Starting point `sp` is fixed to ~0 by
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
check_model.rdm <- function(model, data = NULL, formula = NULL) {
  positive_pars <- setdiff(names(model$links), "mu")
  bad_positive <- positive_pars[vapply(
    positive_pars,
    function(par) !identical(model$links[[par]], "log"),
    logical(1)
  )]

  stopif(
    length(bad_positive) > 0,
    "RDM parameters {collapse_comma(bad_positive)} only support the 'log' link."
  )

  NextMethod("check_model")
}

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

    bad_names <- intersect(tolower(cat_pars), .rdm_stan_reserved)
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
    for (i in seq_along(cat_names)) {
      data[[paste0(".rdm_n", i)]] <- as.integer(num_alt[cat_names[i]])
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
      data[[paste0(".rdm_n", i)]] <- as.integer(data[, num_alt[cat_names[i]]])
    }
  }

  model$other_vars$n_choices <- n_cats
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
  use_start_var <- if (has_sp) 1 else 0

  cat_args <- paste(paste0("vector ", cat_names), collapse = ", ")
  n_args <- paste(paste0("array[] int n", seq_len(n_cats)), collapse = ", ")
  drift_array <- paste0(
    "      array[", n_cats, "] real drift_i = {",
    paste0(cat_names, "[i]", collapse = ", "), "};\n"
  )
  n_array <- paste0(
    "      array[", n_cats, "] int n_i = {",
    paste(paste0("n", seq_len(n_cats), "[i]"), collapse = ", "), "};\n"
  )

  glue(
    "real {family_name}_lpdf(vector rt, vector mu, {cat_args}, ",
    "vector gap, vector ndt, vector s, vector sp, array[] int response, {n_args}) {{\n",
    "  int N = rows(rt);\n",
    "  real log_lik = 0;\n",
    "  for (i in 1:N) {{\n",
    "{drift_array}",
    "{n_array}",
    "    log_lik += rdm_log_lik_one(\n",
    "      rt[i], drift_i, gap[i], ndt[i], s[i], sp[i], response[i], n_i, {use_start_var});\n",
    "  }}\n",
    "  return log_lik;\n",
    "}}"
  )
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
    dpars = c("mu", cat_names, "gap", "ndt", "s", "sp"),
    links = c("identity", model$links$driftc, model$links$drifte,
              model$links$gap, model$links$ndt, model$links$s,
              model$links$sp),
    ub = rep(NA, 7),
    lb = c(NA, 0, 0, 0, 0, 0, 0),
    type = "real",
    vars = c("vint1", "vint2", "vint3"),
    loop = FALSE,
    log_lik = log_lik_rdm_simple,
    posterior_predict = posterior_predict_rdm_simple,
    posterior_epred = posterior_epred_rdm_simple
  )
  formula$family$rdm_has_sp <- has_sp

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
    dpars = c("mu", cat_names, "gap", "ndt", "s", "sp"),
    links = c(
      "identity",
      vapply(cat_names, function(p) model$links[[p]], character(1)),
      model$links$gap, model$links$ndt, model$links$s, model$links$sp
    ),
    ub = rep(NA, n_dpars),
    lb = c(NA, rep(0, n_cats), 0, 0, 0, 0),
    type = "real",
    vars = paste0("vint", seq_len(n_cats + 1)),
    loop = FALSE,
    log_lik = log_lik_rdm_custom,
    posterior_predict = posterior_predict_rdm_custom,
    posterior_epred = posterior_epred_rdm_custom
  )
  formula$family$rdm_has_sp <- has_sp

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
  gap <- brms::get_dpar(prep, "gap", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sp_val <- brms::get_dpar(prep, "sp", i = i)
  has_sp <- isTRUE(prep$family$rdm_has_sp)

  t <- rt - ndt
  t[t <= 0] <- NA
  b <- if (has_sp) gap + sp_val else gap
  A <- if (has_sp) sp_val else 0

  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )

  drift_win <- brms::get_dpar(prep, cat_names[response], i = i)

  if (!has_sp) {
    log_lik <- log(n_cat[response]) +
      .dwald(t, drift = drift_win, bound = b, s = s, log = TRUE)

    if (n_cat[response] > 1) {
      log_lik <- log_lik + (n_cat[response] - 1) *
        .pwald(t, drift = drift_win, bound = b, s = s,
               lower.tail = FALSE, log.p = TRUE)
    }

    for (j in seq_len(n_cats)) {
      if (j == response || n_cat[j] == 0) next
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
      if (j == response || n_cat[j] == 0) next
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
  gap <- brms::get_dpar(prep, "gap", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sp <- brms::get_dpar(prep, "sp", i = i)
  n_draws <- length(ndt)
  has_sp <- isTRUE(prep$family$rdm_has_sp)
  drift <- lapply(cat_names, function(p) brms::get_dpar(prep, p, i = i))
  n_cat <- vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )

  b <- if (has_sp) gap + sp else gap
  min_ft <- rep(Inf, n_draws)
  for (j in seq_len(n_cats)) {
    if (n_cat[j] == 0) next
    for (k in seq_len(n_cat[j])) {
      bound <- if (has_sp) b - stats::runif(n_draws, 0, sp) else b
      min_ft <- pmin(
        min_ft, .rwald_ig(n_draws, drift = drift[[j]], bound = bound, s = s)
      )
    }
  }
  min_ft + ndt
}

.rdm_posterior_epred <- function(prep, cat_names, n_cats, ...) {
  n_obs <- prep$nobs
  n_draws <- prep$ndraws
  n_sim <- 100L
  has_sp <- isTRUE(prep$family$rdm_has_sp)

  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (i in seq_len(n_obs)) {
    gap <- brms::get_dpar(prep, "gap", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    s <- brms::get_dpar(prep, "s", i = i)
    sp <- brms::get_dpar(prep, "sp", i = i)
    drift <- lapply(cat_names, function(p) brms::get_dpar(prep, p, i = i))
    n_cat <- vapply(
      seq_len(n_cats),
      function(j) prep$data[[paste0("vint", j + 1)]][i],
      integer(1)
    )

    b_m <- matrix(if (has_sp) gap + sp else gap, n_draws, n_sim)
    s_m <- matrix(s, n_draws, n_sim)
    min_ft <- matrix(Inf, n_draws, n_sim)
    for (j in seq_len(n_cats)) {
      if (n_cat[j] == 0) next
      dj <- matrix(drift[[j]], n_draws, n_sim)
      for (k in seq_len(n_cat[j])) {
        bound <- if (has_sp) {
          b_m - matrix(stats::runif(n_draws * n_sim), n_draws, n_sim) *
            matrix(sp, n_draws, n_sim)
        } else {
          b_m
        }
        ft <- .rwald_ig(n_draws * n_sim, drift = dj, bound = bound, s = s_m)
        dim(ft) <- c(n_draws, n_sim)
        min_ft <- pmin(min_ft, ft)
      }
    }
    epred[, i] <- rowMeans(min_ft) + ndt
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
