############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.cswald_version_table <- list(
  simple = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary (distance from starting point to correct boundary)",
      ndt = "non-decision time (minimum non-decision time when sndt > 0)",
      s = "diffusion constant",
      sndt = "range of the uniform trial-to-trial variability in the non-decision time"
    ),
    links = list(
      drift = "log",
      bound = "log",
      ndt = "log",
      s = "log",
      sndt = "log"
    ),
    links_fixed = list(sndt = "identity"),
    fixed_parameters = list(
      mu = 0,
      s = 0,
      sndt = 0
    ),
    priors = list(
      drift = list(main = "normal(0,1)", effects = "normal(0,0.3)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)"),
      s = list(main = "normal(0,0.3)", effects = "normal(0,0.2)"),
      sndt = list(main = "normal(-2.5,1)", effects = "normal(0,0.3)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      drift = c(1, 2),
      bound = c(1.5, 2),
      ndt = c(0.025, 0.05),
      s = c(0.95, 1.05),
      sndt = c(0.01, 0.05)
    )
  ),
  crisk = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary separation (total distance between boundaries)",
      ndt = "non-decision time (minimum non-decision time when sndt > 0)",
      zr = "relative starting point",
      s = "diffusion constant",
      sndt = "range of the uniform trial-to-trial variability in the non-decision time"
    ),
    links = list(
      drift = "identity",
      bound = "log",
      ndt = "log",
      zr = "logit",
      s = "log",
      sndt = "log"
    ),
    links_fixed = list(sndt = "identity"),
    fixed_parameters = list(
      mu = 0,
      zr = 0,
      s = 0,
      sndt = 0
    ),
    priors = list(
      drift = list(main = "normal(0,1)", effects = "normal(0,0.5)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)"),
      zr = list(main = "normal(0,0.3)", effects = "normal(0,0.2)"),
      s = list(main = "normal(0,0.5)", effects = "normal(0,0.2)"),
      sndt = list(main = "normal(-2.5,1)", effects = "normal(0,0.3)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      drift = c(-0.5, 0.5),
      bound = c(1.5, 2),
      ndt = c(0.025, 0.05),
      zr = c(0.45, 0.55),
      s = c(0.95, 1.05),
      sndt = c(0.01, 0.05)
    )
  )
)


.model_cswald <- function(
    rt = NULL,
    response = NULL,
    links = NULL,
    version = "simple",
    call = NULL,
    ...) {
  out <- structure(
    list(
      resp_vars = nlist(rt, response),
      other_vars = list(),
      domain = "Decision Making / Response times",
      task = "Choice Reaction Time tasks (with few errors)",
      name = "Censored-Shifted Wald Model",
      citation = "Miller, R., Scherbaum, S., Heck, D. W., Goschke, T., & Enge, S. (2017).
        On the Relation Between the (Censored) Shifted Wald and the Wiener Distribution as Measurement Models
        for Choice Response Times. Applied Psychological Measurement, 42(2), 116-135. https://doi.org/10.1177/0146621617710465",
      version = version,
      requirements = glue(
        "- Reaction times should be passed in seconds", "\n",
        "- The response variable should be passed numerically: 0 = lower response, 1 = upper response"
      ),
      parameters = .cswald_version_table[[version]][["parameters"]],
      links = .cswald_version_table[[version]][["links"]],
      links_fixed = .cswald_version_table[[version]][["links_fixed"]],
      fixed_parameters = .cswald_version_table[[version]][["fixed_parameters"]],
      default_priors = .cswald_version_table[[version]][["priors"]],
      init_ranges = .cswald_version_table[[version]][["init_ranges"]]
    ),
    class = c("bmmodel", "cswald", paste0("cswald_", version)),
    call = call
  )

  out$links[names(links)] <- links
  resolve_fixed_links(out)
}

#' @title `r .model_cswald()$name`
#' @name cswald
#' @details `r model_info(.model_cswald())`
#' @param rt The name of the variable in the dataset containing the response
#'   times. Response times should be coded in seconds (not milliseconds).
#' @param response The name of the variable in the dataset containing the
#'   response/decision. Responses should be coded as 0 (lower boundary) or
#'   1 (upper boundary). Alternatively, character values "lower" and "upper"
#'   or logical values (FALSE/TRUE) are accepted and will be converted
#'   automatically.
#' @param links A named list of link functions for the model parameters.
#'   Available parameters depend on the version: "simple" has `drift`, `bound`,
#'   `ndt`, `s`, and `sndt`; "crisk" additionally has `zr`. Default links are
#'   "log" for most parameters and "logit" for `zr`. For positive parameters, "softplus"
#'   is available as an alternative to "log" that grows linearly for large
#'   values and avoids the numerical blow-up of `exp()`.
#' @param version A character string specifying which version of the cswald
#'   model to use. Options are:
#'   \itemize{
#'     \item `"simple"` (default): The standard censored shifted Wald model,
#'       which treats error responses as censored correct responses. Best suited
#'       for tasks with few errors (<20%). **Note:** The `bound` parameter in
#'       the simple version represents the distance from the starting point to
#'       the correct boundary, which is half the total boundary separation
#'       in the diffusion model (assuming an unbiased starting point). To
#'       convert to the full boundary separation (as in DDM or crisk), multiply
#'       by 2.
#'     \item `"crisk"`: The competing risks version, which models both response
#'       types as arising from racing accumulators toward opposite boundaries.
#'       Better suited for tasks with substantial error rates. The `bound`
#'       parameter represents the total boundary separation, consistent with
#'       the diffusion model parameterization.
#'   }
#'   For more details, see Miller et al. (2017).
#' @param ... Additional arguments passed internally (for testing purposes).
#' @return An object of class `bmmodel`
#' @section Trial-to-trial variability in the non-decision time:
#'
#'   Both versions have an optional parameter `sndt`: the non-decision time is
#'   `Uniform(ndt, ndt + sndt)`, so `ndt` is the *minimum* non-decision time and
#'   `ndt + sndt/2` the mean. This is the `st0` of [rtdists::rdiffusion()];
#'   fast-dm centers the same uniform on its `t0`, which is then the mean.
#'
#'   `sndt` is fixed at 0 by default, which reproduces the standard censored
#'   shifted Wald model. Estimate it with a formula (`bmf(..., sndt ~ 1)`) or
#'   fix it in seconds (`bmf(..., sndt = 0.15)`); fixed values are on the
#'   natural scale, and `sndt` uses the log link only once it is estimated.
#'
#'   Without `sndt`, the fastest responses cap the `ndt` estimate (the
#'   likelihood requires `ndt < min(rt)`), which biases `ndt`, `drift`, and
#'   `bound` when the true non-decision time varies (Miller et al., 2017,
#'   Fig. 3). Estimating `sndt` removes this bias, but `sndt` itself is weakly
#'   identified at typical trial numbers and the prior acts as its regularizer:
#'   estimate it at the population level only (`sndt ~ 1`) and check prior
#'   sensitivity.
#'
#'   In the `"crisk"` version with `sndt > 0`, each accumulator draws its own
#'   non-decision time instead of sharing one draw per trial. Choice
#'   probabilities agree exactly with the shared-draw model at `zr = 0.5`; away
#'   from it they differ by up to ~1 percentage point at `sndt = 0.3` (~4 at
#'   `zr = 0.2`), and densities by up to ~10%. `posterior_predict()` simulates
#'   through [rtdists::rdiffusion()], which shares the draw, so `pp_check()`
#'   compares against a slightly different model whenever `zr` departs from 0.5.
#'
#' @section Estimating `sndt` in the "crisk" version at substantial error rates:
#'
#'   The crisk race only approximates a single diffusion process, and the
#'   approximation degrades as errors become frequent. A free `sndt` can absorb
#'   that approximation error rather than genuine non-decision-time
#'   variability: on diffusion-generated data with 10-17% errors (true `bound`
#'   1.6, `sndt` 0.2), the crisk fit converged to `bound` ~1.2 and `sndt` ~0.4
#'   even with very large samples, while data from its own racing process were
#'   recovered without bias. A shared-draw formulation gives the same estimates,
#'   so this is a property of the crisk approximation itself. Options:
#'   \itemize{
#'     \item Below roughly 10% errors, use the `"simple"` version, where `sndt`
#'       recovery is unbiased.
#'     \item At substantial error rates, prefer the `ddm` model (exact, but
#'       currently without non-decision-time variability) or keep `sndt` fixed
#'       at 0, where the crisk estimates stay close to the diffusion values.
#'     \item Fix `sndt` at a plausible value (e.g. `bmf(..., sndt = 0.1)`) to
#'       bound its influence while letting the mean non-decision time exceed
#'       the fastest response.
#'     \item Treat a free `sndt` that comes out much larger, with a much smaller
#'       `bound`, than in a fit with `sndt` fixed or a `ddm` fit as a sign of
#'       absorbed approximation error, and check with `pp_check()` before
#'       interpreting the parameters.
#'   }
#' @export
#' @keywords bmmodel
#' @seealso [dcswald()] and [rcswald()] for the density and random generation
#'   functions.
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # generate simulated data from the diffusion model
#' dat <- rcswald(n = 500, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5, s = 1)
#'
#' # specify the model
#' model <- cswald(rt = "rt", response = "response", version = "simple")
#'
#' # specify the formula
#' formula <- bmf(
#'   drift ~ 1,
#'   bound ~ 1,
#'   ndt ~ 1
#' )
#'
#' # fit the model
#' fit <- bmm(
#'   formula = formula,
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
cswald <- function(rt, response, links = NULL, version = c("simple", "crisk"), ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  .model_cswald(
    rt = rt,
    response = response,
    links = links,
    version = version,
    call = call,
    ...
  )
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !

#' @export
check_data.cswald <- function(model, data, formula) {
  rt_var <- model$resp_vars$rt
  response_var <- model$resp_vars$response

  stopif(
    not_in(rt_var, colnames(data)),
    "The RT variable '{rt_var}' is not present in the data."
  )

  stopif(
    not_in(response_var, colnames(data)),
    "The response variable '{response_var}' is not present in the data."
  )

  n_na_rt <- sum(is.na(data[, rt_var]))
  n_na_resp <- sum(is.na(data[, response_var]))

  stopif(
    n_na_rt > 0,
    "The RT variable '{rt_var}' contains {n_na_rt} NA values. \\
    Please remove or impute missing values before fitting the model."
  )

  stopif(
    n_na_resp > 0,
    "The response variable '{response_var}' contains {n_na_resp} NA values. \\
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

  resp_type <- typeof(data[, response_var])

  if (is.factor(data[, response_var])) {
    data[, response_var] <- as.character(data[, response_var])
    resp_type <- "character"
  }

  if (resp_type %in% c("integer", "double")) {
    stopif(
      any(!data[, response_var] %in% c(0, 1)),
      "The response variable '{response_var}' contains values other than 0 and 1.\n
      Please pass responses coded as 0 (lower boundary) and 1 (upper boundary)."
    )
  } else if (resp_type == "logical") {
    warning2(
      "The response variable is boolean and will be internally transformed ",
      "to an integer variable with values 0 for FALSE and 1 for TRUE."
    )
    data[, response_var] <- as.integer(data[, response_var])
  } else if (resp_type == "character") {
    data[, response_var] <- tolower(data[, response_var])
    stopif(
      any(!data[, response_var] %in% c("upper", "lower")),
      "The response variable '{response_var}' contains invalid character values.\n
      Please pass only 'upper' or 'lower' as response values, or use \\
      numeric coding (0 = lower, 1 = upper)."
    )
    warning2(
      "The response variable is a character variable and will be internally ",
      "transformed to an integer variable with 0 for 'lower' and 1 for 'upper'."
    )
    data[, response_var] <- ifelse(data[, response_var] == "upper", 1L, 0L)
  } else {
    stop2(
      "The response variable '{response_var}' is of type '{resp_type}'.\n
      Please provide responses as integer (0/1), logical, character \\
      ('upper'/'lower'), or factor."
    )
  }

  if (model$version == "simple") {
    error_rate <- mean(data[, response_var] == 0)
    warnif(
      error_rate > 0.20,
      "Your data has an error rate of {round(error_rate * 100, 1)}%.\n
      The simple censored shifted Wald model assumes few errors. \\
      Consider using version = 'crisk' (competing risks) for data with \\
      substantial error rates."
    )
  }

  NextMethod("check_data")
}

############################################################################# !
# Convert bmmformula to brmsformla methods                               ####
############################################################################# !

#' @export
bmf2bf.cswald <- function(model, formula) {
  rt_var <- model$resp_vars$rt
  response_var <- model$resp_vars$response
  brms::bf(glue(rt_var, " | dec(", response_var, ") ~ 1"))
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

# Only sndt == 0 collapses the likelihood to the closed form the vectorized
# (loop = FALSE) overload evaluates. Any other value, fixed or estimated, needs
# the per-observation convolution, so the test is on the value, not on whether
# sndt is fixed.
#
# brms slices Y per thread but pastes a custom family's `vars` in unsliced, so
# under threading the family must emit "dec[start:end]" itself. start/end only
# exist in threaded Stan code, and threading(force = TRUE) compiles threaded but
# keeps the serial likelihood, so slice only when brms will really thread.
cswald_family_args <- function(model) {
  if (!isTRUE(model$fixed_parameters[["sndt"]] == 0)) {
    return(list(loop = TRUE, vars = "dec[n]"))
  }
  threads <- getOption("brms.threads", NULL)
  # brms also accepts a bare number for this option
  if (is.numeric(threads)) {
    threads <- brms::threading(threads)
  }
  threaded <- is.list(threads) && isTRUE(threads$threads > 0) && !isTRUE(threads$force)
  list(loop = FALSE, vars = if (threaded) "dec[start:end]" else "dec")
}

#' @export
configure_model.cswald_simple <- function(model, data, formula) {
  links <- model$links
  formula <- bmf2bf(model, formula)

  family_args <- cswald_family_args(model)

  formula$family <- brms::custom_family(
    "cswald",
    dpars = c("mu", "drift", "bound", "ndt", "s", "sndt"),
    links = c("identity", links$drift, links$bound, links$ndt, links$s, links$sndt),
    ub = c(NA, NA, NA, NA, NA, NA),
    lb = c(NA, 0, 0, 0, 0, 0),
    type = "real",
    vars = family_args$vars,
    loop = family_args$loop,
    log_lik = log_lik_cswald_simple,
    posterior_predict = posterior_predict_cswald_simple
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_helpers <- read_lines2(paste0(sc_path, "/cswald_helper_functions.stan"))
  stan_functions <- read_lines2(paste0(sc_path, "/cswald_simple_functions.stan"))

  stanvars <- brms::stanvar(scode = stan_helpers, block = "functions") +
    brms::stanvar(scode = stan_functions, block = "functions")

  nlist(formula, data, stanvars)
}

posterior_predict_cswald_simple <- function(i, prep, ...) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sndt <- brms::get_dpar(prep, "sndt", i = i)

  # convert single-boundary bound to total separation for the full DDM generator
  out <- .rcswald(
    n = length(drift),
    drift = drift,
    bound = bound * 2,
    ndt = ndt,
    zr = 0.5,
    s = s,
    sndt = sndt
  )

  dots <- list(...)
  if (!is.null(dots$negative_rt) && dots$negative_rt) {
    out$rt * ifelse(out$response == 1, 1, -1)
  } else {
    out$rt
  }
}

log_lik_cswald_simple <- function(i, prep) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sndt <- brms::get_dpar(prep, "sndt", i = i)

  rt <- rep(prep$data$Y[i], length(drift))
  response <- rep(prep$data$dec[i], length(drift))

  .dcswald(rt, response, drift, bound, ndt,
    zr = 0.5, s = s, sndt = sndt, version = "simple", log = TRUE
  )
}

#' @export
configure_model.cswald_crisk <- function(model, data, formula) {
  links <- model$links
  formula <- bmf2bf(model, formula)

  family_args <- cswald_family_args(model)

  formula$family <- brms::custom_family(
    "cswald_crisk",
    dpars = c("mu", "drift", "bound", "ndt", "zr", "s", "sndt"),
    links = c(
      "identity", links$drift, links$bound, links$ndt, links$zr, links$s,
      links$sndt
    ),
    ub = c(NA, NA, NA, NA, 1, NA, NA),
    lb = c(NA, NA, 0, 0, 0, 0, 0),
    type = "real",
    vars = family_args$vars,
    loop = family_args$loop,
    log_lik = log_lik_cswald_crisk,
    posterior_predict = posterior_predict_cswald_crisk
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_helpers <- read_lines2(paste0(sc_path, "/cswald_helper_functions.stan"))
  stan_functions <- read_lines2(paste0(sc_path, "/cswald_crisk_functions.stan"))

  stanvars <- brms::stanvar(scode = stan_helpers, block = "functions") +
    brms::stanvar(scode = stan_functions, block = "functions")

  nlist(formula, data, stanvars)
}

log_lik_cswald_crisk <- function(i, prep) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sndt <- brms::get_dpar(prep, "sndt", i = i)

  rt <- rep(prep$data$Y[i], length(drift))
  response <- rep(prep$data$dec[i], length(drift))

  .dcswald(rt, response, drift, bound, ndt,
    zr = zr, s = s, sndt = sndt, version = "crisk", log = TRUE
  )
}

posterior_predict_cswald_crisk <- function(i, prep, ...) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sndt <- brms::get_dpar(prep, "sndt", i = i)

  # rtdists shares one non-decision-time draw between the accumulators, the
  # crisk likelihood does not; see the sndt section of ?cswald
  out <- .rcswald(
    n = length(drift),
    drift = drift,
    bound = bound,
    ndt = ndt,
    zr = zr,
    s = s,
    sndt = sndt
  )

  dots <- list(...)
  if (!is.null(dots$negative_rt) && dots$negative_rt) {
    out$rt * ifelse(out$response == 1, 1, -1)
  } else {
    out$rt
  }
}
