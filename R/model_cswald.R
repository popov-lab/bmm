############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.cswald_version_table <- list(
  simple = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary (distance from starting point to correct boundary, half the boundary separation `bound` of version \"crisk\", ddm() and ezdm())",
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
      drift = list(main = "normal(0,1)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      s = list(main = "normal(0,0.3)", effects = "normal(0,0.2)", sd = "exponential(2)"),
      sndt = list(main = "normal(-2.5,1)", effects = "normal(0,0.3)", sd = "exponential(2)")
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
      zr = "relative starting point (proportion of the boundary separation, measured from the lower boundary; 0.5 = unbiased)",
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
      drift = list(main = "normal(0,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      zr = list(main = "normal(0,0.3)", effects = "normal(0,0.2)", sd = "exponential(2)"),
      s = list(main = "normal(0,0.5)", effects = "normal(0,0.2)", sd = "exponential(2)"),
      sndt = list(main = "normal(-2.5,1)", effects = "normal(0,0.3)", sd = "exponential(2)")
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
      domain = "Decision making",
      task = "Choices and response times",
      name = "Censored-Shifted Wald Model",
      citation = glue(
        "Miller, R., Scherbaum, S., Heck, D. W., Goschke, T., & Enge, S. (2018). \\
        On the relation between the (censored) shifted Wald and the Wiener \\
        distribution as measurement models for choice response times. Applied \\
        Psychological Measurement, 42(2), 116-135. \\
        https://doi.org/10.1177/0146621617710465"
      ),
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
      init_ranges = .cswald_version_table[[version]][["init_ranges"]],
      functions = list(density = "dcswald", generator = "rcswald")
    ),
    class = c("bmmodel", "cswald", paste0("cswald_", version)),
    call = call
  )

  resolve_fixed_links(set_links(out, links))
}

#' @title `r .model_cswald()$name`
#' @name cswald
#' @details `r model_docs(.model_cswald())`
#'
#' #### The prior for `ndt` and the fastest response
#' `ndt` has to stay below every response the model fits: all responses in
#' version `"crisk"`, the correct ones in `"simple"`, which censors errors.
#' `bmm()` and `update()` warn when more than half of the prior for `ndt` lies
#' at or above the fastest of these responses. Sampling can then fail to
#' start, and the data pull `ndt` down against the prior. With the default
#' prior, `normal(-2, 0.3)` on the log scale, this happens when the fastest
#' response is faster than 0.135 s. Set a prior with a lower location, as the
#' warning suggests, or filter out anticipations. The check reads `normal()`
#' and `student_t()` priors on the intercept of `ndt`, or on each cell of a
#' formula without an intercept such as `ndt ~ 0 + cond`; other priors and
#' designs are not checked, and neither are truncated priors or the links
#' `"inverse"` and `"logm1"`. In a formula with other fixed effects, such as
#' `ndt ~ 1 + cond`, the prior on the intercept is the prior on the centred
#' intercept (the grand mean), and the check compares it with the fastest
#' response overall. `update()` with `init` or `inits` skips the check, while
#' `bmm()` with `init` still runs it.
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
#'   "log" for `bound`, `ndt`, `s` and `sndt`, and "logit" for `zr`. `drift`
#'   has a "log" link in "simple", where it is positive, and an "identity" link
#'   in "crisk", where it can be negative, as in [ddm()] and [ezdm()]. For
#'   positive parameters, "softplus" is available as an alternative to "log"
#'   that grows linearly for large values and avoids the numerical blow-up of
#'   `exp()`. A name that is not a
#'   parameter of the model is an error, and a link that allows values the
#'   default link excludes (e.g. "identity" for a positive parameter) is a
#'   warning.
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
#'       the diffusion model parameterization. `zr` is fixed at 0 on its
#'       logit scale (an unbiased `zr = 0.5`) unless you give it a formula.
#'   }
#'   For more details, see Miller et al. (2018).
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
#'   natural scale, and `sndt` uses the log link only once it is estimated. A
#'   negative fixed value is an error. `update()` carries a non-zero `sndt`
#'   forward: a formula that omits it leaves the stored value in place, so
#'   pass `sndt = 0` explicitly to return to the model without variability.
#'
#'   Without `sndt`, the fastest responses cap the `ndt` estimate (the
#'   likelihood requires `ndt < min(rt)`), which biases `ndt`, `drift`, and
#'   `bound` when the true non-decision time varies (Miller et al., 2018,
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
#'   through [rtdists::rdiffusion()], so for `"crisk"` `pp_check()` compares
#'   the data against the diffusion the race approximates rather than against
#'   the fitted likelihood -- at drift 1, bound 1.6, `zr = 0.5` the two differ
#'   by 7 percentage points in P(upper) (0.832 vs 0.899), at any `sndt`. That
#'   each accumulator draws its own non-decision time adds a smaller
#'   difference on top.
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
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
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
  .check_constructor_dots(list(...))
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

  data <- check_rt_data(
    data, rt_var, response_var,
    fast_rt_hint = paste(
      "Such fast responses are often anticipations and pull the estimate of",
      "`ndt` down; consider filtering them."
    )
  )
  data[[response_var]] <- code_dec_response(data[[response_var]], response_var)

  if (model$version == "simple") {
    error_rate <- mean(data[[response_var]] == 0)
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
# CHECK_FORMULA S3 METHODS                                               ####
############################################################################# !

# a negative sndt makes every Stan evaluation -inf, which a real backend
# reports as an initial-value storm that never names the parameter
#' @export
check_formula.cswald <- function(model, data, formula) {
  stopif(
    isTRUE(model$fixed_parameters$sndt < 0),
    "The non-decision time variability 'sndt' cannot be fixed to a negative \\
    value. Use sndt = 0 for a model without variability."
  )

  NextMethod("check_formula")
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
# exist inside partial_log_lik, so the decisions are sliced only where brms
# really threads (see brms_slices_likelihood)
cswald_family_args <- function(model) {
  if (!isTRUE(model$fixed_parameters[["sndt"]] == 0)) {
    return(list(loop = TRUE, vars = "dec[n]"))
  }
  list(loop = FALSE, vars = if (brms_slices_likelihood()) "dec[start:end]" else "dec")
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
    posterior_predict = posterior_predict_cswald_simple,
    posterior_epred = posterior_epred_cswald
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
  sndt <- .cswald_sndt(prep, i)

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

# Named after the family, "cswald", rather than the version: restructure()
# finds the function of a fit saved without one by that name. Like
# posterior_predict_cswald_simple(), it reads bound as the distance from the
# start point to either boundary.
posterior_epred_cswald <- function(prep) {
  .epred_matrix(
    with(prep$dpars, .diffusion_mean_rt(drift, bound * 2, ndt, zr = 0.5, s)) +
      .cswald_mean_sndt(prep),
    prep
  )
}

# ndt is the minimum of Uniform(ndt, ndt + sndt), so the mean RT gains sndt / 2.
# Fits saved before sndt existed carry no such dpar.
.cswald_mean_sndt <- function(prep) {
  (prep$dpars$sndt %||% 0) / 2
}

.cswald_sndt <- function(prep, i) {
  if (is.null(prep$dpars$sndt)) 0 else brms::get_dpar(prep, "sndt", i = i)
}

log_lik_cswald_simple <- function(i, prep) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sndt <- .cswald_sndt(prep, i)

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
    posterior_predict = posterior_predict_cswald_crisk,
    posterior_epred = posterior_epred_cswald_crisk
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_helpers <- read_lines2(paste0(sc_path, "/cswald_helper_functions.stan"))
  stan_functions <- read_lines2(paste0(sc_path, "/cswald_crisk_functions.stan"))

  stanvars <- brms::stanvar(scode = stan_helpers, block = "functions") +
    brms::stanvar(scode = stan_functions, block = "functions")

  nlist(formula, data, stanvars)
}

############################################################################# !
# SIMULATE_DATA S3 METHODS                                               ####
############################################################################# !

# the simple version's bound is the distance from a centred start point to
# either boundary, as in pp_simulate.cswald_simple()
#' @export
simulate_data.cswald <- function(model, pars, design = NULL, scale = c("link", "native"), ...) {
  design <- .design_or_n(model, design, character(0), ...)
  sim <- .simulate_pars(model, pars, design, scale)
  p <- sim$native
  out <- if (model$version == "simple") {
    rcswald(nrow(design), p$drift, 2 * p$bound, p$ndt, zr = 0.5, s = p$s)
  } else {
    rcswald(nrow(design), p$drift, p$bound, p$ndt, zr = p$zr, s = p$s)
  }
  .simulate_result(design, .rt_response(model, out), sim)
}

log_lik_cswald_crisk <- function(i, prep) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sndt <- .cswald_sndt(prep, i)

  rt <- rep(prep$data$Y[i], length(drift))
  response <- rep(prep$data$dec[i], length(drift))

  .dcswald(rt, response, drift, bound, ndt,
    zr = zr, s = s, sndt = sndt, version = "crisk", log = TRUE
  )
}

posterior_epred_cswald_crisk <- function(prep) {
  .epred_matrix(
    with(prep$dpars, .diffusion_mean_rt(drift, bound, ndt, zr, s)) +
      .cswald_mean_sndt(prep),
    prep
  )
}

posterior_predict_cswald_crisk <- function(i, prep, ...) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)
  sndt <- .cswald_sndt(prep, i)

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

############################################################################# !
# CREATE_INITFUN S3 method                                               ####
############################################################################# !

# The likelihood is -Inf wherever ndt >= min(rt), so the start for ndt has to
# come from the data; the other ranges do not depend on it. rt = 0 is left out:
# a censored response there never binds, and a response of 1 there cannot be
# fitted with any positive ndt
#' @export
create_initfun.cswald <- function(model, data, formula, prior = NULL, ...) {
  rt <- data[[model$resp_vars$rt]]
  model$init_ranges$ndt <- c(0.4, 0.6) * min(rt[rt > 0])
  # the first step of bmm() and update() that sees both the data and the prior
  # brms will use
  .warn_ndt_prior_conflict(model, data, formula, prior)
  NextMethod()
}

# A prior with most of its mass where the likelihood is -Inf fights the data
# (#548). ndt has to stay below every response it fits: all of them in crisk,
# the correct ones in simple, whose errors are censored. Only the intercept, or
# each cell of a design without one (ndt ~ 0 + cond), has a prior that maps to
# a set of responses, and only normal() and student_t() are read
.warn_ndt_prior_conflict <- function(model, data, formula, prior) {
  bterms <- brms::brmsterms(formula)
  ndt_terms <- bterms$dpars$ndt
  probe <- suppressWarnings(link_transform(c(0.1, 1), model$links$ndt))
  # on a decreasing link, such as inverse, the mass above the fastest response
  # is the mass below its transform; logm1 is undefined at the probe points
  increasing <- all(is.finite(probe)) && probe[2] > probe[1]
  if (NROW(prior) == 0 || is.null(ndt_terms$fe) || !increasing) {
    return(invisible())
  }
  rt <- data[[model$resp_vars$rt]]
  # brms drops every row with a missing value in any variable of the model
  complete <- stats::complete.cases(data[intersect(all.vars(bterms$allvars), names(data))])
  fits <- rt > 0 & complete
  if (model$version == "simple") {
    fits <- fits & data[[model$resp_vars$response]] == 1
  }

  if (has_intercept(ndt_terms$fe)) {
    prior_class <- "Intercept"
    cells <- list(Intercept = fits)
  } else {
    prior_class <- "b"
    # brms adds the Intercept pseudo-column of ndt ~ 0 + Intercept to its own
    # frame only; a diagnostic must not stop a fit it cannot read the design of
    X <- tryCatch(
      get_model_matrix(ndt_terms$fe, transform(data, Intercept = 1)),
      error = function(e) NULL
    )
    if (is.null(X) || !all(X %in% c(0, 1)) || any(rowSums(X) != 1)) {
      return(invisible())
    }
    # model.matrix() drops the rows with a missing predictor, as brms does
    # later, so the cells are placed on the rows it kept
    kept <- match(rownames(X), rownames(data))
    cells <- lapply(seq_cols(X), function(j) {
      in_cell <- logical(nrow(data))
      in_cell[kept] <- X[, j] == 1
      fits & in_cell
    })
    names(cells) <- colnames(X)
  }
  rows <- prior[prior$dpar == "ndt" & prior$class == prior_class & !nzchar(prior$nlpar) &
                  !nzchar(prior$group) & !nzchar(prior$resp), ]
  # a truncated prior has less mass above the fastest response than it shows
  if (any(nzchar(c(rows$lb, rows$ub), keepNA = TRUE), na.rm = TRUE)) {
    return(invisible())
  }

  conflicts <- lapply(names(cells), function(cell) {
    if (!any(cells[[cell]])) {
      return(NULL)
    }
    # as in brms, a row on the coefficient beats the row on its class, and a
    # row without a prior string inherits
    applied <- c(rows$prior[!nzchar(rows$coef)], rows$prior[rows$coef == cell])
    applied <- applied[nzchar(applied)]
    dist <- if (length(applied)) .parse_location_scale_prior(applied[length(applied)])
    if (is.null(dist)) {
      return(NULL)
    }
    fastest <- min(rt[cells[[cell]]])
    limit <- link_transform(fastest, model$links$ndt)
    n <- length(dist$args)
    z <- (limit - dist$args[n - 1]) / dist$args[n]
    mass <- if (dist$family == "normal") {
      stats::pnorm(z, lower.tail = FALSE)
    } else {
      stats::pt(z, dist$args[1], lower.tail = FALSE)
    }
    if (!isTRUE(mass > 0.5)) {
      return(NULL)
    }
    # a mass just above one half would read as "50% ... most of its mass"
    share <- if (round(100 * mass, 1) == 50) {
      "just over 50%"
    } else if (round(100 * mass) == 50) {
      paste0(format(round(100 * mass, 1), nsmall = 1), "%")
    } else {
      paste0(round(100 * mass), "%")
    }
    # two scales below the fastest response leaves about 2.5% of the mass above
    dist$args[n - 1] <- floor((limit - 2 * dist$args[n]) * 100 + 1e-8) / 100
    list(
      limit = limit,
      own_row = any(nzchar(rows$prior[rows$coef == cell])),
      location = dist$args[n - 1],
      line = glue("{cell}: {applied[length(applied)]} puts {share} of \\
                  its mass at or above {signif(fastest, 3)} s"),
      suggestion = paste0(dist$family, "(", paste(dist$args, collapse = ", "), ")"),
      cell = cell
    )
  })
  conflicts <- Filter(Negate(is.null), conflicts)
  if (length(conflicts) == 0) {
    return(invisible())
  }

  kind <- if (model$version == "simple") "correct response" else "response"
  where <- if (prior_class == "b") " in its cell" else ""
  # a long list would push the suggestion past the 1000 characters R prints
  lines <- vapply(conflicts, `[[`, character(1), "line")
  listed <- max(1, min(5, sum(cumsum(nchar(lines)) <= 400)))
  lines <- lines[seq_len(listed)]
  if (length(conflicts) > listed) {
    lines <- c(lines, glue("... and {length(conflicts) - listed} more"))
  }
  lines <- paste0("  ", lines, collapse = "\n")
  # one row on the class serves the conflicting cells only when none has a
  # row of its own, which beats the class row; it is built from the lowest
  # fastest response. Otherwise the suggestion names one cell, one with its
  # own row first, so that applying it makes progress
  own_rows <- vapply(conflicts, `[[`, logical(1), "own_row")
  several <- prior_class == "b" && length(conflicts) > 1 && !any(own_rows)
  lowest <- which.min(vapply(conflicts, `[[`, numeric(1), "limit"))
  first <- conflicts[[if (several) lowest else if (any(own_rows)) which(own_rows)[1] else 1]]
  coef <- if (prior_class == "b" && !several) glue(', coef = "{first$cell}"') else ""
  # a location that maps to a time at or below zero is no prior for a time
  advice <- if (link_transform(first$location, model$links$ndt, inverse = TRUE) <= 0) {
    paste(
      "The fastest response is less than two scales above zero, so the prior",
      "needs a smaller scale as well as a lower location. You can also filter",
      "out anticipations, responses faster than the task allows."
    )
  } else {
    glue(
      "Give `ndt` a prior with a lower location, for example\n",
      "  set_prior(\"{first$suggestion}\", class = \"{prior_class}\"{coef}, dpar = \"ndt\")\n",
      "or filter out anticipations, responses faster than the task allows.",
      .trim = FALSE
    )
  }
  warning2(
    "The prior for `ndt` puts most of its mass at or above the fastest {kind} \\
    time{where}, but `ndt` must stay below it:
    {lines}
    Sampling can fail to start, and the data pull `ndt` below that response \\
    against the prior. {advice}"
  )
}

# The family and numeric arguments of a normal() or student_t() prior string,
# NULL for anything else, including arguments that are not plain numbers
.parse_location_scale_prior <- function(prior) {
  parts <- regmatches(prior, regexec("^\\s*(normal|student_t)\\s*\\(([^()]*)\\)\\s*$", prior))[[1]]
  if (length(parts) == 0) {
    return(NULL)
  }
  args <- suppressWarnings(as.numeric(strsplit(parts[3], ",")[[1]]))
  if (anyNA(args) || length(args) != c(normal = 2, student_t = 3)[[parts[2]]]) {
    return(NULL)
  }
  if (parts[2] == "student_t" && args[1] <= 0) {
    return(NULL)
  }
  list(family = parts[2], args = args)
}

############################################################################# !
# PP_CHECK OBSERVABLES                                                    ####
############################################################################# !

#' @export
pp_observables.cswald <- function(model) {
  .pp_spec_rt_response()
}

#' @export
pp_simulate.cswald_simple <- function(model, prep) {
  # .rcswald() is two-boundary; cswald's bound is the single-boundary distance
  .pp_simulate_joint(prep, .rcswald, c("drift", "ndt", "s"),
                     bound = .pp_dpar_vector(prep, "bound") * 2, zr = 0.5,
                     sndt = .pp_cswald_sndt(prep))
}

#' @export
pp_simulate.cswald_crisk <- function(model, prep) {
  .pp_simulate_joint(prep, .rcswald, c("drift", "bound", "ndt", "zr", "s"),
                     sndt = .pp_cswald_sndt(prep))
}


# .rcswald() defaults sndt = 0, so simply omitting it from the dpar list would
# silently simulate a model the fit did not estimate. Read it when the fit has
# it, and fall back to 0 only for fits saved before sndt existed -- for those
# the absent parameter *is* zero.
.pp_cswald_sndt <- function(prep) {
  if (is.null(prep$dpars$sndt)) {
    return(0)
  }
  .pp_dpar_vector(prep, "sndt")
}
