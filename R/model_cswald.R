############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.cswald_version_table <- list(
  simple = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary (distance from starting point to correct boundary, half the boundary separation `bound` of version \"crisk\", ddm() and ezdm())",
      ndt = "non-decision time",
      s = "diffusion constant"
    ),
    links = list(
      drift = "log",
      bound = "log",
      ndt = "log",
      s = "log"
    ),
    fixed_parameters = list(
      mu = 0,
      s = 0
    ),
    priors = list(
      drift = list(main = "normal(0,1)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      s = list(main = "normal(0,0.3)", effects = "normal(0,0.2)", sd = "exponential(2)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      drift = c(1, 2),
      bound = c(1.5, 2),
      ndt = c(0.025, 0.05),
      s = c(0.95, 1.05)
    )
  ),
  crisk = list(
    parameters = list(
      drift = "drift rate",
      bound = "boundary separation (total distance between boundaries)",
      ndt = "non-decision time",
      zr = "relative starting point (proportion of the boundary separation, measured from the lower boundary; 0.5 = unbiased)",
      s = "diffusion constant"
    ),
    links = list(
      drift = "identity",
      bound = "log",
      ndt = "log",
      zr = "logit",
      s = "log"
    ),
    fixed_parameters = list(
      mu = 0,
      zr = 0,
      s = 0
    ),
    priors = list(
      drift = list(main = "normal(0,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
      bound = list(main = "normal(0,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      ndt = list(main = "normal(-2,0.3)", effects = "normal(0,0.3)", sd = "exponential(2)"),
      zr = list(main = "normal(0,0.3)", effects = "normal(0,0.2)", sd = "exponential(2)"),
      s = list(main = "normal(0,0.5)", effects = "normal(0,0.2)", sd = "exponential(2)")
    ),
    init_ranges = list(
      mu = c(-0.5, 0.5),
      drift = c(-0.5, 0.5),
      bound = c(1.5, 2),
      ndt = c(0.025, 0.05),
      zr = c(0.45, 0.55),
      s = c(0.95, 1.05)
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
      fixed_parameters = .cswald_version_table[[version]][["fixed_parameters"]],
      default_priors = .cswald_version_table[[version]][["priors"]],
      init_ranges = .cswald_version_table[[version]][["init_ranges"]]
    ),
    class = c("bmmodel", "cswald", paste0("cswald_", version)),
    call = call
  )

  out <- set_links(out, links)
  out
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
#' designs are not checked, and neither are truncated priors or a decreasing
#' link such as `"inverse"`. In a formula with other fixed effects, such as
#' `ndt ~ 1 + cond`, the prior on the intercept is the prior on the centred
#' intercept (the grand mean), and the check compares it with the fastest
#' response overall. `update()` with `init` or `inits` skips the check.
#' @param rt The name of the variable in the dataset containing the response
#'   times. Response times should be coded in seconds (not milliseconds).
#' @param response The name of the variable in the dataset containing the
#'   response/decision. Responses should be coded as 0 (lower boundary) or
#'   1 (upper boundary). Alternatively, character values "lower" and "upper"
#'   or logical values (FALSE/TRUE) are accepted and will be converted
#'   automatically.
#' @param links A named list of link functions for the model parameters.
#'   Available parameters depend on the version: "simple" has `drift`, `bound`,
#'   `ndt`, and `s`; "crisk" additionally has `zr`. Default links are "log" for
#'   `bound`, `ndt` and `s`, and "logit" for `zr`. `drift` has a "log" link in
#'   "simple", where it is positive, and an "identity" link in "crisk", where
#'   it can be negative, as in [ddm()] and [ezdm()]. For positive parameters,
#'   "softplus" is available as an alternative to "log" that grows linearly for
#'   large values and avoids the numerical blow-up of `exp()`. A name that is not a
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

# start/end only exist inside partial_log_lik, so the decisions are sliced only
# where brms really threads (see brms_slices_likelihood)
cswald_decision_var <- function() {
  if (brms_slices_likelihood()) "dec[start:end]" else "dec"
}

#' @export
configure_model.cswald_simple <- function(model, data, formula) {
  links <- model$links
  formula <- bmf2bf(model, formula)

  formula$family <- brms::custom_family(
    "cswald",
    dpars = c("mu", "drift", "bound", "ndt", "s"),
    links = c("identity", links$drift, links$bound, links$ndt, links$s),
    ub = c(NA, NA, NA, NA, NA),
    lb = c(NA, 0, 0, 0, 0),
    type = "real",
    vars = cswald_decision_var(),
    loop = FALSE,
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

  # convert single-boundary bound to total separation for the full DDM generator
  out <- .rcswald(
    n = length(drift),
    drift = drift,
    bound = bound * 2,
    ndt = ndt,
    zr = 0.5,
    s = s
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
  .epred_matrix(with(prep$dpars, .diffusion_mean_rt(drift, bound * 2, ndt, zr = 0.5, s)), prep)
}

log_lik_cswald_simple <- function(i, prep) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  rt <- rep(prep$data$Y[i], length(drift))
  response <- rep(prep$data$dec[i], length(drift))

  .dcswald(rt, response, drift, bound, ndt, zr = 0.5, s = s, version = "simple", log = TRUE)
}

#' @export
configure_model.cswald_crisk <- function(model, data, formula) {
  links <- model$links
  formula <- bmf2bf(model, formula)

  formula$family <- brms::custom_family(
    "cswald_crisk",
    dpars = c("mu", "drift", "bound", "ndt", "zr", "s"),
    links = c("identity", links$drift, links$bound, links$ndt, links$zr, links$s),
    ub = c(NA, NA, NA, NA, 1, NA),
    lb = c(NA, NA, 0, 0, 0, 0),
    type = "real",
    vars = cswald_decision_var(),
    loop = FALSE,
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

log_lik_cswald_crisk <- function(i, prep) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  rt <- rep(prep$data$Y[i], length(drift))
  response <- rep(prep$data$dec[i], length(drift))

  .dcswald(rt, response, drift, bound, ndt, zr = zr, s = s, version = "crisk", log = TRUE)
}

posterior_epred_cswald_crisk <- function(prep) {
  .epred_matrix(with(prep$dpars, .diffusion_mean_rt(drift, bound, ndt, zr, s)), prep)
}

posterior_predict_cswald_crisk <- function(i, prep, ...) {
  drift <- brms::get_dpar(prep, "drift", i = i)
  bound <- brms::get_dpar(prep, "bound", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  zr <- brms::get_dpar(prep, "zr", i = i)
  s <- brms::get_dpar(prep, "s", i = i)

  out <- .rcswald(
    n = length(drift),
    drift = drift,
    bound = bound,
    ndt = ndt,
    zr = zr,
    s = s
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
  ndt_terms <- brms::brmsterms(formula)$dpars$ndt
  probe <- suppressWarnings(link_transform(c(0.1, 1), model$links$ndt))
  # on a decreasing link, such as inverse, the mass above the fastest response
  # is the mass below its transform; logm1 is undefined at the probe points
  increasing <- all(is.finite(probe)) && probe[2] > probe[1]
  if (NROW(prior) == 0 || is.null(ndt_terms$fe) || !increasing) {
    return(invisible())
  }
  rt <- data[[model$resp_vars$rt]]
  fits <- rt > 0
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
    mass <- if (dist$family == "normal") stats::pnorm(z, lower.tail = FALSE) else
      stats::pt(z, dist$args[1], lower.tail = FALSE)
    if (!isTRUE(mass > 0.5)) {
      return(NULL)
    }
    # a mass just above one half would read as "50% ... most of its mass"
    percent <- if (round(100 * mass) == 50) {
      format(round(100 * mass, 1), nsmall = 1)
    } else {
      round(100 * mass)
    }
    # two scales below the fastest response leaves about 2.5% of the mass above
    dist$args[n - 1] <- floor((limit - 2 * dist$args[n]) * 100 + 1e-8) / 100
    list(
      location = dist$args[n - 1],
      line = glue("{cell}: {applied[length(applied)]} puts {percent}% of \\
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
  first <- conflicts[[1]]
  coef <- if (prior_class == "b") glue(', coef = "{first$cell}"') else ""
  # a location at or below zero is no prior for a time on its natural scale
  advice <- if (first$location <= 0 && identical(probe, c(0.1, 1))) {
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
                     bound = .pp_dpar_vector(prep, "bound") * 2, zr = 0.5)
}

#' @export
pp_simulate.cswald_crisk <- function(model, prep) {
  .pp_simulate_joint(prep, .rcswald, c("drift", "bound", "ndt", "zr", "s"))
}
