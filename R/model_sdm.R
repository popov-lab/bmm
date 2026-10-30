############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdm <- function(resp_error = NULL, links = NULL, version = "simple",
                       activation = "vonmises", scaling = "density", call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(activation, scaling),
      domain = "Working memory",
      task = "Continuous reproduction",
      name = "Signal Discrimination Model (SDM)",
      citation = glue(
        "Oberauer, K. (2023). Measurement models for visual working \\
        memory\u2014A factorial model comparison. Psychological Review, \\
        130(3), 841-852. https://doi.org/10.1037/rev0000328"
      ),
      version = version,
      requirements = glue(
        "- The response variable should be in radians and represent the angular \\
        error relative to the target"
      ),
      parameters = list(
        mu = glue("Location parameter of the SDM distribution (in radians; \\
                  by default fixed internally to 0)"),
        c = .sdm_c_description[[scaling]],
        kappa = if (activation == "laplace") {
          "Precision parameter: the rate at which the Laplace activation falls with distance from the target"
        } else {
          "Precision parameter of the SDM distribution"
        }
      ),
      links = list(
        mu = "tan_half",
        c = "log",
        kappa = "log"
      ),
      fixed_parameters = list(mu = 0),
      default_priors = list(
        mu = list(main = "normal(0, 0.5)", effects = "normal(0, 0.25)", sd = "exponential(4)"),
        kappa = list(main = "student_t(5, 1.75, 0.75)", effects = "normal(0, 1)", sd = "exponential(1)"),
        # the same for every scaling: with kappa's default prior, a prior
        # predictive check of the circular SD of the response error matched the
        # default cell best at this location for the von Mises peak and range
        # scalings too (#558, see ?sdm)
        c = list(main = "student_t(5, 2, 0.75)", effects = "normal(0, 1)", sd = "exponential(1)")
      ),
      init_ranges = list(
        mu = c(-0.5,0.5),
        kappa = c(2.5,3.5),
        c = c(4,6)
      )
    ),
    class = c("bmmodel", "circular", "sdm", paste0("sdm_", version)),
    call = call
  )
  out <- set_links(out, links)
  out
}

.sdm_c_description <- c(
  density = "Memory strength parameter of the SDM distribution",
  peak = "Memory strength: the activation of the target (scaling = \"peak\")",
  range = "Memory strength: the activation of the target minus that of the location opposite it (scaling = \"range\")"
)

#' @export
print_model_details.sdm <- function(model, ...) {
  meaning <- switch(model$other_vars$scaling,
    peak = "the activation of the target",
    range = "the target-opposite contrast in activation"
  )
  if (!is.null(meaning)) {
    cat(glue("Scaling:    memory strength `c` is {meaning}."), "\n")
  }
  invisible(NULL)
}

# configure_model.sdm declares the family links itself, and the log link of `c`
# is written into the Stan chunk as exp(c), so none of the three can be set
#' @exportS3Method
settable_links.sdm <- function(model) {
  character(0)
}

# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_sdm_simple(NA)$info

#' @title `r .model_sdm()$name`
#' @name sdm
#' @details see [the online article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for a detailed description of the model
#'   and how to use it. `r model_docs(.model_sdm())`
#' @param resp_error The name of the variable in the dataset containing the
#'   response error. The response error should code the response relative to the
#'   to-be-recalled target in radians. You can transform the response error in
#'   degrees to radians using the `deg2rad` function.
#' @param version Character. The version of the model to use. Currently only
#'   "simple" is supported.
#' @param activation Character. The activation kernel of the SDM,
#'   `"vonmises"` (the default, as in Oberauer, 2023) or `"laplace"`. See the
#'   section "Activation functions".
#' @param scaling Character. What the memory strength `c` measures:
#'   `"density"`, `"peak"` or `"range"`. `NULL` (the default) means `"density"`
#'   for the von Mises and `"peak"` for the Laplace activation. See the
#'   sections "Scaling of the activation" and "What `c` means".
#' @param ... used internally for testing, ignore it
#' @section Activation functions: Each response \eqn{x} is activated by
#'   \eqn{c^* s(x - \mu)}, where the activation function \eqn{s} is 1 at the
#'   target and \eqn{c^*} is the activation of the target; the response density
#'   is \eqn{\exp(c^* s(x - \mu)) / Z}, with \eqn{Z} its integral over the
#'   circle.
#'
#'   - `activation = "vonmises"`: \eqn{s(x) = \exp(\kappa (\cos x - 1))}, the
#'   von Mises activation kernel of the signal discrimination model of Oberauer
#'   (2023).
#'   - `activation = "laplace"`: \eqn{s(x) = \exp(-\kappa |x|)}, with
#'   \eqn{|x|} the angular distance and \eqn{\kappa} the rate of the Laplace
#'   activation in Oberauer (2023, eq. 2). This is the SDM with the Laplace
#'   activation kernel of Oberauer's model space (Table 2): the same Gumbel
#'   noise and decision rule, with \eqn{\kappa} estimated freely. The density
#'   has a kink at \eqn{\mu}, so estimating `mu` puts a kink at every
#'   response, which gradient-based sampling handles poorly; `bmm()` warns
#'   when the formula estimates `mu`.
#'
#'   With the Laplace activation \eqn{\kappa} is a rate per radian, and the
#'   default priors, set for the von Mises activation, put 21% of the prior
#'   predictive response distributions below a circular SD of 0.05 radians with
#'   `scaling = "peak"` and 63% with `scaling = "density"`, against 7% for the
#'   default cell. Set priors on `kappa` and `c` where that matters.
#' @section Scaling of the activation: `scaling` decides how `c` translates
#'   into the activation of the target \eqn{c^*}; the table in "What `c`
#'   means" gives the formulas. `"density"` divides `c` by the integral of
#'   \eqn{s} over the circle, so `c` weighs an activation of unit area. For the
#'   von Mises activation that integral is approximated by
#'   \eqn{\sqrt{2\pi/\kappa}}, the form `bmm` has always used (the "sqrtexp"
#'   parametrization; [c_sqrtexp2bessel()] converts to Oberauer's exact
#'   Bessel form), and fits from before `scaling` existed are this cell. For the
#'   Laplace activation the integral is exact, \eqn{2(1 - e^{-\kappa \pi}) /
#'   \kappa}; Oberauer (2023, eq. 2) multiplies by \eqn{\kappa / 2}, which
#'   divides by the integral over the real line, \eqn{2 / \kappa}. The two
#'   differ by the factor \eqn{1 - e^{-\kappa \pi}} (0.998 at
#'   \eqn{\kappa = 2}). `"peak"` takes `c` as the activation of the target, and
#'   `"range"` as the difference between the activation of the target and that
#'   of the location opposite it.
#'
#'   The default prior on `c` is `student_t(5, 2, 0.75)` on the log scale for
#'   every scaling. A prior predictive check of the circular SD of the
#'   response error, with `kappa`'s default prior, matched that of the default
#'   cell best at this location for the peak and range scalings of the von
#'   Mises activation too (Kolmogorov distance 0.06).
#' @section What `c` means:
#'
#'   | `scaling` | \eqn{c^*}, von Mises | \eqn{c^*}, Laplace | `c` is |
#'   |---|---|---|---|
#'   | `"density"` | \eqn{c \sqrt{\kappa / 2\pi}} | \eqn{c \kappa / (2 (1 - e^{-\kappa \pi}))} | the weight of an activation of unit area |
#'   | `"peak"` | \eqn{c} | \eqn{c} | the activation of the target |
#'   | `"range"` | \eqn{c / (1 - e^{-2\kappa})} | \eqn{c / (1 - e^{-\kappa \pi})} | the target-opposite contrast in activation |
#'
#'   For a given data set the three scalings give the same likelihood at
#'   values of `c` related by these formulas; what differs is which quantity
#'   the regression formula for `c` and its prior act on.
#' @return An object of class `bmmodel`
#' @export
#' @keywords bmmodel
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' # simulate data from the model
#' dat <- data.frame(y = rsdm(n = 1000, c = 4, kappa = 3))
#'
#' # specify formula
#' ff <- bmf(
#'   c ~ 1,
#'   kappa ~ 1
#' )
#'
#' # specify the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = sdm(resp_error = "y"),
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
sdm <- function(resp_error, version = "simple", activation = c("vonmises", "laplace"),
                scaling = NULL, ...) {
  call <- match.call()
  .check_constructor_dots(list(...), known = "links")
  stop_missing_args()
  version <- match.arg(version)
  activation <- match.arg(activation)
  scaling <- if (is.null(scaling)) .sdm_default_scaling(activation) else match.arg(scaling, .sdm_scalings)
  .model_sdm(
    resp_error = resp_error, version = version, activation = activation,
    scaling = scaling, call = call, ...
  )
}

#' @rdname sdm
#' @keywords deprecated
#' @export
sdmSimple <- function(resp_error, version = "simple", ...) {
  warning2("The function `sdmSimple()` is deprecated. Please use `sdm()` instead.")
  call <- match.call()
  .check_constructor_dots(list(...), known = "links")
  stop_missing_args()
  .model_sdm(
    resp_error = resp_error, version = version, activation = "vonmises",
    scaling = "density", call = call, ...
  )
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_formula.sdm <- function(model, data, formula) {
  warnif(
    model$other_vars$activation == "laplace" && !is_constant(formula$mu),
    "The formula estimates mu with the Laplace activation, whose density has a \\
    kink at mu: the likelihood has a kink at every response, which gradient-based \\
    sampling handles poorly. Consider keeping mu fixed at 0."
  )
  NextMethod("check_formula")
}

#' @export
check_data.sdm <- function(model, data, formula) {
  # data sorted by predictors is necessary for speedy computation of normalizing constant
  data <- order_data_query(model, data, formula)
  attr(data, "sdm_run_metadata") <- sdm_run_metadata(data, formula, model)
  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.sdm <- function(model, data, formula) {
  activation <- model$other_vars$activation
  scaling <- model$other_vars$scaling
  # the default cell keeps the family, chunks and stanvars of the fits before
  # #558, so their Stan code stays byte-identical (test-model_sdm-stancode.R);
  # the other cells pass the cell to Stan as two integer codes
  family_name <- if (.sdm_is_default_cell(activation, scaling)) "sdm_simple" else "sdm_cont"
  # note - c has a log link, but I've coded it manually for computational efficiency
  sdm_family <- brms::custom_family(
    name = family_name,
    dpars = c("mu", "c", "kappa"),
    links = c("tan_half", "identity", "log"),
    lb = c(NA, NA, NA),
    ub = c(NA, NA, NA),
    type = "real", loop = FALSE,
    vars = if (family_name == "sdm_cont") c("sdm_activation", "sdm_scaling"),
    log_lik = log_lik_sdm_simple,
    posterior_predict = posterior_predict_sdm_simple,
    posterior_epred = posterior_epred_sdm_simple
  )

  # prepare initial stanvars to pass to brms, model formula and priors
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdm_simple_functions.stan"))
  if (family_name == "sdm_cont") {
    stan_funs <- paste(
      read_lines2(paste0(sc_path, "/sdm_kernel_funs.stan")),
      read_lines2(paste0(sc_path, "/sdm_cont_funs.stan")),
      sep = "\n"
    )
  }
  stan_tdata <- read_lines2(paste0(sc_path, "/sdm_simple_tdata.stan"))
  likelihood_file <- if (brms_slices_likelihood()) {
    paste0(family_name, "_likelihood_threaded.stan")
  } else {
    paste0(family_name, "_likelihood.stan")
  }
  stan_likelihood <- read_lines2(paste0(sc_path, "/", likelihood_file))
  run_metadata <- attr(data, "sdm_run_metadata")
  if (is.null(run_metadata)) {
    # Guard direct configure_model.sdm() calls that bypass check_data.sdm().
    run_metadata <- sdm_run_metadata(data, formula, model)
  }
  # 16 nodes per panel and the 32-node wide rule keep log Z within 4e-12
  # where it is below 1e4 and, above it, within 5e-10 to log c = 12 and
  # kappa = 3000 and within 3 ulps beyond; where sdm_simple_ldenom()
  # switches between them they differ by 4e-14 (#615). The Laplace rule
  # reuses the 16 nodes (#558). test-model_sdm-stan.R holds the bars
  gl <- .gl_rule(16L)
  wide <- .sdm_wide_rule()
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(x = gl$nodes, name = "sdm_gl_x") +
    brms::stanvar(x = gl$weights, name = "sdm_gl_w") +
    brms::stanvar(x = wide$cosm1, name = "sdm_wide_cosm1") +
    brms::stanvar(x = wide$weights, name = "sdm_wide_w") +
    brms::stanvar(x = run_metadata$G_sdm_runs, name = "G_sdm_runs") +
    sdm_stanvar_int_array(run_metadata$sdm_run_start, "sdm_run_start", "G_sdm_runs") +
    sdm_stanvar_int_array(run_metadata$sdm_run_count, "sdm_run_count", "G_sdm_runs") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end")
  if (family_name == "sdm_cont") {
    stanvars <- stanvars +
      sdm_stanvar_int(match(activation, .sdm_activations), "sdm_activation") +
      sdm_stanvar_int(match(scaling, .sdm_scalings), "sdm_scaling")
  }

  # construct main brms formula from the bmm formula
  formula <- bmf2bf(model, formula)
  formula$family <- sdm_family

  nlist(formula, data, stanvars)
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

#' @export
postprocess_brm.sdm <- function(model, fit, ...) {
  # manually set link_c to "log" since I coded it manually
  fit$family$link_c <- "log"
  fit$formula$family$link_c <- "log"
  fit
}

#' @export
revert_postprocess_brm.sdm <- function(model, fit, ...) {
  fit$family$link_c <- "identity"
  fit$formula$family$link_c <- "identity"
  fit
}

# serve both families: the codes of the cell are in the Stan data of an
# sdm_cont fit, and an sdm_simple fit is the von Mises density cell
log_lik_sdm_simple <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  c <- brms::get_dpar(prep, "c", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  y <- prep$data$Y[i]
  cell <- .sdm_prep_cell(prep)
  dsdm(y, mu, c, kappa, log = TRUE, activation = cell$activation, scaling = cell$scaling)
}

.sdm_prep_cell <- function(prep) {
  if (is.null(prep$data$sdm_activation)) {
    return(list(activation = "vonmises", scaling = "density"))
  }
  list(
    activation = .sdm_activations[prep$data$sdm_activation],
    scaling = .sdm_scalings[prep$data$sdm_scaling]
  )
}

posterior_epred_sdm_simple <- posterior_epred_undefined("sdm")

posterior_predict_sdm_simple <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  c <- brms::get_dpar(prep, "c", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  cell <- .sdm_prep_cell(prep)
  rsdm(length(mu), mu, c, kappa, activation = cell$activation, scaling = cell$scaling)
}

sdm_run_metadata <- function(data, formula, model) {
  predictors <- data_predictor_vars(data, formula)
  # brms excludes rows with missing values in any model variable before
  # fitting, so run boundaries must be computed on the rows brms will keep
  model_vars <- unique(c(unlist(model$resp_vars), predictors))
  model_vars <- model_vars[model_vars %in% colnames(data)]
  if (length(model_vars) > 0L) {
    data <- data[stats::complete.cases(data[model_vars]), , drop = FALSE]
  }
  if (length(predictors) == 0L) {
    return(list(
      G_sdm_runs = 1L,
      sdm_run_start = 1L,
      sdm_run_count = nrow(data)
    ))
  }

  run_id <- interaction(data[predictors], drop = TRUE, lex.order = TRUE)
  run_start <- c(1L, which(run_id[-1] != run_id[-length(run_id)]) + 1L)
  run_count <- diff(c(run_start, nrow(data) + 1L))
  list(
    G_sdm_runs = length(run_start),
    sdm_run_start = as.integer(run_start),
    sdm_run_count = as.integer(run_count)
  )
}

# The midpoint rule on [0, pi] after mapping t to t - 0.55 sin(t) three times.
# The map fixes 0 and pi and keeps the integrand smooth and periodic, so the
# rule keeps converging geometrically, and it crowds the nodes toward the
# spike at u = 0, where they are 11 times denser than without it. cos(u) - 1
# is formed as -2 sin(u / 2)^2, which does not cancel near the spike.
.sdm_wide_rule <- function() {
  n <- 32L
  u <- (seq_len(n) - 0.5) * pi / n
  weights <- rep(pi / n, n)
  for (i in 1:3) {
    weights <- weights * (1 - 0.55 * cos(u))
    u <- u - 0.55 * sin(u)
  }
  list(cosm1 = -2 * sin(u / 2)^2, weights = weights)
}

sdm_stanvar_int <- function(x, name) {
  out <- brms::stanvar(x = as.integer(x), name = name, scode = paste0("int ", name, ";"))
  out[[name]]$pll_args <- paste("data int", name)
  out
}

sdm_stanvar_int_array <- function(x, name, size) {
  x <- array(as.integer(x), dim = length(x))
  out <- brms::stanvar(x = x, name = name)
  out[[name]]$scode <- paste0("array[", size, "] int ", name, ";")
  out[[name]]$pll_args <- paste("data array[] int", name)
  out
}
