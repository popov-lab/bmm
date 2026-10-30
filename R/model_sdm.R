############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdm <- function(resp_error = NULL, links = NULL, version = "simple",
                       activation = "vonmises", scaling = "density", n_options = NULL,
                       units = "radians", options = NULL, n_alternatives = NULL, call = NULL, ...) {
  # the option arguments enter other_vars only when given, so the continuous
  # model keeps the other_vars of fits from before them
  option_vars <- if (!is.null(n_options)) {
    Filter(Negate(is.null), nlist(n_options, units, options, n_alternatives))
  }
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = c(nlist(activation, scaling), option_vars),
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
        kappa = list(main = .sdm_main_priors[[activation]]$kappa, effects = "normal(0, 1)", sd = "exponential(1)"),
        c = list(main = .sdm_main_priors[[activation]]$c[[scaling]], effects = "normal(0, 1)", sd = "exponential(1)")
      ),
      init_ranges = list(
        mu = c(-0.5,0.5),
        kappa = c(2.5,3.5),
        c = c(4,6)
      )
    ),
    class = c("bmmodel", "circular", if (!is.null(n_options)) "option_set", "sdm",
              paste0("sdm_", version)),
    call = call
  )
  out <- set_links(out, links)
  out
}

# kappa and c mean different things per cell, so the same prior would imply very
# different response errors; these make the implied circular SD of the response
# error comparable across cells, not equal (#558, see ?sdm)
.sdm_main_priors <- list(
  vonmises = list(
    kappa = "student_t(5, 1.75, 0.75)",
    c = c(density = "student_t(5, 2, 0.75)", peak = "student_t(5, 2, 0.75)", range = "student_t(5, 2, 0.75)")
  ),
  laplace = list(
    kappa = "student_t(5, -0.5, 0.75)",
    c = c(density = "student_t(5, 3.25, 0.75)", peak = "student_t(5, 2.5, 0.75)", range = "student_t(5, 2.25, 0.75)")
  )
)

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
  n_options <- model$other_vars$n_options
  if (!is.null(n_options)) {
    n_alternatives <- model$other_vars$n_alternatives
    offered <- if (!is.null(model$other_vars$options)) {
      "the options each trial offers"
    } else if (is.numeric(n_alternatives)) {
      glue("{n_alternatives} equally spaced options")
    } else {
      "all of them"
    }
    cat(glue("Options:    {n_options} positions on the wheel; each response is a choice among {offered}."), "\n")
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
#' @param n_options `NULL` (the default) for a continuous response, or the
#'   number of equally spaced response options on the wheel, an even number of
#'   at least 4: 360 for a colour wheel in 1-degree steps. Each response is
#'   then a choice among these options. See the section "Discrete response
#'   options".
#' @param options The names of the columns holding the options offered on
#'   each trial, one column per option, in the units of `resp_error` and
#'   relative to the target; rows that offer fewer options hold `NA` in the
#'   remaining columns. Needs `n_options`, the grid the options sit on. See the
#'   section "Forced-choice option sets".
#' @param n_alternatives With `options`, optionally the name of a column
#'   holding the number of options offered on each trial, which is checked
#'   against the option columns. Without `options`, a single number `m`: each
#'   trial offers `m` equally spaced options, the target among them, as in the
#'   m-alternative tasks of Robinson et al. (2023).
#' @param units With `n_options`, the units of `resp_error` and of the option
#'   columns: `"radians"` (the default), binned to the nearest option, or
#'   `"bins"`, whole option steps from the target (any integer, taken modulo
#'   `n_options`).
#' @param ... used internally for testing, ignore it
#' @section Discrete response options: With `n_options = K`, the response is a
#'   choice among \eqn{K} equally spaced options, and the probability of the
#'   option at offset \eqn{x} from the target is the softmax of the activation
#'   over the options, \eqn{\exp(c^* s(x - \mu)) / \sum_k \exp(c^* s(x_k -
#'   \mu))}: the SDM with Gumbel noise on a wheel of \eqn{K} options. Use it
#'   when participants chose among discrete options, as on a colour wheel of
#'   360 colours, and especially when the options are coarse.
#'
#'   The log-likelihood is reported on the density scale, \eqn{\log P -
#'   \log(2\pi / K)}, so `loo()` of a discrete fit can be compared with that
#'   of a continuous `sdm()` fit of the same responses. For the von Mises
#'   activation the two models agree wherever the activation spans a few
#'   steps: with \eqn{w = 1 / \sqrt{\max(c^*, 1) \kappa}}, the log-likelihood
#'   of a trial differs by less than 1e-10 once \eqn{w} is at least 2.25
#'   steps of \eqn{2\pi / K} (measured for \eqn{K} from 90 to 720; on coarser
#'   wheels the difference is larger at the same \eqn{w}), by about 1e-5 at
#'   one step, and by more below.
#'   `bmm()` warns after the fit when a posterior draw has \eqn{w} below one
#'   step; the discrete likelihood is then still exact, but a comparison with
#'   a continuous fit is not. For the Laplace activation the two never agree
#'   to that precision, because its kink at the target makes the sum over the
#'   options differ from the integral by a term of the order of the step.
#'
#'   The likelihood needs one sum over the options per combination of
#'   parameter values and option set. As in the continuous model, rows share
#'   parameter values when they share the values of all predictors, which
#'   is fastest when they are adjacent, so bmm offers to sort the data by the
#'   predictors. With a
#'   continuous predictor every trial has its own sum; for the von Mises
#'   activation on a wheel of at least 90 options the sum over the full wheel
#'   is then replaced by the continuous normaliser wherever \eqn{w} is at least
#'   2.25 steps, which changes the log-likelihood by less than 1e-10 and costs
#'   about as much as the
#'   continuous model (1.1 times as much per gradient on the 15,200 trials of
#'   [oberauer_lin_2017], against 7.6 times without the replacement). The
#'   Laplace activation has no such replacement, and there the full wheel of
#'   360 options cost 9.7 times as much as the continuous model per gradient.
#'   `log_lik()`, and so `loo()`, sums over all options of a row for every
#'   posterior draw, which with a continuous predictor takes a while on large
#'   data sets. `posterior_predict()` returns option offsets from the
#'   target in steps; multiply by `2 * pi / n_options` for radians.
#'
#'   Fits from before `n_options` existed are continuous, and `n_options =
#'   NULL` keeps that model unchanged.
#' @section Forced-choice option sets: With `options`, each trial offers the
#'   options in its row of the option columns, which must include the target
#'   (an offset of 0) and the response, and must sit on the grid of
#'   `n_options` positions (values within a tenth of a step are rounded onto
#'   it, with a message). The probability of the response is the softmax over
#'   the offered options only. Sets may differ in size and spacing from trial
#'   to trial, and trials that offer the same options, in any order, share one
#'   normalising sum within each combination of parameter values.
#'   `conditional_effects()` shows the parameters. A prediction of the
#'   response needs the options a trial offers, so predict with
#'   `posterior_predict()` or `predict()` on `newdata` that carries the option
#'   columns; `conditional_effects(method = "posterior_predict")` is not
#'   meaningful for option sets, because its grid of conditions holds no
#'   option set.
#'
#'   A data set could look like this, for trials with two and four options in
#'   radians:
#'
#'   | `error` | `opt1` | `opt2` | `opt3` | `opt4` |
#'   |---|---|---|---|---|
#'   | 0 | 0 | 3.142 | NA | NA |
#'   | 1.571 | 0 | 1.571 | 3.142 | -1.571 |
#'
#'   fitted with `sdm("error", n_options = 360, options = paste0("opt", 1:4))`.
#'   Robinson et al. (2023) compared such m-alternative tasks with 2, 4, 8 and
#'   16 equally spaced colours; for that design `n_alternatives = m` builds the
#'   option columns. They found that the softmax parameter varied with the
#'   number of alternatives while the d' of a signal detection model with
#'   normal noise did not; with `options`, the SDM, whose Gumbel noise makes it
#'   a softmax model, can be fitted across set sizes to test that for its
#'   memory strength `c`.
#'
#'   Robinson, M. M., DeStefano, I. C., Vul, E., & Brady, T. F. (2023). How do
#'   people build up visual memory representations from sensory evidence?
#'   Revisiting two classic models of choice. Journal of Mathematical
#'   Psychology, 117, 102805. https://doi.org/10.1016/j.jmp.2023.102805
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
#'   divides by the integral over the real line, \eqn{2 / \kappa}. `bmm`
#'   keeps the exact integral over the circle, so its `c` is that of eq. 2
#'   times \eqn{1 - e^{-\kappa \pi}} (0.998 at \eqn{\kappa = 2}, 0.79 at
#'   \eqn{\kappa = 0.5}); [c_laplace_circle2line()] converts it to the `c`
#'   of eq. 2 and [c_laplace_line2circle()] back. `"peak"` takes `c` as the
#'   activation of the target, and `"range"` as the difference between the
#'   activation of the target and that of the location opposite it.
#' @section Default priors: `kappa` and `c` mean different things in each
#'   cell, so each cell has its own `main` priors (log scale; the `effects` and
#'   `sd` priors are the same everywhere). They were chosen in a prior
#'   predictive check so that the implied distributions of the circular SD of
#'   the response error are comparable across cells, not exactly matched. The
#'   last column is the share of prior predictive distributions with a circular
#'   SD below 0.05 radians; the Laplace cells still produce more of these
#'   spikes. The von Mises peak and range cells keep the priors of the
#'   default cell, which already give them comparable error distributions
#'   (Kolmogorov-Smirnov distance 0.06 to the default cell); the Laplace cells
#'   needed their own. Set your own priors where that matters.
#'
#'   | `activation` | `scaling` | `c` | `kappa` | circular SD < 0.05 |
#'   |---|---|---|---|---|
#'   | `"vonmises"` | `"density"` | `student_t(5, 2, 0.75)` | `student_t(5, 1.75, 0.75)` | 0.07 |
#'   | `"vonmises"` | `"peak"` | `student_t(5, 2, 0.75)` | `student_t(5, 1.75, 0.75)` | 0.03 |
#'   | `"vonmises"` | `"range"` | `student_t(5, 2, 0.75)` | `student_t(5, 1.75, 0.75)` | 0.03 |
#'   | `"laplace"` | `"density"` | `student_t(5, 3.25, 0.75)` | `student_t(5, -0.5, 0.75)` | 0.16 |
#'   | `"laplace"` | `"peak"` | `student_t(5, 2.5, 0.75)` | `student_t(5, -0.5, 0.75)` | 0.11 |
#'   | `"laplace"` | `"range"` | `student_t(5, 2.25, 0.75)` | `student_t(5, -0.5, 0.75)` | 0.09 |
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
                scaling = NULL, n_options = NULL, options = NULL, n_alternatives = NULL,
                units = c("radians", "bins"), ...) {
  call <- match.call()
  .check_constructor_dots(list(...), known = "links")
  stop_missing_args()
  version <- match.arg(version)
  activation <- match.arg(activation)
  scaling <- if (is.null(scaling)) .sdm_default_scaling(activation) else match.arg(scaling, .sdm_scalings)
  units <- match.arg(units)
  option_args <- .check_option_args(n_options, options, n_alternatives, units)
  if (activation == "laplace") {
    message2(
      'sdm(activation = "laplace") is not the signal discrimination model of \\
      Oberauer (2023), which uses a von Mises activation function; it is the \\
      Laplace x Gumbel cell of his model space (Table 2). The fit is labelled \\
      accordingly.'
    )
  }
  .model_sdm(
    resp_error = resp_error, version = version, activation = activation,
    scaling = scaling, n_options = option_args$n_options, units = units, options = options,
    n_alternatives = option_args$n_alternatives, call = call, ...
  )
}

#' @rdname sdm
#' @keywords deprecated
#' @export
sdmSimple <- function(resp_error, version = "simple", ...) {
  warning2("The function `sdmSimple()` is deprecated. Please use `sdm()` instead.")
  call <- match.call()
  dots <- .check_constructor_dots(list(...), known = "links")
  stop_missing_args()
  # the deprecated alias builds the continuous model only; the dots it warned
  # about must not reach .model_sdm(), whose arguments they can name, so only
  # `links` is read -- from the checked dots, where `link =` is already `links`
  .model_sdm(
    resp_error = resp_error, version = version, activation = "vonmises",
    scaling = "density", n_options = NULL, units = "radians", options = NULL,
    n_alternatives = NULL, call = call, links = dots$links
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
  runs <- sdm_run_metadata(data, formula, model)
  attr(data, "sdm_run_metadata") <- runs
  if (inherits(model, "option_set")) {
    attr(data, "sdm_grid_metadata") <- sdm_grid_metadata(data, formula, model, runs)
  }
  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.sdm <- function(model, data, formula) {
  if (inherits(model, "option_set")) {
    return(.sdm_grid_configure(model, data, formula))
  }
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
      read_lines2(paste0(sc_path, "/sdm_kernel_functions.stan")),
      read_lines2(paste0(sc_path, "/sdm_cont_functions.stan")),
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
  # a stanfit can carry sampler settings and no parameter draws (the stub
  # fits of the update() tests)
  if (inherits(model, "option_set") && has_draws(fit) && length(fit$fit@sim$samples) &&
      length(fit$fit@sim$samples[[1]])) {
    .sdm_check_grid(model, brms::prepare_predictions(fit))
  }
  fit
}

# The discrete model is exact whatever the grid, but it equals the continuous
# one only where the activation spans a few steps (?sdm, "Discrete response
# options"); w = 1 / sqrt(max(c*, 1) kappa) is the width that decides it
.sdm_check_grid <- function(model, prep) {
  if (!inherits(model, "option_set")) {
    return(invisible())
  }
  cell <- .sdm_prep_cell(prep)
  step <- 2 * pi / prep$data$sdm_K
  narrowest <- rep(Inf, prep$ndraws)
  for (idx in split(seq_len(prep$nobs), ceiling(seq_len(prep$nobs) / 1000))) {
    c <- matrix(brms::get_dpar(prep, "c", i = idx), nrow = prep$ndraws)
    kappa <- matrix(brms::get_dpar(prep, "kappa", i = idx), nrow = prep$ndraws)
    log_cstar <- .sdm_log_cstar(log(c), kappa, cell$activation, cell$scaling)
    w <- 1 / sqrt(pmax(exp(log_cstar), 1) * kappa)
    narrowest <- pmin(narrowest, matrixStats::rowMins(matrix(w, nrow = prep$ndraws)))
  }
  narrow <- narrowest < step
  warnif(
    any(narrow),
    "In {sum(narrow)} of {length(narrow)} posterior draws, the activation is narrower \
    than one step of the response grid (2 * pi / {prep$data$sdm_K}): the response \
    grid is coarser than the fitted precision. The discrete likelihood is exact, but \
    comparisons with continuous fits are not."
  )
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
  data <- data[sdm_kept_rows(data, formula, model), , drop = FALSE]
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

############################################################################# !
# DISCRETE RESPONSE OPTIONS (#558, K2)                                   ####
############################################################################# !

# the rows brms keeps: it drops those with a missing value in any model
# variable before fitting, so runs and cells are computed on the rest
sdm_kept_rows <- function(data, formula, model) {
  model_vars <- unique(c(unlist(model$resp_vars), data_predictor_vars(data, formula)))
  model_vars <- model_vars[model_vars %in% colnames(data)]
  if (length(model_vars) == 0L) {
    return(rep(TRUE, nrow(data)))
  }
  stats::complete.cases(data[model_vars])
}

# The cells of a model with response options (see .option_cells()) and the
# option sets they use. Set 0 is the full wheel
sdm_grid_metadata <- function(data, formula, model, runs) {
  kept <- data[sdm_kept_rows(data, formula, model), , drop = FALSE]
  run_id <- rep(seq_along(runs$sdm_run_start), runs$sdm_run_count)
  width <- .option_width(model)
  sets <- if (width > 0) {
    .option_set_table(as.matrix(kept[.option_columns(width)]), model$other_vars$n_options)
  } else {
    list(set_id = rep(0L, nrow(kept)), size = integer(0), off = matrix(0L, 0, 1))
  }
  cells <- .option_cells(run_id, sets$set_id)
  list(
    G_sdm_cells = length(cells$first),
    sdm_cell_first = as.integer(cells$first),
    sdm_cell_count = as.integer(cells$count),
    sdm_cell_set = as.integer(cells$set),
    S_sdm_sets = length(sets$size),
    M_sdm_sets = ncol(sets$off),
    sdm_set_size = as.integer(sets$size),
    sdm_set_off = sets$off
  )
}

.sdm_grid_configure <- function(model, data, formula) {
  n_options <- model$other_vars$n_options
  grid <- attr(data, "sdm_grid_metadata") %||%
    sdm_grid_metadata(data, formula, model, sdm_run_metadata(data, formula, model))
  n_rows <- sum(grid$sdm_cell_count)
  # a cell of an offered subset sums over its few options, so only cells of
  # the full wheel (set 0) are worth reporting
  n_full_cells <- sum(grid$sdm_cell_set == 0)
  if (n_full_cells > n_rows / 2) {
    message2(
      "{n_full_cells} of the {n_rows} rows have a normaliser of their own, because \\
      no other row shares their parameter values. Each normaliser sums over the \\
      options of its row in every gradient, so the fit takes longer than when rows \\
      share them; see ?sdm, section \"Discrete response options\"."
    )
  }
  sdm_family <- brms::custom_family(
    name = "sdm_grid",
    dpars = c("mu", "c", "kappa"),
    links = c("tan_half", "identity", "log"),
    lb = c(NA, NA, NA),
    ub = c(NA, NA, NA),
    type = "int", loop = FALSE,
    vars = c("sdm_K", "sdm_activation", "sdm_scaling"),
    log_lik = log_lik_sdm_grid,
    posterior_predict = posterior_predict_sdm_grid,
    posterior_epred = posterior_epred_sdm_grid
  )
  sc_path <- system.file("stan_chunks", package = "bmm")
  chunk <- function(file) read_lines2(file.path(sc_path, file))
  likelihood_file <- if (brms_slices_likelihood()) "sdm_grid_likelihood_threaded.stan" else "sdm_grid_likelihood.stan"
  step <- 2 * pi / n_options
  half <- n_options / 2
  gl <- .gl_rule(16L)
  wide <- .sdm_wide_rule()
  stanvars <- brms::stanvar(scode = paste(chunk("sdm_kernel_funs.stan"), chunk("sdm_grid_funs.stan"), sep = "\n"),
                            block = "functions") +
    brms::stanvar(scode = chunk("sdm_grid_tdata.stan"), block = "tdata") +
    sdm_stanvar_int(n_options, "sdm_K") +
    sdm_stanvar_int(match(model$other_vars$activation, .sdm_activations), "sdm_activation") +
    sdm_stanvar_int(match(model$other_vars$scaling, .sdm_scalings), "sdm_scaling") +
    brms::stanvar(x = (0:half) * step, name = "sdm_grid_half") +
    brms::stanvar(x = log(c(1, rep(2, half - 1), 1)), name = "sdm_grid_logmult") +
    brms::stanvar(x = ((1 - half):half) * step, name = "sdm_grid_full") +
    brms::stanvar(x = gl$nodes, name = "sdm_gl_x") +
    brms::stanvar(x = gl$weights, name = "sdm_gl_w") +
    brms::stanvar(x = wide$cosm1, name = "sdm_wide_cosm1") +
    brms::stanvar(x = wide$weights, name = "sdm_wide_w") +
    brms::stanvar(x = grid$G_sdm_cells, name = "G_sdm_cells") +
    sdm_stanvar_int_array(grid$sdm_cell_first, "sdm_cell_first", "G_sdm_cells") +
    sdm_stanvar_int_array(grid$sdm_cell_count, "sdm_cell_count", "G_sdm_cells") +
    sdm_stanvar_int_array(grid$sdm_cell_set, "sdm_cell_set", "G_sdm_cells") +
    brms::stanvar(x = grid$S_sdm_sets, name = "S_sdm_sets") +
    brms::stanvar(x = grid$M_sdm_sets, name = "M_sdm_sets") +
    sdm_stanvar_int_array(grid$sdm_set_size, "sdm_set_size", "S_sdm_sets") +
    sdm_stanvar_int_matrix(grid$sdm_set_off, "sdm_set_off", "S_sdm_sets", "M_sdm_sets") +
    brms::stanvar(scode = chunk(likelihood_file), block = "likelihood", position = "end")
  formula <- bmf2bf(model, formula)
  formula$family <- sdm_family
  nlist(formula, data, stanvars)
}

sdm_stanvar_int_matrix <- function(x, name, rows, cols) {
  x <- array(as.integer(x), dim = dim(x))
  out <- brms::stanvar(x = x, name = name)
  out[[name]]$scode <- paste0("array[", rows, ", ", cols, "] int ", name, ";")
  out[[name]]$pll_args <- paste("data array[,] int", name)
  out
}

# The offered options of observation i as offsets in steps, NULL for the full
# wheel: the vint columns of the brms data, without the padding
.sdm_prep_options <- function(prep, i) {
  columns <- grep("^vint[0-9]+$", names(prep$data), value = TRUE)
  if (length(columns) == 0) {
    return(NULL)
  }
  columns <- columns[order(as.integer(sub("vint", "", columns)))]
  offsets <- vapply(columns, function(column) prep$data[[column]][i], integer(1), USE.NAMES = FALSE)
  offsets[offsets != .option_absent(prep$data$sdm_K)]
}

# the log of the sum over a set of options of exp(c* s(x - mu)), one value per
# element of the parameter vectors; NULL options is the full wheel. In the
# form peak + log sum exp(peak expm1(log s)), which loses no digits at a
# large peak; with mu = 0 the full wheel takes each distance once with its
# multiplicity
.sdm_grid_log_sum <- function(options, mu, log_cstar, kappa, activation, n_options) {
  n <- max(length(mu), length(log_cstar), length(kappa))
  mu <- rep_len(mu, n)
  peak <- exp(rep_len(log_cstar, n))
  kappa <- rep_len(kappa, n)
  step <- 2 * pi / n_options
  half <- n_options / 2
  folded <- is.null(options) && all(mu == 0, na.rm = TRUE)
  offsets <- if (folded) 0:half else if (is.null(options)) (1 - half):half else options
  log_mult <- if (folded) log(c(1, rep(2, half - 1), 1)) else rep(0, length(offsets))
  out <- rep(NA_real_, n)
  for (rows in split(seq_len(n), ceiling(seq_len(n) / 2000))) {
    delta <- outer(-mu[rows], offsets * step, `+`)
    terms <- peak[rows] * expm1(.sdm_log_s(delta, kappa[rows], activation)) +
      rep(log_mult, each = length(rows))
    out[rows] <- peak[rows] + matrixStats::rowLogSumExps(terms)
  }
  out
}

# log P of the option x (radians) is taken at, among the offered options,
# minus log(2 pi / n_options): the density scale of the continuous model
.sdm_grid_log_density <- function(x, mu, log_cstar, kappa, activation, n_options, options = NULL) {
  bins <- .bin_options(x, n_options, "radians")
  numerator <- exp(log_cstar + .sdm_log_s(bins * 2 * pi / n_options - mu, kappa, activation))
  out <- numerator - .sdm_grid_log_sum(options, mu, log_cstar, kappa, activation, n_options) -
    log(2 * pi / n_options)
  if (!is.null(options)) {
    out[!is.na(bins) & !bins %in% options] <- -Inf
  }
  out
}

# brms calls log_lik and posterior_predict once per observation, and the rows
# of one cell share their draws and options, so the last sums are kept
.sdm_grid_cache <- new.env(parent = emptyenv())

.sdm_grid_prep_log_sum <- function(prep, i, cell, options, mu, log_cstar, kappa) {
  key <- list(mu, log_cstar, kappa, options, cell, prep$data$sdm_K)
  if (!identical(.sdm_grid_cache$key, key)) {
    .sdm_grid_cache$value <- .sdm_grid_log_sum(options, mu, log_cstar, kappa, cell$activation,
                                               prep$data$sdm_K)
    .sdm_grid_cache$key <- key
  }
  .sdm_grid_cache$value
}

.sdm_grid_prep_draws <- function(i, prep) {
  cell <- .sdm_prep_cell(prep)
  mu <- brms::get_dpar(prep, "mu", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  log_cstar <- .sdm_log_cstar(log(brms::get_dpar(prep, "c", i = i)), kappa, cell$activation, cell$scaling)
  options <- .sdm_prep_options(prep, i)
  log_sum <- .sdm_grid_prep_log_sum(prep, i, cell, options, mu, log_cstar, kappa)
  nlist(cell, mu, kappa, log_cstar, options, log_sum)
}

log_lik_sdm_grid <- function(i, prep) {
  d <- .sdm_grid_prep_draws(i, prep)
  n_options <- prep$data$sdm_K
  y <- prep$data$Y[i]
  # newdata that carries the built columns is not checked again, so a
  # response the trial did not offer can reach here; it has no probability,
  # as in dsdm(options = )
  if (!is.null(d$options) && !is.na(y) && !y %in% d$options) {
    return(rep(-Inf, prep$ndraws))
  }
  exp(d$log_cstar + .sdm_log_s(y * 2 * pi / n_options - d$mu, d$kappa, d$cell$activation)) -
    d$log_sum - log(2 * pi / n_options)
}

# offsets in steps, as the brms response; multiply by 2 * pi / n_options for
# radians
posterior_predict_sdm_grid <- function(i, prep, ...) {
  d <- .sdm_grid_prep_draws(i, prep)
  n_options <- prep$data$sdm_K
  offsets <- d$options %||% seq(1 - n_options / 2, n_options / 2)
  .sdm_draw_options(offsets, d$mu, d$log_cstar, d$kappa, d$cell$activation, n_options)
}

posterior_epred_sdm_grid <- posterior_epred_undefined("sdm")

# one option per element of the parameter vectors, by inverse CDF over the
# offsets; the parameters are recycled to n
.sdm_draw_options <- function(offsets, mu, log_cstar, kappa, activation, n_options, n = NULL) {
  n <- n %||% max(length(mu), length(log_cstar), length(kappa))
  mu <- rep_len(mu, n)
  log_cstar <- rep_len(log_cstar, n)
  kappa <- rep_len(kappa, n)
  out <- integer(n)
  for (rows in split(seq_len(n), ceiling(seq_len(n) / 2000))) {
    delta <- outer(-mu[rows], offsets * 2 * pi / n_options, `+`)
    peak <- exp(log_cstar[rows])
    weights <- exp(peak * expm1(.sdm_log_s(delta, kappa[rows], activation)))
    cum <- matrixStats::rowCumsums(weights)
    draw <- stats::runif(length(rows)) * cum[, ncol(cum)]
    out[rows] <- as.integer(offsets[rowSums(cum < draw) + 1])
  }
  out
}
