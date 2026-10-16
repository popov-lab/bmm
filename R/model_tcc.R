############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_tcc <- function(resp_error = NULL, similarity = NULL, n_options = 360,
                       units = "radians", noise = "gaussian", d_max = 30, links = NULL,
                       call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      # table = "r" passes the d table from R instead of building it in Stan;
      # an internal switch for table transforms (correlated or motor noise)
      other_vars = nlist(similarity, n_options, units, noise, d_max, table = "stan"),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Target Confusability Competition (TCC) model by Schurgin, Wixted & Brady (2020)",
      citation = glue(
        "Schurgin, M. W., Wixted, J. T., & Brady, T. F. (2020). Psychophysical \\
        scaling reveals a unified theory of visual memory strength. Nature Human \\
        Behaviour, 4(11), 1156-1172. https://doi.org/10.1038/s41562-020-00938-0"
      ),
      version = "NA",
      requirements = glue(
        "- The response is the error relative to the target, in radians (binned \\
        to the nearest of `n_options` options) or, with units = 'bins', in whole \\
        option steps\n",
        "- `similarity` is 'laplace' or 'vonmises', with its width estimated, or \\
        gives the psychophysical similarity f of an option to the target at each \\
        distance 0, 1, ..., n_options / 2 options, scaled so that f = 1 at \\
        distance 0"
      ),
      parameters = list(
        d = if (noise == "gumbel") {
          glue(
            "Memory strength d: the mean memory-match signal of the target in units \\
            of the scale of the standard Gumbel noise (whose SD is pi / sqrt(6)); an \\
            option at distance x has mean d * f(x)"
          )
        } else {
          glue(
            "Memory strength d (the d' of the uncorrelated-noise model): the mean \\
            memory-match signal of the target in units of the noise SD; an option \\
            at distance x has mean d * f(x)"
          )
        }
      ),
      # log, not identity as in the sdt_* models: d f(x) >= 0 in TCC, and a
      # subject intercept on set-size cells of oberauer_lin_2017 mixed to Rhat
      # 1.006 on the log link against 1.018 on the identity link (#566)
      links = list(d = "log"),
      fixed_parameters = list(),
      # d of the uncorrelated model runs from 0.6 to 2.0 over set sizes 1-8 in
      # oberauer_lin_2017 with the similarity of the examples, so the prior is
      # centred at exp(0.5) = 1.65; d = 12 lies two SDs above it and the
      # default d_max = 30 almost three. Gumbel noise puts d higher: 1.4 to 6.3
      # for set sizes 8 to 1 with a von Mises similarity, within 1.4 SDs of the
      # prior's centre, so both noises share the prior (#566)
      default_priors = list(
        d = list(main = "normal(0.5, 1)", effects = "normal(0, 1)", sd = "exponential(1)")
      ),
      init_ranges = list(d = c(1, 2))
    ),
    class = c("bmmodel", "tcc"),
    call = call
  )
  shape <- .tcc_shape(similarity)
  if (!is.null(shape)) {
    out$parameters[[shape$name]] <- shape$description
    out$links[[shape$name]] <- "log"
    out$default_priors[[shape$name]] <- shape$priors
    out$init_ranges[[shape$name]] <- shape$init
  }
  set_links(out, links)
}

# The named similarities, each with the parameter it adds and the range of its
# table. The table covers what the 2-D table spike measured (TCC, 2026-10):
# from a similarity so sharp that only the target stands out, to one that falls
# almost linearly (Laplace) or as a raised cosine (von Mises) with distance. On
# oberauer_lin_2017 the spike estimated lambda at 0.6 to 1.2 radians and kappa
# at 3.7 to 4.7 over set sizes 1-8, where the priors are centred
.tcc_shapes <- list(
  laplace = list(
    name = "lambda", code = 1L, range = c(0.004, 20), init = c(0.3, 1),
    description = glue(
      "Width of the Laplace similarity, in radians of the response circle: an \\
      option x radians from the target has similarity exp(-x / lambda), rescaled \\
      to 1 at the target and 0 opposite it"
    ),
    priors = list(main = "normal(-0.5, 1)", effects = "normal(0, 1)", sd = "exponential(1)")
  ),
  vonmises = list(
    name = "kappa", code = 2L, range = c(0.02, 5e4), init = c(2, 8),
    description = glue(
      "Concentration of the von Mises similarity: an option x radians from the \\
      target has similarity exp(-kappa * (1 - cos(x))), rescaled to 1 at the \\
      target and 0 opposite it"
    ),
    priors = list(main = "normal(1.5, 1)", effects = "normal(0, 1)", sd = "exponential(1)")
  )
)

# NULL for a similarity vector
.tcc_shape <- function(similarity) {
  if (is.character(similarity)) .tcc_shapes[[similarity]]
}

#' @title `r .model_tcc()$name`
#' @name tcc
#' @description Continuous report as an `n_options`-alternative forced choice:
#'   every option on the response wheel produces a Gaussian memory-match
#'   signal with unit variance and mean `d * f(x)`, where `f(x)` is the
#'   psychophysical similarity of an option `x` options away from the target,
#'   and the response is the option with the strongest signal. This is the
#'   uncorrelated-noise model of Schurgin et al. (2020). The d' they report
#'   comes from a correlated-noise version and is larger: across the
#'   subjects and set sizes of their Figure 3 data, the median ratio of their
#'   d' to this `d` was 1.60. With `noise = "gumbel"` the signals carry
#'   standard Gumbel noise instead, and the response probabilities become a
#'   softmax (see the section Noise).
#' @details `r model_docs(.model_tcc())`
#'
#' ## Similarity
#'
#'   `similarity` says how similar an option is to the target as a function of
#'   their distance on the wheel. The model assumes the similarity depends only
#'   on the distance, the same for every target and on both sides of it. It is
#'   1 at the target and 0 at the option opposite it. There are three choices:
#'
#'   - A numeric vector holds fixed values, one per distance in option steps,
#'   in the order 0, 1, ..., `n_options / 2`: the first value belongs to the
#'   target itself and must be 1, the last to the option opposite the target.
#'   Schurgin et al. (2020) measured them in a separate perceptual task and
#'   scaled them to range from 0 to 1. With 360 options in 1-degree steps,
#'   `similarity[k + 1]` is the similarity at `k` degrees.
#'   - `"laplace"` estimates a similarity that falls exponentially with the
#'   distance `x` in radians, `exp(-x / lambda)`, and adds the parameter
#'   `lambda`, its width in radians.
#'   - `"vonmises"` estimates a similarity of the shape `exp(-kappa * (1 -
#'   cos(x)))` and adds the parameter `kappa`, its concentration.
#'
#'   The two shapes are the Laplace and von Mises activation functions of
#'   Oberauer (2021, eq. 2); his Laplace rate is `1 / lambda`. Both are scaled
#'   to run from 1 at the target to 0 opposite it, like a measured similarity.
#'   Scaling changes the value of `d` that describes the data, not the shape.
#'   The distance is measured on the response circle: with 180 options for
#'   orientation, `x = pi` is the option opposite the target, 90 degrees of
#'   orientation away.
#'
#'   `lambda` and `kappa` are modelled on the log scale and take formulas like
#'   `d`, including continuous predictors and random effects. Their default
#'   priors are `normal(-0.5, 1)` on `log(lambda)`, centred at 0.61 radians (35
#'   degrees), and `normal(1.5, 1)` on `log(kappa)`, centred at 4.5. The
#'   formula never chooses the similarity: naming `lambda` in it requires
#'   `similarity = "laplace"`, and `kappa` requires `"vonmises"`. A value fixed
#'   in the formula, such as `lambda = -1`, is on the log scale.
#'
#'   Oberauer, K. (2021). Measurement models for visual working memory: A
#'   factorial model comparison. Psychological Review.
#'   https://doi.org/10.1037/rev0000328
#'
#' ## Noise
#'
#'   `noise` sets the distribution of the noise on each option's signal.
#'   `"gaussian"`, the default, is the model of Schurgin et al. (2020), with
#'   normal noise of unit variance. Its response probabilities have no closed
#'   form, so the likelihood reads them from a table (see the grid limits
#'   below). With `"gumbel"` the noise follows the standard Gumbel (type I
#'   extreme value) distribution. The largest of independent Gumbel variables
#'   is again Gumbel distributed, which turns the probability of choosing an
#'   option into the softmax of `d * f(x)` over all options:
#'   `exp(d * f(x)) / sum(exp(d * f(x')))` (Oberauer, 2021, eq. 4). This needs
#'   no table, so there are no grid limits on `d`, `lambda` or `kappa`, and
#'   `d_max` is not used.
#'
#'   Oberauer (2021) crossed three choices: the activation function (von
#'   Mises or Laplace), the response-selection rule (Luce's choice rule, or
#'   signal detection with Gumbel or with normal noise) and whether the
#'   precision `kappa` is the same in all conditions. `tcc()` fits the cells
#'   with a signal-detection rule: `noise` selects the noise, `similarity`
#'   the activation function, and the formula for `kappa` or `lambda` holds it
#'   constant (`kappa ~ 1`) or lets it vary with the conditions. The von Mises
#'   similarity with Gumbel noise is the cell that fit best in his comparison
#'   and that he named the signal discrimination model, which [sdm()] fits on
#'   a continuous response scale. On a wheel of `n_options` options,
#'   `tcc(similarity = "vonmises", noise = "gumbel")` gives the [sdm()] density
#'   at each option's centre times the option width `2 * pi / n_options`, with
#'   the same `kappa` and `c = d * sqrt(2 * pi / kappa) / (1 - exp(-2 *
#'   kappa))` in the default parametrization of [sdm()].
#'
#'   Oberauer found the fits with Gumbel and with normal noise always very
#'   close (p. 5). Gumbel noise needs no table, so it skips the build before
#'   sampling and has no grid limits, and it is the noise to compare with
#'   [sdm()]. Per gradient it is not cheaper: on [oberauer_lin_2017], with `d`
#'   and `kappa` per subject and set size, it cost 1.4 times the normal noise
#'   with a von Mises similarity (about as much as [sdm()]) and 1.5 times with
#'   a fixed similarity; with a continuous predictor of `kappa`, which makes
#'   every trial its own cell, it cost 6.2 times (174 against 28 ms). Keep
#'   normal noise to stay with the model of Schurgin et al. and to compare `d`
#'   with their d'. The two `d` are not on one scale.
#'   The standard Gumbel distribution has an SD of `pi / sqrt(6)` = 1.28, but
#'   matching the response probabilities takes a larger factor. On the
#'   probabilities of the normal-noise model, the Gumbel `d` with the smallest
#'   Kullback-Leibler divergence was 2.6 to 2.9 times the normal `d` for
#'   values of the latter from 0.5 to 2, and 1.9 to 2.4 times at 8, for the
#'   similarities of the examples below.
#'
#' ## Response
#'
#'   With `units = "radians"`, `resp_error` is the response minus the target
#'   in radians, as for [sdm()] and [mixture2p()], and each value is assigned to
#'   the nearest of the `n_options` options (+-pi fall into the same option).
#'   With `units = "bins"`, it is the chosen option's offset from the target in
#'   option steps. The model treats the response as a choice among
#'   `n_options` discrete options, so `n_options` should be the number of
#'   options the participant could choose from.
#'
#' ## Grid limit on d
#'
#'   `d` is modelled on the log scale, so it is positive, effects on it are
#'   multiplicative, and priors on its coefficients are priors on `log(d)`.
#'   Keep those priors proper: as `d` approaches 0 the likelihood levels off
#'   at that of guessing uniformly, so under a flat prior on `log(d)` the
#'   posterior cannot be normalised, and a cell with little information
#'   drifts towards `d = 0` without end. The default prior is proper. This
#'   holds for either noise; the rest of this section and the next apply to
#'   normal noise only. The likelihood reads the response probabilities from a table on a grid
#'   of `d` values up to `d_max` in steps of 0.02, built once per chain
#'   by quadrature and interpolated with a cubic spline (error below 1e-6 in
#'   log probability). Beyond `d_max` the table is extended
#'   linearly, with the log probability of each response allowed to fall but
#'   not to rise. This keeps the sampler moving during warmup but is not the
#'   model: past `d_max` the likelihood no longer favours a larger `d`, so a
#'   posterior that reaches the end of the grid follows the prior beyond it,
#'   with estimates that mean nothing. `bmm()` warns after the fit when any
#'   posterior draw of `d` passes `d_max`; then refit with a larger `d_max`.
#'   The default of 30 lies far above the values of `d` from 0.6 to 2.0 that
#'   `tcc()` estimates for set sizes 1 to 8 in [oberauer_lin_2017] with the
#'   similarity of the examples below (`d` depends on the similarity). Each
#'   extra unit adds 50 table rows, which costs build time and memory but
#'   nothing per gradient. A large `d` is weakly identified
#'   when the similarity is broad, because the probability of reporting the
#'   target then grows slowly with `d`.
#'
#' ## Grid limits with an estimated similarity
#'
#'   With `similarity = "laplace"` or `"vonmises"` the table has two
#'   dimensions: `d` in steps of 0.1 up to `d_max`, and `log(lambda)` or
#'   `log(kappa)` in steps of 0.1. `lambda` runs from 0.004 to 21.7 radians
#'   (0.23 to 1245 degrees) and `kappa` from 0.02 to 53,500, which reaches
#'   from a similarity so sharp that no option but the target resembles it,
#'   to one that falls almost as a straight line (Laplace) or a raised cosine
#'   (von Mises) with distance. Between the grid nodes the interpolated log
#'   probability errs by at most 2e-4 on responses with a probability above
#'   1e-4 up to `d = 12` and by 6e-4 up to `d = 30`; on any response, by
#'   1.5e-3 up to `d = 12` and 1e-2 up to `d = 30`.
#'   Beyond either end of either axis the table is extended linearly, and
#'   above `d_max`, as with a fixed similarity, the log probability of a
#'   response may fall with `d` but never rise. When no response lies far
#'   from the target, a large `d` with a broad similarity describes the data
#'   almost as well as a small `d` with a sharp one, and along that ridge a
#'   rising extension would let `d` run off without bound. `bmm()` warns after
#'   the fit when any posterior draw of `d` passes `d_max`, or of `lambda` or
#'   `kappa` leaves its grid. There is no argument to widen the grid of
#'   `lambda` or `kappa`: use a prior that keeps the posterior inside.
#'
#'   Each chain builds this table once before sampling. With the default
#'   `d_max = 30` that took 10 s (Laplace) and 18 s (von Mises) per chain on
#'   the machine bmm is developed on, and each chain process held 61 MB and
#'   91 MB.
#'
#' ## Posterior predictive checks
#'
#'   `pp_check(fit)` compares observed and predicted offsets in option steps.
#'   `pp_check(fit, resp_var = "error")` shows the same on the scale of
#'   `resp_error`, which for radians is the view [sdm()] and [mixture2p()] fits
#'   give; add `type = "bars_binned"` and name a condition column in `group`
#'   for histograms per condition. `resp_var = "abs_error"` shows the distance from the target,
#'   the quantity the likelihood depends on. See [pp_check.bmmfit()].
#'   `posterior_predict()` and `predict()` return offsets in option steps
#'   whatever `units` is; multiply by `2 * pi / n_options` for radians.
#'
#' ## Comparing with continuous models
#'
#'   The log-likelihood of `tcc()` is a log probability mass over options,
#'   whereas [sdm()] and [mixture2p()] return a log density in radians. The two
#'   differ by the log of the option width, `log(2 * pi / n_options)` per
#'   observation, so summed over `N` observations by `N * log(2 * pi /
#'   n_options)`; for 360 options, `log(2 * pi / 360)` is -4.048. Before
#'   comparing LOO estimates, convert the TCC log-likelihood into a density:
#'
#'   ```
#'   ll <- brms::log_lik(fit_tcc) + log(n_options / (2 * pi))
#'   loo_tcc <- loo::loo(ll, r_eff = loo::relative_eff(exp(ll), chain_id =
#'     rep(seq_len(brms::nchains(fit_tcc)), each = brms::ndraws(fit_tcc) /
#'     brms::nchains(fit_tcc))))
#'   loo::loo_compare(loo_tcc, brms::loo(fit_sdm))
#'   ```
#'
#'   `loo_compare()` warns that the models do not share the same y variable,
#'   because `tcc()` stores the binned response; the warning is expected here.
#'   The comparison is fair only if the responses were themselves recorded on
#'   `n_options` discrete options, as on a 360-degree colour wheel.
#'
#' ## Speed
#'
#'   The fit is fastest when rows that share a value of `d` are adjacent, so
#'   bmm offers to sort the data by the predictors of `d`, as for [sdm()].
#'   Continuous predictors of `d` are supported at the per-trial cost.
#'
#'   With `lambda` or `kappa` estimated, rows are grouped by their values of
#'   both parameters. On [oberauer_lin_2017] (15,200 trials in 152 cells of
#'   subject and set size) one gradient took 2.4 ms with both parameters
#'   predicted by set size and subject, and 28 ms with a continuous predictor
#'   of `lambda` as well, which makes every trial its own cell. Estimating
#'   the similarity also about doubles the number of gradients per iteration.
#' @param resp_error The name of the variable in the data holding the response
#'   error. With `units = "radians"` (the default) it is the error relative to
#'   the target in radians, as for the other circular models, and is binned to
#'   the nearest of the `n_options` options; with `units = "bins"` it is the
#'   offset of the chosen option from the target in whole option steps (any
#'   integer, taken modulo `n_options`).
#' @param similarity Either a numeric vector of `n_options / 2 + 1` values in
#'   \[0, 1\], the fixed similarity of an option to the target at each
#'   distance 0, 1, ..., `n_options / 2` options, starting with 1 at distance
#'   0; or `"laplace"` or `"vonmises"` to estimate the similarity with the
#'   parameter `lambda` or `kappa`. See the section Similarity.
#' @param n_options The number of response options on the wheel, an even
#'   number: 360 for a colour wheel in 1-degree steps, 180 for orientation.
#' @param units Either `"radians"` or `"bins"`, see `resp_error`.
#' @param noise The noise on each option's signal: `"gaussian"` (normal, as in
#'   Schurgin et al., 2020) or `"gumbel"` (standard Gumbel, which makes the
#'   response probabilities a softmax). See the section Noise.
#' @param d_max The upper end of the `d` grid of the likelihood table; draws
#'   of `d` beyond it are reported after the fit. Not used with `noise =
#'   "gumbel"`, which needs no table.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' f <- exp(-(0:180) / 20)
#' f <- (f - min(f)) / (1 - min(f))
#' dat <- data.frame(err = rtcc(1000, d = 2, similarity = f))
#' fit <- bmm(bmf(d ~ 1), dat, tcc("err", similarity = f), cores = 4)
#'
#' # the same data, with the width of a Laplace similarity estimated
#' fit <- bmm(bmf(d ~ 1, lambda ~ 1), dat, tcc("err", similarity = "laplace"), cores = 4)
#'
#' # Gumbel noise with a von Mises similarity: Oberauer's (2021) signal
#' # discrimination model on 360 discrete options
#' fit <- bmm(bmf(d ~ 1, kappa ~ 1), dat,
#'            tcc("err", similarity = "vonmises", noise = "gumbel"), cores = 4)
tcc <- function(resp_error, similarity, n_options = 360, units = c("radians", "bins"),
                noise = c("gaussian", "gumbel"), d_max = 30, ...) {
  call <- match.call()
  stop_missing_args()
  units <- match.arg(units)
  noise <- match.arg(noise)
  .tcc_check_similarity(similarity, n_options)
  stopif(!is.numeric(d_max) || length(d_max) != 1 || !isTRUE(d_max > 0),
         "`d_max` must be a single positive number.")
  .model_tcc(resp_error = resp_error, similarity = similarity, n_options = n_options,
             units = units, noise = noise, d_max = d_max, call = call, ...)
}

.tcc_check_similarity <- function(similarity, n_options) {
  stopif(!is.numeric(n_options) || length(n_options) != 1 || !isTRUE(n_options >= 4) ||
           n_options %% 2 != 0,
         "`n_options` must be an even number of at least 4.")
  n_values <- n_options / 2 + 1
  if (is.character(similarity) && length(similarity) == 1 && similarity %in% names(.tcc_shapes)) {
    return(invisible())
  }
  stopif(!is.numeric(similarity),
         "`similarity` must be 'laplace', 'vonmises' or a numeric vector of similarities.")
  stopif(length(similarity) != n_values,
         "`similarity` must hold {n_values} values, one for each distance from 0 to \\
         {n_options / 2} options, but it holds {length(similarity)}.")
  stopif(anyNA(similarity) || any(similarity < 0 | similarity > 1),
         "`similarity` must lie between 0 and 1.")
  stopif(similarity[1] != 1,
         "`similarity` must be 1 at distance 0 (the target), as d is the target's signal.")
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.tcc <- function(model, data, formula) {
  resp <- model$resp_vars$resp_error
  units <- model$other_vars$units
  stopif(not_in(resp, colnames(data)),
         "The response variable '{resp}' is not present in the data.")
  stopif(!is.numeric(data[[resp]]),
         "The response variable '{resp}' must be numeric.")
  stopif(units == "radians" && isTRUE(max(abs(data[[resp]]), na.rm = TRUE) > 2 * pi),
         "The response variable '{resp}' holds values beyond 2*pi, so it is not in \\
         radians. Convert it with deg2rad(), or set units = 'bins' if it counts options.")
  stopif(units == "bins" && any(data[[resp]] %% 1 != 0, na.rm = TRUE),
         "With units = 'bins', the response variable '{resp}' must hold whole numbers.")
  warnif("tcc_bin" %in% colnames(data),
         "Column 'tcc_bin' in your data is reserved by the TCC model and will be overwritten.")
  data$tcc_bin <- .tcc_bin(data[[resp]], model$other_vars$n_options, units)
  data <- order_data_query(model, data, formula)
  NextMethod("check_data")
}

# signed offset of the chosen option from the target, in (-n_options / 2, n_options / 2]
.tcc_bin <- function(x, n_options, units) {
  steps <- if (units == "radians") floor(x * n_options / (2 * pi) + 0.5) else x
  half <- n_options / 2
  as.integer(half - (half - steps) %% n_options)
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.tcc <- function(model, formula) {
  brms::bf("tcc_bin ~ 0")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.tcc <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)
  shape <- .tcc_shape(model$other_vars$similarity)
  gumbel <- .tcc_gumbel(model)
  formula$family <- brms::custom_family(
    paste0("tcc", if (gumbel) "_gumbel", if (!is.null(shape)) "_parametric"),
    dpars = c("mu", "d", shape$name),
    links = c("identity", unlist(model$links[c("d", shape$name)], use.names = FALSE)),
    type = "int",
    loop = FALSE,
    log_lik = if (gumbel) log_lik_tcc_gumbel else log_lik_tcc,
    posterior_predict = if (gumbel) posterior_predict_tcc_gumbel else posterior_predict_tcc,
    posterior_epred = posterior_epred_tcc,
    vars = .tcc_family_vars(gumbel, shape)
  )
  stanvars <- .tcc_stanvars(model, model$other_vars$table)
  nlist(formula, data, stanvars)
}

.tcc_gumbel <- function(model) {
  identical(model$other_vars$noise, "gumbel")
}

.tcc_family_vars <- function(gumbel, shape) {
  if (gumbel) {
    return(if (is.null(shape)) "tcc_similarity" else c("tcc_delta", "tcc_shape"))
  }
  if (is.null(shape)) {
    return(c("tcc_tab", "tcc_grid_lo", "tcc_grid_step"))
  }
  c("tcc_tab2", "tcc_grid_lo", "tcc_grid_step", "tcc_u_lo", "tcc_u_step")
}

# One kernel, two sources for its table: built in transformed data from the
# similarity, or computed in R and passed as data under the same name. A named
# similarity has a 2-D table, always built in transformed data: from R it would
# be 40-70 MB of data per fit
.tcc_stanvars <- function(model, table) {
  if (.tcc_gumbel(model)) {
    return(.tcc_gumbel_stanvars(model))
  }
  grid <- .tcc_grid(model)
  sc_path <- system.file("stan_chunks", package = "bmm")
  out <- brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_funs.stan")),
                       block = "functions") +
    brms::stanvar(x = grid$lo, name = "tcc_grid_lo") +
    brms::stanvar(x = grid$step, name = "tcc_grid_step")
  shape <- .tcc_shape(model$other_vars$similarity)
  if (!is.null(shape)) {
    return(out + .tcc_shape_stanvars(model, shape, grid, sc_path))
  }
  similarity <- .tcc_similarity(model)
  if (table == "r") {
    return(out + brms::stanvar(x = .tcc_table(similarity, grid), name = "tcc_tab"))
  }
  rule <- .gl_rule(.tcc_gl_nodes)
  out +
    brms::stanvar(x = grid$rows, name = "tcc_grid_rows") +
    brms::stanvar(x = similarity, name = "tcc_similarity") +
    brms::stanvar(x = rule$nodes, name = "tcc_gl_nodes") +
    brms::stanvar(x = log(rule$weights), name = "tcc_gl_logw") +
    brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_tdata.stan")),
                  block = "tdata", pll_args = "data matrix tcc_tab")
}

.tcc_shape_stanvars <- function(model, shape, grid, sc_path) {
  ugrid <- .tcc_ugrid(shape)
  rule <- .gl_rule(.tcc_gl_nodes)
  n_options <- model$other_vars$n_options
  brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_parametric_funs.stan")),
                block = "functions") +
    brms::stanvar(x = grid$rows, name = "tcc_grid_rows") +
    brms::stanvar(x = ugrid$lo, name = "tcc_u_lo") +
    brms::stanvar(x = ugrid$step, name = "tcc_u_step") +
    brms::stanvar(x = ugrid$rows, name = "tcc_u_rows") +
    brms::stanvar(x = shape$code, name = "tcc_shape") +
    brms::stanvar(x = 2 * pi * (0:(n_options / 2)) / n_options, name = "tcc_delta") +
    brms::stanvar(x = rule$nodes, name = "tcc_gl_nodes") +
    brms::stanvar(x = log(rule$weights), name = "tcc_gl_logw") +
    brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_parametric_tdata.stan")),
                  block = "tdata", pll_args = "data array[] matrix tcc_tab2")
}

# Gumbel noise has a closed form and needs no table (#566)
.tcc_gumbel_stanvars <- function(model) {
  sc_path <- system.file("stan_chunks", package = "bmm")
  funs <- function(file) {
    brms::stanvar(scode = read_lines2(file.path(sc_path, file)), block = "functions")
  }
  shape <- .tcc_shape(model$other_vars$similarity)
  if (is.null(shape)) {
    return(funs("tcc_gumbel_funs.stan") +
             brms::stanvar(x = .tcc_similarity(model), name = "tcc_similarity"))
  }
  n_options <- model$other_vars$n_options
  funs("tcc_gumbel_funs.stan") + funs("tcc_gumbel_parametric_funs.stan") +
    brms::stanvar(x = shape$code, name = "tcc_shape") +
    brms::stanvar(x = 2 * pi * (0:(n_options / 2)) / n_options, name = "tcc_delta")
}

############################################################################# !
# THE LIKELIHOOD TABLE                                                   ####
############################################################################# !

# The similarity vector of a model with a fixed similarity. The similarity
# argument selects fixed or named, never the formula (#566)
.tcc_similarity <- function(model) {
  model$other_vars$similarity
}

# A named similarity at distances `delta` in radians, scaled to 1 at the target
# and 0 opposite it, as Schurgin et al. (2020) scale theirs. The shapes are
# Oberauer's (2021, eq. 2 and footnote 2): exp(-delta / lambda) and
# exp(-kappa * (1 - cos(delta))). Scaling is an affine map of the signals, and
# the part that is the same for every option does not change which one is
# largest, so it only rescales d. expm1() and 1 - cos(x) = 2 sin(x / 2)^2 keep
# the forms exact at both ends of the grids, where the scaled shapes tend to a
# straight line (Laplace) or a raised cosine (von Mises).
# tcc_similarity() in tcc_parametric_funs.stan is the Stan twin
.tcc_shape_similarity <- function(similarity, par, delta) {
  if (similarity == "laplace") {
    return(exp(-delta / par) * expm1(-(pi - delta) / par) / expm1(-pi / par))
  }
  s <- sin(delta / 2)^2
  exp(-2 * par * s) * expm1(-2 * par * (1 - s)) / expm1(-2 * par)
}

# The quadrature error peaks at small |d|, where the winning signal is a
# narrow bump: the probabilities sum to 1 within 6e-7 at 100 nodes, 1e-8 at
# 120 and 2e-10 at 140 (d from -2 to 12, n_options = 360). 120 nodes cost
# 1.2 times the table build and nothing per gradient. A d step of 0.02 adds
# about 1e-7 of interpolation error (TCC cost spike, 2026-10). Row 1 and the
# last row pad the grid by one step, so that every node inside it has a
# central-difference slope
.tcc_gl_nodes <- 120L

# The grid starts at 0 under a link that keeps d positive (the default log):
# rows below 0 are never read and would hide the lower end from the warning on
# draws beyond the table. Beyond the grid the table is extended along its end
# slopes. A rising log P extended in a straight line outgrows any prior on
# log(d) and makes the posterior improper, so the padding row above the grid
# caps each response's end slope at 0 (.tcc_table(), tcc_tdata.stan), and one
# extra node keeps the last cell below d_max exact (#566)
.tcc_grid <- function(model) {
  lo <- if (.link_ranges[[model$links$d]][1] >= 0) 0 else -2
  d_max <- model$other_vars$d_max
  if (is.null(.tcc_shape(model$other_vars$similarity))) {
    return(.tcc_axis(lo, d_max + 0.02, 0.02))
  }
  # In the 2-D table, d steps of 0.1 add 1.5e-6 of error against 1.2e-3 from
  # the steps of 0.1 in log lambda or log kappa (TCC 2-D table spike, 2026-10).
  # With an estimated similarity, d and lambda or kappa can trade off along a
  # ridge that leaves the grid when no response is far from the target: the
  # log-likelihood of such data was flat within 3 nats from (d, kappa) = (3, 20)
  # to (100, 0.15), and before the slope cap d ran to 1e5 in recovery (#566)
  .tcc_axis(lo, d_max + 0.1, 0.1)
}

# The grid of u = log(lambda) or log(kappa). Up to d = 12 a step of 0.1 errs by
# up to 1.2e-3 in log P on any response and 1.7e-4 on responses with P > 1e-4,
# measured at random points: at cell midpoints, where the spike measured, the
# leading error of the cubic vanishes (#566). 0.05 is 16 times more accurate
# for twice the build time and memory and the same cost per gradient (TCC 2-D
# table spike, 2026-10)
.tcc_ugrid <- function(shape) {
  .tcc_axis(log(shape$range[1]), log(shape$range[2]), 0.1)
}

# value of row g is lo + (g - 2) * step; rows 1 and `rows` pad the grid
.tcc_axis <- function(lo, hi, step) {
  rows <- as.integer(ceiling(round((hi - lo) / step, 6)) + 3)
  nlist(lo, step, rows, top = lo + (rows - 3) * step)
}

# log P of responding one particular option at each distance 0, ..., n / 2
# from the target. With z_i = x - d f_i,
#   P(r) = int phi(z_r) prod_{i != r} Phi(z_i) dx = int phi(z_r) / Phi(z_r) G(x) dx,
# where G = prod_i Phi(z_i) is shared by all options, and log G is a sum over
# distances weighted by the number of options at each (1, 2, ..., 2, 1).
# phi(z_r) kills every integrand above d f_r + 8 <= m + 8, with m = max(d f).
# Below, the Phi factor of the largest mean kills it from about m - 8, except
# for an option far from the target under a sharp similarity: its integrand
# peaks halfway between its own mean and m, so the window reaches down to
# (min(d f) + m) / 2 - 6 as well (1.5e-4 off in log P at d = 12 with a 1-degree
# Laplace similarity before, 1e-13 after; TCC 2-D table spike, 2026-10).
# tcc_logp_by_distance() in tcc_funs.stan is the Stan twin
.tcc_logp <- function(d, similarity) {
  n_dist <- length(similarity)
  mult <- c(1, rep(2, n_dist - 2), 1)
  means <- d * similarity
  top <- max(means)
  half <- (top + 8 - min(top - 8, (min(means) + top) / 2 - 6)) / 2
  rule <- .gl_rule(.tcc_gl_nodes)
  z <- outer(top + 8 - half + half * rule$nodes, means, "-")
  lphi <- stats::pnorm(z, log.p = TRUE)
  terms <- log(half * rule$weights) + drop(lphi %*% mult) + stats::dnorm(z, log = TRUE) - lphi
  matrixStats::colLogSumExps(terms)
}

# how many options lie at each distance 0, 1, ..., n / 2 from the target
.tcc_mult <- function(n_dist) {
  c(1, rep(2, n_dist - 2), 1)
}

# log P of one option at each distance under standard Gumbel noise. The largest
# of independent Gumbel(mu_x, 1) signals is Gumbel(log sum exp(mu), 1), and the
# argmax integral reduces to the softmax of the means (Yellott, 1977;
# Oberauer, 2021, eq. 4; derivation in the TCC design note, #566).
# tcc_gumbel_lpmf() in tcc_gumbel_funs.stan is the Stan twin
.tcc_gumbel_logp <- function(d, similarity) {
  a <- d * similarity
  a - matrixStats::logSumExp(log(.tcc_mult(length(similarity))) + a)
}

.tcc_table <- function(similarity, grid) {
  d_values <- grid$lo + (seq_len(grid$rows) - 2) * grid$step
  tab <- t(vapply(d_values, .tcc_logp, numeric(length(similarity)), similarity = similarity))
  tab[grid$rows, ] <- pmin(tab[grid$rows, ], tab[grid$rows - 2, ])
  tab
}

############################################################################# !
# POSTPROCESS: DRAWS BEYOND THE TABLE                                    ####
############################################################################# !

# Past the end of a table the likelihood is extrapolated, so a posterior that
# reaches it can run away: with the linear extension before the slope cap, a
# true d of 11.5 came out at 236 with d_max = 12 and at 13.2 [11.8, 14.7] with
# d_max = 30 (#566)
#' @export
postprocess_brm.tcc <- function(model, fit, ...) {
  fit <- NextMethod("postprocess_brm")
  # a stanfit can carry sampler settings and no parameter draws (the stub
  # fits of the update() tests)
  if (has_draws(fit) && length(fit$fit@sim$samples) && length(fit$fit@sim$samples[[1]])) {
    .tcc_check_grid(model, brms::prepare_predictions(fit))
  }
  fit
}

.tcc_check_grid <- function(model, prep) {
  if (.tcc_gumbel(model)) {
    return(invisible())
  }
  d_max <- model$other_vars$d_max
  shape <- .tcc_shape(model$other_vars$similarity)
  # with an estimated similarity the table no longer rises past d_max, so a
  # posterior there is the d-lambda/kappa ridge, which a larger grid does not end
  remedy <- if (is.null(shape)) {
    glue("Refit with d_max = {2 * d_max}.")
  } else {
    glue("With no responses far from the target, d and {shape$name} trade off: a \\
         large d with a broad similarity fits almost as well as a small d with a \\
         sharp one. Use priors on log(d) and log({shape$name}) that keep the \\
         posterior inside the table.")
  }
  .tcc_warn_beyond_grid(.tcc_dpar_extremes(prep, "d"), .tcc_grid(model)$lo, d_max, "d", remedy)
  if (!is.null(shape)) {
    grid <- .tcc_ugrid(shape)
    .tcc_warn_beyond_grid(
      .tcc_dpar_extremes(prep, shape$name), exp(grid$lo), exp(grid$top), shape$name,
      glue("Outside this range the similarity hardly changes with {shape$name}: at \\
           one end no option but the target resembles it, at the other the \\
           similarity falls almost as a straight line (Laplace) or a raised cosine \\
           (von Mises) with distance. Use a prior on log({shape$name}) that keeps \\
           the posterior inside the range.")
    )
  }
}

# The smallest and largest value of a parameter over the observations, per
# draw, read in blocks that keep the draws x observations matrices small
.tcc_dpar_extremes <- function(prep, dpar) {
  lo <- rep(Inf, prep$ndraws)
  hi <- rep(-Inf, prep$ndraws)
  for (idx in split(seq_len(prep$nobs), ceiling(seq_len(prep$nobs) / 1000))) {
    values <- matrix(brms::get_dpar(prep, dpar, i = idx), nrow = prep$ndraws)
    lo <- pmin(lo, matrixStats::rowMins(values))
    hi <- pmax(hi, matrixStats::rowMaxs(values))
  }
  nlist(lo, hi)
}

.tcc_warn_beyond_grid <- function(extremes, lower, upper, name, remedy) {
  beyond <- extremes$lo < lower | extremes$hi > upper
  warnif(
    any(beyond),
    "In {sum(beyond)} of {length(beyond)} posterior draws, {name} passes the range \\
    from {signif(lower, 3)} to {signif(upper, 3)} that the likelihood table covers. \\
    Beyond it the table is extrapolated, which keeps the sampler moving but is \\
    not the model, so the posterior can run far past the table and \\
    its estimates mean nothing. {remedy}"
  )
}

############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT AND PP_CHECK                                ####
############################################################################# !

# brms calls these once per observation, so the table a fit built in
# transformed data is rebuilt in R from the standata once and kept for the next
# call. One table node takes about 1 ms, so the 1-D table costs 1.6 s and the
# 2-D one 30 to 50 s in full; only the nodes the draws reach are computed, as
# they are reached, by the quadrature the fit used
.tcc_cache <- new.env(parent = emptyenv())

# The table the fit read, from the standata in the prep: passed from R, or
# rebuilt from the similarity (fixed, or named with the u grid) and the d grid.
# The R table is an array rows of d x distances x rows of u, with one row of u
# for a fixed similarity
.tcc_prep_table <- function(prep) {
  data <- prep$data
  table <- list(lo = data$tcc_grid_lo, step = data$tcc_grid_step, rows = data$tcc_grid_rows,
                u_rows = 1L)
  if (!is.null(data$tcc_tab)) {
    table$rows <- nrow(data$tcc_tab)
    cache <- list2env(list(tab = data$tcc_tab, done = matrix(TRUE, table$rows, 1)))
    return(c(table, list(cache = cache)))
  }
  if (is.null(data$tcc_shape)) {
    similarity <- function(j) as.numeric(data$tcc_similarity)
  } else {
    table$u_lo <- data$tcc_u_lo
    table$u_step <- data$tcc_u_step
    table$u_rows <- data$tcc_u_rows
    shape <- names(.tcc_shapes)[data$tcc_shape]
    similarity <- function(j) {
      .tcc_shape_similarity(shape, exp(table$u_lo + (j - 2) * table$u_step), as.numeric(data$tcc_delta))
    }
  }
  key <- rlang::hash(list(table, data$tcc_shape, as.numeric(data$tcc_delta),
                          as.numeric(data$tcc_similarity)))
  if (!identical(.tcc_cache$key, key)) {
    .tcc_cache$tab <- array(NA_real_, c(table$rows, length(similarity(2)), table$u_rows))
    .tcc_cache$done <- matrix(FALSE, table$rows, table$u_rows)
    .tcc_cache$similarity <- similarity
    .tcc_cache$key <- key
  }
  c(table, list(cache = .tcc_cache))
}

# Catmull-Rom rows and weights at each x on one axis, as tcc_first_row() and
# tcc_row_weights() compute them in Stan, linear beyond the grid
.tcc_weights <- function(x, lo, step, n_rows) {
  g <- (x - lo) / step + 2
  first <- pmin(pmax(floor(g) - 1, 1), n_rows - 3)
  t <- g - first - 1
  w <- 0.5 * cbind(-t + 2 * t^2 - t^3, 2 - 5 * t^2 + 3 * t^3, t + 4 * t^2 - 3 * t^3, t^3 - t^2)
  below <- t < 0
  above <- t > 1
  w[below, ] <- cbind(-0.5 * t, 1, 0.5 * t, 0)[below, ]
  w[above, ] <- cbind(0, -0.5 * (t - 1), 1, 0.5 * (t - 1))[above, ]
  nlist(first, w)
}

# The table entries that interpolate each (d, u) and their weights, 4 for a
# fixed similarity and 4 x 4 for a named one; `entry` indexes the table at the
# first distance, and the distance k adds (k - 1) * rows. Nodes not yet in the
# table are computed first
.tcc_lookup <- function(table, d, u = NULL) {
  at_d <- .tcc_weights(d, table$lo, table$step, table$rows)
  a <- if (table$u_rows > 1) rep(0:3, 4) else 0:3
  rows <- outer(at_d$first, a, "+")
  w <- at_d$w[, a + 1, drop = FALSE]
  urows <- 1
  if (table$u_rows > 1) {
    b <- rep(0:3, each = 4)
    at_u <- .tcc_weights(u, table$u_lo, table$u_step, table$u_rows)
    urows <- outer(at_u$first, b, "+")
    w <- w * at_u$w[, b + 1, drop = FALSE]
  }
  node <- rows + (urows - 1) * table$rows
  .tcc_fill(table, unique(as.vector(node)))
  tab <- table$cache$tab
  list(entry = rows + (urows - 1) * table$rows * dim(tab)[2], w = w, tab = tab, rows = table$rows)
}

.tcc_fill <- function(table, node) {
  cache <- table$cache
  node <- node[!cache$done[node]]
  if (!length(node)) {
    return(invisible())
  }
  r <- (node - 1) %% table$rows + 1
  j <- (node - 1) %/% table$rows + 1
  values <- vapply(seq_along(node), function(m) {
    at <- function(g) .tcc_logp(table$lo + (g - 2) * table$step, cache$similarity(j[m]))
    # the capped padding row, as in .tcc_table() and the Stan tables (.tcc_grid())
    if (r[m] == table$rows) pmin(at(r[m]), at(r[m] - 2)) else at(r[m])
  }, numeric(dim(cache$tab)[2]))
  n_dist <- nrow(values)
  cache$tab[cbind(rep(r, each = n_dist), seq_len(n_dist), rep(j, each = n_dist))] <- values
  cache$done[node] <- TRUE
}

# log P of distance index `col` (1 = the target) at each (d, u), the R twin of
# tcc_lpmf() and tcc_parametric_lpmf()
.tcc_interp <- function(look, col) {
  rowSums(look$w * matrix(look$tab[look$entry + (col - 1) * look$rows], nrow(look$w)))
}

# signed option offsets drawn at each (d, u) from the interpolated table, by
# inverse CDF over all options at once, in chunks to bound the memory
.tcc_rbins <- function(d, u, table, n_options, chunk = 5000) {
  out <- integer(length(d))
  for (start in seq(1, length(d), by = chunk)) {
    idx <- start:min(length(d), start + chunk - 1)
    look <- .tcc_lookup(table, d[idx], u[idx])
    shift <- (seq_len(n_options / 2 + 1) - 1) * look$rows
    logp <- Reduce(`+`, lapply(seq_len(ncol(look$w)), function(m) {
      # c(): a 3-column index matrix (n_options = 4) would index the 3-D table by dimension
      look$w[, m] * matrix(look$tab[c(outer(look$entry[, m], shift, "+"))], length(idx))
    }))
    out[idx] <- .tcc_draw_offsets(logp, n_options)
  }
  out
}

# one signed offset per row of logp, which holds log P of one option at each
# distance, by inverse CDF over all options at once
.tcc_draw_offsets <- function(logp, n_options) {
  offsets <- seq(-n_options / 2 + 1, n_options / 2)
  cum <- matrixStats::rowCumsums(exp(logp[, abs(offsets) + 1, drop = FALSE]))
  draw <- stats::runif(nrow(logp)) * cum[, n_options]
  as.integer(offsets[rowSums(cum < draw) + 1])
}

# log P of one option at each distance under Gumbel noise, one row per value of
# d (and of u = log lambda or log kappa), from the similarity in the standata
.tcc_gumbel_logp_rows <- function(data, d, u) {
  f <- if (is.null(data$tcc_shape)) {
    matrix(as.numeric(data$tcc_similarity), length(d), length(data$tcc_similarity), byrow = TRUE)
  } else {
    shape <- names(.tcc_shapes)[data$tcc_shape]
    outer(exp(u), as.numeric(data$tcc_delta), function(p, x) .tcc_shape_similarity(shape, p, x))
  }
  a <- d * f
  a - matrixStats::rowLogSumExps(a + rep(log(.tcc_mult(ncol(a))), each = nrow(a)))
}

# brms calls log_lik and posterior_predict once per observation, and the
# observations of one cell share their draws, so the last rows are kept
.tcc_gumbel_logp_draws <- function(prep, i) {
  d <- brms::get_dpar(prep, "d", i = i)
  u <- .tcc_prep_u(prep, i)
  key <- list(d, u, prep$data$tcc_similarity, prep$data$tcc_delta, prep$data$tcc_shape)
  if (!identical(.tcc_cache$gumbel_key, key)) {
    .tcc_cache$gumbel_logp <- .tcc_gumbel_logp_rows(prep$data, d, u)
    .tcc_cache$gumbel_key <- key
  }
  .tcc_cache$gumbel_logp
}

.tcc_gumbel_rbins <- function(data, d, u, n_options, chunk = 5000) {
  out <- integer(length(d))
  for (start in seq(1, length(d), by = chunk)) {
    idx <- start:min(length(d), start + chunk - 1)
    out[idx] <- .tcc_draw_offsets(.tcc_gumbel_logp_rows(data, d[idx], u[idx]), n_options)
  }
  out
}

# log of the similarity parameter at observation i, NULL for a fixed similarity
.tcc_prep_u <- function(prep, i = NULL) {
  shape <- prep$data$tcc_shape
  if (is.null(shape)) {
    return(NULL)
  }
  name <- .tcc_shapes[[shape]]$name
  log(if (is.null(i)) .pp_dpar_vector(prep, name) else brms::get_dpar(prep, name, i = i))
}

log_lik_tcc <- function(i, prep) {
  look <- .tcc_lookup(.tcc_prep_table(prep), brms::get_dpar(prep, "d", i = i), .tcc_prep_u(prep, i))
  .tcc_interp(look, abs(prep$data$Y[i]) + 1)
}

posterior_predict_tcc <- function(i, prep, ...) {
  table <- .tcc_prep_table(prep)
  .tcc_rbins(brms::get_dpar(prep, "d", i = i), .tcc_prep_u(prep, i), table,
             2 * (dim(table$cache$tab)[2] - 1))
}

log_lik_tcc_gumbel <- function(i, prep) {
  .tcc_gumbel_logp_draws(prep, i)[, abs(prep$data$Y[i]) + 1]
}

posterior_predict_tcc_gumbel <- function(i, prep, ...) {
  logp <- .tcc_gumbel_logp_draws(prep, i)
  .tcc_draw_offsets(logp, 2 * (ncol(logp) - 1))
}

posterior_epred_tcc <- posterior_epred_undefined("tcc")

# pp_check(fit) without resp_var plots the signed offsets in options (the brms
# response); "error" shows them on the user's scale, which for radians is the
# view pp_check() gives for sdm() and mixture2p()
#' @export
pp_observables.tcc <- function(model) {
  radians <- model$other_vars$units == "radians"
  scale <- if (radians) 2 * pi / model$other_vars$n_options else 1
  unit <- if (radians) "radians" else "options"
  list(
    observed = c(error = "Y"),
    checks = list(
      error = .pp_observable(function(x) x$error * scale,
                             label = glue("Response error ({unit})")),
      abs_error = .pp_observable(function(x) abs(x$error * scale),
                                 label = glue("Absolute response error ({unit})"),
                                 type = "bars_binned")
    )
  )
}

#' @export
pp_simulate.tcc <- function(model, prep) {
  d <- .pp_dpar_vector(prep, "d")
  n_options <- model$other_vars$n_options
  bins <- if (.tcc_gumbel(model)) {
    .tcc_gumbel_rbins(prep$data, d, .tcc_prep_u(prep), n_options)
  } else {
    .tcc_rbins(d, .tcc_prep_u(prep), .tcc_prep_table(prep), n_options)
  }
  list(error = matrix(bins, nrow = prep$ndraws))
}

############################################################################# !
# STORED FRAME                                                           ####
############################################################################# !

# brms keeps the binned response but not the user's column, which is rebuilt
# at the bin centre
#' @exportS3Method
revert_check_data.tcc <- function(model, data) {
  resp <- model$resp_vars$resp_error
  if (!resp %in% colnames(data)) {
    data[[resp]] <- if (model$other_vars$units == "radians") {
      data$tcc_bin * 2 * pi / model$other_vars$n_options
    } else {
      data$tcc_bin
    }
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), resp)
  }
  data$tcc_bin <- NULL
  NextMethod("revert_check_data")
}
