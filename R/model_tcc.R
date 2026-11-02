############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_tcc <- function(resp_error = NULL, similarity = NULL, n_options = 360,
                       units = "radians", noise = "gaussian", d_max = 30, options = NULL,
                       n_alternatives = NULL, links = NULL, call = NULL, ...) {
  gumbel <- noise == "gumbel"
  # the option arguments enter other_vars only when given, so that a model
  # stored before they existed is the full-wheel model it was
  option_vars <- Filter(Negate(is.null), nlist(options, n_alternatives))
  schurgin <- glue(
    "Schurgin, M. W., Wixted, J. T., & Brady, T. F. (2020). Psychophysical \\
    scaling reveals a unified theory of visual memory strength. Nature Human \\
    Behaviour, 4(11), 1156-1172. https://doi.org/10.1038/s41562-020-00938-0"
  )
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      # table = "r" passes the d table from R instead of building it in Stan;
      # an internal switch for table transforms (correlated or motor noise)
      other_vars = c(nlist(similarity, n_options, units, noise, d_max, table = "stan"), option_vars),
      domain = "Working memory",
      task = "Continuous reproduction",
      # Gumbel noise is the first step from the TCC towards the SDM in Oberauer
      # (2023, eq. 4): the same model with another noise, so both are cited
      name = if (gumbel) {
        "TCC with Gumbel noise"
      } else {
        "Target Confusability Competition (TCC) model"
      },
      citation = if (gumbel) {
        c(
          schurgin,
          glue(
            "Oberauer, K. (2023). Measurement models for visual working \\
            memory\u2014A factorial model comparison. Psychological Review, \\
            130(3), 841-852. https://doi.org/10.1037/rev0000328"
          )
        )
      } else {
        schurgin
      },
      version = "NA",
      requirements = glue(
        "- The response is the error relative to the target, in radians (binned \\
        to the nearest of `n_options` options) or, with units = 'bins', in whole \\
        option steps\n",
        "- `similarity` gives the psychophysical similarity f of an option to the \\
        target at each distance 0, 1, ..., n_options / 2 options, scaled so that \\
        f = 1 at distance 0"
      ),
      parameters = list(
        d = if (gumbel) {
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
    # option_set after tcc: check_data.tcc() runs first and leaves the binning
    # to the shared contract, and bmf2bf.tcc() keeps mu out of the model
    class = c("bmmodel", "tcc", if (length(option_vars)) "option_set"),
    call = call
  )
  set_links(out, links)
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
#'   `similarity` holds the psychophysical similarity of an option to the
#'   target as a function of their distance on the wheel, one value per
#'   distance in option steps, in the order 0, 1, ..., `n_options / 2`: the
#'   first value belongs to the target itself and must be 1, the last to the
#'   option opposite the target. The model assumes the similarity depends only
#'   on the distance, the same for every target and on both sides of it. The
#'   values are fixed, not estimated; Schurgin et al. (2020) measured them in
#'   a separate perceptual matching task. With 360 options in 1-degree steps,
#'   `similarity[k + 1]` is the similarity at `k` degrees. Without a measured
#'   similarity, [tcc_similarity()] returns a von Mises or Laplace shape for a
#'   chosen `kappa`; the fit is then a TCC with that shape fixed. A model that
#'   estimates the width of the gradient from the responses is the [sdm()].
#'
#' ## Noise
#'
#'   `noise` sets the distribution of the noise on each option's signal.
#'   `"gaussian"`, the default, is the model of Schurgin et al. (2020), with
#'   normal noise of unit variance. Its response probabilities have no closed
#'   form, so on the full wheel the likelihood reads them from a table (see
#'   the section Grid limit on d); with `options` it integrates them per
#'   cell. With `"gumbel"` the noise follows the standard Gumbel (type I
#'   extreme value) distribution. The largest of independent Gumbel
#'   variables is again Gumbel distributed, which turns the probability of
#'   choosing an option into the softmax of `d * f(x)` over all options,
#'   `exp(d * f(x)) / sum(exp(d * f(x')))` (Oberauer, 2023, eq. 4). This needs
#'   no table, so `d` has no grid limit and `d_max` is not used.
#'
#'   Gumbel noise is the first step on Oberauer's (2023) path from the TCC to
#'   his signal discrimination model, and it keeps the similarity measured and
#'   fixed. Estimating the width of the similarity gradient is the next step,
#'   and that model is the [sdm()]. The two meet at a von Mises similarity:
#'   with `similarity = tcc_similarity("vonmises", kappa)` and `noise =
#'   "gumbel"`, [dtcc()] is the [sdm()] density at each option's centre times
#'   the option width `2 * pi / n_options`, with the same `kappa` and `c = d *
#'   sqrt(2 * pi / kappa) / (1 - exp(-2 * kappa))` in the default
#'   parametrization of [sdm()]. The sum over options stands in for the
#'   integral of the [sdm()], so the two agree to 1e-8 in log probability
#'   while the gradient's SD, `1 / sqrt(kappa)` radians, spans about 4 options
#'   or more (`kappa` up to 200 on 360 options), and differ by 0.013 at `kappa
#'   = 1000`.
#'
#'   Per gradient, Gumbel noise is not cheaper than normal noise: with `d` per
#'   subject and set size on [oberauer_lin_2017] it cost 1.5 times as much.
#'   It skips the table build before sampling and has no grid limit. Keep
#'   normal noise to stay with the model of Schurgin et al. (2020) and to
#'   compare `d` with their d'. The two `d` are not on one scale. The
#'   standard Gumbel distribution has an SD of `pi / sqrt(6)` = 1.28, but
#'   matching the response probabilities takes a larger factor: the Gumbel
#'   `d` with the smallest Kullback-Leibler divergence from the probabilities
#'   of the normal-noise model was 2.6 to 2.9 times the normal `d` for values
#'   of the latter from 0.5 to 2, and 1.9 to 2.4 times at 8, for Laplace
#'   similarities with `kappa` 2.9 and 1.7 and a von Mises similarity with
#'   `kappa` 4.
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
#'   `n_options` is part of the model, not a resolution setting. Each option
#'   draws its own noise, and the more options lie near the target, the more
#'   often one of them wins: at `d = 2`, with the Laplace similarity of the
#'   examples built for each wheel (`tcc_similarity(..., n_options = n)`), the
#'   probability of choosing an option at most 10 degrees from the target is
#'   0.67 on 360 options, 0.82 on 3,600 and 0.89 on 36,000. The model has no
#'   continuous limit, and values of `d` fitted with different `n_options` are
#'   not on one scale.
#'
#' ## Grid limit on d
#'
#'   `d` is modelled on the log scale, so it is positive, effects on it are
#'   multiplicative, and priors on its coefficients are priors on `log(d)`.
#'   Keep those priors proper: as `d` approaches 0 the likelihood levels off
#'   at that of guessing uniformly, so under a flat prior on `log(d)` the
#'   posterior cannot be normalised, and a cell with little information
#'   drifts towards `d = 0` without end. The default prior is proper. This
#'   holds for either noise; the rest of this section applies to normal noise
#'   on the full wheel only, not with `options`. The likelihood reads the
#'   response probabilities from a table on a grid of `d` values up to
#'   `d_max` in steps of 0.02, built once per chain
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
#' ## Forced-choice option sets
#'
#'   With `options`, each trial offers only some of the `n_options` positions,
#'   as in an m-alternative forced choice, and the response is the offered
#'   option with the strongest signal: the probabilities above are taken over
#'   the offered options only. The columns named in `options` hold the offered
#'   options of each row, relative to the target and in the units of
#'   `resp_error`, with `NA` where a row offers fewer options than there are
#'   columns. They must include the target (an offset of 0) and the response,
#'   and must sit on the grid of `n_options` positions (values within a tenth
#'   of a step are rounded onto it). Sets may differ in size and spacing from
#'   trial to trial; rows that offer the same options, in any order, share
#'   one set. A data set with two- and four-option trials on a 360-option
#'   wheel, in degrees for readability (convert them with [deg2rad()]):
#'
#'   ```
#'   err  opt1 opt2 opt3 opt4
#'     0     0  180   NA   NA
#'    90     0   90  180  -90
#'   ```
#'
#'   fitted with `tcc("err", similarity, options = paste0("opt", 1:4))`.
#'   When every trial offers `m` options spaced equally around the wheel from
#'   the target, `n_alternatives = m` builds these columns; with `options`,
#'   `n_alternatives` may instead name a column holding the number of options
#'   of each row, which is checked against the option columns. The data
#'   contract is the one of `sdm(options = )`, so the same data fit both models.
#'
#'   The design is that of Robinson et al. (2023), who offered 2, 4, 8 or 16
#'   equally spaced colours (16 options lie 22.5 degrees apart, so they need
#'   `n_options = 720`). Their critical test fixes the parameter across m and
#'   compares the two noises: with one value for all m, the normal-noise TCC
#'   fitted better than the softmax (Gumbel) model, while both fitted alike
#'   with the parameter free per m. In bmm the test compares `d ~ 1` under
#'   `noise = "gaussian"` with `d ~ 1` under `noise = "gumbel"`, on the same
#'   similarity and `options`. The Gumbel fit is their softmax model, and both
#'   log-likelihoods are log probabilities over the offered options, so
#'   `loo_compare()` reads them as they are (see the examples). Comparing
#'   `d ~ 1` with `d ~ 0 + m_condition` within one noise checks whether that
#'   model needs a parameter per m. [sdm()] with the same `options` gives this
#'   check for `c` (`c ~ 1` against `c ~ 0 + m_condition`); with `kappa`
#'   estimated it is not their softmax model, which keeps the measured
#'   similarity. `d` under normal noise and the Gumbel parameters (`d` with
#'   `noise = "gumbel"`, or `c` of [sdm()]) are not on one scale (see the
#'   section Noise), so compare the fits, or how each parameter varies across
#'   m, not their values.
#'   For either noise, the log-likelihood is a log probability over the
#'   offered options, as for the full wheel (see the section Comparing with
#'   continuous models); `sdm(options = )` reports it on the density scale,
#'   so add `log(n_options / (2 * pi))` per observation before comparing a
#'   `tcc()` fit with an [sdm()] fit.
#'
#'   Under normal noise the likelihood integrates over the offered options
#'   once per cell, a combination of a value of `d` and an option set, and
#'   needs no table, so `d_max` is not used. Measured on one chain and one
#'   thread with fixed-step NUTS on an Apple M3 Max laptop, a cell cost about
#'   0.07 ms per gradient with 2 options and 0.18 ms with 16. The design of
#'   Robinson et al. with `d` per participant and m (30 participants, 15,000
#'   trials, 120 cells) took 14 ms per gradient, and 500 trials with random
#'   foils, each trial its own cell, took 41 ms with 4 options and 87 ms
#'   with 16. Gumbel noise needs no integral and took 1 ms and 0.2 to 0.3 ms
#'   for the same designs. `bmm()` says so when there are more cells than half
#'   the rows, as with a continuous predictor of `d` or foils drawn anew on
#'   each trial. `log_lik()`, `posterior_predict()` and `pp_check()` read each
#'   row's options from the data; to predict, pass `newdata` that carries the
#'   option columns and a response that is one of the offered options (any
#'   one will do, as `posterior_predict()` ignores it). `log_lik()` and `loo()` evaluate each
#'   row on its own under normal noise, unless consecutive rows share their
#'   draws and options; for 4,000 draws on the same laptop this took 0.13 s
#'   per row with 4 options and 0.38 s with 16.
#' @param resp_error The name of the variable in the data holding the response
#'   error. With `units = "radians"` (the default) it is the error relative to
#'   the target in radians, as for the other circular models, and is binned to
#'   the nearest of the `n_options` options; with `units = "bins"` it is the
#'   offset of the chosen option from the target in whole option steps (any
#'   integer, taken modulo `n_options`).
#' @param similarity A numeric vector of `n_options / 2 + 1` values in \[0, 1\]:
#'   the similarity of an option to the target at each distance 0, 1, ...,
#'   `n_options / 2` options, starting with 1 at distance 0.
#' @param n_options The number of response options on the wheel, an even
#'   number: 360 for a colour wheel in 1-degree steps, 180 for orientation.
#' @param units Either `"radians"` or `"bins"`, see `resp_error`.
#' @param noise The noise on each option's signal: `"gaussian"` (normal, as in
#'   Schurgin et al., 2020) or `"gumbel"` (standard Gumbel, which makes the
#'   response probabilities a softmax). See the section Noise.
#' @param d_max The upper end of the `d` grid of the likelihood table; draws
#'   of `d` beyond it are reported after the fit. Not used with `noise =
#'   "gumbel"` or with `options`, which need no table.
#' @param options The names of the columns holding the options offered on
#'   each trial, one column per option, in the units of `resp_error` and
#'   relative to the target; rows that offer fewer options hold `NA` in the
#'   remaining columns. See the section "Forced-choice option sets".
#' @param n_alternatives With `options`, optionally the name of a column
#'   holding the number of options offered on each trial, which is checked
#'   against the option columns. Without `options`, a single number `m`: each
#'   trial offers `m` options spaced equally around the wheel, the target
#'   among them.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references Robinson, M. M., DeStefano, I. C., Vul, E., & Brady, T. F.
#'   (2023). How do people build up visual memory representations from sensory
#'   evidence? Revisiting two classic models of choice. *Journal of
#'   Mathematical Psychology, 117*, 102805.
#'   https://doi.org/10.1016/j.jmp.2023.102805
#' @keywords bmmodel
#' @export
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' f <- tcc_similarity("laplace", kappa = 1 / deg2rad(20))
#' dat <- data.frame(err = rtcc(1000, d = 2, similarity = f))
#' fit <- bmm(bmf(d ~ 1), dat, tcc("err", similarity = f), cores = 4)
#'
#' # Gumbel noise with a fixed von Mises similarity
#' f_vm <- tcc_similarity("vonmises", kappa = 4)
#' fit <- bmm(bmf(d ~ 1), dat, tcc("err", similarity = f_vm, noise = "gumbel"), cores = 4)
#'
#' # Forced choice among 2, 4, 8 or 16 equally spaced options (Robinson et al.,
#' # 2023), 125 trials each, on a wheel of 720 half-degree positions
#' f720 <- tcc_similarity("vonmises", kappa = 4, n_options = 720)
#' m <- c(2, 4, 8, 16)
#' afc <- data.frame(m_condition = factor(rep(m, each = 125)))
#' offered <- lapply(m, function(k) (seq_len(k) - 1) * 2 * pi / k)
#' afc[paste0("opt", 1:16)] <- t(vapply(rep(offered, each = 125), function(x) {
#'   c(x, rep(NA, 16 - length(x)))
#' }, numeric(16)))
#' afc$err <- unlist(lapply(offered, function(x) {
#'   rtcc(125, d = 1.8, similarity = f720, n_options = 720, options = x)
#' }))
#' tcc_afc <- tcc("err", similarity = f720, n_options = 720, options = paste0("opt", 1:16))
#' fit_d1 <- bmm(bmf(d ~ 1), afc, tcc_afc, cores = 4)
#' fit_dm <- bmm(bmf(d ~ 0 + m_condition), afc, tcc_afc, cores = 4)
#' tcc_gumbel <- tcc("err", similarity = f720, n_options = 720, options = paste0("opt", 1:16),
#'                   noise = "gumbel")
#' fit_g1 <- bmm(bmf(d ~ 1), afc, tcc_gumbel, cores = 4)
#' # Robinson et al.'s critical test: one d for all m, normal against Gumbel noise
#' loo::loo_compare(brms::loo(fit_d1), brms::loo(fit_g1))
#' sdm_afc <- sdm("err", n_options = 720, options = paste0("opt", 1:16))
#' fit_c1 <- bmm(bmf(c ~ 1, kappa ~ 1), afc, sdm_afc, cores = 4)
#' fit_cm <- bmm(bmf(c ~ 0 + m_condition, kappa ~ 1), afc, sdm_afc, cores = 4)
#' # within one model, does one parameter hold across m? Here the data come
#' # from one d, so the normal-noise model should need no d per m
#' loo::loo_compare(brms::loo(fit_d1), brms::loo(fit_dm))
#' loo::loo_compare(brms::loo(fit_c1), brms::loo(fit_cm))
tcc <- function(resp_error, similarity, n_options = 360, units = c("radians", "bins"),
                noise = c("gaussian", "gumbel"), d_max = 30, options = NULL,
                n_alternatives = NULL, ...) {
  call <- match.call()
  .check_constructor_dots(list(...), known = "links")
  stop_missing_args()
  units <- match.arg(units)
  noise <- match.arg(noise)
  .tcc_check_similarity(similarity, n_options)
  stopif(!is.numeric(d_max) || length(d_max) != 1 || !isTRUE(d_max > 0),
         "`d_max` must be a single positive number.")
  option_args <- if (!is.null(options) || !is.null(n_alternatives)) {
    .check_option_args(n_options, options, n_alternatives, units)
  }
  .model_tcc(resp_error = resp_error, similarity = similarity, n_options = n_options,
             units = units, noise = noise, d_max = d_max, options = options,
             n_alternatives = option_args$n_alternatives, call = call, ...)
}

# the full wheel prints nothing here, so that print() of tcc() without option
# sets is unchanged; the line is the one of sdm(options = )
#' @export
print_model_details.tcc <- function(model, ...) {
  options <- model$other_vars$options
  n_alternatives <- model$other_vars$n_alternatives
  if (is.null(options) && is.null(n_alternatives)) {
    return(invisible(NULL))
  }
  offered <- if (is.null(options)) {
    glue("{n_alternatives} equally spaced options")
  } else {
    "the options each trial offers"
  }
  cat(glue("Options:    {model$other_vars$n_options} positions on the wheel; \\
           each response is a choice among {offered}."), "\n")
  invisible(NULL)
}

.tcc_check_n_options <- function(n_options) {
  stopif(!is.numeric(n_options) || length(n_options) != 1 || !isTRUE(is.finite(n_options)) ||
           !isTRUE(n_options >= 4) || n_options %% 2 != 0,
         "`n_options` must be an even number of at least 4.")
}

.tcc_check_similarity <- function(similarity, n_options) {
  .tcc_check_n_options(n_options)
  n_values <- n_options / 2 + 1
  stopif(!is.numeric(similarity), "`similarity` must be a numeric vector.")
  stopif(length(similarity) != n_values,
         "`similarity` must hold {n_values} values, one for each distance from 0 to \\
         {n_options / 2} options, but it holds {length(similarity)}.")
  stopif(anyNA(similarity) || any(similarity < 0 | similarity > 1),
         "`similarity` must lie between 0 and 1.")
  stopif(similarity[1] != 1,
         "`similarity` must be 1 at distance 0 (the target), as d is the target's signal.")
}

#' @title A parametric similarity for [tcc()]
#' @description Returns the similarity vector that [tcc()] takes, from a von
#'   Mises or Laplace shape with a given `kappa`, for when no measured
#'   similarity is at hand.
#' @details The shapes are `exp(kappa * (cos(x) - 1))` for `"vonmises"`, the
#'   activation function of [sdm()], and `exp(-kappa * |x|)` for `"laplace"`,
#'   with `x` the distance from the target in radians and `kappa` the
#'   concentration or the rate (Oberauer, 2023). Each is scaled to 1
#'   at the target and 0 at the opposite option, as Schurgin et al. (2020)
#'   scaled their measured similarity.
#'
#'   `kappa` is fixed here, so a fit with this vector is still the TCC: the
#'   similarity enters the model as known. Schurgin et al. (2020) measured
#'   theirs in a separate perceptual task, and a measured similarity is the
#'   model as its authors intended it. To estimate `kappa` from the responses,
#'   fit [sdm()] instead, which estimates it for the von Mises shape.
#' @param shape `"vonmises"` or `"laplace"`.
#' @param kappa The concentration (von Mises) or rate (Laplace) of the shape, a
#'   single positive finite number, for distances in radians: a Laplace
#'   similarity that falls by a factor e over 20 degrees has `kappa = 1 /
#'   deg2rad(20)`.
#' @param n_options The number of response options on the wheel, an even
#'   number, as in [tcc()].
#' @return A numeric vector of `n_options / 2 + 1` similarities, at the
#'   distances 0, 1, ..., `n_options / 2` options from the target.
#' @references Oberauer, K. (2023). Measurement models for visual working
#'   memory - A factorial model comparison. *Psychological Review, 130*(3),
#'   841-852. https://doi.org/10.1037/rev0000328
#'
#'   Schurgin, M. W., Wixted, J. T., & Brady, T. F. (2020). Psychophysical
#'   scaling reveals a unified theory of visual memory strength. *Nature Human
#'   Behaviour, 4*(11), 1156-1172. https://doi.org/10.1038/s41562-020-00938-0
#' @keywords transform
#' @export
#' @examples
#' f <- tcc_similarity("laplace", kappa = 1 / deg2rad(20))
#' plot(0:180, f, type = "l", xlab = "Distance from the target (degrees)",
#'      ylab = "Similarity")
#' lines(0:180, tcc_similarity("vonmises", kappa = 8), lty = 2)
tcc_similarity <- function(shape = c("vonmises", "laplace"), kappa, n_options = 360) {
  stop_missing_args()
  shape <- match.arg(shape)
  stopif(!is.numeric(kappa) || length(kappa) != 1 || !isTRUE(is.finite(kappa) && kappa > 0),
         "`kappa` must be a single positive finite number.")
  .tcc_check_n_options(n_options)
  .tcc_shape_similarity(shape, kappa, seq(0, pi, length.out = n_options / 2 + 1))
}

# expm1() and 1 - cos(x) = 2 sin(x / 2)^2 keep the scaled shapes exact where
# they tend to a straight line (Laplace) or a raised cosine (von Mises), at a
# small kappa
.tcc_shape_similarity <- function(shape, kappa, delta) {
  if (shape == "laplace") {
    return(exp(-kappa * delta) * expm1(-kappa * (pi - delta)) / expm1(-kappa * pi))
  }
  s <- sin(delta / 2)^2
  exp(-2 * kappa * s) * expm1(-2 * kappa * (1 - s)) / expm1(-2 * kappa)
}

# a model stored before `noise` existed is the Gaussian cell
cell_args.tcc <- function(model) {
  if (is.null(model$other_vars$noise)) list() else list(noise = model$other_vars$noise)
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @exportS3Method
built_data_columns.tcc <- function(model) {
  c("tcc_bin", NextMethod("built_data_columns"))
}

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
  # with option sets, check_data.option_set() bins the response next
  if (!inherits(model, "option_set")) {
    warnif("tcc_bin" %in% colnames(data),
           "Column 'tcc_bin' in your data is reserved by the TCC model and will be overwritten.")
    data$tcc_bin <- .tcc_bin(data[[resp]], model$other_vars$n_options, units)
  }
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
  if (!inherits(model, "option_set")) {
    return(brms::bf("tcc_bin ~ 0"))
  }
  # the option offsets ride along as vint terms, so that log_lik() and
  # posterior_predict() read each trial's options from the data, newdata included
  brms::bf(glue("tcc_bin | vint({paste(.option_columns(.option_width(model)), collapse = ', ')}) ~ 0"))
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.tcc <- function(model, data, formula) {
  if (inherits(model, "option_set")) {
    return(.tcc_options_configure(model, data, formula))
  }
  formula <- bmf2bf(model, formula)
  gumbel <- .tcc_gumbel(model)
  formula$family <- brms::custom_family(
    if (gumbel) "tcc_gumbel" else "tcc",
    dpars = c("mu", "d"),
    links = c("identity", model$links$d),
    type = "int",
    loop = FALSE,
    log_lik = if (gumbel) log_lik_tcc_gumbel else log_lik_tcc,
    posterior_predict = if (gumbel) posterior_predict_tcc_gumbel else posterior_predict_tcc,
    posterior_epred = posterior_epred_tcc,
    vars = if (gumbel) "tcc_similarity" else c("tcc_tab", "tcc_grid_lo", "tcc_grid_step")
  )
  stanvars <- .tcc_stanvars(model, model$other_vars$table)
  nlist(formula, data, stanvars)
}

.tcc_gumbel <- function(model) {
  identical(model$other_vars$noise, "gumbel")
}

# One kernel, two sources for its table: built in transformed data from the
# similarity, or computed in R and passed as data under the same name
.tcc_stanvars <- function(model, table) {
  similarity <- .tcc_similarity(model)
  sc_path <- system.file("stan_chunks", package = "bmm")
  # Gumbel noise has a closed form and needs no table (#566)
  if (.tcc_gumbel(model)) {
    return(brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_gumbel_functions.stan")),
                         block = "functions") +
             brms::stanvar(x = similarity, name = "tcc_similarity"))
  }
  grid <- .tcc_grid(model)
  out <- brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_functions.stan")),
                       block = "functions") +
    brms::stanvar(x = grid$lo, name = "tcc_grid_lo") +
    brms::stanvar(x = grid$step, name = "tcc_grid_step")
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

############################################################################# !
# THE LIKELIHOOD TABLE                                                   ####
############################################################################# !

# The one place the model's similarity by distance is read, so that a
# parametric similarity function changes only this function. The similarity
# argument selects that function, never the formula (#566)
.tcc_similarity <- function(model) {
  model$other_vars$similarity
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
  step <- 0.02
  rows <- as.integer(ceiling(round((model$other_vars$d_max + step - lo) / step, 6)) + 3)
  nlist(lo, step, rows)
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
# tcc_logp_by_distance() in tcc_functions.stan is the Stan twin
.tcc_logp <- function(d, similarity, mult = .tcc_mult(length(similarity))) {
  means <- d * similarity
  top <- max(means)
  half <- (top + 8 - min(top - 8, (min(means) + top) / 2 - 6)) / 2
  rule <- .gl_rule(.tcc_gl_nodes)
  z <- outer(top + 8 - half + half * rule$nodes, means, "-")
  lphi <- stats::pnorm(z, log.p = TRUE)
  terms <- log(half * rule$weights) + drop(lphi %*% mult) + stats::dnorm(z, log = TRUE) - lphi
  matrixStats::colLogSumExps(terms)
}

# each distance but 0 and n / 2 has an option on either side of the target
.tcc_mult <- function(n_dist) {
  c(1, rep(2, n_dist - 2), 1)
}

# log P of one option at each distance under standard Gumbel noise. The largest
# of independent Gumbel(mu_x, 1) signals is Gumbel(log sum exp(mu), 1), and the
# argmax integral reduces to the softmax of the means (Yellott, 1977;
# Oberauer, 2023, eq. 4; derivation in the TCC design note, #566).
# tcc_gumbel_lpmf() in tcc_gumbel_functions.stan is the Stan twin
.tcc_gumbel_logp <- function(d, similarity, mult = .tcc_mult(length(similarity))) {
  a <- d * similarity
  a - matrixStats::logSumExp(log(mult) + a)
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
  if (.tcc_gumbel(model) || inherits(model, "option_set")) {
    return(invisible())
  }
  d_max <- model$other_vars$d_max
  .tcc_warn_beyond_grid(.tcc_dpar_extremes(prep, "d"), .tcc_grid(model)$lo, d_max, "d",
                        glue("Refit with d_max = {2 * d_max}."))
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
    Beyond it the table is extrapolated, which keeps the sampler \\
    moving but is not the model, so the posterior can run far past the table and \\
    its estimates mean nothing. {remedy}"
  )
}

############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT AND PP_CHECK                                ####
############################################################################# !

# brms calls these once per observation, so the table that a fit built in
# transformed data is rebuilt in R (about 1 s) once and kept for the next call
.tcc_cache <- new.env(parent = emptyenv())

# The table the fit read, from the standata in the prep: passed from R, or
# rebuilt from the similarity and grid the fit built it from
.tcc_prep_table <- function(prep) {
  data <- prep$data
  grid <- list(lo = data$tcc_grid_lo, step = data$tcc_grid_step, rows = data$tcc_grid_rows)
  if (!is.null(data$tcc_tab)) {
    return(c(grid, list(tab = data$tcc_tab)))
  }
  key <- rlang::hash(list(as.numeric(data$tcc_similarity), grid))
  if (!identical(.tcc_cache$key, key)) {
    .tcc_cache$tab <- .tcc_table(as.numeric(data$tcc_similarity), grid)
    .tcc_cache$key <- key
  }
  c(grid, list(tab = .tcc_cache$tab))
}

# Catmull-Rom rows and weights at each d, as tcc_first_row() and
# tcc_row_weights() compute them in Stan, linear beyond the grid
.tcc_weights <- function(d, table) {
  u <- (d - table$lo) / table$step + 2
  first <- pmin(pmax(floor(u) - 1, 1), nrow(table$tab) - 3)
  t <- u - first - 1
  w <- 0.5 * cbind(-t + 2 * t^2 - t^3, 2 - 5 * t^2 + 3 * t^3, t + 4 * t^2 - 3 * t^3, t^3 - t^2)
  below <- t < 0
  above <- t > 1
  w[below, ] <- cbind(-0.5 * t, 1, 0.5 * t, 0)[below, ]
  w[above, ] <- cbind(0, -0.5 * (t - 1), 1, 0.5 * (t - 1))[above, ]
  nlist(first, w)
}

# log P of distance index `col` (1 = the target) at each d, the R twin of tcc_lpmf()
.tcc_interp <- function(d, col, table) {
  at <- .tcc_weights(d, table)
  col <- rep_len(col, length(d))
  Reduce(`+`, lapply(1:4, function(k) at$w[, k] * table$tab[cbind(at$first + k - 1, col)]))
}

# signed option offsets drawn at each d from the interpolated table, by
# inverse CDF over all options at once, in chunks to bound the memory
.tcc_rbins <- function(d, table, n_options, chunk = 5000) {
  out <- integer(length(d))
  for (start in seq(1, length(d), by = chunk)) {
    idx <- start:min(length(d), start + chunk - 1)
    at <- .tcc_weights(d[idx], table)
    logp <- Reduce(`+`, lapply(1:4, function(k) at$w[, k] * table$tab[at$first + k - 1, , drop = FALSE]))
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

# log P of one option at each distance under Gumbel noise, one row per value of d
.tcc_gumbel_logp_rows <- function(d, similarity) {
  a <- outer(d, similarity)
  a - matrixStats::rowLogSumExps(a + rep(log(.tcc_mult(ncol(a))), each = nrow(a)))
}

# brms calls log_lik and posterior_predict once per observation, and the
# observations of one cell share their draws, so the last rows are kept
.tcc_gumbel_logp_draws <- function(prep, i) {
  d <- brms::get_dpar(prep, "d", i = i)
  similarity <- as.numeric(prep$data$tcc_similarity)
  key <- list(d, similarity)
  if (!identical(.tcc_cache$gumbel_key, key)) {
    .tcc_cache$gumbel_logp <- .tcc_gumbel_logp_rows(d, similarity)
    .tcc_cache$gumbel_key <- key
  }
  .tcc_cache$gumbel_logp
}

.tcc_gumbel_rbins <- function(d, similarity, n_options, chunk = 5000) {
  out <- integer(length(d))
  for (start in seq(1, length(d), by = chunk)) {
    idx <- start:min(length(d), start + chunk - 1)
    out[idx] <- .tcc_draw_offsets(.tcc_gumbel_logp_rows(d[idx], similarity), n_options)
  }
  out
}

log_lik_tcc <- function(i, prep) {
  d <- brms::get_dpar(prep, "d", i = i)
  .tcc_interp(d, abs(prep$data$Y[i]) + 1, .tcc_prep_table(prep))
}

posterior_predict_tcc <- function(i, prep, ...) {
  table <- .tcc_prep_table(prep)
  .tcc_rbins(brms::get_dpar(prep, "d", i = i), table, 2 * (ncol(table$tab) - 1))
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
  if (inherits(model, "option_set")) {
    draws <- vapply(seq_len(prep$nobs), posterior_predict_tcc_options, integer(prep$ndraws), prep = prep)
    return(list(error = matrix(draws, nrow = prep$ndraws)))
  }
  d <- .pp_dpar_vector(prep, "d")
  n_options <- model$other_vars$n_options
  bins <- if (.tcc_gumbel(model)) {
    .tcc_gumbel_rbins(d, as.numeric(prep$data$tcc_similarity), n_options)
  } else {
    .tcc_rbins(d, .tcc_prep_table(prep), n_options)
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

############################################################################# !
# FORCED-CHOICE OPTION SETS (#566, D4)                                   ####
############################################################################# !

# Each trial offers some of the n_options positions (the option-set contract,
# R/helpers-options.R), and the response is the offered option with the
# strongest signal. The probabilities depend on d and on the offered set, so
# the rows that share both share them: a cell is a run of equal predictors
# times an option set. Cells are numbered by first appearance, and each row
# carries its cell, from which transformed data counts the choices per offered
# option of each cell; the likelihood then evaluates each cell once
.tcc_option_cells <- function(model, data, formula) {
  predictors <- data_predictor_vars(data, formula)
  # brms drops the rows with a missing value in a model variable before fitting
  used <- intersect(c(model$resp_vars$resp_error, predictors), colnames(data))
  kept <- data[stats::complete.cases(data[used]), , drop = FALSE]
  # a new run wherever any predictor changes value; interaction() would compare
  # pasted labels, which coincide for ("a.b", "c") and ("a", "b.c")
  changed <- Reduce(`|`, lapply(kept[predictors], function(x) x[-1] != x[-length(x)]),
                    logical(max(nrow(kept) - 1, 0)))
  run_id <- cumsum(c(TRUE, changed))[seq_len(nrow(kept))]
  sets <- .option_set_table(as.matrix(kept[.option_columns(.option_width(model))]),
                            model$other_vars$n_options)
  key <- paste(run_id, sets$set_id)
  list(row_cell = match(key, unique(key)), cell_set = sets$set_id[!duplicated(key)],
       size = sets$size, off = sets$off)
}

# data for the sliced likelihood of a threaded fit needs its own argument
.tcc_data_stanvar <- function(x, name, arg, scode = NULL) {
  out <- brms::stanvar(x = x, name = name, scode = scode)
  out[[name]]$pll_args <- paste("data", arg, name)
  out
}

# brms declares an integer array of length 1 as a scalar and an integer matrix
# as a one-dimensional array, so the cell tables get their declarations here
.tcc_int_array_stanvar <- function(x, name, dims, arg = NULL) {
  x <- array(as.integer(x), dim = if (is.matrix(x)) dim(x) else length(x))
  scode <- paste0("array[", paste(dims, collapse = ", "), "] int ", name, ";")
  if (is.null(arg)) {
    return(brms::stanvar(x = x, name = name, scode = scode))
  }
  .tcc_data_stanvar(x, name, arg, scode = scode)
}

.tcc_options_configure <- function(model, data, formula) {
  cells <- .tcc_option_cells(model, data, formula)
  n_rows <- length(cells$row_cell)
  n_cells <- length(cells$cell_set)
  # one integral over its options per cell and gradient: 0.07 ms with 2
  # options to 0.18 ms with 16 (#566, D4; ?tcc, "Forced-choice option sets")
  if (!.tcc_gumbel(model) && n_cells > n_rows / 2) {
    message2(
      "The {n_rows} rows fall into {n_cells} cells of equal d and offered options. \\
      Under normal noise each cell takes one integral over its options per gradient, \\
      so the fit takes longer than when rows share cells; see ?tcc, section \\
      \"Forced-choice option sets\"."
    )
  }
  formula <- bmf2bf(model, formula)
  formula$family <- brms::custom_family(
    "tcc_options",
    dpars = c("mu", "d"),
    links = c("identity", model$links$d),
    type = "int",
    loop = FALSE,
    log_lik = log_lik_tcc_options,
    posterior_predict = posterior_predict_tcc_options,
    posterior_epred = posterior_epred_tcc
  )
  sc_path <- system.file("stan_chunks", package = "bmm")
  chunk <- function(file) read_lines2(file.path(sc_path, file))
  rule <- .gl_rule(.tcc_gl_nodes)
  stanvars <- brms::stanvar(scode = chunk("tcc_options_functions.stan"), block = "functions") +
    .tcc_data_stanvar(match(model$other_vars$noise, c("gaussian", "gumbel")), "tcc_noise", "int") +
    .tcc_data_stanvar(.tcc_similarity(model), "tcc_similarity", "vector") +
    .tcc_data_stanvar(rule$nodes, "tcc_gl_nodes", "vector") +
    .tcc_data_stanvar(log(rule$weights), "tcc_gl_logw", "vector") +
    brms::stanvar(x = n_cells, name = "G_tcc_cells") +
    .tcc_int_array_stanvar(cells$row_cell, "tcc_row_cell", "N") +
    .tcc_int_array_stanvar(cells$cell_set, "tcc_cell_set", "G_tcc_cells", "array[] int") +
    brms::stanvar(x = nrow(cells$off), name = "S_tcc_sets") +
    brms::stanvar(x = ncol(cells$off), name = "M_tcc_sets") +
    .tcc_int_array_stanvar(cells$size, "tcc_set_size", "S_tcc_sets", "array[] int") +
    .tcc_int_array_stanvar(cells$off, "tcc_set_off", c("S_tcc_sets", "M_tcc_sets"), "array[,] int") +
    brms::stanvar(scode = chunk("tcc_options_tdata.stan"), block = "tdata",
                  pll_args = "data array[] int tcc_cell_first, data array[,] int tcc_cell_n") +
    brms::stanvar(
      scode = chunk(if (brms_slices_likelihood()) "tcc_options_likelihood_threaded.stan"
                    else "tcc_options_likelihood.stan"),
      block = "likelihood", position = "end"
    )
  nlist(formula, data, stanvars)
}

# The offered options of observation i as offsets in steps: the vint columns
# of the brms data, without the padding
.tcc_prep_options <- function(prep, i) {
  columns <- grep("^vint[0-9]+$", names(prep$data), value = TRUE)
  columns <- columns[order(as.integer(sub("vint", "", columns)))]
  offsets <- vapply(columns, function(column) prep$data[[column]][i], integer(1), USE.NAMES = FALSE)
  offsets[offsets != .option_absent(2 * (length(prep$data$tcc_similarity) - 1))]
}

# log P of each offered option, whose similarities to the target are f, at one
# value of d; the Stan twin is tcc_set_logp() in tcc_options_functions.stan
.tcc_set_logp <- function(d, f, noise) {
  if (noise == "gumbel") {
    return(.tcc_gumbel_logp(d, f, mult = rep(1, length(f))))
  }
  .tcc_logp(d, f, mult = rep(1, length(f)))
}

# one row per value of d, one column per offered option
.tcc_set_logp_rows <- function(d, f, noise) {
  values <- unique(d)
  rows <- matrix(vapply(values, .tcc_set_logp, numeric(length(f)), f = f, noise = noise),
                 ncol = length(f), byrow = TRUE)
  rows[match(d, values), , drop = FALSE]
}

# brms calls log_lik and posterior_predict once per observation, and the rows
# of one cell share their draws and options, so the last rows are kept
.tcc_options_logp_draws <- function(prep, i, offsets) {
  d <- brms::get_dpar(prep, "d", i = i)
  f <- as.numeric(prep$data$tcc_similarity)[abs(offsets) + 1]
  noise <- c("gaussian", "gumbel")[prep$data$tcc_noise]
  key <- list(d, f, noise)
  if (!identical(.tcc_cache$options_key, key)) {
    .tcc_cache$options_logp <- .tcc_set_logp_rows(d, f, noise)
    .tcc_cache$options_key <- key
  }
  .tcc_cache$options_logp
}

log_lik_tcc_options <- function(i, prep) {
  offsets <- .tcc_prep_options(prep, i)
  y <- prep$data$Y[i]
  # newdata that carries the built columns is not checked again, so a
  # response the trial did not offer can reach here; it has no probability
  if (!is.na(y) && !y %in% offsets) {
    return(rep(-Inf, prep$ndraws))
  }
  .tcc_options_logp_draws(prep, i, offsets)[, match(y, offsets)]
}

posterior_predict_tcc_options <- function(i, prep, ...) {
  offsets <- .tcc_prep_options(prep, i)
  .tcc_draw_from_set(.tcc_options_logp_draws(prep, i, offsets), offsets)
}

# one offset per row of logp, which holds log P of each offered option
.tcc_draw_from_set <- function(logp, offsets) {
  cum <- matrixStats::rowCumsums(exp(logp))
  draw <- stats::runif(nrow(logp)) * cum[, ncol(cum)]
  as.integer(offsets[rowSums(cum < draw) + 1])
}
