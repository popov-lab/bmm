############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_yn <- function(response = NULL, stimulus = NULL,
                          n_trials = NULL, dist = "normal",
                          links = NULL, call = NULL, ...) {
  parameters <- list(
    d = paste0(
      "d_a sensitivity (= d' when sdratio is fixed): the distance between ",
      "the signal and noise distributions in units of their root-mean-square SD"
    ),
    criterion = "Response bias: location of decision boundary",
    sdratio = paste0(
      "SD ratio signal/noise, log link (0 = equal SDs): summary() prints the ",
      "log of this ratio, not the ratio itself, so exp() a posterior value ",
      "to read it"
    )
  )
  # d is d_a, and the noise-standardized separation is d * sqrt((1 + r^2) / 2).
  # Over plausible SD ratios, r in [0.56, 1.80], that factor runs from
  # 0.81 to 1.46, which normal(1, 1) is wide enough to absorb on either scale.
  # sd rates track between-subject SDs fitted on broeder_schuetz_2009_e3 (40
  # subjects) with the wide brms sd default in force: d ~0.6, criterion ~0.15,
  # sdratio ~0.2, against the rule that the prior median covers the typical SD
  # and its 95% quantile twice the largest. Rate 2's median (0.35) sits below
  # d's SD, so d takes rate 1. criterion's 0.15 is deflated by one subject
  # intercept shared across five conditions -- a single condition alone gives
  # 0.24, above rate 4's median (0.17) -- so criterion takes rate 2. sdratio's
  # 0.2 clears that same median on its own, so it takes rate 2 too.
  default_priors <- list(
    d = list(main = "normal(1, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
    criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)", sd = "exponential(2)"),
    # sdratio is on the log scale, so normal(0, 0.3) puts 95% of the prior on
    # ratios in [0.56, 1.80] and next to nothing above 2, which recognition
    # does not produce. Both group-level estimates in hand sit inside that:
    # Mickes et al. (2007) Table 1 averages sd(lure)/sd(target) = 0.79 over
    # their 13 retained subjects, i.e. 1.26 signal over noise, and the
    # broeder_schuetz_2009_e3 fit puts the group intercept at 1.46 [1.25, 1.71].
    # Their subject-to-subject spread is the sd prior's business, not this
    # one's. 0.3 is 3.7 times that fit's posterior SD of 0.081, so the prior
    # regularizes without standing in for the data. effects keeps the
    # intercept's scale because nothing read here bounds how far a condition
    # moves the ratio -- that fit estimated sdratio ~ 1.
    sdratio = list(main = "normal(0, 0.3)", effects = "normal(0, 0.3)", sd = "exponential(2)")
  )
  requirements <- glue(
    "Provide pre-aggregated data with the following columns:", "\n\n",
    "  - Response counts (response): number of 'old'/'signal' responses", "\n",
    "  - Stimulus type (stimulus): 0 = noise, 1 = signal", "\n",
    "  - Number of trials (n_trials): total trials per cell"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(stimulus, n_trials, dist),
      domain = "Perception & Recognition Memory",
      task = "Yes/No Detection or Old/New Recognition",
      name = "Signal Detection Theory (Yes/No)",
      citation = glue(
        "Green, D. M., & Swets, J. A. (1966). Signal detection theory ",
        "and psychophysics. Wiley."
      ),
      version = "NA",
      requirements = requirements,
      parameters = parameters,
      links = list(d = "identity", criterion = "identity",
                   sdratio = "log"),
      # on the link scale, so exp(0) = 1 = equal variance
      fixed_parameters = list(sdratio = 0),
      default_priors = default_priors,
      # on the natural scale; create_initfun() applies the forward link
      init_ranges = list(
        d         = c(0.5, 1.5),
        criterion = c(-0.5, 0.5),
        sdratio   = c(0.75, 1.35)
      )
    ),
    class = c("bmmodel", "sdt", "sdt_yn"),
    call = call
  )
  out <- set_links(out, links)
  out
}

# `sdratio` is fixed at 0, and that 0 is read on whatever link the model
# carries: only log maps it to the equal-variance ratio of 1. Every other link
# bmm offers turns it into a different ratio -- identity, sqrt, log1p and
# tan_half into 0, which makes eta infinite on every signal row; softplus,
# logm1, logit, probit, cloglog into a finite ratio the user never asked for,
# which samples to completion and is wrong; and inverse into an infinite
# ratio, which makes eta NaN (Inf/Inf) rather than Inf -- worse than the other
# two failure modes, because a NaN poisons the adjoint of every shared
# parameter instead of just killing one row's likelihood. `d` and `criterion`
# fix nothing, so their links stay settable.
#' @exportS3Method
settable_links.sdt_yn <- function(model) {
  c("d", "criterion")
}


#' @title Yes/No Signal Detection Theory Model
#' @name sdt_yn
#' @description Estimates sensitivity (`d`, which is \eqn{d'} under equal
#'   variance and \eqn{d_a} once `sdratio` is estimated) and response bias
#'   (`criterion`) from yes/no detection or old/new recognition counts.
#' @details `r model_info(.model_sdt_yn())`
#' @param response The name of the variable in the dataset containing the
#'   count of "old"/"signal" responses for each cell.
#' @param stimulus The name of the variable in the dataset coding the stimulus
#'   type: 0 (noise/new) and 1 (signal/old). Logical, and factor or character
#'   columns holding "0"/"1", are coerced automatically; anything else (e.g.
#'   "noise"/"signal" labels) must be recoded by hand, since bmm cannot guess
#'   which level is the signal. [dsdt_yn()] and [rsdt_yn()] take the same
#'   column but are stricter, accepting only numeric or logical input.
#' @param n_trials The name of the variable in the dataset containing the
#'   total number of trials for each cell. It may differ from cell to cell.
#'   Trial-level data also works: keep one row per trial, with the column
#'   named by `response` holding 0 or 1 and the column named by `n_trials`
#'   holding `1`s.
#' @param dist The noise distribution assumed for the latent evidence variable,
#'   given here by its cumulative distribution function. One of:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT, \eqn{\Phi(x)}
#'     \item "gumbel_min": smallest-extreme-value SDT,
#'       \eqn{1 - \exp(-\exp(x))}{1 - exp(-exp(x))} (complementary log-log)
#'     \item "gumbel_max": largest-extreme-value SDT,
#'       \eqn{\exp(-\exp(-x))}{exp(-exp(-x))}
#'       (log-log, as in \code{evd::pgumbel})
#'     \item "logistic": logistic SDT,
#'       \eqn{1 / (1 + \exp(-x))}{1 / (1 + exp(-x))}
#'   }
#' @param links A named list of link functions for the parameters, one entry
#'   per parameter you want to change, e.g. `links = list(d = "log")`. Only `d`
#'   and `criterion` can be set. `sdratio` keeps its log link, because the
#'   model's default of an equal-variance SD ratio is stored as the 0 that the
#'   log link maps to 1; read on any other link that same 0 is a different
#'   ratio, and on an identity link it is a ratio of zero that makes the signal
#'   trials' evidence infinite.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#'
#' @section Which sensitivity measure `d` is:
#' **`d` is the familiar \eqn{d'}** whenever `sdratio` keeps its default of 0
#' (an SD ratio of 1, equal variance), which is every fit that does not give
#' `sdratio` a formula. The rest of this section only matters once you estimate
#' `sdratio`.
#'
#' When the signal and noise distributions have different widths, their
#' separation only becomes dimensionless after choosing an SD to divide by.
#' `bmm` then reports \eqn{d_a}, the separation divided by the root-mean-square
#' of the two SDs:
#' \deqn{d_a = \sqrt{2}\,\delta / \sqrt{1 + r^2},}{d_a = sqrt(2) * delta / sqrt(1 + r^2),}
#' where \eqn{\delta} is the separation in noise-SD units and
#' \eqn{r} is the SD ratio, `exp(sdratio)`. This weights the two
#' distributions equally, and it is the measure Simpson and Fitter (1973),
#' Macmillan and Creelman (2005), and Mickes et al. (2007) recommend under
#' unequal variance. The classical noise-standardized index is
#' \eqn{d_N = \delta = d_a \sqrt{(1 + r^2)/2}}{d_N = delta = d_a * sqrt((1 + r^2)/2)}, so a published \eqn{d'} from an
#' unequal-variance analysis is larger than `d` when \eqn{r > 1}: by 13% at
#' \eqn{r = 1.25} and by 33% at \eqn{r = 1.6}.
#'
#' `d` is \eqn{d_a} rather than \eqn{d_N} because only \eqn{d_a} is comparable
#' across conditions or subjects that differ in `sdratio`: two conditions that
#' are equally discriminable can show a large, confidently estimated difference
#' in \eqn{d_N}. When `sdratio` is estimated but constant across the conditions
#' you compare, the two indices differ by one common factor and give the same
#' contrasts up to scale.
#'
#' **Units of the other parameters.** `criterion` is *not* rescaled. It is the
#' location of the decision boundary relative to the midpoint between the two
#' distributions, in noise-SD units, so under unequal variance `d` (in
#' root-mean-square SD units) and `criterion` (in noise-SD units) are on
#' different scales, and a ratio such as `criterion / d` mixes them.
#'
#' **Extreme-value distributions.** For `dist = "normal"`, \eqn{d_a} is also the
#' AUC-equivalent index, \eqn{d_a = \sqrt{2}\,\Phi^{-1}(\mathrm{AUC})}{d_a = sqrt(2) * Phi^-1(AUC)}, so it
#' carries the same information as 2AFC accuracy. For `"gumbel_min"` and
#' `"gumbel_max"` that identity holds only under equal variance. With `sdratio`
#' estimated, \eqn{d_a} keeps its balanced geometry but drifts away from the
#' AUC-equivalent index as the SD ratio moves away from 1 (for `"gumbel_min"`
#' at \eqn{\delta = 1.5} and \eqn{r = 2}, \eqn{d_a} is 31% larger), so compare
#' such fits on the AUC rather than on `d`.
#'
#' Because `d` is a short name, a column called `d` in your data that is also
#' used as a predictor will collide with this parameter; `bmm()` warns when that
#' happens.
#'
#' @section Identifying `sdratio`:
#' `sdratio` needs a design that supplies more than one operating point. A single
#' (hit, false-alarm) pair is two numbers for three unknowns, so when every
#' parameter is intercept-only with no random effects, `sdratio ~ 1` returns its
#' prior and `d` is pulled along the resulting ridge: sampling converges, `Rhat`
#' is fine, and the profile likelihood over `sdratio` is flat to 1e-12. `bmm()`
#' warns for that one design, which is the only shape that provably cannot work.
#'
#' Any linear predictor that moves the operating point along the ROC supplies
#' what is missing, and it need not sit on `criterion`: a sensitivity
#' manipulation (`d ~ 0 + condition`, a study-time or strength manipulation with
#' bias held constant) identifies `sdratio` just as a criterion manipulation
#' does, and so does between-subject variation entering through a random effect
#' such as `criterion ~ 1 + (1 | id)` — though, like any predictor, only in
#' proportion to how far it actually moves the operating point: a random
#' effect with little between-subject spread carries little information and
#' will not trigger the warning above, because that warning counts formula
#' terms, not how much they move the design.
#'
#' A criterion manipulation is still the cleanest design, because it traces the
#' ROC at fixed sensitivity: give `criterion` a predictor that shifts the
#' decision boundary — a base-rate, payoff, or confidence manipulation — as in
#' `criterion ~ 0 + condition`; see [broeder_schuetz_2009_e3]. Leaving `sdratio`
#' at its default is always identified.
#'
#' @section Reading `sdratio` and carrying it to [dsdt_yn()]/[rsdt_yn()]:
#' As in every `bmm` model, the parameters the model *estimates* are on their
#' link scale, while the distribution functions take their arguments on the
#' *natural* scale. `sdratio` has a log link, so `summary()` reports
#' \eqn{\log r}{log(r)} whereas [dsdt_yn()] and [rsdt_yn()] expect the ratio
#' \eqn{r} itself (their default is 1, equal variance). Exponentiate before
#' carrying a posterior value across:
#'
#' ```
#' r <- exp(as_draws_matrix(fit)[, "b_sdratio_Intercept"])
#' ```
#'
#' A posterior mean of `sdratio = 0.22` is a ratio of `exp(0.22) = 1.25`.
#' Passing `0.22` straight to `rsdt_yn()` instead asks for a signal
#' distribution 4.5 times *narrower* than the noise — a legal value that
#' raises no error, and the one mistake worth checking for in a posterior
#' predictive check written by hand. `d` and `criterion` have identity links,
#' so they carry across unchanged.
#'
#' The **zROC slope** reported in the recognition-memory literature is the
#' reciprocal of that ratio, `1 / exp(sdratio)`, so a `sdratio` posterior mean
#' of 0.375 is a zROC slope of 0.69.
#'
#' The same log link applies going the other way: a constant you supply
#' yourself, whether as `bmf(sdratio = )` or through a hand-written
#' `brms::set_prior(..., dpar = "sdratio")`, is read on it too.
#' `bmf(sdratio = 1)` does not fix a ratio of 1 — it fixes `exp(1) = 2.72`,
#' and `bmm()` raises no warning; a fixed ratio of 1.25 needs
#' `bmf(sdratio = log(1.25))`. Both `sdratio` defaults `default_prior()`
#' reports — `normal(0, 0.3)` on the intercept and `exponential(2)` on the
#' random-effect SDs — are on that same scale, unannotated.
#'
#' @section Terms used on this page:
#' \itemize{
#'   \item **noise-standardized axis** — latent evidence expressed in units of
#'     the noise distribution's SD. `criterion` always lives on this axis.
#'   \item \eqn{\delta} (**separation**) — the distance between the signal and
#'     noise means, on that axis.
#'   \item \eqn{d_N} — the classical \eqn{d'}: \eqn{\delta} itself, i.e. the
#'     separation in noise-SD units.
#'   \item \eqn{d_a} — the separation divided by the root-mean-square of the
#'     two SDs, which is what `d` reports. \eqn{d_a = d_N} under equal
#'     variance.
#'   \item **operating point** — one (false-alarm rate, hit rate) pair, i.e.
#'     one point of an ROC curve. One condition gives one operating point.
#'   \item **AUC** — the area under that ROC curve, equivalently the
#'     probability that a random signal trial yields more evidence than a
#'     random noise trial. Obtain it from the posterior with
#'     `pnorm(d / sqrt(2))` for `dist = "normal"`.
#'   \item `main` / `effects` / `sd` — the keys of the default priors shown in
#'     the model description above: `main` is the prior on the intercept,
#'     `effects` the prior on regression coefficients, and `sd` the prior on
#'     the standard deviations of the parameter's random effects.
#' }
#'
#' @seealso [sdt_d()] and [sdt_criterion()] compute the `d` and `criterion` of
#'   this model in closed form from a single pair of observed hit and false-alarm
#'   rates, without fitting: use them for a quick check of a fitted value, and
#'   this model when you need a hierarchical or condition-wise estimate, or
#'   unequal variance.
#'
#' @references
#' Green, D. M., & Swets, J. A. (1966). \emph{Signal detection theory and
#'   psychophysics}. Wiley.
#'
#' DeCarlo, L. T. (1998). Signal detection theory and generalized linear
#'   models. \emph{Psychological Methods}, \emph{3}(2), 186--205.
#'   \doi{10.1037/1082-989X.3.2.186}
#'
#' Simpson, A. J., & Fitter, M. J. (1973). What is the best index of
#'   detectability? \emph{Psychological Bulletin}, \emph{80}(6), 481--488.
#'   \doi{10.1037/h0035203}
#'
#' Macmillan, N. A., & Creelman, C. D. (2005). \emph{Detection theory: A user's
#'   guide} (2nd ed.). Erlbaum.
#'
#' Mickes, L., Wixted, J. T., & Wais, P. E. (2007). A direct test of the
#'   unequal-variance signal detection model of recognition memory.
#'   \emph{Psychonomic Bulletin & Review}, \emph{14}(5), 858--865.
#'   \doi{10.3758/BF03194112}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
#' dat$n_trials <- 100L
#' dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus,
#'                      d = 1.5, criterion = 0.2)
#'
#' model <- sdt_yn(
#'   response = "n_old",
#'   stimulus = "stimulus",
#'   n_trials = "n_trials"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(d ~ 1, criterion ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Sensitivity and bias per participant
#' fit_re <- bmm(
#'   formula = bmf(d ~ 1 + (1 | id), criterion ~ 1 + (1 | id)),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Unequal-variance yes/no SDT. sdratio needs more than one operating point:
#' # on the single-condition `dat` above it would not be identified.
#' # `model` already names the columns this dataset uses.
#' fit_uv <- bmm(
#'   formula = bmf(d ~ 1, criterion ~ 0 + condition, sdratio ~ 1),
#'   data = broeder_schuetz_2009_e3,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Simulating from a fitted unequal-variance model: the model reports
#' # sdratio on its log link, rsdt_yn() takes the ratio itself.
#' sdratio_posterior <- 0.375
#' rsdt_yn(2, 100L, c(0L, 1L),
#'         d = 1.3, criterion = 0.1, sdratio = exp(sdratio_posterior))
#' }
sdt_yn <- function(response, stimulus, n_trials,
                   dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
                   links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)

  .model_sdt_yn(response = response, stimulus = stimulus,
                n_trials = n_trials, dist = dist,
                links = links, call = call, ...)
}


############################################################################# !
# CHECK_FORMULA S3 METHODS                                               ####
############################################################################# !

# Per-format, not shared: sdt_rating and sdt_cdp carry n_ratings - 1 thresholds,
# so one condition already gives 2(n_ratings - 1) cumulative statistics against
# n_ratings + 1 parameters and sdratio is identified within it -- the
# zROC-slope-from-confidence-ratings logic those formats exist for. A shared rule
# would warn on every legitimate single-condition rating fit.
#' @export
check_formula.sdt_yn <- function(model, data, formula) {
  # rhs_vars() reports random-effect grouping variables and the `Intercept`
  # pseudo-predictor, so `1 + (1 | id)` and `0 + Intercept` are both visible
  flat_design <- all(vapply(formula, function(x) {
    all(rhs_vars(x) %in% "Intercept")
  }, logical(1)))
  warnif(
    flat_design && !is_constant(formula)[["sdratio"]],
    "'sdratio' is estimated, but every parameter formula is intercept-only with \\
    no random effects. Such a design has two sufficient statistics (one hit \\
    rate, one false-alarm rate) for three parameters, so 'sdratio' will return \\
    its prior and 'd' will slide along the resulting ridge while sampling \\
    converges and Rhat stays fine. Give a parameter a predictor or a random \\
    effect, or leave 'sdratio' at its default. See the 'Identifying sdratio' \\
    section of ?sdt_yn"
  )
  NextMethod("check_formula")
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_yn <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  data[[stim_var]] <- .validate_sdt_stimulus(data, stim_var)

  .validate_sdt_counts(data, model$resp_vars$response,
                       model$other_vars$n_trials)

  warnif("dist_type" %in% colnames(data),
         "Column 'dist_type' in your data is reserved by {model$name} and \\
         will be overwritten with the code for dist = '{model$other_vars$dist}'")
  data$dist_type <- .sdt_dist_id(model$other_vars$dist)

  NextMethod("check_data")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt_yn <- function(model, formula) {
  resp_var <- model$resp_vars$response
  stim_var <- model$other_vars$stimulus
  n_trials_var <- model$other_vars$n_trials

  brms::bf(paste0(
    resp_var, " | vint(", stim_var, ", dist_type) + trials(",
    n_trials_var, ") ~ 0"
  ))
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_yn <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  formula$family <- brms::custom_family(
    "sdt_yn",
    dpars = c("mu", "d", "criterion", "sdratio"),
    links = c("identity", model$links$d, model$links$criterion,
              model$links$sdratio),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_sdt_yn,
    posterior_predict = posterior_predict_sdt_yn,
    vars = c("vint1[n]", "vint2[n]", "trials[n]")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- paste(
    read_lines2(paste0(sc_path, "/sdt_dist_funs.stan")),
    read_lines2(paste0(sc_path, "/sdt_yn_funs.stan")),
    sep = "\n"
  )
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

log_lik_sdt_yn <- function(i, prep) {
  d <- brms::get_dpar(prep, "d", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  sdratio <- brms::get_dpar(prep, "sdratio", i = i)
  dist <- .sdt_dist_names[prep$data$vint2[i]]

  dsdt_yn(prep$data$Y[i], prep$data$trials[i], prep$data$vint1[i],
              d, criterion, sdratio = sdratio, dist = dist, log = TRUE)
}

posterior_predict_sdt_yn <- function(i, prep, ...) {
  d <- brms::get_dpar(prep, "d", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  sdratio <- brms::get_dpar(prep, "sdratio", i = i)
  dist <- .sdt_dist_names[prep$data$vint2[i]]

  rsdt_yn(length(d), prep$data$trials[i], prep$data$vint1[i],
              d, criterion, sdratio = sdratio, dist = dist)
}
