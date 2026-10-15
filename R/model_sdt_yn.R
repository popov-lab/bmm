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
    criterion = paste0(
      "Decision criterion, measured from the mean of the noise distribution in ",
      "noise-SD units; the false-alarm rate depends on it alone"
    ),
    bias = paste0(
      "Alternative to criterion: the decision criterion measured from the ",
      "midpoint between the noise and signal means (0 = neutral, positive = ",
      "conservative). Use one of the two"
    ),
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
  # bias is the midpoint-relative criterion those fits estimated. criterion
  # adds half the separation, so it inherits half of d's spread: closed-form
  # per-subject values on the same data spread 0.42 against bias's 0.17, which
  # puts criterion at rate 1. Its intercept centres on 0.5, the unbiased
  # observer at d's prior mean of 1; condition means there run 0.28 to 1.47.
  default_priors <- list(
    d = list(main = "normal(1, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
    criterion = list(main = "normal(0.5, 1.5)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
    bias = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)", sd = "exponential(2)"),
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
      links = list(d = "identity", criterion = "identity", bias = "identity",
                   sdratio = "log"),
      # on the link scale, so exp(0) = 1 = equal variance
      fixed_parameters = list(sdratio = 0),
      alternative_parameters = list(c("criterion", "bias")),
      default_priors = default_priors,
      # on the natural scale; create_initfun() applies the forward link
      init_ranges = list(
        d         = c(0.5, 1.5),
        criterion = c(0, 1),
        bias      = c(-0.5, 0.5),
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
# parameter instead of just killing one row's likelihood. `d`, `criterion` and
# `bias` fix nothing, so their links stay settable.
#' @exportS3Method
settable_links.sdt_yn <- function(model) {
  c("d", "criterion", "bias")
}


#' @title Yes/No Signal Detection Theory Model
#' @name sdt_yn
#' @description Estimates sensitivity (`d`, which is \eqn{d'} under equal
#'   variance and \eqn{d_a} once `sdratio` is estimated) and the decision
#'   criterion from yes/no detection or old/new recognition counts. The
#'   criterion is estimated either as `criterion`, measured from the mean of the
#'   noise distribution, or as `bias`, measured from the midpoint between the
#'   noise and signal means; the formula names the one to estimate.
#' @details `r model_docs(.model_sdt_yn())`
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
#'   per parameter you want to change, e.g. `links = list(d = "log")`. Only `d`,
#'   `criterion` and `bias` can be set. `sdratio` keeps its log link, because the
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
#' \eqn{r = 1.25} and by 33% at \eqn{r = 1.6}. [sdt_sensitivity()] returns
#' \eqn{d_N} and the signal-standardized \eqn{d_S} from a fitted model, and
#' `summary()` notes when `d` is \eqn{d_a}.
#'
#' `d` is \eqn{d_a} rather than \eqn{d_N} because only \eqn{d_a} is comparable
#' across conditions or subjects that differ in `sdratio`: two conditions that
#' are equally discriminable can show a large, confidently estimated difference
#' in \eqn{d_N}. When `sdratio` is estimated but constant across the conditions
#' you compare, the two indices differ by one common factor and give the same
#' contrasts up to scale.
#'
#' **Units of the other parameters.** `criterion` and `bias` are *not*
#' rescaled: both are in noise-SD units, so under unequal variance `d` (in
#' root-mean-square SD units) and the criterion are on different scales, and a
#' ratio such as `bias / d` mixes them.
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
#' @section Choosing between `criterion` and `bias`:
#' Both parameters place the same decision boundary on the noise-standardized
#' evidence axis, and they differ only in the point they measure it from:
#' \itemize{
#'   \item `criterion` is the distance from the mean of the noise distribution.
#'     The false-alarm rate depends on it alone: \eqn{FA = 1 - F(criterion)}.
#'   \item `bias` is the distance from the midpoint between the noise and signal
#'     means: 0 is the neutral point, positive values a conservative bias toward
#'     "new"/"noise". Under equal variance it is the index usually written
#'     \eqn{c}.
#' }
#' They are related by \eqn{bias = criterion - \delta / 2}{bias = criterion -
#' delta / 2}, where \eqn{\delta = d \sqrt{(1 + r^2)/2}}{delta = d * sqrt((1 +
#' r^2) / 2)} is the separation in noise-SD units. Write one of them in the
#' formula; a formula that names neither estimates `criterion ~ 1`, and one
#' that names both is refused. [sdt_criterion()] and [sdt_bias()] compute the
#' two from a single hit and false-alarm rate.
#'
#' **The choice matters whenever `d` or `sdratio` differs between rows that
#' share a criterion term**, because a linear formula holds constant the
#' parameter you write. `criterion ~ 1` with `d ~ 1 + strength` says the
#' observer keeps one boundary while targets get stronger; `bias ~ 1` with the
#' same `d` formula says the boundary moves up by half of every change in
#' separation. When the criterion and `d` share condition-wise (factor)
#' predictors and `sdratio` does not vary within those cells, or share
#' predictors while `d` keeps its identity link and `sdratio` is the same in
#' every row, the two give the same likelihood and their coefficients differ
#' by \eqn{\delta / 2}.
#' Otherwise they are different models. Even with the same likelihood, the
#' posteriors differ through the default priors: `criterion` takes
#' `normal(0.5, 1.5)` for its intercept, or for each cell under
#' `~ 0 + condition`, and `exponential(1)` for its random-effect SDs, while
#' `bias` takes `normal(0, 1.5)` and `exponential(2)`.
#'
#' Use `criterion` when the lures, and the observer's boundary, do not depend
#' on the manipulation of `d`: the classic recognition design in which weak and
#' strong targets are tested against one shared set of lures under the same
#' criteria. Its lure rows have no strength level, and with `criterion` they
#' need none, because the false-alarm rate does not involve `d`. Use `bias`
#' when you want to hold, or compare, the observer's position relative to the
#' two distributions, e.g. a bias contrast between groups that differ in
#' sensitivity.
#'
#' `bias` is measured from the midpoint, which is the point where the two
#' densities are equal -- the optimal boundary for equal base rates and
#' payoffs -- only for the symmetric distributions (`"normal"`, `"logistic"`)
#' under equal variance. Under unequal variance, or for the Gumbel
#' distributions, `bias = 0` is the midpoint, not the optimum.
#'
#' @section Identifying `sdratio`:
#' `sdratio` needs a design that supplies more than one operating point, and
#' whether it gets one depends on which criterion parameter the formula holds
#' constant.
#'
#' **With `criterion`**, the boundary must move among rows that share the
#' same `d`: give `criterion` a predictor that shifts it while `d` stays put,
#' such as a base-rate, payoff, or confidence manipulation (`d ~ 1`,
#' `criterion ~ 0 + condition`; see [broeder_schuetz_2009_e3]). A criterion
#' that is the same in every row is fixed by the false-alarm rate, and each
#' hit rate then fixes only \eqn{(\delta - criterion) / r}{(delta -
#' criterion) / r}, so `sdratio` trades against `d`. A change in `sdratio` is
#' matched by shifting the separation of every cell by one common amount and
#' rescaling it by one common factor, and whenever the formula for `d` can
#' make that change, `sdratio` returns its prior. A factor can, under the
#' identity and the log link alike (`d ~ 0 + condition` or `d ~ 1 +
#' condition`), because it gives every cell a value of its own. An intercept
#' with a numeric predictor can under the identity link (`d ~ 1 + strength`),
#' because the intercept shifts and the slope rescales. `sdratio` is
#' identified only when the formula constrains `d` across cells beyond that.
#' A numeric predictor without an intercept (`d ~ 0 + strength`) holds `d`
#' proportional to `strength`. It identifies `sdratio` when `strength` takes
#' at least two nonzero values and the criterion lies clearly away from the
#' noise mean (for `"normal"`, a false-alarm rate clearly away from .5); at a
#' criterion of 0 it does not. The estimate then rests on the assumption that
#' `d` is zero where `strength` is zero, not on a second operating point in
#' the data. A numeric predictor under a log link on `d` barely identifies
#' `sdratio`, and in practice `sdratio` returns its prior. `bmm()` warns for
#' every criterion that is the same in every row, so it also warns for
#' `d ~ 0 + strength`. Giving `d` and `criterion` the same condition-wise
#' predictors (`d ~ 0 + condition`, `criterion ~ 0 + condition`) fails too,
#' because each condition then has two rates for two parameters of its own,
#' and no warning catches it. A random effect such as
#' `criterion ~ 1 + (1 | id)` moves the boundary between subjects and does not
#' trigger the warning, but it informs `sdratio`
#' only weakly, through the assumption that the subject effects on `d` and on
#' the criterion are independent.
#'
#' **With `bias`**, any predictor that moves the operating point supplies what
#' is missing, including a sensitivity manipulation such as `d ~ 0 + condition`
#' with `bias ~ 1`: holding the bias constant moves the criterion with the
#' separation, and that movement traces the ROC. The identification then rests
#' on that assumption about the observer. The design in which every parameter
#' is intercept-only with no random effects cannot work: a single (hit,
#' false-alarm) pair is two numbers for three unknowns, and `bmm()` warns for
#' it. Neither can `bias` and `d` with the same condition-wise predictors, for
#' the reason given above, and that design raises no warning.
#'
#' Sampling converges and `Rhat` is fine in every unidentified case, so the
#' warnings are the only signal. They count formula terms, not how far the
#' terms move the operating point, so a predictor with a small effect carries
#' little information without a warning. Leaving `sdratio` at its default is
#' always identified.
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
#' predictive check written by hand. `d`, `criterion` and `bias` have identity
#' links, so they carry across unchanged.
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
#'     the noise distribution's SD, with its origin at the noise mean.
#'     `criterion` is a position on this axis; `bias` is a distance on it,
#'     measured from the midpoint between the means.
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
#' @seealso [sdt_d()], [sdt_criterion()] and [sdt_bias()] compute the `d`,
#'   `criterion` and `bias` of this model in closed form from a single pair of
#'   observed hit and false-alarm rates, without fitting: use them for a quick check of a fitted value, and
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
#'                      d = 1.5, criterion = 0.9)
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
#' # The same observer, with the boundary measured from the midpoint
#' fit_bias <- bmm(
#'   formula = bmf(d ~ 1 + (1 | id), bias ~ 1 + (1 | id)),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Weak and strong targets against one shared set of lures, one criterion per
#' # subject: `strong` is 1 for strong targets and 0 for weak targets and lures
#' # (the lure rows' value does not matter, their false-alarm rate depends on
#' # the criterion alone)
#' lure_dat <- expand.grid(id = 1:20, cell = c("lure", "weak", "strong"))
#' lure_dat$stimulus <- as.integer(lure_dat$cell != "lure")
#' lure_dat$strong <- as.integer(lure_dat$cell == "strong")
#' lure_dat$n_trials <- 100L
#' lure_dat$n_old <- rsdt_yn(nrow(lure_dat), lure_dat$n_trials,
#'                           lure_dat$stimulus, d = 1 + 0.8 * lure_dat$strong,
#'                           criterion = 0.7)
#' fit_strength <- bmm(
#'   formula = bmf(d ~ 1 + strong + (1 | id), criterion ~ 1 + (1 | id)),
#'   data = lure_dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Simulating from a fitted unequal-variance model: the model reports
#' # sdratio on its log link, rsdt_yn() takes the ratio itself.
#' sdratio_posterior <- 0.375
#' rsdt_yn(2, 100L, c(0L, 1L),
#'         d = 1.3, criterion = 0.75, sdratio = exp(sdratio_posterior))
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
  is_flat <- function(x) all(rhs_vars(x) %in% "Intercept")
  free_sdratio <- !is_constant(formula)[["sdratio"]]
  # Measured from the noise mean, a criterion that is the same in every row is
  # fixed by the false-alarm rate, and each hit rate then fixes only
  # (d * s - criterion) / sdratio: a change in sdratio is matched by shifting
  # every cell's separation by one amount, proportional to the criterion, and
  # rescaling it by one factor. A factor on d absorbs that under any link, an
  # intercept plus a numeric slope under the identity link, and under the log
  # link a numeric slope leaves only a trace. Only a numeric slope without an
  # intercept rules out the shift, and only while the criterion is away from
  # 0. Terms do not reveal that, so the warning errs on the side of firing.
  # A bias held constant moves the criterion with d and sdratio instead, which
  # is why only the all-flat design is caught there.
  if (.sdt_anchor_parameter(formula) == "criterion") {
    warnif(
      free_sdratio && (is.null(formula$criterion) || is_flat(formula$criterion)),
      "'sdratio' is estimated, but 'criterion' is the same in every row (no \\
      predictor, no random effect). Measured from the noise mean, such a \\
      criterion is fixed by the false-alarm rate, so every hit rate can be \\
      matched by trading 'sdratio' against 'd'. This leaves 'sdratio' at its \\
      prior, while sampling converges and Rhat stays fine, unless the formula \\
      for 'd' rules out shifting every cell's separation by one common amount \\
      and rescaling it by one common factor (e.g. a numeric predictor without \\
      an intercept and an identity link, with the criterion away from the \\
      noise mean). Give \\
      'criterion' a predictor that moves it (a base-rate, payoff, or \\
      confidence manipulation), or leave 'sdratio' at its default. See the \\
      'Identifying sdratio' section of ?sdt_yn"
    )
    return(NextMethod("check_formula"))
  }
  warnif(
    all(vapply(formula, is_flat, logical(1))) && free_sdratio,
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
  anchor <- .sdt_anchor_parameter(formula)
  formula <- bmf2bf(model, formula)

  # the dpar names must match the parameter the formula estimates, so each
  # parameterization is its own family; the Stan functions share one kernel
  formula$family <- brms::custom_family(
    if (anchor == "bias") "sdt_yn_bias" else "sdt_yn",
    dpars = c("mu", "d", anchor, "sdratio"),
    links = c("identity", model$links$d, model$links[[anchor]],
              model$links$sdratio),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_sdt_yn,
    posterior_predict = posterior_predict_sdt_yn,
    posterior_epred = posterior_epred_sdt_yn,
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
  sdratio <- brms::get_dpar(prep, "sdratio", i = i)
  criterion <- .sdt_criterion_dpar(prep, d, sdratio, i = i)
  dist <- .sdt_dist_names[prep$data$vint2[i]]

  dsdt_yn(prep$data$Y[i], prep$data$trials[i], prep$data$vint1[i],
              d, criterion, sdratio = sdratio, dist = dist, log = TRUE)
}

posterior_predict_sdt_yn <- function(i, prep, ...) {
  d <- brms::get_dpar(prep, "d", i = i)
  sdratio <- brms::get_dpar(prep, "sdratio", i = i)
  criterion <- .sdt_criterion_dpar(prep, d, sdratio, i = i)
  dist <- .sdt_dist_names[prep$data$vint2[i]]

  rsdt_yn(length(d), prep$data$trials[i], prep$data$vint1[i],
              d, criterion, sdratio = sdratio, dist = dist)
}

posterior_epred_sdt_yn <- function(prep) {
  stimulus <- .epred_data(prep$data$vint1, prep)
  dist <- .sdt_dist_names[prep$data$vint2[1]]
  criterion <- .sdt_criterion_dpar(prep, prep$dpars$d, prep$dpars$sdratio)
  p_old <- with(prep$dpars, exp(.sdt_log_p_old(.sdt_eta(d, criterion, stimulus, sdratio), dist)))
  .epred_matrix(.epred_data(prep$data$trials, prep) * p_old, prep)
}
