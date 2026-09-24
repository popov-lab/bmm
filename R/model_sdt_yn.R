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
      "SD ratio signal/noise, log link (0 = equal SDs): the parameter is the ",
      "ratio itself, so exp() a value reported on this scale to read it"
    )
  )
  # d is d_a, and the noise-standardized separation is d * sqrt((1 + r^2) / 2).
  # Over plausible SD ratios, r in [0.56, 1.80], that factor runs from
  # 0.81 to 1.46, which normal(1, 1) is wide enough to absorb on either scale.
  default_priors <- list(
    d = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
    criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)"),
    sdratio = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)")
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
  out$links[names(links)] <- links
  out
}


#' @title Yes/No Signal Detection Theory Model
#' @name sdt_yn
#' @details `r model_info(.model_sdt_yn())`
#' @param response The name of the variable in the dataset containing the
#'   count of "old"/"signal" responses for each cell.
#' @param stimulus The name of the variable in the dataset coding the stimulus
#'   type. Stimuli should be coded as 0 (noise/new) and 1 (signal/old).
#' @param n_trials The name of the variable in the dataset containing the
#'   total number of trials for each cell. It may differ from cell to cell.
#'   Trial-level data also works: keep one row per trial, with `response`
#'   0 or 1 and a column of `1`s for `n_trials`.
#' @param dist The noise distribution assumed for the latent evidence variable,
#'   given here by its cumulative distribution function. One of:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT, \eqn{\Phi(x)}
#'     \item "gumbel_min": smallest-extreme-value SDT,
#'       \eqn{1 - \exp(-\exp(x))} (complementary log-log)
#'     \item "gumbel_max": largest-extreme-value SDT, \eqn{\exp(-\exp(-x))}
#'       (log-log, as in \code{evd::pgumbel})
#'     \item "logistic": logistic SDT, \eqn{1 / (1 + \exp(-x))}
#'   }
#' @param links A named list of link functions for the parameters.
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
#' \deqn{d_a = \sqrt{2}\,\delta / \sqrt{1 + r^2},}
#' where \eqn{\delta} is the separation in noise-SD units and
#' \eqn{r} is the SD ratio, `exp(sdratio)`. This weights the two
#' distributions equally, and it is the measure Simpson and Fitter (1973),
#' Macmillan and Creelman (2005), and Mickes et al. (2007) recommend under
#' unequal variance. The classical noise-standardized index is
#' \eqn{d_N = \delta = d_a \sqrt{(1 + r^2)/2}}, so a published \eqn{d'} from an
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
#' AUC-equivalent index, \eqn{d_a = \sqrt{2}\,\Phi^{-1}(\mathrm{AUC})}, so it
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
#' is fine, and the profile likelihood over `sdratio` is flat to 1e-12.
#'
#' Any linear predictor that moves the operating point along the ROC supplies
#' what is missing, and it need not sit on `criterion`: a sensitivity
#' manipulation (`d ~ 0 + condition`, a study-time or strength manipulation with
#' bias held constant) identifies `sdratio` just as a criterion manipulation
#' does, and so does between-subject variation entering through a random effect
#' such as `criterion ~ 1 + (1 | id)`.
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
#'   \item `main` / `effects` — the keys of the default priors shown in the
#'     model description above: `main` is the prior on the intercept,
#'     `effects` the prior on regression coefficients.
#' }
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

  # an unmatched name is otherwise appended, so print() advertises a parameter
  # that does not exist and the intended link is never applied
  valid_links <- names(.model_sdt_yn()$links)
  stopif(!all(names(links) %in% valid_links),
         "Unrecognized link target(s): \\
         {collapse_comma(setdiff(names(links), valid_links))}. \\
         sdt_yn() takes links for {collapse_comma(valid_links)}")

  .model_sdt_yn(response = response, stimulus = stimulus,
                n_trials = n_trials, dist = dist,
                links = links, call = call, ...)
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
