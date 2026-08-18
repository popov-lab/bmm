############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_yn <- function(response = NULL, stimulus = NULL,
                          n_trials = NULL, dist = "normal",
                          links = NULL, call = NULL, ...) {
  parameters <- list(
    d = paste0(
      "Sensitivity: the balanced discriminability index d_a, which measures ",
      "the distance between the signal and noise distributions in units of ",
      "their root-mean-square SD, so it equals d' when sdratio is 1"
    ),
    criterion = "Response bias: location of decision boundary",
    sdratio = paste0(
      "Log SD ratio: the log of the signal-to-noise standard deviation ",
      "ratio, so 0 means equal variance (an SD ratio of 1)"
    )
  )
  # d_a and the noise-standardized separation differ by at most ~16% over the
  # plausible sdratio range, so the sensitivity prior needs no recalibration.
  default_priors <- list(
    d = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
    criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)"),
    # sdratio is on the log scale, so normal(0, 0.3) puts 95% of the prior on
    # ratios in [0.56, 1.80], i.e. zROC slopes in [0.56, 1.80]. That spans the
    # empirical range with room -- Mickes et al. (2007) report a mean ratio of
    # 1.25 with subjects between 0.85 and 1.82, and the broeder_schuetz_2009_e3
    # fit gives 1.46 [1.25, 1.72] -- while excluding ratios above 2, which are
    # not observed in recognition and which no quadrature-based SDT model in
    # bmm is calibrated for.
    sdratio = list(main = "normal(0, 0.3)", effects = "normal(0, 0.15)")
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
#'   total number of trials for each cell.
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
#' When the signal and noise distributions have different widths there is no
#' single natural way to express their separation, because the separation only
#' becomes dimensionless after choosing a scale to divide by. `bmm` reports
#' \eqn{d_a}, the separation divided by the root-mean-square of the two SDs:
#' \deqn{d_a = \sqrt{2}\,\delta / \sqrt{1 + r^2},}
#' where \eqn{\delta} is the separation in noise-SD units and \eqn{r} is
#' `sdratio`. This weights the two distributions equally, and it is the measure
#' Macmillan and Creelman (2005) and Mickes et al. (2007) recommend under
#' unequal variance. `sdratio` is sampled on the log scale and fixed at 0 by
#' default (an SD ratio of 1), where \eqn{d_a} equals the familiar \eqn{d'}, so
#' the choice only matters for unequal-variance fits.
#'
#' The alternative — dividing by the noise SD alone — is not comparable across
#' conditions that differ in `sdratio`: two conditions that are equally
#' discriminable can then show a large, confidently estimated difference in
#' sensitivity. Note that `criterion` is **not** rescaled and stays on the
#' noise-standardized axis, following the same convention.
#'
#' Because `d` is a short name, a column called `d` in your data that is also
#' used as a predictor will collide with this parameter; `bmm()` warns when that
#' happens.
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
#' # Unequal-variance yes/no SDT
#' fit_uv <- bmm(
#'   formula = bmf(d ~ 1, criterion ~ 1, sdratio ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
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
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_yn <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stim_vals <- data[[stim_var]]
  stopif(!is.numeric(stim_vals) || !all(stim_vals %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (noise) and 1 (signal)")

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
