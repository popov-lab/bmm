############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_mafc <- function(response = NULL, n_trials = NULL, m = NULL,
                            dist = "normal",
                            links = NULL, call = NULL, ...) {
  parameters <- list(
    dprime = "Sensitivity: distance between signal and noise distributions"
  )
  default_priors <- list(
    dprime = list(main = "normal(1, 1)", effects = "normal(0, 0.5)")
  )
  requirements <- glue(
    "Provide pre-aggregated accuracy data with the following columns:", "\n\n",
    "  - Response counts (n_correct): number of correct responses", "\n",
    "  - Number of trials (n_trials): total trials per cell", "\n",
    "  No stimulus column needed (each trial has exactly one signal alternative)"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(n_trials, dist, m),
      domain = "Perception & Recognition Memory",
      task = "m-Alternative Forced Choice",
      name = "Signal Detection Theory (m-AFC)",
      citation = glue(
        "Green, D. M., & Swets, J. A. (1966). Signal detection theory ",
        "and psychophysics. Wiley.\n",
        "DeCarlo, L. T. (2012). On a signal detection approach to ",
        "m-alternative forced choice with bias, with maximum likelihood ",
        "and Bayesian approaches to estimation. JMASM, 11(1), 257-282."
      ),
      version = "mafc",
      requirements = requirements,
      parameters = parameters,
      links = list(dprime = "identity"),
      fixed_parameters = list(),
      default_priors = default_priors,
      init_ranges = list(dprime = c(0.5, 1.5)),
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdt", "sdt_mafc"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title m-Alternative Forced Choice Signal Detection Theory Model
#' @name sdt_mafc
#' @details `r model_info(.model_sdt_mafc())`
#'
#' Models accuracy in m-AFC tasks where each trial presents one signal among
#' `m` alternatives and the observer chooses the strongest one. Only `dprime`
#' is estimated (m-AFC has no response bias). The probability correct,
#' \eqn{P_c = \int f(x - d')\, F(x)^{m-1}\, dx}, is computed per noise
#' distribution: a closed-form softmax for `gumbel_min`, a closed-form Gamma
#' ratio for `gumbel_max`, Gauss-Hermite quadrature for `normal`
#' (\eqn{\Phi(d'/\sqrt{2})} at m = 2), and Gauss-Legendre quadrature for
#' `logistic`.
#' @param response A single string naming the column with counts of correct
#'   responses.
#' @param n_trials The name of the variable containing the total number of
#'   trials per cell.
#' @param m Either a single integer >= 2 giving the number of alternatives
#'   (constant across all rows), or a single string naming a data column that
#'   gives the number of alternatives per row. A column lets trials with
#'   different set sizes be fit jointly.
#' @param dist The noise distribution assumed for the latent evidence variable.
#'   One of:
#'   \itemize{
#'     \item "normal" (default): Gaussian m-AFC
#'     \item "logistic": Logistic m-AFC
#'     \item "gumbel_min": Extreme-value (Gumbel minimum) m-AFC
#'     \item "gumbel_max": Gumbel maximum m-AFC
#'   }
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Green, D. M., & Swets, J. A. (1966). \emph{Signal detection theory and
#'   psychophysics}. Wiley.
#'
#' DeCarlo, L. T. (2012). On a signal detection approach to m-alternative
#'   forced choice with bias, with maximum likelihood and Bayesian approaches
#'   to estimation. \emph{Journal of Mathematical and Statistical Psychology},
#'   \emph{11}(1), 257--282.
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' dat <- data.frame(id = 1:20, n_trials = 200L)
#' dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4,
#'                            dprime = rnorm(20, 1.5, 0.4))
#'
#' model <- sdt_mafc(
#'   response = "n_correct",
#'   n_trials = "n_trials",
#'   m = 4
#' )
#'
#' fit <- bmm(
#'   formula = bmf(dprime ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' }
sdt_mafc <- function(response, n_trials, m,
                     dist = c("normal", "logistic", "gumbel_min", "gumbel_max"),
                     links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)

  stopif(!((is.numeric(m) && length(m) == 1 && m >= 2) ||
           (is.character(m) && length(m) == 1)),
         "m must be a single integer >= 2, or the name of a set-size column in the data")
  if (is.numeric(m)) m <- as.integer(m)

  .model_sdt_mafc(response = response, n_trials = n_trials,
                  m = m, dist = dist,
                  links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_mafc <- function(model, data, formula) {
  resp_var <- model$resp_vars$response
  n_trials_var <- model$other_vars$n_trials

  .validate_sdt_counts(data, resp_var, n_trials_var)

  data$m_afc <- .sdt_resolve_set_size(model$other_vars$m, data)
  data$dist_type <- .sdt_dist_id(model$other_vars$dist)

  NextMethod("check_data")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt_mafc <- function(model, formula) {
  resp_var <- model$resp_vars$response
  n_trials_var <- model$other_vars$n_trials

  brms::bf(paste0(
    resp_var, " | vint(m_afc, dist_type) + trials(",
    n_trials_var, ") ~ 0"
  ))
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_mafc <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  formula$family <- brms::custom_family(
    "sdt_mafc",
    dpars = c("mu", "dprime"),
    links = c("identity", model$links$dprime),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_sdt_mafc,
    posterior_predict = posterior_predict_sdt_mafc,
    vars = c("vint1[n]", "vint2[n]", "trials[n]")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- paste(
    read_lines2(paste0(sc_path, "/sdt_dist_funs.stan")),
    read_lines2(paste0(sc_path, "/sdt_mafc_funs.stan")),
    sep = "\n"
  )
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

log_lik_sdt_mafc <- function(i, prep) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  dist <- .sdt_dist_names[prep$data$vint2[i]]

  dsdt_mafc(prep$data$Y[i], prep$data$trials[i], prep$data$vint1[i],
            dprime, dist = dist, log = TRUE)
}

posterior_predict_sdt_mafc <- function(i, prep, ...) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  dist <- .sdt_dist_names[prep$data$vint2[i]]

  rsdt_mafc(length(dprime), prep$data$trials[i], prep$data$vint1[i],
            dprime, dist = dist)
}
