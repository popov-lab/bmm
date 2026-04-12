############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_mafc <- function(response = NULL, n_trials = NULL, m = NULL,
                            dist = "normal",
                            links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)

  parameters <- list(
    mu = glue("Internal parameter (fixed to 0)"),
    dprime = glue("Sensitivity: distance between signal and noise distributions")
  )
  default_priors <- list(
    dprime = list(main = "normal(1, 1)", effects = "normal(0, 0.5)")
  )
  param_links <- list(mu = "identity", dprime = "identity")

  requirements <- glue(
    "Provide pre-aggregated accuracy data with the following columns:", "\n\n",
    "  - Correct responses ({response}): number of correct responses", "\n",
    "  - Number of trials ({n_trials}): total trials per cell", "\n",
    "  No stimulus column needed (each trial has exactly one signal alternative)"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(n_trials, dist, dist_int, m),
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
      links = param_links,
      fixed_parameters = list(mu = 0),
      default_priors = default_priors,
      init_ranges = list(mu = c(0, 0), dprime = c(0.5, 1.5)),
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
#' Models accuracy in m-AFC tasks. Only the `dprime` parameter is estimated
#' (no criterion). For 2-AFC, uses the closed-form `Phi(d'/sqrt(2))`. For
#' m >= 3, uses 20-point Gauss-Hermite quadrature.
#' @param response A single string naming the column with counts of correct
#'   responses.
#' @param n_trials The name of the variable containing the total number of
#'   trials per cell.
#' @param m Integer. Number of alternatives (must be >= 2).
#' @param dist Character. Currently only "normal" is supported.
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' DeCarlo, L. T. (2012). On a signal detection approach to m-alternative
#'   forced choice with bias, with maximum likelihood and Bayesian approaches
#'   to estimation. \emph{Journal of Mathematical and Statistical Psychology},
#'   \emph{11}(1), 257--282.
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' dat <- rsdt_mafc(n_per_cell = 200, n_subjects = 20,
#'                  dprime = 1.5, m = 4)
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

  stopif(!is.numeric(m) || length(m) != 1 || m < 2,
         "m must be a single integer >= 2 for m-AFC models")
  m <- as.integer(m)
  stopif(dist != "normal",
         "m-AFC currently only supports dist = 'normal'")

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

  data$m_afc <- model$other_vars$m

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
    resp_var, " | vint(m_afc) + trials(",
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
    vars = c("vint1[n]", "trials[n]")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdt_mafc_funs.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

# R-side Gauss-Hermite computation of P(correct) for m-AFC
# Mirrors the Stan mafc_pc function
.mafc_pc_r <- function(dprime, m) {
  if (m == 2L) return(pnorm(dprime / sqrt(2)))

  # 20-point Gauss-Hermite nodes and weights for standard normal integration
  gh_nodes <- c(
    -7.6190485416797546e+00, -6.5105901570136488e+00,
    -5.5787388058932059e+00, -4.7345813340460463e+00,
    -3.9439673506573110e+00, -3.1890148165533843e+00,
    -2.4586636111723603e+00, -1.7452473208141255e+00,
    -1.0429453488027509e+00, -3.4696415708135458e-01,
     3.4696415708135830e-01,  1.0429453488027574e+00,
     1.7452473208141317e+00,  2.4586636111723683e+00,
     3.1890148165533900e+00,  3.9439673506573163e+00,
     4.7345813340460552e+00,  5.5787388058932033e+00,
     6.5105901570136551e+00,  7.6190485416797591e+00
  )
  gh_weights <- c(
    1.2578006724378954e-13, 2.4820623623151972e-10,
    6.1274902599825256e-08, 4.4021210902309806e-06,
    1.2882627996193093e-04, 1.8301031310804826e-03,
    1.3997837447100857e-02, 6.1506372063977507e-02,
    1.6173933398399959e-01, 2.6079306344955683e-01,
    2.6079306344955305e-01, 1.6173933398399776e-01,
    6.1506372063977438e-02, 1.3997837447101162e-02,
    1.8301031310805052e-03, 1.2882627996193072e-04,
    4.4021210902309052e-06, 6.1274902599829068e-08,
    2.4820623623151936e-10, 1.2578006724379269e-13
  )

  log_terms <- log(gh_weights) + (m - 1) * pnorm(gh_nodes + dprime, log.p = TRUE)
  exp(matrixStats::logSumExp(log_terms))
}

log_lik_sdt_mafc <- function(i, prep) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  m <- prep$data$vint1[i]
  n_trials <- prep$data$trials[i]
  y <- prep$data$Y[i]

  pc <- vapply(dprime, .mafc_pc_r, numeric(1), m = m)
  dbinom(y, n_trials, pc, log = TRUE)
}

posterior_predict_sdt_mafc <- function(i, prep, ...) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  m <- prep$data$vint1[i]
  n_trials <- prep$data$trials[i]

  pc <- vapply(dprime, .mafc_pc_r, numeric(1), m = m)
  rbinom(length(pc), n_trials, pc)
}
