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
#' dat <- rsdt_mafc(dprime = rnorm(20, 1.5, 0.4), m = 4, n_trials = 200)
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
  data$dist_type <- model$other_vars$dist_int

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

# R-side probability correct for m-AFC, mirroring the Stan mafc_pc function.
# gumbel_min has a closed-form softmax, gumbel_max a closed-form Gamma ratio,
# normal uses Gauss-Hermite quadrature, logistic uses Gauss-Legendre on [0, 1].
.mafc_pc_r <- function(dprime, m, dist = "normal") {
  if (dist == "gumbel_min") {
    return(1 / (1 + (m - 1) * exp(-dprime)))
  }
  if (dist == "gumbel_max") {
    return(exp(lgamma(1 + exp(-dprime)) + lgamma(m) - lgamma(m + exp(-dprime))))
  }
  if (dist == "normal") {
    if (m == 2L) return(stats::pnorm(dprime / sqrt(2)))
    gh_nodes <- c(
      -1.14533778415487379e+01, -1.04815605346742640e+01, -9.67355636693402765e+00, -8.94950454385556249e+00,
      -8.27894062365948535e+00, -7.64616376454146440e+00, -7.04173840645382576e+00, -6.45942337758375906e+00,
      -5.89480567537201416e+00, -5.34460544572008622e+00, -4.80628719209386723e+00, -4.27782615636274777e+00,
      -3.75755977616898207e+00, -3.24408873299986844e+00, -2.73620834046542960e+00, -2.23285921863486791e+00,
      -1.73309059063171489e+00, -1.23603200479915287e+00, -7.40870725285924792e-01, -2.46832896022723958e-01,
       2.46832896022727510e-01,  7.40870725285931897e-01,  1.23603200479916175e+00,  1.73309059063172377e+00,
       2.23285921863487502e+00,  2.73620834046543315e+00,  3.24408873299987022e+00,  3.75755977616898384e+00,
       4.27782615636274954e+00,  4.80628719209387523e+00,  5.34460544572008622e+00,  5.89480567537201683e+00,
       6.45942337758376706e+00,  7.04173840645382842e+00,  7.64616376454145907e+00,  8.27894062365947647e+00,
       8.94950454385555538e+00,  9.67355636693403120e+00,  1.04815605346742657e+01,  1.14533778415487308e+01
    )
    gh_weights <- c(
       1.46183987386930516e-29,  4.82046794020072741e-25,  1.44860943155167746e-21,  1.12227520682703716e-18,
       3.38985344324777725e-16,  4.96808852919722085e-14,  4.03763858169491567e-12,  1.98911852602780986e-10,
       6.32589718854883025e-09,  1.36034242157482606e-07,  2.04889743608149897e-06,  2.22117714324753619e-05,
       1.77072928799239520e-04,  1.05587901690180051e-03,  4.77354488182319455e-03,  1.65378441425691192e-02,
       4.42745552022761890e-02,  9.21765791700618065e-02,  1.49921111763569481e-01,  1.91059009661991935e-01,
       1.91059009661987633e-01,  1.49921111763571979e-01,  9.21765791700600856e-02,  4.42745552022768551e-02,
       1.65378441425699553e-02,  4.77354488182340705e-03,  1.05587901690182349e-03,  1.77072928799244128e-04,
       2.22117714324759446e-05,  2.04889743608150575e-06,  1.36034242157490811e-07,  6.32589718854897914e-09,
       1.98911852602780831e-10,  4.03763858169524929e-12,  4.96808852919782859e-14,  3.38985344324820570e-16,
       1.12227520682709321e-18,  1.44860943155158925e-21,  4.82046794020079904e-25,  1.46183987386941726e-29
    )
    log_terms <- log(gh_weights) + (m - 1) * stats::pnorm(gh_nodes + dprime, log.p = TRUE)
    return(exp(matrixStats::logSumExp(log_terms)))
  }

  gl_nodes <- c(
    3.47479132114081324e-04, 1.82994161402261213e-03, 4.49331426162824510e-03, 8.33187305768723352e-03,
    1.33365861050445123e-02, 1.94956001739736706e-02, 2.67943125707985619e-02, 3.52154139340299377e-02,
    4.47389314607484767e-02, 5.53422770024430966e-02, 6.70003009229536151e-02, 7.96853518737098421e-02,
    9.33673424386013417e-02, 1.08013820528329307e-01, 1.23590046369734252e-01, 1.40059074914194670e-01,
    1.57381843472883531e-01, 1.75517264372671455e-01, 1.94422322413803195e-01, 2.14052176898682944e-01,
    2.34360267990052940e-01, 2.55298427146473550e-01, 2.76816991373267984e-01, 2.98864921018004326e-01,
    3.21389920831166132e-01, 3.44338564004894487e-01, 3.67656418895616288e-01, 3.91288178129996389e-01,
    4.15177789788003682e-01, 4.39268590351939658e-01, 4.63503439106100479e-01, 4.87824853668287650e-01,
    5.12175146331712128e-01, 5.36496560893899632e-01, 5.60731409648060231e-01, 5.84822210211996207e-01,
    6.08711821870003611e-01, 6.32343581104383823e-01, 6.55661435995105624e-01, 6.78610079168834091e-01,
    7.01135078981995896e-01, 7.23183008626732016e-01, 7.44701572853526450e-01, 7.65639732009947283e-01,
    7.85947823101317056e-01, 8.05577677586196583e-01, 8.24482735627328656e-01, 8.42618156527116580e-01,
    8.59940925085805441e-01, 8.76409953630265970e-01, 8.91986179471670804e-01, 9.06632657561398769e-01,
    9.20314648126290269e-01, 9.32999699077046385e-01, 9.44657722997557014e-01, 9.55261068539251412e-01,
    9.64784586065969840e-01, 9.73205687429201438e-01, 9.80504399826026884e-01, 9.86663413894955488e-01,
    9.91668126942312989e-01, 9.95506685738372088e-01, 9.98170058385977610e-01, 9.99652520867886141e-01
  )
  gl_weights <- c(
    8.91640360848292403e-04, 2.07351663028077990e-03, 3.25222898448934228e-03, 4.42337991318164214e-03,
    5.58406973006490056e-03, 6.73152394835961187e-03, 7.86301523801327153e-03, 8.97585788784879823e-03,
    1.00674115767644488e-02, 1.11350869041918748e-02, 1.21763512843553277e-02, 1.31887348575272928e-02,
    1.41698363071303297e-02, 1.51173285362008816e-02, 1.60289641774256156e-02, 1.69025809185709350e-02,
    1.77361066284417446e-02, 1.85275642701207242e-02, 1.92750765893085693e-02, 1.99768705663597619e-02,
    2.06312816213114793e-02, 2.12367575618269550e-02, 2.17918622646615483e-02, 2.22952790818779396e-02,
    2.27458139637090571e-02, 2.31423982906572012e-02, 2.34840914081049928e-02, 2.37700828574146442e-02,
    2.39996942982292662e-02, 2.41723811174017686e-02, 2.42877337207512492e-02, 2.43454785045694837e-02,
    2.43454785045695184e-02, 2.42877337207515927e-02, 2.41723811174016194e-02, 2.39996942982290164e-02,
    2.37700828574152270e-02, 2.34840914081047535e-02, 2.31423982906576141e-02, 2.27458139637086477e-02,
    2.22952790818782449e-02, 2.17918622646617877e-02, 2.12367575618269445e-02, 2.06312816213118297e-02,
    1.99768705663600533e-02, 1.92750765893081426e-02, 1.85275642701202004e-02, 1.77361066284411686e-02,
    1.69025809185707407e-02, 1.60289641774257301e-02, 1.51173285362012598e-02, 1.41698363071297347e-02,
    1.31887348575274056e-02, 1.21763512843555601e-02, 1.11350869041916285e-02, 1.00674115767651357e-02,
    8.97585788784880864e-03, 7.86301523801245968e-03, 6.73152394835926579e-03, 5.58406973006549817e-03,
    4.42337991318194831e-03, 3.25222898448917618e-03, 2.07351663028124350e-03, 8.91640360848207835e-04
  )
  cdf <- .SDT_DISTS[[dist]]$cdf
  qf <- .SDT_DISTS[[dist]]$qf
  sum(gl_weights * cdf(qf(gl_nodes) + dprime)^(m - 1))
}

log_lik_sdt_mafc <- function(i, prep) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  m <- prep$data$vint1[i]
  dist <- .sdt_dist_names[prep$data$vint2[i]]
  n_trials <- prep$data$trials[i]
  y <- prep$data$Y[i]

  pc <- vapply(dprime, .mafc_pc_r, numeric(1), m = m, dist = dist)
  stats::dbinom(y, n_trials, pc, log = TRUE)
}

posterior_predict_sdt_mafc <- function(i, prep, ...) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  m <- prep$data$vint1[i]
  dist <- .sdt_dist_names[prep$data$vint2[i]]
  n_trials <- prep$data$trials[i]

  pc <- vapply(dprime, .mafc_pc_r, numeric(1), m = m, dist = dist)
  stats::rbinom(length(pc), n_trials, pc)
}
