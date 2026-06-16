############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_ranking <- function(response = NULL, rank = NULL, m = NULL,
                               dist = "gumbel_min",
                               links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)

  parameters <- list(
    mu = glue("Internal parameter (fixed to 0)"),
    dprime = glue("Sensitivity: ranking discrimination parameter")
  )
  default_priors <- list(
    dprime = list(main = "normal(1, 1)", effects = "normal(0, 0.5)")
  )
  param_links <- list(mu = "identity", dprime = "identity")

  fixed_pars <- list(mu = 0)

  # Gaussian ranking supports sdratio as overridable fixed parameter
  if (dist == "normal") {
    parameters$sdratio <- glue(
      "SD ratio: log ratio of signal to noise standard deviations ",
      "(exp(sdratio) ensures positivity, 0 = equal variance)"
    )
    default_priors$sdratio <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
    )
    param_links$sdratio <- "identity"
    fixed_pars$sdratio <- 0
  }

  requirements <- glue(
    "Provide pre-aggregated ranking data in long format:", "\n\n",
    "  - Response counts ({response}): number of times target received this rank", "\n",
    "  - Rank position ({rank}): rank position (1 = most likely target, m = least)", "\n",
    "  No stimulus column needed (all items include exactly one target)"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(rank, dist, dist_int, m),
      domain = "Perception & Recognition Memory",
      task = "Ranking Task",
      name = "Signal Detection Theory (Ranking)",
      citation = glue(
        "Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H. ",
        "(2025). Extreme-value signal detection theory for recognition memory: ",
        "The parametric road not taken. PsyArXiv. ",
        "https://doi.org/10.31234/osf.io/qhrfj"
      ),
      version = "ranking",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = fixed_pars,
      default_priors = default_priors,
      init_ranges = list(mu = c(0, 0), dprime = c(0.5, 1.5),
                         sdratio = c(-0.3, 0.3)),
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdt", "sdt_ranking"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Ranking Signal Detection Theory Model
#' @name sdt_ranking
#' @details `r model_info(.model_sdt_ranking())`
#'
#' Models rank ordering of m items by perceived strength. Only `dprime`
#' is estimated (no criterion). Supports `dist = "gumbel_min"` (closed-form
#' via lgamma ratios) and `dist = "normal"` (numerical integration).
#'
#' For Gaussian ranking (`dist = "normal"`), sdratio is fixed to 0 by
#' default. Add `sdratio ~ 1` to the formula for unequal-variance ranking.
#'
#' The set size `m` may be a constant or the name of a data column. Supply a
#' column to fit trials with different set sizes in a single model.
#' @param response A single string naming the column with rank frequency
#'   counts.
#' @param rank The name of the variable coding the rank position (1 = most
#'   likely target, m = least likely).
#' @param m Either a single integer >= 2 giving the number of ranked items
#'   (constant across all rows), or a single string naming a data column that
#'   gives the number of ranked items per row.
#' @param dist Character. The noise distribution:
#'   \itemize{
#'     \item "gumbel_min" (default): Closed-form Gumbel-min ranking
#'     \item "normal": Gaussian ranking (supports unequal variance via
#'       sdratio)
#'   }
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H. (2025).
#'   \emph{Extreme-value signal detection theory for recognition memory: The
#'   parametric road not taken}. PsyArXiv preprint.
#'   \doi{10.31234/osf.io/qhrfj}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' dat <- rsdt_ranking(n_per_cell = 200, n_subjects = 20,
#'                     dprime = 1.5, m = 4)
#'
#' model <- sdt_ranking(
#'   response = "observed",
#'   rank = "rank",
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
sdt_ranking <- function(response, rank, m,
                        dist = c("gumbel_min", "normal"),
                        links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)

  stopif(!((is.numeric(m) && length(m) == 1 && m >= 2) ||
           (is.character(m) && length(m) == 1)),
         "m must be a single integer >= 2, or the name of a set-size column in the data")
  if (is.numeric(m)) m <- as.integer(m)

  .model_sdt_ranking(response = response, rank = rank,
                     m = m, dist = dist,
                     links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_ranking <- function(model, data, formula) {
  resp_var <- model$resp_vars$response
  rank_var <- model$other_vars$rank

  stopif(!resp_var %in% colnames(data),
         "Response variable '{resp_var}' missing in the data")
  resp_vals <- data[[resp_var]]
  stopif(any(resp_vals < 0, na.rm = TRUE),
         "Response variable '{resp_var}' must contain non-negative counts")
  warnif(any(resp_vals != round(resp_vals), na.rm = TRUE),
         "Response variable '{resp_var}' should contain integer counts")

  stopif(!rank_var %in% colnames(data),
         "Rank variable '{rank_var}' missing in the data")
  max_rank <- .sdt_resolve_set_size(model$other_vars$m, data)
  rank_vals <- as.integer(data[[rank_var]])
  stopif(any(rank_vals < 1 | rank_vals > max_rank, na.rm = TRUE),
         "Rank variable '{rank_var}' must be between 1 and the set size in every row")
  warnif(any(data[[rank_var]] != round(data[[rank_var]]), na.rm = TRUE),
         "Rank variable '{rank_var}' should contain integer values")

  data$rank_pos <- rank_vals
  data$max_rank <- max_rank

  NextMethod("check_data")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt_ranking <- function(model, formula) {
  resp_var <- model$resp_vars$response
  brms::bf(paste0(resp_var, " | vint(rank_pos, max_rank) ~ 0"))
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_ranking <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdt_ranking_funs.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  # Gaussian ranking always carries sdratio as a dpar (fixed to 0 for the
  # equal-variance default, sampled when predicted), mirroring sdt_binary.
  if (model$other_vars$dist == "normal") {
    formula$family <- brms::custom_family(
      "sdt_ranking_uv",
      dpars = c("mu", "dprime", "sdratio"),
      links = c("identity", model$links$dprime, model$links$sdratio),
      type = "int",
      loop = TRUE,
      log_lik = log_lik_sdt_ranking_uv,
      posterior_predict = posterior_predict_sdt_ranking,
      vars = c("vint1[n]", "vint2[n]")
    )
  } else {
    formula$family <- brms::custom_family(
      "sdt_ranking",
      dpars = c("mu", "dprime"),
      links = c("identity", model$links$dprime),
      type = "int",
      loop = TRUE,
      log_lik = log_lik_sdt_ranking,
      posterior_predict = posterior_predict_sdt_ranking,
      vars = c("vint1[n]", "vint2[n]")
    )
  }

  nlist(formula, data, stanvars)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

# R-side computation of rank probability
# Mirrors the Stan sdt_ranking_lpmf kernel
.ranking_prob_r <- function(dprime, rank_pos, m, dist = "gumbel_min",
                            sdratio = 0) {
  if (dist == "gumbel_min") {
    g <- dprime
    e_neg_g <- exp(-g)
    log_p <- -g + lgamma(m) + lgamma(rank_pos - 1 + e_neg_g) -
             lgamma(rank_pos) - lgamma(m + e_neg_g)
    exp(log_p)
  } else {
    # Gaussian UV-SDT: fixed Gauss-Hermite quadrature over the target distribution
    sigma <- exp(sdratio)
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
    eta <- dprime + sigma * gh_nodes
    log_terms <- log(gh_weights) +
      (m - rank_pos) * pnorm(eta, log.p = TRUE) +
      (rank_pos - 1) * pnorm(eta, lower.tail = FALSE, log.p = TRUE)
    exp(lchoose(m - 1, rank_pos - 1) + matrixStats::logSumExp(log_terms))
  }
}

# R-side computation of all rank probabilities (vectorized over ranks 1..m)
.ranking_all_probs_r <- function(dprime, m, dist = "gumbel_min",
                                 sdratio = 0) {
  probs <- vapply(seq_len(m), function(r) {
    .ranking_prob_r(dprime, r, m, dist, sdratio)
  }, numeric(1))
  probs / sum(probs)
}

log_lik_sdt_ranking <- function(i, prep) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  rank_pos <- prep$data$vint1[i]
  max_rank <- prep$data$vint2[i]
  y <- prep$data$Y[i]

  log_p <- vapply(dprime, function(d) {
    p <- .ranking_prob_r(d, rank_pos, max_rank, dist = "gumbel_min")
    log(p)
  }, numeric(1))
  y * log_p
}

log_lik_sdt_ranking_uv <- function(i, prep) {
  dprime   <- brms::get_dpar(prep, "dprime",   i = i)
  sdratio  <- brms::get_dpar(prep, "sdratio",  i = i)
  rank_pos <- prep$data$vint1[i]
  max_rank <- prep$data$vint2[i]
  y        <- prep$data$Y[i]

  log_p <- mapply(function(d, s) {
    p <- .ranking_prob_r(d, rank_pos, max_rank, dist = "normal", sdratio = s)
    log(p)
  }, dprime, sdratio)
  y * log_p
}

posterior_predict_sdt_ranking <- function(i, prep, ...) {
  dprime   <- brms::get_dpar(prep, "dprime", i = i)
  rank_pos <- prep$data$vint1[i]
  max_rank <- prep$data$vint2[i]

  if ("sdratio" %in% names(prep$dpars)) {
    dist    <- "normal"
    sdratio <- brms::get_dpar(prep, "sdratio", i = i)
  } else {
    dist    <- "gumbel_min"
    sdratio <- rep(0, length(dprime))
  }

  vapply(seq_along(dprime), function(s) {
    probs <- .ranking_all_probs_r(dprime[s], max_rank, dist, sdratio[s])
    probs[rank_pos]
  }, numeric(1))
}
