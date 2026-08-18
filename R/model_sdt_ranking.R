############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_ranking <- function(response = NULL, m = NULL,
                               dist = "gumbel_min",
                               links = NULL, call = NULL, ...) {
  parameters <- list(
    d = paste0(
      "Sensitivity: the balanced discriminability index d_a, the distance ",
      "between the target and lure distributions in units of their ",
      "root-mean-square SD; equals d' (g' for gumbel_min) under equal variance"
    )
  )
  default_priors <- list(
    d = list(main = "normal(1, 1)", effects = "normal(0, 0.5)")
  )
  param_links <- list(d = "identity")
  fixed_pars <- list()
  init_ranges <- list(d = c(0.5, 1.5))

  # Gaussian ranking carries sdratio as an overridable fixed parameter (fixed to
  # 0 = equal variance, sampled when the user adds sdratio ~ ...), mirroring the
  # rating and binary SDT models.
  if (dist == "normal") {
    parameters$sdratio <- paste0(
      "Log SD ratio: log ratio of signal to noise standard deviations ",
      "(exp(sdratio) is the ratio, 0 = equal variance)"
    )
    # Matches sdt_yn: on the log scale, normal(0, 0.3) covers ratios in
    # [0.56, 1.80] at 95%, which spans the empirical recognition range and keeps
    # the prior inside the interval the Gauss-Hermite ladder is calibrated over
    # (see .ranking_gh_n).
    default_priors$sdratio <- list(
      main = "normal(0, 0.3)", effects = "normal(0, 0.15)"
    )
    param_links$sdratio <- "identity"
    fixed_pars$sdratio <- 0
    init_ranges$sdratio <- c(-0.3, 0.3)
  }

  requirements <- glue(
    "Provide pre-aggregated ranking counts in wide format:", "\n\n",
    "  - Rank-count columns (response): one column per rank position, each ",
    "giving the number of trials in which the target received that rank ",
    "(column 1 = most likely target, column m = least)", "\n",
    "  - Set size (m): a constant or a column giving the number of ranked ",
    "items per row; rows with fewer ranks leave the surplus columns at 0", "\n",
    "  No stimulus column needed (all trials include exactly one target)"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(m, dist),
      domain = "Perception & Recognition Memory",
      task = "Ranking Task",
      name = "Signal Detection Theory (Ranking)",
      citation = glue(
        "Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H. ",
        "(2025). Extreme-value signal detection theory for recognition memory: ",
        "The parametric road not taken. PsyArXiv. ",
        "https://doi.org/10.31234/osf.io/qhrfj"
      ),
      version = "NA",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = fixed_pars,
      default_priors = default_priors,
      init_ranges = init_ranges
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
#' Models the rank ordering of `m` items by perceived strength. Only `d`
#' is estimated (no criterion). Supports `dist = "gumbel_min"` (closed-form via
#' lgamma ratios) and `dist = "normal"` (Gauss-Hermite quadrature).
#'
#' The model uses the native brms multinomial family: each rank position is a
#' multinomial category whose logit is set to `log p(rank)`, so `softmax`
#' recovers the rank distribution exactly. This means `log_lik`,
#' `posterior_predict`, `posterior_epred`, and `pp_check` come from brms as
#' proper joint multinomial draws.
#'
#' @section Sensitivity is on the same scale as [sdt_yn()]:
#' `d` is the balanced index \eqn{d_a} that the rest of the SDT family reports:
#' the separation between the target and lure distributions divided by the
#' root-mean-square of their SDs. It equals \eqn{d'} whenever the two share a
#' scale, which is always the case for `dist = "gumbel_min"` (there it is the
#' \eqn{g'} of Meyer-Grant et al.) and is the default for `dist = "normal"`.
#'
#' Ranking is the one SDT design that identifies the variance ratio from a
#' single condition. The rank distribution supplies `m - 1` free probabilities
#' per set size, so with `m >= 3` there is enough information to separate `d`
#' from `sdratio` without the criterion sweep that [sdt_yn()] needs -- the
#' *shape* of the rank distribution, not just its mean, carries the ratio.
#'
#' At `m = 2` the model reduces to 2AFC and the two parameters are no longer
#' separable: the probability of ranking the target first is the area under the
#' yes/no ROC, which for Gaussian noise is \eqn{\Phi(d/\sqrt{2})} whatever
#' `sdratio` is. Keep `sdratio` fixed for two-item designs.
#'
#' For Gaussian ranking (`dist = "normal"`), `sdratio` is fixed to 0 by default.
#' Add `sdratio ~ 1` to the formula for unequal-variance ranking.
#'
#' The set size `m` may be a constant or the name of a data column. Supply a
#' column to fit trials with different set sizes in a single model: the response
#' has `max(m)` columns, and rows with a smaller set size switch off the surplus
#' rank categories (the multinomial then renormalizes over the valid ranks).
#' @param response A character vector of column names with the target rank-count
#'   columns, ordered from rank 1 (most likely target) to rank `m` (least).
#'   With a constant `m`, supply exactly `m` columns; with a varying set size,
#'   supply `max(m)` columns.
#' @param m Either a single integer >= 2 giving the number of ranked items
#'   (constant across all rows), or a single string naming a data column that
#'   gives the number of ranked items per row.
#' @param dist Character. The distribution assumed for the latent evidence,
#'   given here by its cumulative distribution function:
#'   \itemize{
#'     \item "gumbel_min" (default): smallest extreme value,
#'       \eqn{1 - \exp(-\exp(x))} (complementary log-log). Closed form via
#'       gamma-function ratios.
#'     \item "normal": Gaussian, \eqn{\Phi(x)} (supports unequal variance via
#'       `sdratio`)
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
#' dat <- data.frame(id = 1:20)
#' dat <- cbind(dat, rsdt_ranking(20, 200, m = 4, d = 1.5))
#'
#' model <- sdt_ranking(
#'   response = c("rank1", "rank2", "rank3", "rank4"),
#'   m = 4
#' )
#'
#' fit <- bmm(
#'   formula = bmf(d ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Mixed set sizes in one model: response has max(m) columns, m is a column
#' model_mixed <- sdt_ranking(
#'   response = c("rank1", "rank2", "rank3", "rank4", "rank5"),
#'   m = "set_size"
#' )
#' }
sdt_ranking <- function(response, m,
                        dist = c("gumbel_min", "normal"),
                        links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)

  stopif(!((is.numeric(m) && length(m) == 1 && m >= 2) ||
           (is.character(m) && length(m) == 1)),
         "m must be a single integer >= 2, or the name of a set-size column in the data")
  if (is.numeric(m)) m <- as.integer(m)

  stopif(length(response) < 2,
         "response must name at least 2 rank-count columns")
  stopif(is.numeric(m) && length(response) != m,
         "With a constant m ({m}), response must name exactly m columns ({length(response)} supplied)")

  .model_sdt_ranking(response = response, m = m, dist = dist,
                     links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_ranking <- function(model, data, formula) {
  resp_cols <- model$resp_vars$response
  .validate_sdt_count_cols(data, resp_cols)

  max_rank <- .sdt_resolve_set_size(model$other_vars$m, data)
  n_ranks <- length(resp_cols)
  stopif(any(max_rank > n_ranks, na.rm = TRUE),
         "Set size must not exceed the number of rank columns ({n_ranks})")

  Y <- as.matrix(data[resp_cols])
  # NAs are only allowed as structural blanks beyond a row's set size; an NA
  # within the set size would silently shrink nTrials if zeroed.
  stopif(anyNA(Y[col(Y) <= max_rank]),
         "Response columns must not contain NA counts within the row's set \\
         size (m); only rank columns beyond the set size may be NA")
  Y[is.na(Y)] <- 0
  stopif(any(rowSums(Y) <= 0, na.rm = TRUE),
         "Row sums of response columns must be positive (no empty rows)")

  # Counts in rank columns beyond a row's set size must be structural zeros: a
  # nonzero count there means the set size is mislabelled for that row.
  surplus <- col(Y) > max_rank
  stopif(any(Y[surplus] != 0),
         "Rank columns beyond the row's set size (m) must be 0")

  reserved <- intersect(c("Y", "nTrials", "max_rank"), colnames(data))
  warnif(length(reserved) > 0,
         "Column(s) {collapse_comma(reserved)} in your data are reserved by \\
         {model$name} and will be overwritten")
  data <- data[!colnames(data) %in% resp_cols]
  data$Y <- Y
  data$nTrials <- rowSums(Y)
  data$max_rank <- as.numeric(max_rank)

  free_sdratio <- "sdratio" %in% names(model$parameters) &&
    !"sdratio" %in% names(model$fixed_parameters)
  warnif(free_sdratio && max(max_rank, na.rm = TRUE) > 8,
         "Estimating `sdratio` with set sizes above 8 pushes the Gaussian \\
         quadrature past the node count bmm calibrates for (128), so the rank \\
         probabilities may carry errors above 1e-6. Consider fixing `sdratio`, \\
         or `dist = \"gumbel_min\"`, which is closed form")
  warnif(free_sdratio && max(max_rank, na.rm = TRUE) < 3,
         "`sdratio` cannot be identified from two-item rankings: the model \\
         reduces to 2AFC, where the rank probabilities depend only on `d`. \\
         Remove `sdratio` from the formula so it stays fixed at equal variance")

  NextMethod("check_data")
}


############################################################################# !
# MULTINOMIAL FORMULA & FAMILY CONSTRUCTION                              ####
############################################################################# !

# sdt_ranking uses brms' native multinomial family: each rank position's logit
# is set to log(p(rank)) so softmax recovers the ranking probabilities. The
# ranking math is computed by the static Stan function sdt_ranking_logmu for
# fitting and by the exported R companion sdt_ranking_logmu() for
# posterior_predict/epred (brms evaluates the non-linear formula in R for
# prediction). log_lik and posterior_predict therefore come from brms — proper
# joint multinomial draws — while both noise distributions stay supported.

#' @export
bmf2bf.sdt_ranking <- function(model, formula) {
  resp_cats <- model$resp_vars$response

  # max_rank (per-row set size) travels as a covariate; sdratio is the
  # parameter for the normal distribution and the literal 0 for gumbel_min;
  # the distribution id selects the kernel in both the Stan function and the
  # R companion.
  sdratio_arg <- if ("sdratio" %in% names(model$parameters)) "sdratio" else "0"
  args <- paste("max_rank", "d", sdratio_arg,
                .sdt_dist_id(model$other_vars$dist), sep = ", ")

  bform <- brms::bf(
    glue("Y | trials(nTrials) ~ sdt_ranking_logmu(1, {args})"),
    nl = TRUE
  )
  for (k in seq_along(resp_cats)[-1]) {
    bform <- bform + brms::nlf(stats::as.formula(
      glue("mu{resp_cats[k]} ~ sdt_ranking_logmu({k}, {args})")
    ))
  }
  bform
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_ranking <- function(model, data, formula) {
  resp_cats <- model$resp_vars$response

  formula <- bmf2bf(model, formula)
  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- resp_cats
  formula$family$dpars <- paste0("mu", resp_cats)

  # sdt_dist_funs.stan supplies sdt_rms_scale(), which the Gaussian branch needs
  # to turn d_a into the noise-standardized separation
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- paste(
    read_lines2(paste0(sc_path, "/sdt_dist_funs.stan")),
    .ranking_fill_quadrature(
      read_lines2(paste0(sc_path, "/sdt_ranking_funs.stan")),
      max_m = max(data$max_rank),
      free_sdratio = "sdratio" %in% names(model$parameters) &&
        !"sdratio" %in% names(model$fixed_parameters)
    ),
    sep = "\n"
  )
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}


# Substitute the quadrature tokens in sdt_ranking_funs.stan with the rule that
# .ranking_gh_n() selects, so the compiled Stan code and .ranking_prob_r() use
# the same nodes. 17 significant digits round-trips a double exactly.
.ranking_fill_quadrature <- function(scode, max_m, free_sdratio) {
  gh <- .gh_rule(.ranking_gh_n(max_m, free_sdratio))
  brace <- function(x) paste0("({", paste(sprintf("%.17g", x), collapse = ", "), "})")
  scode <- gsub("{{N_GH}}", length(gh$nodes), scode, fixed = TRUE)
  scode <- gsub("({{GH_NODES}})", brace(gh$nodes), scode, fixed = TRUE)
  gsub("({{GH_WEIGHTS}})", brace(gh$weights), scode, fixed = TRUE)
}


#' @title Ranking SDT rank log-probability (multinomial logit)
#' @description R companion to the Stan `sdt_ranking_logmu` function. It returns
#'   `log(p(rank))` for rank position `cat`, which the [sdt_ranking()]
#'   multinomial formula uses as the category logit (so `softmax` recovers the
#'   ranking probabilities). `brms` evaluates the non-linear formula in R for
#'   `posterior_predict()` and `posterior_epred()`, so this function must be on
#'   the search path; it is exported for that reason and is not called directly.
#' @param cat Integer rank position index.
#' @param max_rank Set size (number of ranked items) for the observation; ranks
#'   above it return the finite `-100` sentinel so the category switches off.
#' @param d Sensitivity (draws-by-observation matrix supplied by brms).
#' @param sdratio Log SD ratio; `0` for the gumbel_min distribution.
#' @param dist Integer noise-distribution id (2 = gumbel_min, 1 = normal).
#' @return `log(p(cat))`, matching the shape of `d`.
#' @keywords internal
#' @export
sdt_ranking_logmu <- function(cat, max_rank, d, sdratio = 0, dist = 2L) {
  dist_name <- .sdt_dist_names[dist]

  shape <- dim(d)
  d <- as.vector(d)
  n <- length(d)
  max_rank <- rep_len(as.vector(max_rank), n)
  sdratio <- rep_len(as.vector(sdratio), n)

  # grouping by the (few) unique set sizes vectorizes the quadrature over all
  # draws-by-observations elements sharing a set size
  out <- rep(-100, n)
  for (m_j in unique(max_rank[max_rank >= cat])) {
    idx <- which(max_rank == m_j)
    probs <- rbind(.ranking_all_probs_r(d[idx], m_j, dist_name,
                                        sdratio[idx]))
    out[idx] <- log(probs[, cat])
  }

  if (!is.null(shape)) dim(out) <- shape
  out
}
