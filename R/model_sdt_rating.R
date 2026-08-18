############################################################################# !
# CONFIDENCE-THRESHOLD PARAMETERIZATIONS                                  ####
############################################################################# !

# Threshold-parameterization registry: the position defines the integer
# thresh_type code passed to Stan -- reordering entries changes the R <-> Stan
# contract (see sdt_make_thresholds_rating in sdt_rating_funs.stan)
.sdt_threshold_types <- c("parsimonious", "equidistant", "log_distance",
                          "log_ratio", "softmax")

# Parameter spec for one threshold parameterization. parsimonious/equidistant
# need a single spacing; log_distance/log_ratio need K-2 distance deltas;
# softmax needs spacing plus K-3 allocation deltas. `anchor` is the threshold
# index that carries no delta (it sits at the criterion); it defaults to the
# same middle threshold the builders skip ((K - 1) %/% 2 + 1, so delta labels
# name the threshold they control for odd K); sdt_cdp passes its old/new
# boundary (n_new) instead.
.sdt_threshold_parameter_parts <- function(n_ratings, threshold_type,
                                           anchor = (n_ratings - 1L) %/% 2L + 1L) {
  parameters <- list()
  default_priors <- list()
  param_links <- list()

  if (threshold_type %in% c("equidistant", "parsimonious")) {
    parameters$spacing <- paste0(
      "Threshold spacing: controls distance between adjacent thresholds ",
      "(exp(spacing) ensures positive spacing)"
    )
    default_priors$spacing <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
    )
    param_links$spacing <- "identity"
  } else if (threshold_type %in% c("log_distance", "log_ratio")) {
    n_deltas <- n_ratings - 2L
    for (i in seq_len(n_deltas)) {
      idx <- if (i < anchor) i else i + 1L
      pname <- paste0("delta", idx)
      parameters[[pname]] <- glue("Threshold parameter for threshold {idx}")
      default_priors[[pname]] <- list(
        main = "normal(0, 1)", effects = "normal(0, 0.5)"
      )
      param_links[[pname]] <- "identity"
    }
  } else if (threshold_type == "softmax") {
    parameters$spacing <- paste0(
      "Average log spacing across threshold intervals ",
      "(exp(spacing) is the mean interval size)"
    )
    default_priors$spacing <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
    )
    param_links$spacing <- "identity"

    n_deltas <- max(0L, n_ratings - 3L)
    for (i in seq_len(n_deltas)) {
      pname <- paste0("delta", i)
      parameters[[pname]] <- glue(
        "Softmax threshold allocation parameter for interval {i}"
      )
      default_priors[[pname]] <- list(
        main = "normal(0, 1)", effects = "normal(0, 0.5)"
      )
      param_links[[pname]] <- "identity"
    }
  }

  nlist(parameters, default_priors, param_links)
}

.sdt_threshold_delta_names <- function(object) {
  grep("^delta", names(object$parameters), value = TRUE)
}


############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# Per-version specification for the rating SDT lattice. `standard` is a single
# familiarity process; `dpsdt` adds a recollection threshold process (Yonelinas,
# 1994) whose Ro/Rn are fixed near zero by default (recovering standard SDT) and
# freed through the formula; `metad` adds a type-2 metacognitive process
# (Maniscalco & Lau, 2012) parameterized by the log M-ratio. Each entry carries
# the extra parameters plus the Stan
# logmu function name and call the shared multinomial pipeline assembles, so the
# formula, family, and prediction code stay version-agnostic.
.sdt_rating_variants <- list(
  standard = list(
    parameters = list(),
    default_priors = list(),
    links = list(),
    init_ranges = list(),
    fixed_parameters = list(),
    extra_params = character(0),
    logmu_fun = "sdt_rating_logmu",
    logmu_cat_call = "sdt_rating_logmu_cat(cat, thr, d, sdratio, stimulus, dist_type)",
    stan_chunk = "sdt_rating_funs.stan"
  ),
  dpsdt = list(
    parameters = list(
      Ro = glue("Recollection of old items: probability inv_logit(Ro) that an ",
                "old item is recollected as old, loading the most-confident ",
                "'signal' category"),
      Rn = glue("Recollection of new items: probability inv_logit(Rn) that a ",
                "new item is recall-rejected, loading the most-confident ",
                "'noise' category")
    ),
    default_priors = list(
      Ro = list(main = "normal(0, 1)", effects = "normal(0, 0.5)"),
      Rn = list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
    ),
    links = list(Ro = "identity", Rn = "identity"),
    init_ranges = list(Ro = c(-0.5, 0.5), Rn = c(-0.5, 0.5)),
    # inv_logit(-100) is numerically 0: recollection off unless freed by formula
    fixed_parameters = list(Ro = -100, Rn = -100),
    extra_params = c("Ro", "Rn"),
    logmu_fun = "sdt_dpsdt_logmu",
    logmu_cat_call = "sdt_dpsdt_logmu_cat(cat, thr, d, sdratio, stimulus, dist_type, Ro, Rn)",
    stan_chunk = "sdt_dpsdt_funs.stan"
  ),
  metad = list(
    parameters = list(
      logmratio = glue("Log M-ratio, log(meta-d/d): metacognitive efficiency. ",
                       "0 is ideal metacognition (meta-d = d), negative is ",
                       "inefficiency, positive is hyper-efficiency. meta-d is ",
                       "recovered as exp(logmratio) * d, on the same d_a scale")
    ),
    # Estimating log(meta-d/d) rather than meta-d directly keeps meta-d
    # positive, regularizes it toward d, and anchors the field-standard
    # metacognitive-efficiency measure (M-ratio) at the ideal point of 0
    # (Maniscalco & Lau, 2014; Fleming, 2017). Both sensitivities are d_a
    # indices, so the ratio is invariant to sdratio.
    default_priors = list(
      logmratio = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)")
    ),
    links = list(logmratio = "identity"),
    init_ranges = list(logmratio = c(-0.3, 0.3)),
    fixed_parameters = list(),
    extra_params = "logmratio",
    logmu_fun = "sdt_metad_logmu",
    logmu_cat_call = "sdt_metad_logmu_cat(cat, thr, d, exp(logmratio) * d, sdratio, stimulus, dist_type)",
    stan_chunk = "sdt_metad_funs.stan"
  )
)

.model_sdt_rating <- function(response = NULL, stimulus = NULL,
                              dist = "normal", n_ratings = NULL,
                              threshold_type = "parsimonious",
                              version = "standard",
                              links = NULL, call = NULL, ...) {
  if (is.null(n_ratings) && length(response) > 1) {
    n_ratings <- length(response)
  }

  parameters <- list(
    d = paste0(
      "Sensitivity: the balanced discriminability index d_a, which measures ",
      "the distance between the signal and noise distributions in units of ",
      "their root-mean-square SD, so it equals d' when sdratio is 1"
    ),
    criterion = paste0(
      "Response bias: location of the decision boundary on the ",
      "noise-standardized axis"
    )
  )
  default_priors <- list(
    d = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
    criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)")
  )
  param_links <- list(d = "identity", criterion = "identity")

  threshold_parts <- .sdt_threshold_parameter_parts(n_ratings, threshold_type)
  parameters <- c(parameters, threshold_parts$parameters)
  default_priors <- c(default_priors, threshold_parts$default_priors)
  param_links <- c(param_links, threshold_parts$param_links)

  parameters$sdratio <- paste0(
    "Log SD ratio: the log of the signal-to-noise standard deviation ratio, ",
    "so exp(sdratio) is the ratio itself and 0 means equal variance"
  )
  # Matches sdt_yn and sdt_ranking: on the log scale, normal(0, 0.3) covers
  # ratios in [0.56, 1.80] at 95%, spanning the empirical recognition range.
  default_priors$sdratio <- list(
    main = "normal(0, 0.3)", effects = "normal(0, 0.15)"
  )
  param_links$sdratio <- "identity"

  # Bounded starting values for every estimated parameter. Without them the
  # flexible threshold types start from brms' wide random init, where adjacent
  # thresholds can land in the same distribution tail and the category mass
  # underflows to log_diff_exp(-Inf, -Inf) = NaN, so every chain rejects its
  # initial value.
  init_ranges <- list(
    d = c(0.5, 1.5), criterion = c(-0.5, 0.5), sdratio = c(-0.3, 0.3)
  )
  for (p in names(threshold_parts$parameters)) {
    init_ranges[[p]] <- if (p == "spacing") c(-0.7, -0.2) else c(-0.5, 0.2)
  }

  variant <- .sdt_rating_variants[[version]]
  parameters <- c(parameters, variant$parameters)
  default_priors <- c(default_priors, variant$default_priors)
  param_links <- c(param_links, variant$links)
  init_ranges <- c(init_ranges, variant$init_ranges)
  fixed_parameters <- c(list(sdratio = 0), variant$fixed_parameters)

  requirements <- glue(
    "Provide pre-aggregated data with the following columns:", "\n\n",
    "  - Response counts: one column per rating category (K columns)", "\n",
    "  - Stimulus type (stimulus): 0 = noise, 1 = signal", "\n",
    "  Categories should be ordered: 1 = 'definitely noise' to ",
    "K = 'definitely signal'"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(stimulus, dist, n_ratings, threshold_type),
      domain = "Perception & Recognition Memory",
      task = "Signal/Noise or Old/New Recognition",
      name = "Signal Detection Theory (Confidence Rating)",
      citation = glue(
        "Green, D. M., & Swets, J. A. (1966). Signal detection theory ",
        "and psychophysics. Wiley."
      ),
      version = version,
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = fixed_parameters,
      default_priors = default_priors,
      init_ranges = init_ranges
    ),
    class = c("bmmodel", "sdt", "sdt_rating", paste0("sdt_rating_", version)),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Confidence Rating Signal Detection Theory Model
#' @name sdt_rating
#' @details Three versions share the confidence-rating response interface,
#'   selected with `version`:
#'
#' #### Version: `standard` (default)
#' `r model_info(.model_sdt_rating(version = "standard"))`
#'
#' #### Version: `dpsdt`
#' Dual-process SDT (Yonelinas, 1994): a familiarity SDT process plus a
#' recollection threshold process. `Ro` is recollection of old items (loads the
#' most-confident "signal" category) and `Rn` is recall-to-reject of new items
#' (loads the most-confident "noise" category); `inv_logit(Ro)`/`inv_logit(Rn)`
#' are the recollection probabilities. Both are fixed off by default
#' (recovering `standard`); add `Ro ~ 1` for the one-sided model and
#' `Ro ~ 1, Rn ~ 1` for the two-sided model. `d` describes the familiarity
#' distributions only -- the observed ROC is a mixture of familiarity and
#' recollection, so it is not the discriminability of that mixture.
#' `r model_info(.model_sdt_rating(version = "dpsdt"))`
#'
#' #### Version: `metad`
#' Meta-d' (Maniscalco & Lau, 2012): a type-2 metacognitive sensitivity governs
#' confidence-threshold placement, with the total "old"/"new" response rates
#' held to what type-1 `d` implies. Rather than estimating meta-d'
#' directly, the model estimates `logmratio`, the log M-ratio
#' \eqn{\log(\mathrm{meta\text{-}d'}/d')}, and recovers meta-d' as
#' `exp(logmratio) * d`. The M-ratio is the field-standard measure of
#' metacognitive efficiency (Maniscalco & Lau, 2014; Fleming, 2017): estimating
#' it on the log scale keeps meta-d' positive, regularizes it toward `d`,
#' and anchors the ideal point (meta-d' = `d`, perfect metacognition) at
#' `logmratio = 0`, which recovers `standard`. Type-1 and type-2 sensitivity are
#' both on the \eqn{d_a} scale, so the M-ratio is unaffected by `sdratio`.
#' Extract the M-ratio posterior with [mratio()].
#' `r model_info(.model_sdt_rating(version = "metad"))`
#'
#' By default, the model assumes equal variance (sdratio fixed to 0). To
#' estimate unequal variance, add `sdratio ~ 1` (or `sdratio ~ predictors`)
#' to the formula.
#'
#' @section Sensitivity is on the same scale as [sdt_yn()]:
#' `d` is the balanced index \eqn{d_a} that [sdt_yn()] reports: the separation
#' between the signal and noise distributions divided by the root-mean-square
#' of their SDs. It equals \eqn{d'} whenever `sdratio` stays fixed at 0, and
#' unlike \eqn{d'} it remains comparable across conditions that differ in
#' `sdratio` -- see the sensitivity section of [sdt_yn()] for the reasoning.
#' The `criterion` and the confidence thresholds are **not** rescaled and stay
#' on the noise-standardized axis. The `dpsdt` and `metad` versions inherit the
#' same convention: `d` there is the \eqn{d_a} of the familiarity (type-1)
#' distributions, and meta-d' is scaled the same way, which leaves the M-ratio
#' invariant to `sdratio`.
#' @param response A character vector of K column names containing response
#'   counts per rating category, ordered from "definitely noise" to
#'   "definitely signal".
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (noise/new) and 1 (signal/old).
#' @param n_ratings Integer. Number of response categories. When `response`
#'   is a character vector, defaults to its length. Must be > 2.
#' @param dist The distribution assumed for the latent evidence, given here by
#'   its cumulative distribution function. One of:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT, \eqn{\Phi(x)}
#'     \item "gumbel_min": smallest-extreme-value SDT,
#'       \eqn{1 - \exp(-\exp(x))} (complementary log-log)
#'     \item "gumbel_max": largest-extreme-value SDT, \eqn{\exp(-\exp(-x))}
#'       (log-log, as in \code{evd::pgumbel})
#'     \item "logistic": logistic SDT, \eqn{1 / (1 + \exp(-x))}
#'   }
#' @param threshold_type Character. Threshold parameterization:
#'   \itemize{
#'     \item "parsimonious" (default): 2 parameters (criterion + spacing).
#'       Thresholds follow logit-spaced canonical positions (Selker et al.,
#'       2019).
#'     \item "equidistant": 2 parameters (criterion + spacing). Thresholds
#'       are equally spaced around criterion.
#'     \item "log_distance": K-2 parameters. Each threshold distance is
#'       exp(delta), guaranteeing ordering.
#'     \item "log_ratio": K-2 parameters. Threshold distances as ratios
#'       (Paulewicz & Blaut, 2020).
#'     \item "softmax": K-2 parameters. A shared spacing parameter controls
#'       the average interval size, while softmax-transformed delta parameters
#'       allocate interval widths smoothly across the scale.
#'   }
#' @param version Character. The latent-process version of the rating model.
#'   One of `"standard"` (default, a single familiarity SDT process),
#'   `"dpsdt"` (dual-process SDT with recollection parameters `Ro`/`Rn`), or
#'   `"metad"` (meta-d' with a metacognitive efficiency parameter `logmratio`,
#'   the log M-ratio). All three use the same confidence-rating response
#'   interface. See Details.
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Yonelinas, A. P. (1994). Receiver-operating characteristics in recognition
#'   memory: Evidence for a dual-process model. \emph{Journal of Experimental
#'   Psychology: Learning, Memory, and Cognition}, \emph{20}(6), 1341--1354.
#'
#' Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for
#'   estimating metacognitive sensitivity from confidence ratings.
#'   \emph{Consciousness and Cognition}, \emph{21}(1), 422--430.
#'
#' Maniscalco, B., & Lau, H. (2014). Signal detection theory analysis of type 1
#'   and type 2 data: meta-d', response-specific meta-d', and the unequal
#'   variance SDT model. In S. M. Fleming & C. D. Frith (Eds.), \emph{The
#'   cognitive neuroscience of metacognition} (pp. 25--66). Springer.
#'
#' Fleming, S. M. (2017). HMeta-d: hierarchical Bayesian estimation of
#'   metacognitive efficiency from confidence ratings. \emph{Neuroscience of
#'   Consciousness}, \emph{2017}(1), nix007. \doi{10.1093/nc/nix007}
#'
#' Green, D. M., & Swets, J. A. (1966). \emph{Signal detection theory and
#'   psychophysics}. Wiley.
#'
#' Selker, R., van den Bergh, D., Criss, A. H., & Wagenmakers, E.-J. (2019).
#'   Parsimonious estimation of signal detection models from confidence ratings.
#'   \emph{Behavior Research Methods}, \emph{51}(5), 1953--1967.
#'   \doi{10.3758/s13428-019-01231-3}
#'
#' Paulewicz, B., & Blaut, A. (2020). The bhsdtr package: A general-purpose
#'   method of Bayesian inference for signal detection theory models.
#'   \emph{Behavior Research Methods}, \emph{52}, 2122--2141.
#'   \doi{10.3758/s13428-020-01370-y}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # EV-SDT rating model
#' dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
#' dat <- cbind(dat, rsdt_rating(nrow(dat), 200, dat$stimulus,
#'                               d = 1.5, thresholds = c(-0.5, 0, 0.5)))
#'
#' model <- sdt_rating(
#'   response = c("r1", "r2", "r3", "r4"),
#'   stimulus = "stimulus"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # UV-SDT: add sdratio to the formula
#' fit_uv <- bmm(
#'   formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # Dual-process SDT: free recollection of old items (one-sided) or both
#' # old and new items (two-sided) via the formula
#' model_dp <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
#' fit_dp <- bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1),
#'               data = dat, model = model_dp, backend = "cmdstanr")
#'
#' # Meta-d': estimate metacognitive efficiency (log M-ratio)
#' model_md <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
#' fit_md <- bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1),
#'               data = dat, model = model_md, backend = "cmdstanr")
#' mratio(fit_md) # posterior M-ratio (meta-d'/d')
#' }
sdt_rating <- function(response, stimulus,
                       n_ratings = NULL,
                       dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
                       threshold_type = c("parsimonious", "equidistant",
                                          "log_distance", "log_ratio",
                                          "softmax"),
                       version = c("standard", "dpsdt", "metad"),
                       links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)
  threshold_type <- match.arg(threshold_type)
  version <- match.arg(version)

  if (length(response) > 1 && is.null(n_ratings)) {
    n_ratings <- length(response)
  }

  stopif(is.null(n_ratings) || n_ratings <= 2,
         "Rating models require n_ratings > 2 (or pass a response vector with > 2 columns)")

  stopif(threshold_type == "log_ratio" && n_ratings < 4,
         "log_ratio thresholds require n_ratings >= 4 (the anchor ratio needs an interval above the criterion)")

  if (length(response) > 1) {
    stopif(n_ratings != length(response),
           "n_ratings ({n_ratings}) must match the number of response columns ({length(response)})")
  }

  .model_sdt_rating(response = response, stimulus = stimulus,
                    dist = dist, n_ratings = n_ratings,
                    threshold_type = threshold_type, version = version,
                    links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_rating <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stim_vals <- data[[stim_var]]
  stopif(!is.numeric(stim_vals) || !all(stim_vals %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (noise) and 1 (signal)")

  resp_cols <- model$resp_vars$response
  .validate_sdt_count_cols(data, resp_cols)

  Y <- as.matrix(data[resp_cols])
  stopif(anyNA(Y), "Response columns must not contain NA counts")
  stopif(any(rowSums(Y) <= 0),
         "Row sums of response columns must be positive (no empty rows)")

  reserved <- intersect(c("Y", "nTrials"), colnames(data))
  warnif(length(reserved) > 0,
         "Column(s) {collapse_comma(reserved)} in your data are reserved by \\
         {model$name} and will be overwritten")
  data <- data[!colnames(data) %in% resp_cols]
  data$Y <- Y
  data$nTrials <- rowSums(Y)

  NextMethod("check_data")
}


############################################################################# !
# MULTINOMIAL FORMULA & FAMILY CONSTRUCTION                              ####
############################################################################# !

# sdt_rating uses brms' native multinomial family: each category's logit is set
# to log(p_k^SDT) so softmax recovers the SDT probabilities exactly. The SDT math
# (thresholds + noise CDF) is computed by a per-model Stan function for fitting
# and by the exported R companion sdt_rating_logmu() for posterior_predict/epred
# (brms evaluates the non-linear formula in R for prediction, looking the
# function up on the search path). log_lik and posterior_predict therefore come
# from brms — proper joint multinomial draws — while the four distributions and
# five threshold parameterizations stay supported.

# The dispatch pieces (kernel names, extra parameters, Stan chunk) are pure
# derived state: look them up in the variant registry instead of storing them
# in the model object. Legacy model objects from before the version machinery
# carry version = "rating" or "NA" and resolve to the standard variant.
.sdt_rating_variant <- function(model) {
  .sdt_rating_variants[[model$version]] %||% .sdt_rating_variants$standard
}

# cat, K, dist, and threshold type travel as integer literals so the same
# non-linear call resolves against both the generated Stan function and the R
# companion. spacing is the parameter when present, otherwise the literal 0
# (threshold types without spacing ignore it).
.sdt_rating_logmu_args <- function(model) {
  has_spacing <- "spacing" %in% names(model$parameters)
  c(model$other_vars$n_ratings, .sdt_dist_id(model$other_vars$dist),
    match(model$other_vars$threshold_type, .sdt_threshold_types),
    "d", "criterion",
    if (has_spacing) "spacing" else "0", "sdratio",
    .sdt_rating_variant(model)$extra_params,
    model$other_vars$stimulus, .sdt_threshold_delta_names(model))
}

# Stan function generated per model: the delta arity is fixed here, while cat,
# K, dist, and threshold type arrive as the integer arguments described above.
.sdt_rating_logmu_stan <- function(model) {
  variant <- .sdt_rating_variant(model)
  delta_names <- .sdt_threshold_delta_names(model)
  nd <- length(delta_names)
  extra <- variant$extra_params

  signature <- paste(
    c("int cat", "int K", "int dist_type", "int thresh_type",
      "real d", "real criterion", "real spacing", "real sdratio",
      if (length(extra)) paste("real", extra),
      "real stimulus", if (nd) paste("real", delta_names)),
    collapse = ", "
  )
  delta_decl <- if (nd == 0L) {
    "  array[0] real deltas;\n"
  } else {
    paste0("  array[", nd, "] real deltas;\n",
           paste(sprintf("  deltas[%d] = %s;", seq_len(nd), delta_names),
                 collapse = "\n"), "\n")
  }

  paste0(
    "real ", variant$logmu_fun, "(", signature, ") {\n",
    delta_decl,
    "  vector[K - 1] thr = sdt_make_thresholds_rating(criterion, spacing, deltas, K, thresh_type);\n",
    "  return ", variant$logmu_cat_call, ";\n",
    "}\n"
  )
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

# Base multinomial brmsformula: Y | trials(nTrials) carries the counts and each
# category gets a non-linear mu calling the version's logmu function. The user's
# parameter formulas (d, criterion, spacing/deltas, sdratio) are added afterwards
# by bmf2bf.bmmodel, so they are deliberately not added here.
#' @export
bmf2bf.sdt_rating <- function(model, formula) {
  resp_cats <- model$resp_vars$response
  fun <- .sdt_rating_variant(model)$logmu_fun
  args <- paste(.sdt_rating_logmu_args(model), collapse = ", ")

  bform <- brms::bf(
    glue("Y | trials(nTrials) ~ {fun}(1, {args})"),
    nl = TRUE
  )
  for (k in seq_along(resp_cats)[-1]) {
    bform <- bform + brms::nlf(stats::as.formula(
      glue("mu{resp_cats[k]} ~ {fun}({k}, {args})")
    ))
  }
  bform
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_rating <- function(model, data, formula) {
  resp_cats <- model$resp_vars$response

  formula <- bmf2bf(model, formula)
  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- resp_cats
  formula$family$dpars <- paste0("mu", resp_cats)

  sc_path <- system.file("stan_chunks", package = "bmm")
  chunks <- unique(c("sdt_dist_funs.stan", "sdt_rating_funs.stan",
                     .sdt_rating_variant(model)$stan_chunk))
  stan_funs <- paste(
    paste(vapply(chunks, function(f) read_lines2(paste0(sc_path, "/", f)),
                 character(1)), collapse = "\n"),
    .sdt_rating_logmu_stan(model),
    sep = "\n"
  )
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}


#' @title Rating SDT category log-probability (multinomial logit)
#' @description R companion to the Stan `sdt_rating_logmu` function. It returns
#'   `log(p_k)` for rating category `cat`, which the [sdt_rating()] multinomial
#'   formula uses as the category logit (so `softmax` recovers the SDT category
#'   probabilities). `brms` evaluates the non-linear formula in R for
#'   `posterior_predict()` and `posterior_epred()`, so this function must be on
#'   the search path; it is exported for that reason and is not called directly.
#' @param cat Integer rating category index.
#' @param K Integer number of rating categories.
#' @param dist Integer noise-distribution id (see the `.sdt_dists` registry).
#' @param thresh Integer threshold-parameterization id.
#' @param d,criterion,spacing,sdratio Model parameters (draws-by-observation
#'   matrices supplied by brms). `d` is the balanced sensitivity index d_a;
#'   `spacing` is `0` for threshold types without it.
#' @param Ro,Rn Linear-scale recollection parameters for the `dpsdt` version;
#'   `inv_logit(Ro)`/`inv_logit(Rn)` are the recollection probabilities.
#' @param logmratio Log M-ratio for the `metad` version; meta-d' is recovered
#'   as `exp(logmratio) * d`, on the same d_a scale as `d`.
#' @param stimulus Stimulus covariate (0 = noise, 1 = signal).
#' @param ... Threshold `delta` parameters, when the threshold type uses them.
#' @return `log(p_cat)`, matching the shape of `d`.
#' @keywords internal
#' @export
sdt_rating_logmu <- function(cat, K, dist, thresh, d, criterion, spacing,
                             sdratio, stimulus, ...) {
  dist_name <- .sdt_dist_names[dist]
  thresh_name <- .sdt_threshold_types[thresh]

  shape <- dim(d)
  d <- as.vector(d)
  n <- length(d)
  criterion <- rep_len(as.vector(criterion), n)
  spacing <- rep_len(as.vector(spacing), n)
  sdratio <- rep_len(as.vector(sdratio), n)
  stimulus <- rep_len(as.vector(stimulus), n)
  deltas <- if (...length() > 0L) {
    do.call(cbind, lapply(list(...), function(x) rep_len(as.vector(x), n)))
  }

  thr <- .sdt_make_thresholds(criterion, K, thresh_name, spacing, deltas)
  probs <- rbind(.sdt_category_probs(rbind(thr), d, exp(sdratio),
                                     stimulus, dist_name))
  out <- log(probs[, cat])

  if (!is.null(shape)) dim(out) <- shape
  out
}


#' @rdname sdt_rating_logmu
#' @export
sdt_dpsdt_logmu <- function(cat, K, dist, thresh, d, criterion, spacing,
                            sdratio, Ro, Rn, stimulus, ...) {
  dist_name <- .sdt_dist_names[dist]
  thresh_name <- .sdt_threshold_types[thresh]

  shape <- dim(d)
  d <- as.vector(d)
  n <- length(d)
  criterion <- rep_len(as.vector(criterion), n)
  spacing <- rep_len(as.vector(spacing), n)
  sdratio <- rep_len(as.vector(sdratio), n)
  Ro <- rep_len(as.vector(Ro), n)
  Rn <- rep_len(as.vector(Rn), n)
  stimulus <- rep_len(as.vector(stimulus), n)
  deltas <- if (...length() > 0L) {
    do.call(cbind, lapply(list(...), function(x) rep_len(as.vector(x), n)))
  }

  thr <- .sdt_make_thresholds(criterion, K, thresh_name, spacing, deltas)
  probs <- rbind(.sdt_dpsdt_category_probs(rbind(thr), d, exp(sdratio),
                                           stimulus, dist_name,
                                           stats::plogis(Ro),
                                           stats::plogis(Rn)))
  out <- log(probs[, cat])

  if (!is.null(shape)) dim(out) <- shape
  out
}


#' @rdname sdt_rating_logmu
#' @export
sdt_metad_logmu <- function(cat, K, dist, thresh, d, criterion, spacing,
                            sdratio, logmratio, stimulus, ...) {
  dist_name <- .sdt_dist_names[dist]
  thresh_name <- .sdt_threshold_types[thresh]

  shape <- dim(d)
  d <- as.vector(d)
  n <- length(d)
  criterion <- rep_len(as.vector(criterion), n)
  spacing <- rep_len(as.vector(spacing), n)
  sdratio <- rep_len(as.vector(sdratio), n)
  # meta-d is derived from the estimated log M-ratio, mirroring the Stan call
  metad <- exp(rep_len(as.vector(logmratio), n)) * d
  stimulus <- rep_len(as.vector(stimulus), n)
  deltas <- if (...length() > 0L) {
    do.call(cbind, lapply(list(...), function(x) rep_len(as.vector(x), n)))
  }

  thr <- .sdt_make_thresholds(criterion, K, thresh_name, spacing, deltas)
  probs <- rbind(.sdt_metad_category_probs(rbind(thr), d, metad, stimulus,
                                           exp(sdratio), dist_name))
  out <- log(probs[, cat])

  if (!is.null(shape)) dim(out) <- shape
  out
}
