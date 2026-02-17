############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# Internal helper to check if model is a rating model
.is_sdt_rating <- function(model) {
  n_ratings <- model$other_vars$n_ratings
  !is.null(n_ratings) && n_ratings > 2
}

.model_sdt <- function(response = NULL, stimulus = NULL, n_trials = NULL,
                       version = "evsdt", dist = "normal",
                       variances = "equal", n_ratings = NULL,
                       threshold_type = "equidistant",
                       m = NULL, links = NULL, call = NULL, ...) {
  # Map distribution to integer for Stan dispatch (from registry)
  dist_int <- .sdt_dist_id(dist)

  # Determine if this is a rating model
  is_rating <- !is.null(n_ratings) && n_ratings > 2

  # m-AFC has a simpler parameter structure
  if (version == "mafc") {
    return(.build_sdt_mafc(response, n_trials, m, dist, dist_int,
                           links, call))
  }

  # Build parameters and priors based on model type
  parameters <- list(
    dprime = glue("Sensitivity: distance between signal and noise distributions"),
    criterion = glue("Response bias: location of decision boundary")
  )
  default_priors <- list(
    dprime = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
    criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)")
  )
  param_links <- list(dprime = "identity", criterion = "identity")

  if (is_rating) {
    if (threshold_type %in% c("equidistant", "parsimonious")) {
      parameters$spacing <- glue(
        "Threshold spacing: controls distance between adjacent thresholds ",
        "(exp(spacing) ensures positive spacing)"
      )
      default_priors$spacing <- list(
        main = "normal(0, 1)", effects = "normal(0, 0.5)"
      )
      param_links$spacing <- "identity"
    } else if (threshold_type == "log_distance") {
      n_deltas <- n_ratings - 2L
      mid <- n_ratings %/% 2L
      for (i in seq_len(n_deltas)) {
        # Skip the middle index (that's criterion)
        idx <- if (i < mid) i else i + 1L
        pname <- paste0("delta", idx)
        parameters[[pname]] <- glue(
          "Log-distance threshold parameter for threshold {idx}"
        )
        default_priors[[pname]] <- list(
          main = "normal(0, 1)", effects = "normal(0, 0.5)"
        )
        param_links[[pname]] <- "identity"
      }
    }
  }

  # Unequal variance adds sdratio parameter (requires ratings for identification)
  if (variances == "unequal") {
    parameters$sdratio <- glue(
      "SD ratio: ratio of signal to noise standard deviations ",
      "(exp(sdratio) ensures positivity, 0 = equal variance)"
    )
    default_priors$sdratio <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
    )
    param_links$sdratio <- "identity"
  }

  # DPSDT adds Ro (recollection) parameter
  if (version == "dpsdt") {
    parameters$Ro <- glue(
      "Recollection: probability of threshold-based recollection ",
      "for old items (inv_logit(Ro) maps to [0,1])"
    )
    default_priors$Ro <- list(
      main = "normal(0, 1)", effects = "normal(0, 0.5)"
    )
    param_links$Ro <- "identity"
  }

  # Meta-d' adds metad (metacognitive sensitivity) parameter
  if (version == "metad") {
    parameters$metad <- glue(
      "Metacognitive sensitivity: type-2 d' controlling confidence ",
      "threshold placement. If metad = dprime, metacognition is ideal."
    )
    default_priors$metad <- list(
      main = "normal(1, 1)", effects = "normal(0, 0.5)"
    )
    param_links$metad <- "identity"
  }

  # Binary models need mu as internal fixed parameter
  if (!is_rating) {
    parameters <- c(
      list(mu = glue("Internal parameter (fixed to 0)")),
      parameters
    )
    param_links <- c(list(mu = "identity"), param_links)
  }

  # Determine class
  if (is_rating) {
    model_class <- c("bmmodel", "sdt", paste0("sdt_", version, "_rating"))
  } else {
    model_class <- c("bmmodel", "sdt", paste0("sdt_", version))
  }

  # Build requirements text
  if (is_rating) {
    requirements <- glue(
      "Provide pre-aggregated data with the following columns:", "\n\n",
      "  - Response counts: one column per rating category ({n_ratings} columns)", "\n",
      "  - Stimulus type (stimulus): 0 = noise, 1 = signal", "\n",
      "  Categories should be ordered: 1 = 'definitely noise' to ",
      "{n_ratings} = 'definitely signal'"
    )
  } else {
    requirements <- glue(
      "Provide pre-aggregated data with the following columns:", "\n\n",
      "  - Response counts (n_old): number of 'old'/'signal' responses", "\n",
      "  - Stimulus type (stimulus): 0 = noise, 1 = signal", "\n",
      "  - Number of trials (n_trials): total trials per cell"
    )
  }

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(stimulus, n_trials, dist, dist_int, variances,
                         n_ratings, threshold_type),
      domain = "Perception & Recognition Memory",
      task = "Signal/Noise or Old/New Recognition",
      name = "Signal Detection Theory (SDT)",
      citation = glue(
        "Green, D. M., & Swets, J. A. (1966). Signal detection theory ",
        "and psychophysics. Wiley."
      ),
      version = version,
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = if (!is_rating) list(mu = 0) else list(),
      default_priors = default_priors,
      void_mu = is_rating
    ),
    class = model_class,
    call = call
  )
  out$links[names(links)] <- links
  out
}

# Internal: build m-AFC model object (no criterion, no stimulus)
.build_sdt_mafc <- function(response, n_trials, m, dist, dist_int,
                            links, call) {
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
      other_vars = nlist(n_trials, dist, dist_int, m,
                         stimulus = NULL, variances = "equal",
                         n_ratings = NULL, threshold_type = NULL),
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
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdt", "sdt_mafc"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Signal Detection Theory Model
#' @name sdt
#' @details `r model_info(.model_sdt())`
#' @param response For binary SDT: a single string naming the column with
#'   "old"/"signal" response counts. For confidence ratings: a character vector
#'   of K column names containing response counts per category, ordered from
#'   "definitely noise" to "definitely signal".
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (noise/new) and 1 (signal/old).
#' @param n_trials The name of the variable containing the total number of
#'   trials per cell. Required for binary SDT. For ratings, if omitted,
#'   computed as row sums of the response count columns.
#' @param version Character. The version of the model to use:
#'   \itemize{
#'     \item "evsdt" (default): Equal-variance SDT (standard model)
#'     \item "dpsdt": Dual Process SDT (Yonelinas, 1994). Adds a recollection
#'       parameter `Ro` representing threshold-based recollection for old items.
#'       Requires confidence rating data (`n_ratings` > 2).
#'     \item "metad": Meta-d' model (Maniscalco & Lau, 2012). Adds a `metad`
#'       parameter for metacognitive sensitivity. If metad = dprime,
#'       metacognition is ideal. Requires confidence rating data.
#'     \item "mafc": m-Alternative Forced Choice (DeCarlo, 2012). Models
#'       accuracy in m-AFC tasks. Only `dprime` parameter (no criterion).
#'       Requires `m` argument specifying number of alternatives.
#'   }
#' @param dist Character. The noise distribution to use:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT. The sensitivity parameter
#'       `dprime` equals d' = Phi^{-1}(H) - Phi^{-1}(FA), where Phi is the
#'       standard normal CDF.
#'     \item "logistic": Logistic SDT.
#'     \item "gumbel_min": Extreme-value (Gumbel minimum) SDT, theoretically
#'       motivated by recognition memory as event minima (Meyer-Grant et al.,
#'       2025). The sensitivity parameter `dprime` estimates g' =
#'       log(FA) - log(H), which is invariant under uniform expansion of the
#'       choice set — a testable property violated by Gaussian SDT.
#'     \item "gumbel_max": Gumbel maximum SDT (event maxima processes).
#'   }
#'   The parameter `dprime` is always used regardless of distribution, as it
#'   indexes sensitivity in the scale natural to each distribution family.
#' @param variances Character. Variance assumption for signal and noise
#'   distributions:
#'   \itemize{
#'     \item "equal" (default): Equal-variance SDT
#'     \item "unequal": Unequal-variance SDT, adds an `sd_ratio` parameter.
#'       Requires confidence rating data (`n_ratings` > 2) for identification.
#'   }
#' @param n_ratings Integer. Number of response categories for confidence
#'   rating data. When `response` is a character vector, defaults to its length.
#'   Values > 2 trigger the multinomial rating model. Use NULL or omit for
#'   binary SDT.
#' @param threshold_type Character. Threshold parameterization for rating
#'   models:
#'   \itemize{
#'     \item "equidistant" (default): 2 parameters (criterion + spacing).
#'       Thresholds are equally spaced around criterion in decision-variable
#'       space: `theta[k] = criterion + (k - mid) * exp(spacing)`.
#'     \item "parsimonious": 2 parameters (criterion + spacing). Thresholds
#'       follow logit-spaced canonical positions (Selker et al., 2019):
#'       `theta[k] = criterion + exp(spacing) * logit(k/K)`. For an unbiased
#'       observer with d'=0, all K rating categories are used equally often.
#'     \item "log_distance": K-2 parameters (one per non-criterion threshold).
#'       Each threshold distance is exp(delta), guaranteeing ordering.
#'   }
#'   Ignored for binary SDT.
#' @param m Integer. Number of alternatives for m-AFC models. Required when
#'   `version = "mafc"`. Must be >= 2.
#' @param links A named list of link functions for the parameters. The default
#'   link for dprime is "identity" and for criterion is "identity". The dprime
#'   link can be set to "log" to constrain sensitivity to positive values, which
#'   may be useful when including random effects with near-chance performance.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Green, D. M., & Swets, J. A. (1966). \emph{Signal detection theory and
#'   psychophysics}. Wiley.
#'
#' DeCarlo, L. T. (1998). Signal detection theory and generalized linear
#'   models. \emph{Psychological Methods}, \emph{3}(2), 186--205.
#'   \doi{10.1037/1082-989X.3.2.186}
#'
#' Yonelinas, A. P. (1994). Receiver-operating characteristics in recognition
#'   memory: Evidence for a dual-process model. \emph{Journal of Experimental
#'   Psychology: Learning, Memory, and Cognition}, \emph{20}(6), 1341--1354.
#'   \doi{10.1037/0278-7393.20.6.1341}
#'
#' Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach
#'   for estimating metacognitive sensitivity from confidence ratings.
#'   \emph{Consciousness and Cognition}, \emph{21}(1), 422--430.
#'   \doi{10.1016/j.concog.2011.09.021}
#'
#' DeCarlo, L. T. (2012). On a signal detection approach to m-alternative
#'   forced choice with bias, with maximum likelihood and Bayesian approaches
#'   to estimation. \emph{Journal of Mathematical and Statistical Psychology},
#'   \emph{11}(1), 257--282.
#'
#' Selker, R., van den Bergh, D., Criss, A. H., & Wagenmakers, E.-J. (2019).
#'   Parsimonious estimation of signal detection models from confidence ratings.
#'   \emph{Behavior Research Methods}, \emph{51}(5), 1953--1967.
#'   \doi{10.3758/s13428-019-01231-3}
#'
#' Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H. (2025).
#'   \emph{Extreme-value signal detection theory for recognition memory: The
#'   parametric road not taken}. PsyArXiv preprint.
#'   \doi{10.31234/osf.io/qhrfj}
#'
#' Paulewicz, B., & Blaut, A. (2020). The bhsdtr package: A general-purpose
#'   method of Bayesian inference for signal detection theory models.
#'   \emph{Behavior Research Methods}, \emph{52}, 2122--2141.
#'   \doi{10.3758/s13428-020-01370-y}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # --- Binary SDT ---
#' set.seed(42)
#' dat <- rsdt(n_per_cell = 100, n_subjects = 20,
#'             dprime = 1.5, criterion = 0.2)
#'
#' model <- sdt(
#'   response = "n_old",
#'   stimulus = "stimulus",
#'   n_trials = "n_trials",
#'   dist = "normal"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # --- Confidence Rating SDT ---
#' dat_rating <- rsdt(n_per_cell = 200, n_subjects = 20,
#'                    dprime = 1.5, criterion = 0, n_ratings = 4,
#'                    spacing = 0.5, version = "rating")
#'
#' model_rating <- sdt(
#'   response = c("r1", "r2", "r3", "r4"),
#'   stimulus = "stimulus",
#'   dist = "normal",
#'   threshold_type = "equidistant"
#' )
#'
#' fit_rating <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1),
#'   data = dat_rating,
#'   model = model_rating,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # --- Unequal-Variance Rating SDT ---
#' dat_uv <- rsdt(n_per_cell = 200, n_subjects = 20,
#'                dprime = 1.5, criterion = 0, sd_ratio = 1.3,
#'                n_ratings = 4, spacing = 0.5, version = "rating")
#'
#' model_uv <- sdt(
#'   response = c("r1", "r2", "r3", "r4"),
#'   stimulus = "stimulus",
#'   variances = "unequal"
#' )
#'
#' fit_uv <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
#'   data = dat_uv,
#'   model = model_uv,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # --- Dual Process SDT (DPSDT) ---
#' dat_dp <- rsdt(n_per_cell = 200, n_subjects = 20,
#'                dprime = 1.0, criterion = 0, Ro = 0.3,
#'                n_ratings = 4, spacing = 0.5, version = "dpsdt")
#'
#' model_dp <- sdt(
#'   response = c("r1", "r2", "r3", "r4"),
#'   stimulus = "stimulus",
#'   version = "dpsdt"
#' )
#'
#' fit_dp <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1),
#'   data = dat_dp,
#'   model = model_dp,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # --- m-AFC SDT ---
#' dat_mafc <- rsdt(n_per_cell = 200, n_subjects = 20,
#'                  dprime = 1.5, version = "mafc", m = 4)
#'
#' model_mafc <- sdt(
#'   response = "n_correct",
#'   n_trials = "n_trials",
#'   version = "mafc",
#'   m = 4
#' )
#'
#' fit_mafc <- bmm(
#'   formula = bmf(dprime ~ 1),
#'   data = dat_mafc,
#'   model = model_mafc,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' }
sdt <- function(response, stimulus = NULL, n_trials = NULL,
                version = c("evsdt", "dpsdt", "metad", "mafc"),
                dist = c("normal", "logistic", "gumbel_min", "gumbel_max"),
                variances = c("equal", "unequal"),
                n_ratings = NULL,
                threshold_type = c("equidistant", "parsimonious", "log_distance"),
                m = NULL, links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  dist <- match.arg(dist)
  variances <- match.arg(variances)
  threshold_type <- match.arg(threshold_type)

  # m-AFC has different requirements
  if (version == "mafc") {
    stopif(is.null(m) || !is.numeric(m) || length(m) != 1 || m < 2,
           "m must be a single integer >= 2 for m-AFC models")
    m <- as.integer(m)
    stopif(is.null(n_trials),
           "n_trials is required for m-AFC models")
    stopif(dist != "normal",
           "m-AFC currently only supports dist = 'normal'")
    return(.model_sdt(response = response, n_trials = n_trials,
                      version = version, dist = dist, m = m,
                      links = links, call = call, ...))
  }

  # Non-mafc versions require stimulus
  stopif(is.null(stimulus),
         "stimulus is required for version = '{version}'")

  # Infer n_ratings from response vector length
  if (length(response) > 1) {
    if (is.null(n_ratings)) {
      n_ratings <- length(response)
    }
    stopif(
      n_ratings != length(response),
      "n_ratings ({n_ratings}) must match the number of response columns ({length(response)})"
    )
  }

  # Binary SDT requires n_trials; UV-SDT and DPSDT require ratings
  is_rating <- !is.null(n_ratings) && n_ratings > 2
  if (!is_rating) {
    stopif(is.null(n_trials),
           "n_trials is required for binary SDT models")
    stopif(variances == "unequal",
           "Unequal-variance SDT requires confidence rating data (n_ratings > 2)")
    stopif(version == "dpsdt",
           "DPSDT requires confidence rating data (n_ratings > 2)")
    stopif(version == "metad",
           "Meta-d' requires confidence rating data (n_ratings > 2)")
  }

  .model_sdt(response = response, stimulus = stimulus,
             n_trials = n_trials, version = version,
             dist = dist, variances = variances,
             n_ratings = n_ratings, threshold_type = threshold_type,
             links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt <- function(model, data, formula) {
  # m-AFC has no stimulus column
  if (model$version == "mafc") {
    data <- .check_data_sdt_mafc(model, data)
    return(NextMethod("check_data"))
  }

  stim_var <- model$other_vars$stimulus

  # Validate stimulus column exists
  stopif(
    !stim_var %in% colnames(data),
    "Stimulus variable '{stim_var}' missing in the data"
  )

  # Validate stimulus is 0/1
  stim_vals <- unique(data[[stim_var]])
  stopif(
    !all(stim_vals %in% c(0, 1)),
    "Stimulus variable '{stim_var}' must be coded as 0 (noise) and 1 (signal)"
  )

  if (.is_sdt_rating(model)) {
    data <- .check_data_sdt_rating(model, data)
  } else {
    data <- .check_data_sdt_binary(model, data)
  }

  NextMethod("check_data")
}

# Internal: shared validation for count + trial columns
# Checks existence, non-negativity, integerishness, and count <= trials
.validate_sdt_counts <- function(data, resp_var, n_trials_var) {
  required <- c(resp_var, n_trials_var)
  missing <- setdiff(required, colnames(data))
  stopif(length(missing) > 0,
    "Variables {collapse_comma(missing)} missing in the data")

  resp_vals <- data[[resp_var]]
  stopif(any(resp_vals < 0, na.rm = TRUE),
    "Response variable '{resp_var}' must contain non-negative counts")
  warnif(any(resp_vals != round(resp_vals), na.rm = TRUE),
    "Response variable '{resp_var}' should contain integer counts")

  trial_vals <- data[[n_trials_var]]
  stopif(any(trial_vals <= 0, na.rm = TRUE),
    "Variable '{n_trials_var}' must contain positive values")

  stopif(any(resp_vals > trial_vals, na.rm = TRUE),
    "Response counts in '{resp_var}' must not exceed '{n_trials_var}'")
}

# Internal: binary-specific data checks
.check_data_sdt_binary <- function(model, data) {
  resp_var <- model$resp_vars$response
  n_trials_var <- model$other_vars$n_trials

  .validate_sdt_counts(data, resp_var, n_trials_var)

  # Binary-specific: warn on non-integer trial counts
  warnif(any(data[[n_trials_var]] != round(data[[n_trials_var]]), na.rm = TRUE),
    "Variable '{n_trials_var}' should contain integer counts")

  # Add dist_type column for Stan dispatch
  data$dist_type <- model$other_vars$dist_int
  data
}

# Internal: rating-specific data checks
.check_data_sdt_rating <- function(model, data) {
  resp_cols <- model$resp_vars$response

  # Validate all response columns exist
  missing <- setdiff(resp_cols, colnames(data))
  stopif(length(missing) > 0,
    "Response columns {collapse_comma(missing)} missing in the data")

  # Validate all response columns are non-negative
  for (col in resp_cols) {
    vals <- data[[col]]
    stopif(any(vals < 0, na.rm = TRUE),
      "Response column '{col}' must contain non-negative counts")
    warnif(any(vals != round(vals), na.rm = TRUE),
      "Response column '{col}' should contain integer counts")
  }

  # Create response matrix Y and compute nTrials
  Y <- as.matrix(data[resp_cols])
  data$nTrials <- rowSums(Y)
  data$Y <- Y

  stopif(any(data$nTrials <= 0, na.rm = TRUE),
    "Row sums of response columns must be positive (no empty rows)")

  data
}

# Internal: m-AFC-specific data checks
.check_data_sdt_mafc <- function(model, data) {
  resp_var <- model$resp_vars$response
  n_trials_var <- model$other_vars$n_trials

  .validate_sdt_counts(data, resp_var, n_trials_var)

  # Add m as data column for Stan dispatch via vint
  data$m_afc <- model$other_vars$m
  data
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt <- function(model, formula) {
  if (model$version == "mafc") {
    return(.bmf2bf_sdt_mafc(model))
  }
  if (.is_sdt_rating(model)) {
    if (model$version == "dpsdt") {
      return(.bmf2bf_sdt_dpsdt(model))
    }
    if (model$version == "metad") {
      return(.bmf2bf_sdt_metad(model))
    }
    return(.bmf2bf_sdt_rating(model))
  }
  .bmf2bf_sdt_binary(model)
}

# Internal: binary SDT base formula
.bmf2bf_sdt_binary <- function(model) {
  resp_var <- model$resp_vars$response
  stim_var <- model$other_vars$stimulus
  n_trials_var <- model$other_vars$n_trials

  brms::bf(paste0(
    resp_var, " | vint(", stim_var, ", dist_type) + trials(",
    n_trials_var, ") ~ 1"
  ))
}

# Internal: m-AFC SDT base formula (binomial on accuracy)
.bmf2bf_sdt_mafc <- function(model) {
  resp_var <- model$resp_vars$response
  n_trials_var <- model$other_vars$n_trials

  # Pass m via vint()
  brms::bf(paste0(
    resp_var, " | vint(m_afc) + trials(",
    n_trials_var, ") ~ 1"
  ))
}

# Internal: shared infrastructure for rating SDT formula construction
# Extracts model properties and builds the base formula, CDF expression builder,
# threshold expressions, and cumulative probability function.
# shift_expr: the Stan expression for the decision variable shift
# Returns a list with: base_formula, cum_prob, cdf_fn, theta_expr,
#   n_ratings, mid, resp_cols
.sdt_rating_formula_parts <- function(model,
                                      shift_expr = "dprime / 2 * (2 * stimulus - 1)") {
  resp_cols <- model$resp_vars$response
  n_ratings <- model$other_vars$n_ratings
  threshold_type <- model$other_vars$threshold_type
  variances <- model$other_vars$variances
  dist <- model$other_vars$dist
  mid <- n_ratings %/% 2L

  resp_str <- paste(resp_cols, collapse = ", ")
  base_formula <- brms::bf(
    paste0("cbind(", resp_str, ") | trials(nTrials) ~ 1"),
    nl = TRUE
  )

  cdf_fn <- .sdt_cdf_expr(dist)

  theta_expr <- if (threshold_type == "equidistant") {
    .sdt_thresholds_equidistant(n_ratings, mid)
  } else if (threshold_type == "parsimonious") {
    .sdt_thresholds_parsimonious(n_ratings)
  } else {
    .sdt_thresholds_log_distance(n_ratings, mid)
  }

  # Build cumulative probability function: F(theta[k] - shift)
  # UV-SDT: F((theta[k] - shift) / scale)
  cum_prob <- if (variances == "unequal") {
    scale <- "(stimulus * exp(sdratio) + (1 - stimulus))"
    function(k) cdf_fn(paste0("(", theta_expr(k), " - ", shift_expr, ") / ", scale))
  } else {
    function(k) cdf_fn(paste0(theta_expr(k), " - ", shift_expr))
  }

  list(base_formula = base_formula, cum_prob = cum_prob, cdf_fn = cdf_fn,
       theta_expr = theta_expr, n_ratings = n_ratings, mid = mid,
       resp_cols = resp_cols)
}

# Internal: equidistant threshold expressions
# theta[k] = criterion + (k - mid) * exp(spacing)
.sdt_thresholds_equidistant <- function(n_ratings, mid) {
  function(k) {
    offset <- k - mid
    if (offset == 0L) return("criterion")
    if (offset > 0L) {
      return(paste0("criterion + ", offset, " * exp(spacing)"))
    }
    paste0("criterion - ", abs(offset), " * exp(spacing)")
  }
}

# Internal: log-distance threshold expressions
# theta[mid] = criterion
# theta[mid+j] = theta[mid+j-1] + exp(delta_{mid+j})
# theta[mid-j] = theta[mid-j+1] - exp(delta_{mid-j})
.sdt_thresholds_log_distance <- function(n_ratings, mid) {
  K1 <- n_ratings - 1L
  theta_exprs <- character(K1)
  theta_exprs[mid] <- "criterion"

  if (mid < K1) {
    for (k in (mid + 1L):K1) {
      delta_name <- paste0("delta", k)
      theta_exprs[k] <- paste0(theta_exprs[k - 1L], " + exp(", delta_name, ")")
    }
  }

  if (mid > 1L) {
    for (k in (mid - 1L):1L) {
      delta_name <- paste0("delta", k)
      theta_exprs[k] <- paste0(theta_exprs[k + 1L], " - exp(", delta_name, ")")
    }
  }

  function(k) theta_exprs[k]
}

# Internal: parsimonious threshold expressions (Selker et al., 2019)
# theta[k] = criterion + exp(spacing) * logit(k/K)
# where logit(k/K) = log(k / (K - k)) are fixed canonical positions.
# For K even, theta[K/2] = criterion exactly (logit(0.5) = 0).
# For an unbiased observer (d'=0), all K categories are used equally often.
.sdt_thresholds_parsimonious <- function(n_ratings) {
  K <- n_ratings
  K1 <- K - 1L
  canonical <- log(seq_len(K1) / (K - seq_len(K1)))  # logit(k/K)
  function(k) {
    gk <- canonical[k]
    if (gk == 0) return("criterion")
    if (gk > 0) {
      sprintf("criterion + exp(spacing) * %.10f", gk)
    } else {
      sprintf("criterion - exp(spacing) * %.10f", abs(gk))
    }
  }
}

# Internal: standard SDT nlf formulas (log of raw category probabilities)
.sdt_build_nlf_standard <- function(parts) {
  nlf_formulas <- list()
  for (k in seq_len(parts$n_ratings)) {
    dpar <- paste0("mu", parts$resp_cols[k])
    if (k == 1L) {
      expr <- paste0("log(", parts$cum_prob(1), ")")
    } else if (k == parts$n_ratings) {
      expr <- paste0("log(1 - ", parts$cum_prob(parts$n_ratings - 1L), ")")
    } else {
      expr <- paste0("log(", parts$cum_prob(k), " - ", parts$cum_prob(k - 1L), ")")
    }
    nlf_formulas[[k]] <- glue_nlf("{dpar} ~ {expr}")
  }
  nlf_formulas
}

# Internal: DPSDT nlf formulas (familiarity-scaled probs + recollection)
# For old items (stimulus=1): p[k] = (1-Ro)*sdt_prob[k], p[K] += Ro
# For new items (stimulus=0): standard SDT probabilities (unmodified)
# Using: (1 - inv_logit(Ro) * stimulus) handles both cases in one expression
.sdt_build_nlf_dpsdt <- function(parts) {
  fam_scale <- "(1 - inv_logit(Ro) * stimulus)"
  rec_add <- "inv_logit(Ro) * stimulus"

  nlf_formulas <- list()
  for (k in seq_len(parts$n_ratings)) {
    dpar <- paste0("mu", parts$resp_cols[k])
    if (k == 1L) {
      expr <- paste0("log(", fam_scale, " * ", parts$cum_prob(1), ")")
    } else if (k == parts$n_ratings) {
      expr <- paste0(
        "log(", fam_scale, " * (1 - ", parts$cum_prob(parts$n_ratings - 1L), ")",
        " + ", rec_add, ")"
      )
    } else {
      expr <- paste0(
        "log(", fam_scale, " * (",
        parts$cum_prob(k), " - ", parts$cum_prob(k - 1L), "))"
      )
    }
    nlf_formulas[[k]] <- glue_nlf("{dpar} ~ {expr}")
  }
  nlf_formulas
}

# Internal: Meta-d' nlf formulas (normalized meta-d' category probabilities)
# Uses meta-d' for confidence threshold placement, but normalizes so that the
# total "old"/"new" response rates match what type-1 d' predicts.
# For categories below criterion (k <= mid):
#   p[k] = raw_metad_prob[k] * F(crit - d_shift) / F(crit - metad_shift)
# For categories above criterion (k > mid):
#   p[k] = raw_metad_prob[k] * (1-F(crit-d_shift)) / (1-F(crit-metad_shift))
.sdt_build_nlf_metad <- function(parts_metad, parts_dprime) {
  crit_expr <- parts_metad$theta_expr(parts_metad$mid)
  d_shift <- "dprime / 2 * (2 * stimulus - 1)"
  metad_shift <- "metad / 2 * (2 * stimulus - 1)"
  cdf_fn <- parts_metad$cdf_fn

  norm_below <- paste0(
    cdf_fn(paste0(crit_expr, " - ", d_shift)),
    " / ",
    cdf_fn(paste0(crit_expr, " - ", metad_shift))
  )
  norm_above <- paste0(
    "(1 - ", cdf_fn(paste0(crit_expr, " - ", d_shift)), ")",
    " / ",
    "(1 - ", cdf_fn(paste0(crit_expr, " - ", metad_shift)), ")"
  )

  nlf_formulas <- list()
  for (k in seq_len(parts_metad$n_ratings)) {
    dpar <- paste0("mu", parts_metad$resp_cols[k])
    norm <- if (k <= parts_metad$mid) norm_below else norm_above

    if (k == 1L) {
      expr <- paste0("log(", parts_metad$cum_prob(1), " * ", norm, ")")
    } else if (k == parts_metad$n_ratings) {
      expr <- paste0("log((1 - ", parts_metad$cum_prob(parts_metad$n_ratings - 1L),
                      ") * ", norm, ")")
    } else {
      expr <- paste0("log((", parts_metad$cum_prob(k), " - ",
                      parts_metad$cum_prob(k - 1L), ") * ", norm, ")")
    }
    nlf_formulas[[k]] <- glue_nlf("{dpar} ~ {expr}")
  }
  nlf_formulas
}

# Internal: rating SDT formula (standard EV/UV-SDT)
.bmf2bf_sdt_rating <- function(model) {
  parts <- .sdt_rating_formula_parts(model)
  nlf_formulas <- .sdt_build_nlf_standard(parts)
  Reduce(`+`, nlf_formulas, init = parts$base_formula)
}

# Internal: DPSDT rating formula (recollection mixture)
.bmf2bf_sdt_dpsdt <- function(model) {
  parts <- .sdt_rating_formula_parts(model)
  nlf_formulas <- .sdt_build_nlf_dpsdt(parts)
  Reduce(`+`, nlf_formulas, init = parts$base_formula)
}

# Internal: Meta-d' rating formula (normalized meta-d' category probabilities)
.bmf2bf_sdt_metad <- function(model) {
  parts_metad <- .sdt_rating_formula_parts(
    model, shift_expr = "metad / 2 * (2 * stimulus - 1)"
  )
  parts_dprime <- .sdt_rating_formula_parts(model)
  nlf_formulas <- .sdt_build_nlf_metad(parts_metad, parts_dprime)
  Reduce(`+`, nlf_formulas, init = parts_metad$base_formula)
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_evsdt <- function(model, data, formula) {
  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)

  # construct the family & add to formula object
  formula$family <- brms::custom_family(
    "sdt_binary",
    dpars = c("mu", "dprime", "criterion"),
    links = c("identity", model$links$dprime, model$links$criterion),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_sdt_binary,
    posterior_predict = posterior_predict_sdt_binary,
    vars = c("vint1[n]", "vint2[n]", "trials[n]")
  )

  # prepare initial stanvars to pass to brms
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdt_binary_funs.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  # return the list
  nlist(formula, data, stanvars)
}

# Internal: shared rating model configuration (multinomial + nlf)
.configure_sdt_rating <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  resp_cols <- model$resp_vars$response
  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- resp_cols
  formula$family$dpars <- paste0("mu", resp_cols)

  nlist(formula, data)
}

#' @export
configure_model.sdt_evsdt_rating <- function(model, data, formula) {
  .configure_sdt_rating(model, data, formula)
}

#' @export
configure_model.sdt_dpsdt_rating <- function(model, data, formula) {
  .configure_sdt_rating(model, data, formula)
}

#' @export
configure_model.sdt_metad_rating <- function(model, data, formula) {
  .configure_sdt_rating(model, data, formula)
}

#' @export
configure_model.sdt_mafc <- function(model, data, formula) {
  # construct brms formula
  formula <- bmf2bf(model, formula)

  # Define custom family: only dprime (+ mu as internal)
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

  # Load Stan functions
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdt_mafc_funs.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT (Binary SDT)                               ####
############################################################################# !

# Map dist_type integer (from Stan vint) back to distribution name
# Derived from registry: position i corresponds to dist with id == i
.sdt_dist_names <- names(.SDT_DISTS)[order(vapply(.SDT_DISTS, `[[`, 0L, "id"))]

log_lik_sdt_binary <- function(i, prep) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  stimulus <- prep$data$vint1[i]
  dist <- .sdt_dist_names[prep$data$vint2[i]]
  n_trials <- prep$data$trials[i]
  y <- prep$data$Y[i]

  eta <- .sdt_eta(dprime, criterion, stimulus)
  p <- .sdt_cdf(eta, dist)
  dbinom(y, n_trials, p, log = TRUE)
}

posterior_predict_sdt_binary <- function(i, prep, ...) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  stimulus <- prep$data$vint1[i]
  dist <- .sdt_dist_names[prep$data$vint2[i]]
  n_trials <- prep$data$trials[i]

  eta <- .sdt_eta(dprime, criterion, stimulus)
  p <- .sdt_cdf(eta, dist)
  rbinom(length(p), n_trials, p)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT (m-AFC SDT)                                ####
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

  sum(gh_weights * pnorm(gh_nodes + dprime)^(m - 1))
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


############################################################################# !
# HELPER FUNCTIONS                                                       ####
############################################################################# !

#' Combine stimulus type and confidence into SDT response categories
#'
#' Creates a combined response variable from separate stimulus and confidence
#' columns, suitable for use with SDT rating models. The combined response
#' maps "noise" trials to categories 1..K/2 (from highest to lowest confidence)
#' and "signal" trials to categories K/2+1..K (from lowest to highest).
#'
#' @param stimulus Integer vector (0/1). Stimulus type: 0 = noise, 1 = signal.
#' @param confidence Integer vector. Confidence rating (1 = lowest confidence,
#'   n_levels = highest confidence).
#' @param n_levels Integer. Number of confidence levels per stimulus type.
#'
#' @return Integer vector of combined response categories (1 to 2*n_levels).
#' @export
#' @examples
#' # 3 confidence levels per stimulus type -> K=6 combined categories
#' stim <- c(0, 0, 0, 1, 1, 1)
#' conf <- c(3, 2, 1, 1, 2, 3)
#' combine_sdt_response(stim, conf, n_levels = 3)
#' # Returns: 1, 2, 3, 4, 5, 6
combine_sdt_response <- function(stimulus, confidence, n_levels) {
  stopif(any(!stimulus %in% c(0, 1)),
         "stimulus must be 0 (noise) or 1 (signal)")
  stopif(any(confidence < 1 | confidence > n_levels),
         "confidence must be between 1 and {n_levels}")

  # Noise: high conf = cat 1, low conf = cat n_levels
  # Signal: low conf = cat n_levels+1, high conf = cat 2*n_levels
  ifelse(stimulus == 0,
         n_levels - confidence + 1L,
         n_levels + confidence)
}
