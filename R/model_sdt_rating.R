############################################################################# !
# CONFIDENCE-THRESHOLD PARAMETERIZATIONS                                  ####
############################################################################# !

.sdt_threshold_type_ids <- c(
  parsimonious = 1L,
  equidistant = 2L,
  log_distance = 3L,
  log_ratio = 4L,
  softmax = 5L
)

.sdt_threshold_type_id <- function(threshold_type) {
  .sdt_threshold_type_ids[[threshold_type]]
}

.sdt_threshold_type_name <- function(thresh_type) {
  names(.sdt_threshold_type_ids)[match(thresh_type, .sdt_threshold_type_ids)]
}

# Parameter spec for one threshold parameterization. parsimonious/equidistant
# need a single spacing; log_distance/log_ratio need K-2 distance deltas;
# softmax needs spacing plus K-3 allocation deltas.
.sdt_threshold_parameter_parts <- function(n_ratings, threshold_type) {
  parameters <- list()
  default_priors <- list()
  param_links <- list()

  if (threshold_type %in% c("equidistant", "parsimonious")) {
    parameters$spacing <- glue(
      "Threshold spacing: controls distance between adjacent thresholds ",
      "(exp(spacing) ensures positive spacing)"
    )
    default_priors$spacing <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
    )
    param_links$spacing <- "identity"
  } else if (threshold_type %in% c("log_distance", "log_ratio")) {
    n_deltas <- n_ratings - 2L
    # skip the same middle threshold as the builders ((K - 1) %/% 2 + 1), so
    # the delta labels name the threshold they control for odd K
    mid <- (n_ratings - 1L) %/% 2L + 1L
    for (i in seq_len(n_deltas)) {
      idx <- if (i < mid) i else i + 1L
      pname <- paste0("delta", idx)
      parameters[[pname]] <- glue("Threshold parameter for threshold {idx}")
      default_priors[[pname]] <- list(
        main = "normal(0, 1)", effects = "normal(0, 0.5)"
      )
      param_links[[pname]] <- "identity"
    }
  } else if (threshold_type == "softmax") {
    parameters$spacing <- glue(
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

.model_sdt_rating <- function(response = NULL, stimulus = NULL,
                              dist = "normal", n_ratings = NULL,
                              threshold_type = "parsimonious",
                              links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)
  thresh_type_int <- .sdt_threshold_type_id(threshold_type)

  if (is.null(n_ratings) && length(response) > 1) {
    n_ratings <- length(response)
  }

  parameters <- list(
    dprime = glue("Sensitivity: distance between signal and noise distributions"),
    criterion = glue("Response bias: location of decision boundary")
  )
  default_priors <- list(
    dprime = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
    criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)")
  )
  param_links <- list(dprime = "identity", criterion = "identity")

  threshold_parts <- .sdt_threshold_parameter_parts(n_ratings, threshold_type)
  parameters <- c(parameters, threshold_parts$parameters)
  default_priors <- c(default_priors, threshold_parts$default_priors)
  param_links <- c(param_links, threshold_parts$param_links)

  parameters$sdratio <- glue(
    "Log SD ratio: the log of the signal-to-noise standard deviation ratio, ",
    "so exp(sdratio) is the ratio itself and 0 means equal variance"
  )
  default_priors$sdratio <- list(
    main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
  )
  param_links$sdratio <- "identity"

  # Bounded starting values for every estimated parameter. Without them the
  # flexible threshold types start from brms' wide random init, where adjacent
  # thresholds can land in the same distribution tail and the category mass
  # underflows to log_diff_exp(-Inf, -Inf) = NaN, so every chain rejects its
  # initial value.
  init_ranges <- list(
    dprime = c(0.5, 1.5), criterion = c(-0.5, 0.5), sdratio = c(-0.3, 0.3)
  )
  for (p in names(threshold_parts$parameters)) {
    init_ranges[[p]] <- if (p == "spacing") c(-0.7, -0.2) else c(-0.5, 0.2)
  }

  requirements <- glue(
    "Provide pre-aggregated data with the following columns:", "\n\n",
    "  - Response counts: one column per rating category ({n_ratings} columns)", "\n",
    "  - Stimulus type (stimulus): 0 = noise, 1 = signal", "\n",
    "  Categories should be ordered: 1 = 'definitely noise' to ",
    "{n_ratings} = 'definitely signal'"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(stimulus, dist, dist_int, n_ratings,
                         threshold_type, thresh_type_int),
      domain = "Perception & Recognition Memory",
      task = "Signal/Noise or Old/New Recognition",
      name = "Signal Detection Theory (Confidence Rating)",
      citation = glue(
        "Green, D. M., & Swets, J. A. (1966). Signal detection theory ",
        "and psychophysics. Wiley."
      ),
      version = "rating",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = list(sdratio = 0),
      default_priors = default_priors,
      init_ranges = init_ranges,
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdt", "sdt_rating"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Confidence Rating Signal Detection Theory Model
#' @name sdt_rating
#' @details `r model_info(.model_sdt_rating())`
#'
#' By default, the model assumes equal variance (sdratio fixed to 0). To
#' estimate unequal variance, add `sdratio ~ 1` (or `sdratio ~ predictors`)
#' to the formula.
#' @param response A character vector of K column names containing response
#'   counts per rating category, ordered from "definitely noise" to
#'   "definitely signal".
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (noise/new) and 1 (signal/old).
#' @param n_ratings Integer. Number of response categories. When `response`
#'   is a character vector, defaults to its length. Must be > 2.
#' @param dist The noise distribution assumed for the latent evidence variable.
#'   One of:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT
#'     \item "logistic": Logistic SDT
#'     \item "gumbel_min": Extreme-value (Gumbel minimum) SDT
#'     \item "gumbel_max": Gumbel maximum SDT
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
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
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
#' dat <- rsdt_rating(n_per_cell = 200, n_subjects = 20,
#'                    dprime = 1.5, criterion = 0, n_ratings = 4,
#'                    spacing = 0.5)
#'
#' model <- sdt_rating(
#'   response = c("r1", "r2", "r3", "r4"),
#'   stimulus = "stimulus"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # UV-SDT: add sdratio to the formula
#' fit_uv <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' }
sdt_rating <- function(response, stimulus,
                       n_ratings = NULL,
                       dist = c("normal", "logistic", "gumbel_min", "gumbel_max"),
                       threshold_type = c("parsimonious", "equidistant",
                                          "log_distance", "log_ratio",
                                          "softmax"),
                       links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)
  threshold_type <- match.arg(threshold_type)

  if (length(response) > 1 && is.null(n_ratings)) {
    n_ratings <- length(response)
  }

  stopif(is.null(n_ratings) || n_ratings <= 2,
         "Rating models require n_ratings > 2 (or pass a response vector with > 2 columns)")

  if (length(response) > 1) {
    stopif(n_ratings != length(response),
           "n_ratings ({n_ratings}) must match the number of response columns ({length(response)})")
  }

  .model_sdt_rating(response = response, stimulus = stimulus,
                    dist = dist, n_ratings = n_ratings,
                    threshold_type = threshold_type,
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
  stopif(!all(unique(data[[stim_var]]) %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (noise) and 1 (signal)")

  resp_cols <- model$resp_vars$response
  missing <- setdiff(resp_cols, colnames(data))
  stopif(length(missing) > 0,
         "Response columns {collapse_comma(missing)} missing in the data")

  for (col in resp_cols) {
    vals <- data[[col]]
    stopif(any(vals < 0, na.rm = TRUE),
           "Response column '{col}' must contain non-negative counts")
    warnif(any(vals != round(vals), na.rm = TRUE),
           "Response column '{col}' should contain integer counts")
  }

  Y <- as.matrix(data[resp_cols])
  stopif(any(rowSums(Y) <= 0, na.rm = TRUE),
         "Row sums of response columns must be positive (no empty rows)")

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

# cat, K, dist, and threshold type travel as integer literals so the same
# non-linear call resolves against both the generated Stan function and the R
# companion. spacing is the parameter when present, otherwise the literal 0
# (threshold types without spacing ignore it).
.sdt_rating_logmu_args <- function(model) {
  has_spacing <- "spacing" %in% names(model$parameters)
  c(model$other_vars$n_ratings, model$other_vars$dist_int,
    model$other_vars$thresh_type_int, "dprime", "criterion",
    if (has_spacing) "spacing" else "0", "sdratio",
    model$other_vars$stimulus, .sdt_threshold_delta_names(model))
}

# Stan function generated per model: the delta arity is fixed here, while cat,
# K, dist, and threshold type arrive as the integer arguments described above.
.sdt_rating_logmu_stan <- function(model) {
  delta_names <- .sdt_threshold_delta_names(model)
  nd <- length(delta_names)

  signature <- paste(
    c("int cat", "int K", "int dist_type", "int thresh_type",
      "real dprime", "real criterion", "real spacing", "real sdratio",
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
    "real sdt_rating_logmu(", signature, ") {\n",
    delta_decl,
    "  vector[K - 1] thr = sdt_make_thresholds_rating(criterion, spacing, deltas, K, thresh_type);\n",
    "  return sdt_rating_logmu_cat(cat, thr, dprime, sdratio, stimulus, dist_type);\n",
    "}\n"
  )
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

# Base multinomial brmsformula: Y | trials(nTrials) carries the counts and each
# category gets a non-linear mu calling sdt_rating_logmu. The user's parameter
# formulas (dprime, criterion, spacing/deltas, sdratio) are added afterwards by
# bmf2bf.bmmodel, so they are deliberately not added here.
#' @export
bmf2bf.sdt_rating <- function(model, formula) {
  resp_cats <- model$resp_vars$response
  args <- paste(.sdt_rating_logmu_args(model), collapse = ", ")

  bform <- brms::bf(
    glue("Y | trials(nTrials) ~ sdt_rating_logmu(1, {args})"),
    nl = TRUE
  )
  for (k in seq_along(resp_cats)[-1]) {
    bform <- bform + brms::nlf(stats::as.formula(
      glue("mu{resp_cats[k]} ~ sdt_rating_logmu({k}, {args})")
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
  stan_funs <- paste(
    read_lines2(paste0(sc_path, "/sdt_dist_funs.stan")),
    read_lines2(paste0(sc_path, "/sdt_rating_funs.stan")),
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
#' @param dist Integer noise-distribution id (see the `.SDT_DISTS` table).
#' @param thresh Integer threshold-parameterization id.
#' @param dprime,criterion,spacing,sdratio Model parameters (draws-by-observation
#'   matrices supplied by brms). `spacing` is `0` for threshold types without it.
#' @param stimulus Stimulus covariate (0 = noise, 1 = signal).
#' @param ... Threshold `delta` parameters, when the threshold type uses them.
#' @return `log(p_cat)`, matching the shape of `dprime`.
#' @keywords internal
#' @export
sdt_rating_logmu <- function(cat, K, dist, thresh, dprime, criterion, spacing,
                             sdratio, stimulus, ...) {
  dist_name <- .sdt_dist_names[dist]
  thresh_name <- .sdt_threshold_type_name(thresh)

  shape <- dim(dprime)
  dprime <- as.vector(dprime)
  n <- length(dprime)
  criterion <- as.vector(criterion)
  spacing <- rep_len(as.vector(spacing), n)
  sdratio <- as.vector(sdratio)
  stimulus <- rep_len(as.vector(stimulus), n)
  deltas <- lapply(list(...), as.vector)

  out <- vapply(seq_len(n), function(j) {
    deltas_j <- if (length(deltas)) vapply(deltas, `[`, numeric(1), j) else NULL
    thr <- .sdt_make_thresholds(criterion[j], K, thresh_name, spacing[j], deltas_j)
    ratio <- if (stimulus[j] > 0.5) exp(sdratio[j]) else 1
    log(.sdt_category_probs(thr, dprime[j], ratio, stimulus[j], dist_name)[cat])
  }, numeric(1))

  if (!is.null(shape)) dim(out) <- shape
  out
}
