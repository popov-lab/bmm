############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_metad <- function(response = NULL, stimulus = NULL,
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

  parameters$metad <- glue(
    "Metacognitive sensitivity: type-2 d' controlling confidence ",
    "threshold placement. If metad = dprime, metacognition is ideal."
  )
  default_priors$metad <- list(
    main = "normal(1, 1)", effects = "normal(0, 0.5)"
  )
  param_links$metad <- "identity"

  # sdratio as overridable fixed parameter (default = equal variance)
  parameters$sdratio <- glue(
    "SD ratio: log ratio of signal to noise standard deviations ",
    "(exp(sdratio) ensures positivity, 0 = equal variance)"
  )
  default_priors$sdratio <- list(
    main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
  )
  param_links$sdratio <- "identity"

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
      name = "Meta-d' Signal Detection Theory Model",
      citation = glue(
        "Maniscalco, B., & Lau, H. (2012). A signal detection theoretic ",
        "approach for estimating metacognitive sensitivity from confidence ",
        "ratings. Consciousness and Cognition, 21(1), 422-430."
      ),
      version = "metad",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = list(sdratio = 0),
      default_priors = default_priors,
      void_mu = TRUE
    ),
    class = c("bmmodel", "sdt", "sdt_metad"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Meta-d' Signal Detection Theory Model
#' @name sdt_metad
#' @details `r model_info(.model_sdt_metad())`
#'
#' The meta-d' model adds a `metad` parameter for metacognitive sensitivity.
#' If metad equals dprime, metacognition is ideal. Requires confidence
#' rating data.
#'
#' By default, the model assumes equal variance (sdratio fixed to 0). To
#' estimate unequal variance, add `sdratio ~ 1` to the formula.
#' @param response A character vector of K column names containing response
#'   counts per rating category, ordered from "definitely noise" to
#'   "definitely signal". Requires K > 2.
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (noise/new) and 1 (signal/old).
#' @param n_ratings Integer. Number of response categories. Defaults to
#'   `length(response)`. Must be > 2.
#' @param dist Character. The noise distribution to use:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT
#'     \item "logistic": Logistic SDT
#'     \item "gumbel_min": Extreme-value (Gumbel minimum) SDT
#'     \item "gumbel_max": Gumbel maximum SDT
#'   }
#' @param threshold_type Character. Threshold parameterization (see
#'   \code{\link{sdt_rating}} for details).
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach
#'   for estimating metacognitive sensitivity from confidence ratings.
#'   \emph{Consciousness and Cognition}, \emph{21}(1), 422--430.
#'   \doi{10.1016/j.concog.2011.09.021}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' dat <- rsdt_metad(n_per_cell = 200, n_subjects = 20,
#'                   dprime = 1.5, criterion = 0, metad = 1.0,
#'                   n_ratings = 4, spacing = 0.5)
#'
#' model <- sdt_metad(
#'   response = c("r1", "r2", "r3", "r4"),
#'   stimulus = "stimulus"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, metad ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' }
sdt_metad <- function(response, stimulus,
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
         "Meta-d' requires confidence rating data (n_ratings > 2)")

  if (length(response) > 1) {
    stopif(n_ratings != length(response),
           "n_ratings ({n_ratings}) must match the number of response columns ({length(response)})")
  }

  .model_sdt_metad(response = response, stimulus = stimulus,
                   dist = dist, n_ratings = n_ratings,
                   threshold_type = threshold_type,
                   links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_metad <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stim_vals <- unique(data[[stim_var]])
  stopif(!all(stim_vals %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (noise) and 1 (signal)")

  data <- .check_data_sdt_rating(model, data)
  data <- .sdt_rating_long_data(model, data)

  NextMethod("check_data")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt_metad <- function(model, formula) {
  .sdt_rating_base_formula()
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_metad <- function(model, data, formula) {
  formula <- .sdt_remap_formula_dpars(formula, model)
  formula <- bmf2bf(model, formula)
  formula$family <- .sdt_rating_custom_family(
    model,
    family_name = "sdt_metad",
    log_lik = log_lik_sdt_metad,
    posterior_predict = posterior_predict_sdt_metad
  )
  stanvars <- .sdt_rating_stanvars(model, "sdt_metad", "sdt_metad_row_lpmf")

  nlist(formula, data, stanvars)
}

log_lik_sdt_metad <- function(i, prep) {
  .sdt_loglik_rating_common(i, prep, variant = "metad")
}

posterior_predict_sdt_metad <- function(i, prep, ...) {
  .sdt_posterior_predict_rating_common(i, prep, variant = "metad")
}
