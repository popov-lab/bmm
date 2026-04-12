############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_rating <- function(response = NULL, stimulus = NULL,
                              dist = "normal", n_ratings = NULL,
                              threshold_type = "parsimonious",
                              log_scale = FALSE,
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

  # sdratio as overridable fixed parameter (default = equal variance)
  parameters$sdratio <- glue(
    "Log SD ratio: log ratio of signal to noise standard deviations ",
    "(exp(sdratio) is the natural SD ratio, 0 = equal variance)"
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
                         threshold_type, thresh_type_int, log_scale),
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
      void_mu = TRUE
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
#' @param dist Character. The noise distribution to use:
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
#' @param log_scale Logical. Deprecated and ignored. Confidence-rating SDT
#'   likelihoods now use stable log-scale computations internally.
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
                       log_scale = FALSE,
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
  warnif(isTRUE(log_scale),
    paste0(
      "`log_scale` is deprecated and ignored; rating SDT likelihoods now use ",
      "stable log-scale computations internally."
    )
  )

  .model_sdt_rating(response = response, stimulus = stimulus,
                    dist = dist, n_ratings = n_ratings,
                    threshold_type = threshold_type,
                    log_scale = log_scale,
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
bmf2bf.sdt_rating <- function(model, formula) {
  .sdt_rating_base_formula()
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_rating <- function(model, data, formula) {
  formula <- .sdt_remap_formula_dpars(formula, model)
  formula <- bmf2bf(model, formula)
  formula$family <- .sdt_rating_custom_family(
    model,
    family_name = "sdt_rating",
    log_lik = log_lik_sdt_rating,
    posterior_predict = posterior_predict_sdt_rating
  )
  stanvars <- .sdt_rating_stanvars(model, "sdt_rating", "sdt_rating_row_lpmf")

  nlist(formula, data, stanvars)
}

log_lik_sdt_rating <- function(i, prep) {
  .sdt_loglik_rating_common(i, prep, variant = "rating")
}

posterior_predict_sdt_rating <- function(i, prep, ...) {
  .sdt_posterior_predict_rating_common(i, prep, variant = "rating")
}
