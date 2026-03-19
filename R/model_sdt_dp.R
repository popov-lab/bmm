############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_dp <- function(response = NULL, stimulus = NULL,
                          dist = "normal", n_ratings = NULL,
                          threshold_type = "parsimonious",
                          links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)

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
    mid <- n_ratings %/% 2L
    for (i in seq_len(n_deltas)) {
      idx <- if (i < mid) i else i + 1L
      pname <- paste0("delta", idx)
      parameters[[pname]] <- glue("Threshold parameter for threshold {idx}")
      default_priors[[pname]] <- list(
        main = "normal(0, 1)", effects = "normal(0, 0.5)"
      )
      param_links[[pname]] <- "identity"
    }
  }

  parameters$Ro <- glue(
    "Recollection: probability of threshold-based recollection ",
    "for old items (inv_logit(Ro) maps to [0,1])"
  )
  default_priors$Ro <- list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
  param_links$Ro <- "identity"

  parameters$Rn <- glue(
    "Recollection: probability of threshold-based recollection ",
    "for new items (inv_logit(Rn) maps to [0,1])"
  )
  default_priors$Rn <- list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
  param_links$Rn <- "identity"

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
                         threshold_type),
      domain = "Perception & Recognition Memory",
      task = "Signal/Noise or Old/New Recognition",
      name = "Dual Process Signal Detection Theory (DPSDT)",
      citation = glue(
        "Yonelinas, A. P. (1994). Receiver-operating characteristics in ",
        "recognition memory: Evidence for a dual-process model. JEPLMC, ",
        "20(6), 1341-1354."
      ),
      version = "dp",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = list(sdratio = 0),
      default_priors = default_priors,
      void_mu = TRUE
    ),
    class = c("bmmodel", "sdt", "sdt_dp"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Dual Process Signal Detection Theory Model
#' @name sdt_dp
#' @details `r model_info(.model_sdt_dp())`
#'
#' The DPSDT model adds recollection parameters to the standard SDT framework.
#' Recollection for old items (`Ro`) adds mass to the highest confidence "old"
#' category. Recollection for new items (`Rn`) adds mass to the highest
#' confidence "new" category.
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
#' Yonelinas, A. P. (1994). Receiver-operating characteristics in recognition
#'   memory: Evidence for a dual-process model. \emph{Journal of Experimental
#'   Psychology: Learning, Memory, and Cognition}, \emph{20}(6), 1341--1354.
#'   \doi{10.1037/0278-7393.20.6.1341}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' dat <- rsdt_dp(n_per_cell = 200, n_subjects = 20,
#'                dprime = 1.0, criterion = 0, Ro = 0.3, Rn = 0.1,
#'                n_ratings = 4, spacing = 0.5)
#'
#' model <- sdt_dp(
#'   response = c("r1", "r2", "r3", "r4"),
#'   stimulus = "stimulus"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#' }
sdt_dp <- function(response, stimulus,
                   n_ratings = NULL,
                   dist = c("normal", "logistic", "gumbel_min", "gumbel_max"),
                   threshold_type = c("parsimonious", "equidistant",
                                      "log_distance", "log_ratio"),
                   links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)
  threshold_type <- match.arg(threshold_type)

  if (length(response) > 1 && is.null(n_ratings)) {
    n_ratings <- length(response)
  }

  stopif(is.null(n_ratings) || n_ratings <= 2,
         "DPSDT requires confidence rating data (n_ratings > 2)")

  if (length(response) > 1) {
    stopif(n_ratings != length(response),
           "n_ratings ({n_ratings}) must match the number of response columns ({length(response)})")
  }

  .model_sdt_dp(response = response, stimulus = stimulus,
                dist = dist, n_ratings = n_ratings,
                threshold_type = threshold_type,
                links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_dp <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stim_vals <- unique(data[[stim_var]])
  stopif(!all(stim_vals %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (noise) and 1 (signal)")

  data <- .check_data_sdt_rating(model, data)

  NextMethod("check_data")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt_dp <- function(model, formula) {
  parts <- .sdt_rating_formula_parts(model)
  nlf_formulas <- .sdt_build_nlf_dpsdt(parts)
  Reduce(`+`, nlf_formulas, init = parts$base_formula)
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_dp <- function(model, data, formula) {
  .configure_sdt_rating(model, data, formula)
}
