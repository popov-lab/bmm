############################################################################# !
# EZDM SUMMARY STATISTICS                                                 ####
############################################################################# !

.validate_contaminant_bounds <- function(contaminant_bound) {
  stopif(length(contaminant_bound) != 2, "contaminant_bound must be a vector of length 2")
  is_valid_bound <- function(x) {
    tolower(x) %in% c("min", "max") | suppressWarnings(!is.na(as.numeric(x)))
  }
  stopif(
    any(!is_valid_bound(contaminant_bound)),
    "contaminant_bound elements must be numeric or 'min'/'max'"
  )
  stopif(
    all(is.numeric(contaminant_bound)) && contaminant_bound[1] >= contaminant_bound[2],
    "contaminant_bound[1] must be less than contaminant_bound[2]"
  )
}

.validate_contaminant_params <- function(init_contaminant, max_contaminant) {
  stopif(
    !is.numeric(init_contaminant) || init_contaminant <= 0 || init_contaminant >= 1,
    "init_contaminant must be between 0 and 1 (exclusive)"
  )
  stopif(
    !is.numeric(max_contaminant) || max_contaminant <= 0 || max_contaminant > 1,
    "max_contaminant must be between 0 (exclusive) and 1 (inclusive)"
  )
  stopif(
    init_contaminant >= max_contaminant,
    "init_contaminant must be less than max_contaminant"
  )
}

#' Compute Robust Summary Statistics for EZ-Diffusion Model
#'
#' @description Computes robust summary statistics for the EZ-Diffusion Model
#'   by fitting mixture models to raw trial-level RT data, separating
#'   contaminant responses from true responses.
#'
#' @param rt Numeric vector of reaction times in seconds.
#' @param response Vector of response indicators. Accepts multiple formats:
#'   \itemize{
#'     \item Numeric: 1 = upper/correct, 0 = lower/error
#'     \item Logical: TRUE = upper/correct, FALSE = lower/error
#'     \item Character/Factor: "upper"/"lower", "correct"/"error",
#'       "acc"/"err", "hit"/"miss", "yes"/"no" (case-insensitive)
#'   }
#' @param version Character. Either "3par" (default) for pooled RTs or "4par"
#'   for separate upper/lower boundary RTs. Controls the output columns.
#' @param distribution Character. The parametric distribution for the RT
#'   component. One of "exgaussian" (default), "lognormal", or "invgaussian"
#' @param method Character. One of "mixture" (default) for robust estimation
#'   via mixture modeling, "robust" for non-parametric robust estimation using
#'   median and IQR/MAD-based variance, or "simple" for standard moment
#'   calculation. The "robust" method is faster and requires no distributional
#'   assumptions, but its moments are not the mean and variance the EZ
#'   equations were derived for. For right-skewed RTs the median lies below the
#'   mean and (IQR/1.349)^2 below the variance, so the EZ estimates are offset
#'   systematically, not just noisier: a larger drift rate and a smaller
#'   boundary separation than the same data give with "mixture".
#' @param robust_scale Character. Scale estimator for robust method. Either
#'   "iqr" (default) for IQR-based variance estimation (variance = (IQR/1.349)^2)
#'   or "mad" for MAD-based estimation (variance = MAD^2, where MAD is scaled
#'   to be consistent with SD for normal data). Only used when method = "robust".
#' @param contaminant_bound Vector of length 2 specifying the bounds (in
#'   seconds) for the uniform contaminant distribution. Can be numeric values
#'   or the special strings "min" and "max" to use data-driven bounds (default):
#'   \itemize{
#'     \item "min": Use the minimum RT in each group, minus a 50\% buffer
#'     \item "max": Use the maximum RT in each group, plus a 50\% buffer
#'     \item Numeric: Fixed bounds, e.g., c(0.1, 3.0)
#'   }
#'   The buffer extends data-driven bounds to ensure conservative estimates.
#'   Examples: c(0.1, 3.0), c("min", "max"), c(0.1, "max"), c("min", 3.0)
#' @param min_trials Integer. Minimum number of trials required for the RT
#'   summaries, which are NA with fewer. Compared against the total number of
#'   trials for `version = "3par"` and against each boundary's own count for
#'   `version = "4par"`, so the two versions can disagree about whether the
#'   same cell is corrected. For `version = "4par"`, [bmm()] still fits a cell
#'   whose summaries at one or both boundaries are NA, through its response
#'   counts. Default is 10
#' @param init_contaminant Numeric. Initial proportion of contaminants for EM
#'   algorithm. Default is 0.05
#' @param max_contaminant Numeric. Maximum allowed contaminant proportion
#'   (0 < max <= 1). Estimates are clipped to this value to prevent inflated
#'   contaminant proportions. Default is 0.5
#' @param maxit Integer. Maximum number of EM iterations. Default is 500. The
#'   EM converges slowly when the contaminant proportion is close to zero, and
#'   a cell that reaches `maxit` without converging falls back to robust
#'   moments with a warning.
#' @param tol Numeric. Convergence tolerance for EM algorithm. Default is 1e-6
#' @param guess_rate Numeric. Accuracy expected of a contaminant response, used
#'   to split the estimated contaminants of `version = "3par"` across the two
#'   boundaries. Default is 0.5, appropriate for a two-choice task; use 0.25 for
#'   a 4AFC task, or 0 if contaminants are never correct. Ignored for
#'   `version = "4par"`, which estimates contamination separately per boundary
#'   and needs no such assumption.
#'
#' @return A 1-row `data.frame`. For version = "3par": `mean_rt`, `var_rt`,
#'   `n_upper`, `n_trials`, `contaminant_prop`. For version = "4par":
#'   `mean_rt_upper`, `mean_rt_lower`, `var_rt_upper`, `var_rt_lower`,
#'   `n_upper`, `n_trials`, `contaminant_prop_upper`, `contaminant_prop_lower`.
#'
#'   `n_upper` and `n_trials` count the responses the reported moments rest on.
#'   Whenever a contaminant proportion was estimated, the expected number of
#'   contaminants is removed from both counts and the result is rounded; with
#'   `method = "simple"`, `method = "robust"`, or an EM that did not converge,
#'   they are the raw counts.
#'
#' @details RT outliers and contaminant responses (fast guesses, lapses of
#'   attention) can distort the mean and variance estimates used as input to
#'   the EZ-Diffusion equations. This function addresses this by fitting a
#'   mixture model with two components: a uniform distribution for
#'   contaminants and a parametric RT distribution for true responses.
#'   Robust moments are then extracted from the fitted parametric component.
#'
#'   The returned counts describe the same responses as the returned moments.
#'   For `version = "3par"` the estimated contaminants are shared between the
#'   boundaries at `guess_rate`; a cell whose observed accuracy lies outside
#'   `[guess_rate * p, 1 - p * (1 - guess_rate)]` for an estimated proportion
#'   `p` cannot have arisen that way, so the contaminants are removed
#'   proportionally from both boundaries instead, leaving the observed accuracy
#'   unchanged apart from rounding, and a warning is issued. For
#'   `version = "4par"` each boundary is corrected by its own estimate.
#'
#'   This function is designed to work with [dplyr::group_by()] and
#'   [dplyr::reframe()] for grouped operations.
#'
#' @seealso [flag_contaminant_rts()] for trial-level contamination
#'   probabilities, [ezdm()] for fitting the EZ-Diffusion Model
#'
#' @keywords transform
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' rt <- rgamma(100, shape = 5, rate = 10) + 0.3
#' response <- rbinom(100, 1, 0.8)
#'
#' # 3par summary stats
#' ezdm_summary_stats(rt, response)
#'
#' # With dplyr for grouped operations
#' # library(dplyr)
#' # mydata |>
#' #   group_by(subject) |>
#' #   reframe(ezdm_summary_stats(rt, response))
#'
#' # 4par version with separate upper/lower moments
#' ezdm_summary_stats(rt, response, version = "4par")
#'
ezdm_summary_stats <- function(
    rt,
    response,
    version = c("3par", "4par"),
    distribution = c("exgaussian", "lognormal", "invgaussian"),
    method = c("mixture", "simple", "robust"),
    robust_scale = c("iqr", "mad"),
    contaminant_bound = c("min", "max"),
    min_trials = 10,
    init_contaminant = 0.05,
    max_contaminant = 0.5,
    maxit = 500,
    tol = 1e-6,
    guess_rate = 0.5) {
  stop_missing_args()
  version <- match.arg(version)
  distribution <- match.arg(distribution)
  method <- match.arg(method)
  robust_scale <- match.arg(robust_scale)

  stopif(!is.numeric(rt), "Argument 'rt' must be a numeric vector")
  stopif(length(rt) == 0L, "Argument 'rt' has length 0")
  stopif(length(response) != length(rt), "Arguments 'rt' and 'response' must have the same length")
  .validate_contaminant_bounds(contaminant_bound)
  stopif(!is.numeric(min_trials) || min_trials < 1, "min_trials must be a positive integer")
  .validate_contaminant_params(init_contaminant, max_contaminant)
  stopif(
    !is.numeric(guess_rate) || length(guess_rate) != 1L || is.na(guess_rate) ||
      guess_rate < 0 || guess_rate > 1,
    "guess_rate must be a single number between 0 and 1"
  )

  complete <- !is.na(rt)
  rt <- rt[complete]
  response <- response[complete]
  warnif(any(rt > 10), "Some RT values > 10. Ensure RTs are in seconds, not milliseconds.")
  stopif(any(rt <= 0), "Non-positive RT values found.")

  is_upper <- .convert_response_to_upper(response)
  n_upper <- sum(is_upper, na.rm = TRUE)
  n_trials <- length(rt)
  resolved_bounds <- .resolve_contaminant_bounds(contaminant_bound, rt)

  compute_moments <- function(rt_vec, n_obs) {
    if (n_obs < min_trials) {
      list(mean = NA_real_, var = NA_real_, contaminant_prop = NA_real_)
    } else if (method == "simple") {
      c(.simple_aggregation(rt_vec), contaminant_prop = NA_real_)
    } else if (method == "robust") {
      c(.robust_aggregation(rt_vec, scale_method = robust_scale),
        contaminant_prop = NA_real_)
    } else {
      fit <- .fit_rt_mixture(
        rt_vec, distribution, resolved_bounds,
        init_contaminant, max_contaminant, maxit, tol
      )
      if (!fit$converged || is.null(fit$params)) {
        warning2("EM did not converge{.em_maxit_hint(fit, maxit)}. Using robust moments.")
        c(.robust_aggregation(rt_vec, scale_method = robust_scale),
          contaminant_prop = NA_real_)
      } else {
        c(.dist_moments(fit$params, distribution),
          contaminant_prop = fit$contaminant_prop)
      }
    }
  }

  if (version == "3par") {
    moments <- compute_moments(rt, n_trials)
    counts <- .contaminant_free_counts(
      n_upper, n_trials, moments$contaminant_prop, guess_rate
    )
    data.frame(
      mean_rt = moments$mean,
      var_rt = moments$var,
      n_upper = counts$n_upper,
      n_trials = counts$n_trials,
      contaminant_prop = moments$contaminant_prop
    )
  } else {
    rt_upper <- rt[is_upper]
    rt_lower <- rt[!is_upper]
    moments_upper <- compute_moments(rt_upper, length(rt_upper))
    moments_lower <- compute_moments(rt_lower, length(rt_lower))
    counts <- .contaminant_free_counts_per_boundary(
      n_upper, n_trials,
      moments_upper$contaminant_prop, moments_lower$contaminant_prop
    )
    data.frame(
      mean_rt_upper = moments_upper$mean,
      mean_rt_lower = moments_lower$mean,
      var_rt_upper = moments_upper$var,
      var_rt_lower = moments_lower$var,
      n_upper = counts$n_upper,
      n_trials = counts$n_trials,
      contaminant_prop_upper = moments_upper$contaminant_prop,
      contaminant_prop_lower = moments_lower$contaminant_prop
    )
  }
}

# Raising maxit helps only when the EM ran out of iterations, not when it
# stopped early on too few in-range RTs or a numerical failure
.em_maxit_hint <- function(fit, maxit) {
  if (fit$iterations >= maxit) {
    glue(" within maxit = {maxit} iterations (increase `maxit`)")
  } else {
    ""
  }
}

# The EM cleans the RT moments but says nothing about which boundary the
# contaminants went to, so 3par needs the guess rate and 4par does not.
.contaminant_free_counts <- function(n_upper, n_trials, contaminant_prop,
                                     guess_rate) {
  if (is.na(contaminant_prop) || contaminant_prop <= 0) {
    return(nlist(n_upper, n_trials))
  }

  accuracy <- n_upper / n_trials
  upper_limit <- 1 - contaminant_prop * (1 - guess_rate)
  lower_limit <- contaminant_prop * guess_rate
  if (accuracy > upper_limit || accuracy < lower_limit) {
    warning2(
      "Observed accuracy ({round(accuracy, 3)}) lies outside \\
       [{round(lower_limit, 3)}, {round(upper_limit, 3)}], the range a \\
       contaminant proportion of {round(contaminant_prop, 3)} and a guess \\
       rate of {guess_rate} can produce. The correction is degraded for this \\
       cell: the contaminants are removed proportionally from both \\
       boundaries, which preserves the observed accuracy up to rounding. \\
       Check that guess_rate matches the task, or that the contaminant \\
       reaction times are distinguishable from the cognitive ones."
    )
    return(list(
      n_upper = as.integer(round(n_upper * (1 - contaminant_prop))),
      n_trials = as.integer(round(n_trials * (1 - contaminant_prop)))
    ))
  }

  list(
    n_upper = as.integer(round(n_upper - n_trials * contaminant_prop * guess_rate)),
    n_trials = as.integer(round(n_trials * (1 - contaminant_prop)))
  )
}

# Each boundary is rounded on its own because the 4par likelihood reads hits
# and trials - hits, not the total, as the two boundaries' sample sizes.
.contaminant_free_counts_per_boundary <- function(n_upper, n_trials,
                                                  contaminant_prop_upper,
                                                  contaminant_prop_lower) {
  clean_upper <- as.integer(round(
    .clean_boundary_count(n_upper, contaminant_prop_upper)
  ))
  list(
    n_upper = clean_upper,
    n_trials = clean_upper + as.integer(round(
      .clean_boundary_count(n_trials - n_upper, contaminant_prop_lower)
    ))
  )
}

.clean_boundary_count <- function(n, contaminant_prop) {
  if (is.na(contaminant_prop) || contaminant_prop <= 0) n else n * (1 - contaminant_prop)
}

.simple_aggregation <- function(x) {
  list(
    mean = mean(x, na.rm = TRUE),
    var = var(x, na.rm = TRUE)
  )
}

.robust_aggregation <- function(x, scale_method = "iqr") {
  x <- x[!is.na(x)]

  list(
    mean = median(x),
    var = switch(scale_method,
      iqr = (IQR(x) / 1.349)^2,
      mad = mad(x)^2,
      stop2("scale_method must be iqr or mad")
    )
  )
}

.convert_response_to_upper <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(as.logical(x))
  }

  if (is.character(x) || is.factor(x)) {
    x <- tolower(x)
    upper_patterns <- c("upper", "correct", "acc", "1", "true", "yes", "hit")
    lower_patterns <- c(
      "lower", "error", "err", "incorrect", "0", "false",
      "no", "miss", "fa"
    )

    is_upper <- x %in% upper_patterns
    is_lower <- x %in% lower_patterns

    unrecognized <- !is_upper & !is_lower & !is.na(x)
    stopif(
      any(unrecognized),
      "Unrecognized response values: \\
      {collapse_comma(unique(x[unrecognized]))}. Expected values like \\
      'upper', 'lower', 'correct', 'error', 1, 0, TRUE, or FALSE."
    )

    return(is_upper)
  }

  stop2("Response variable must be numeric, logical, character, or factor. \\
        Got class: {class(x)[1]}")
}

.resolve_contaminant_bounds <- function(contaminant_bound, rt_data,
                                        bound_buffer = 0.5) {
  data_min <- min(rt_data, na.rm = TRUE)
  data_max <- max(rt_data, na.rm = TRUE)
  data_range <- data_max - data_min

  resolve_single <- function(bound_val) {
    if (is.numeric(bound_val)) return(bound_val)
    val <- tolower(bound_val)
    if (val == "min") return(data_min)
    if (val == "max") return(data_max)
    # numeric string (e.g., "0.1" from c(0.1, "max"))
    as.numeric(bound_val)
  }

  resolved <- vapply(contaminant_bound, resolve_single, numeric(1))
  lower_is_data <- tolower(contaminant_bound[1]) %in% c("min", "max")
  upper_is_data <- tolower(contaminant_bound[2]) %in% c("min", "max")

  # Buffer data-driven bounds to improve mixture identifiability
  if (lower_is_data && bound_buffer > 0) {
    resolved[1] <- max(0.001, resolved[1] - max(bound_buffer * data_range, 0.1))
  }
  if (upper_is_data && bound_buffer > 0) {
    resolved[2] <- resolved[2] + max(bound_buffer * data_range, 0.1)
  }

  if (resolved[1] >= resolved[2]) {
    warning2("Resolved contaminant bounds are invalid (lower >= upper). \\
             Using data range with buffer.", env.frame = -1)
    buffer_amount <- max(bound_buffer * data_range, 0.1)
    resolved[1] <- max(0.001, data_min - buffer_amount)
    resolved[2] <- data_max + buffer_amount
  }

  warnif(resolved[1] > data_min,
    "Lower contaminant bound ({round(resolved[1], 3)}) is greater than the \\
     minimum observed RT ({round(data_min, 3)}). Observations below the \\
     bound cannot be classified as contaminants.")

  warnif(resolved[2] < data_max,
    "Upper contaminant bound ({round(resolved[2], 3)}) is less than the \\
     maximum observed RT ({round(data_max, 3)}). Observations above the \\
     bound cannot be classified as contaminants.")

  unname(resolved)
}

.adjust_accuracy_counts <- function(n_upper, n_trials, contaminant_prop,
                                    guess_rate) {
  if (is.na(contaminant_prop) || contaminant_prop <= 0) {
    return(list(
      n_upper_adj = as.integer(n_upper),
      n_trials_adj = as.integer(n_trials)
    ))
  }

  n_contam <- rbinom(1, size = n_trials, prob = contaminant_prop)
  # contaminants that happened to be "correct" by chance
  n_contam_upper <- rbinom(1, size = n_contam, prob = guess_rate)

  n_trials_adj <- n_trials - n_contam
  n_upper_adj <- n_upper - n_contam_upper
  n_upper_adj <- max(0L, min(n_upper_adj, n_trials_adj))

  list(
    n_upper_adj = as.integer(n_upper_adj),
    n_trials_adj = as.integer(n_trials_adj)
  )
}

#' Adjust Accuracy Counts for Contamination (deprecated)
#'
#' @description **Deprecated.** [ezdm_summary_stats()] now returns
#'   contaminant-free `n_upper` and `n_trials`, so this second step is no
#'   longer needed; applying it corrects the same cell twice. Set the guess
#'   rate with the `guess_rate` argument of [ezdm_summary_stats()] instead.
#'
#'   Adjusts accuracy counts (`n_upper`, `n_trials`) by removing
#'   estimated contaminant trials using binomial sampling. Contaminant trials
#'   are assumed to produce correct responses at a fixed guess rate (e.g., 0.5
#'   for 2AFC tasks).
#'
#' @param n_upper Numeric. Count of upper boundary (correct) responses, a
#'   single value.
#' @param n_trials Numeric. Total number of trials, a single value.
#' @param contaminant_prop Numeric. Estimated proportion of contaminant trials
#'   (e.g., from the `contaminant_prop` column of [ezdm_summary_stats()]), a
#'   single value.
#' @param guess_rate Numeric. Assumed accuracy rate for contaminant trials
#'   (random guessing). Default is 0.5 (appropriate for 2AFC tasks).
#'
#' @return A 1-row `data.frame` with columns `n_upper_adj` and `n_trials_adj`
#'   (integers). When `contaminant_prop` is `NA` or <= 0, returns the original
#'   counts unchanged.
#'
#' @details The function adjusts one cell at a time. Vectors of length greater
#'   than one are an error; to correct several cells, use
#'   [ezdm_summary_stats()] on each group.
#'
#'   Uses binomial sampling to estimate the number of contaminant trials
#'   and contaminant correct responses, then subtracts these from the raw
#'   counts. Because of the stochastic sampling, results will vary across
#'   calls unless a seed is set by the user.
#'
#' @seealso [ezdm_summary_stats()] for computing the summary statistics and
#'   contamination proportions
#'
#' @keywords transform
#' @export
#'
#' @examples
#' # Deprecated. ezdm_summary_stats() returns contaminant-free counts already:
#' set.seed(42)
#' rt <- c(rnorm(80, 0.55, 0.05), runif(20, 0.1, 4))
#' response <- c(rbinom(80, 1, 0.85), rbinom(20, 1, 0.5))
#' ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4))
#'
adjust_ezdm_accuracy <- function(n_upper, n_trials, contaminant_prop,
                                 guess_rate = 0.5) {
  warning2("The function `adjust_ezdm_accuracy()` is deprecated. \\
            `ezdm_summary_stats()` already returns contaminant-free `n_upper` \\
            and `n_trials`, so applying this correction on top of them counts \\
            the same contaminants twice. Use its `guess_rate` argument to set \\
            the accuracy expected of a contaminant response.")
  stopif(
    length(n_upper) != 1L || length(n_trials) != 1L || length(contaminant_prop) != 1L,
    "`adjust_ezdm_accuracy()` takes one cell at a time. `n_upper`, `n_trials` and \\
     `contaminant_prop` must each be a single value. For a data frame of cells, \\
     use `ezdm_summary_stats()` per group."
  )
  stopif(!is.numeric(n_upper), "n_upper must be numeric")
  stopif(!is.numeric(n_trials), "n_trials must be numeric")
  stopif(!is.numeric(guess_rate) || guess_rate < 0 || guess_rate > 1,
         "guess_rate must be between 0 and 1")
  as.data.frame(.adjust_accuracy_counts(n_upper, n_trials, contaminant_prop,
                                        guess_rate))
}

############################################################################# !
# FLAG_CONTAMINANT_RTS                                                   ####
############################################################################# !

#' Flag contaminant reaction times using mixture modeling
#'
#' @description Identifies contaminant RTs (fast guesses, attention lapses) at
#'   the trial level using mixture modeling. For each trial, it computes the
#'   posterior probability of being a contaminant given a mixture of a uniform
#'   distribution (contaminants) and an RT distribution.
#'
#'   The function takes a numeric vector of RTs and returns a numeric vector of
#'   contamination probabilities, making it compatible with `dplyr::mutate()` and
#'   `dplyr::group_by()` workflows.
#'
#' @param rt Numeric vector. Reaction times in seconds. Must be positive.
#' @param distribution Character. RT distribution for the mixture model:
#'   "exgaussian" (default), "lognormal", or "invgaussian".
#' @param contaminant_bound Vector of length 2. Bounds `[lower, upper]`
#'   for the uniform contaminant distribution. Can be numeric values or
#'   "min"/"max" for data-driven bounds. Default `c("min", "max")`.
#' @param init_contaminant Numeric. Initial contaminant proportion for EM
#'   algorithm. Must be in (0, 1). Default 0.05.
#' @param max_contaminant Numeric. Maximum allowed contaminant proportion. Values
#'   exceeding this are clipped with a warning. Must be in (0, 1]. Default 0.5.
#' @param maxit Integer. Maximum EM iterations. Default 500. A group whose EM
#'   reaches `maxit` without converging gets `NA` probabilities and a warning.
#' @param tol Numeric. Convergence tolerance for log-likelihood. Default 1e-6.
#'
#' @details
#'
#' ## Mixture Model
#'
#' The function fits: `f(RT) = pi_c * Uniform(a,b) + (1-pi_c) * f_RT(RT|theta)`
#'
#' where pi_c is the contaminant proportion, Uniform(a,b) is the contaminant
#' distribution over `contaminant_bound`, and f_RT is the specified RT
#' distribution with parameters theta.
#'
#' ## Grouping
#'
#' To fit separate mixtures by condition or response boundary, use
#' `dplyr::group_by()` before calling this function inside `dplyr::mutate()`.
#'
#' ## Diagnostics
#'
#' Mixture fit diagnostics (parameters, convergence, log-likelihood) are
#' attached as the `"diagnostics"` attribute of the returned vector.
#' Access them with `attr(result, "diagnostics")`.
#'
#' @return Numeric vector of posterior contamination probabilities P(contaminant
#'   | RT), `NA` for every trial of a group whose EM did not converge (a
#'   threshold such as `filter(contam_prob <= 0.5)` drops that whole group, so
#'   check for `NA` first), with a `"diagnostics"` attribute containing a
#'   one-row data.frame with columns: `mixture_params` (list),
#'   `contaminant_prop`, `converged`, `iterations`, `loglik`, `n_trials`,
#'   `distribution`, `method`.
#'
#' @seealso [ezdm_summary_stats()] for aggregated RT statistics with contamination
#'   handling, [validate_fast_guesses()] for testing whether flagged contaminants
#'   show random guessing behavior
#'
#' @keywords transform
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with contaminants
#' library(bmm)
#' set.seed(123)
#' rt_clean <- rgamma(150, shape = 5, rate = 10)
#' rt_contam <- runif(50, 0.1, 0.2)
#'
#' data <- data.frame(
#'   rt = c(rt_clean, rt_contam),
#'   subject = 1,
#'   response = sample(c("upper", "lower"), 200, replace = TRUE)
#' )
#'
#' # Basic usage with mutate
#' library(dplyr)
#' data <- data |>
#'   mutate(contam_prob = flag_contaminant_rts(rt))
#'
#' # Hard threshold: remove trials with P(contaminant) > 0.5
#' data_clean <- data |> filter(contam_prob <= 0.5)
#'
#' # Separate fits by response boundary
#' data <- data |>
#'   group_by(subject, response) |>
#'   mutate(contam_prob = flag_contaminant_rts(rt))
#'
#' # Access diagnostics
#' probs <- flag_contaminant_rts(data$rt)
#' attr(probs, "diagnostics")
#' }
flag_contaminant_rts <- function(
    rt,
    distribution = c("exgaussian", "lognormal", "invgaussian"),
    contaminant_bound = c("min", "max"),
    init_contaminant = 0.05,
    max_contaminant = 0.5,
    maxit = 500,
    tol = 1e-6) {
  distribution <- match.arg(distribution)
  stopif(!is.numeric(rt), "Argument 'rt' must be a numeric vector")
  stopif(length(rt) == 0L, "Argument 'rt' has length 0")

  .validate_contaminant_bounds(contaminant_bound)
  .validate_contaminant_params(init_contaminant, max_contaminant)

  rt_clean <- rt[!is.na(rt)]
  warnif(any(rt_clean > 10), "Some RT values > 10. Ensure RTs are in seconds, not milliseconds.")
  stopif(any(rt_clean <= 0), "Non-positive RT values found.")

  resolved_bounds <- .resolve_contaminant_bounds(contaminant_bound, rt_clean)

  fit <- .fit_rt_mixture(
    rt_clean, distribution, resolved_bounds,
    init_contaminant, max_contaminant, maxit, tol
  )

  if (!fit$converged || is.null(fit$params)) {
    warning2("EM did not converge{.em_maxit_hint(fit, maxit)}. Returning NA.")
    contam_prob_clean <- rep(NA_real_, length(rt_clean))
  } else {
    par <- fit$params
    pi_c <- fit$contaminant_prop
    in_bounds <- rt_clean >= resolved_bounds[1] & rt_clean <= resolved_bounds[2]
    uniform_dens <- ifelse(in_bounds, 1 / (resolved_bounds[2] - resolved_bounds[1]), 0)
    dens_rt <- switch(distribution,
      exgaussian = dexgauss(rt_clean, par["mu"], par["sigma"], par["tau"]),
      lognormal = dlnorm(rt_clean, par["mu"], par["sigma"]),
      invgaussian = dinvgauss(rt_clean, par["mu"], par["lambda"])
    )
    numer_c <- pi_c * uniform_dens
    contam_prob_clean <- numer_c / (numer_c + (1 - pi_c) * pmax(dens_rt, 1e-300))
  }

  result <- rep(NA_real_, length(rt))
  result[!is.na(rt)] <- contam_prob_clean

  attr(result, "diagnostics") <- data.frame(
    mixture_params = I(list(fit$params)),
    contaminant_prop = fit$contaminant_prop,
    converged = fit$converged,
    iterations = fit$iterations,
    loglik = if (fit$converged) fit$loglik else NA_real_,
    n_trials = length(rt_clean),
    distribution = distribution,
    method = "mixture_em",
    stringsAsFactors = FALSE
  )
  result
}

#' Test if fast contaminants show random guessing behavior
#'
#' @description Uses Bayesian Beta-Binomial conjugate analysis to test whether
#'   fast flagged contaminants show random guessing (~50% accuracy for 2AFC).
#'   The test computes the posterior distribution for the proportion of "upper"
#'   responses and uses a Savage-Dickey Bayes Factor to quantify evidence for
#'   or against the guessing hypothesis.
#'
#' @param contam_flag Logical vector indicating which trials were flagged as
#'   contaminants
#' @param rt_data Numeric vector of reaction times (in seconds)
#' @param response Response data in any format accepted by `.convert_response_to_upper()`
#'   (numeric 0/1, logical, character, factor)
#' @param threshold_type Character. How to interpret `rt_threshold`:
#'   - `"quantile"` (default): Use rt_threshold as quantile (0-1)
#'   - `"absolute"`: Use rt_threshold as absolute RT in seconds
#' @param rt_threshold Numeric. Threshold for defining "fast" trials.
#'   Interpretation depends on `threshold_type`:
#'   - If `threshold_type = "quantile"` (default): Quantile of RT distribution
#'     (e.g., 0.25 = 25th percentile). Default 0.25.
#'   - If `threshold_type = "absolute"`: Absolute RT value in seconds
#'     (e.g., 0.25 = 250ms).
#' @param prior_alpha,prior_beta Numeric. Parameters for Beta prior distribution.
#'   Default 1,1 gives uniform prior. Values > 1 express prior belief about
#'   response proportions.
#' @param guess_prob Numeric. Null hypothesis value for guessing probability.
#'   Default 0.5 (equal probability of upper/lower responses).
#' @param credible_mass Numeric. Probability mass for Highest Density Interval.
#'   Default 0.95 for 95% HDI. Common alternatives: 0.90, 0.99.
#'
#' @return List with components:
#'   - `method`: "bayesian"
#'   - `prop_upper`: Observed proportion of upper responses
#'   - `hdi_lower`, `hdi_upper`: 95% Highest Density Interval bounds
#'   - `bf_01`: Bayes Factor for H0 (guessing) vs H1 (non-random)
#'   - `guess_in_hdi`: Logical, whether `guess_prob` is in HDI
#'   - `bf_evidence`: Character, evidence category on Jeffreys scale
#'   - `posterior_alpha`, `posterior_beta`: Posterior Beta parameters
#'   - `n_tested`: Number of fast flagged trials tested
#'   - `rt_threshold`: Actual RT threshold value used (in seconds)
#'   - `threshold_type`: Type of threshold used ("quantile" or "absolute")
#'   - `credible_mass`: Credible mass used for HDI computation
#'   - `mean_rt_tested`: Mean RT of tested trials
#'
#' @details
#' The function performs a Bayesian test using the Beta-Binomial conjugate
#' prior-posterior relationship. With a Beta(alpha, beta) prior and observing n_upper
#' "upper" responses out of n_tested trials, the posterior is:
#'
#' Beta(alpha + n_upper, beta + n_lower)
#'
#' The Savage-Dickey Bayes Factor compares the posterior and prior densities
#' at the null hypothesis value (default 0.5):
#'
#' BF_01 = posterior_density(guess_prob) / prior_density(guess_prob)
#'
#' Evidence categories follow Jeffreys (1961) scale:
#' - BF > 10: Strong evidence for guessing
#' - BF > 3: Moderate evidence for guessing
#' - BF > 1: Anecdotal evidence for guessing
#' - BF < 1/3: Moderate evidence against guessing
#' - BF < 1/10: Strong evidence against guessing
#'
#' **Note**: This function can be used as a standalone validation step after
#' obtaining contamination probabilities from [flag_contaminant_rts()].
#'
#' @references
#' Jeffreys, H. (1961). Theory of Probability (3rd ed.). Oxford University Press.
#'
#' @seealso [flag_contaminant_rts()] for obtaining contamination probabilities
#'
#' @keywords transform
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with random guessing on fast trials
#' set.seed(123)
#' n <- 100
#' rt <- c(runif(20, 0.15, 0.30), rgamma(80, 5, 10))
#' response <- c(rbinom(20, 1, 0.5), rbinom(80, 1, 0.7))
#' contam_flag <- rt < 0.35
#'
#' # Test using quantile threshold (default, adaptive)
#' result1 <- validate_fast_guesses(
#'   contam_flag = contam_flag,
#'   rt_data = rt,
#'   response = response,
#'   threshold_type = "quantile",
#'   rt_threshold = 0.30 # 30th percentile
#' )
#'
#' # Test using absolute threshold (fixed RT)
#' result2 <- validate_fast_guesses(
#'   contam_flag = contam_flag,
#'   rt_data = rt,
#'   response = response,
#'   threshold_type = "absolute",
#'   rt_threshold = 0.30 # 300ms
#' )
#'
#' print(result1$bf_01) # Bayes Factor
#' print(result1$bf_evidence) # Evidence category
#' print(result1$guess_in_hdi) # Is 0.5 in 95% HDI?
#' print(result1$threshold_type) # "quantile"
#' }
validate_fast_guesses <- function(contam_flag, rt_data, response,
                                  threshold_type = c("quantile", "absolute"),
                                  rt_threshold = 0.25,
                                  prior_alpha = 1,
                                  prior_beta = 1,
                                  guess_prob = 0.5,
                                  credible_mass = 0.95) {
  threshold_type <- match.arg(threshold_type)
  stopif(prior_alpha <= 0 || prior_beta <= 0, "Beta prior parameters must be positive")
  stopif(guess_prob <= 0 || guess_prob >= 1, "guess_prob must be between 0 and 1 (exclusive)")

  stopif(
    threshold_type == "quantile" && (rt_threshold <= 0 || rt_threshold >= 1),
    "rt_threshold must be a quantile between 0 and 1 (exclusive) when threshold_type='quantile'"
  )
  stopif(
    threshold_type == "absolute" && rt_threshold <= 0,
    "rt_threshold must be positive when threshold_type='absolute'"
  )

  if (threshold_type == "quantile") {
    rt_threshold <- stats::quantile(rt_data, rt_threshold, na.rm = TRUE)
  } 

  fast_flagged_idx <- contam_flag & (rt_data < rt_threshold) & !is.na(contam_flag) & !is.na(rt_data)
  n_tested <- sum(fast_flagged_idx)

  is_upper <- .convert_response_to_upper(response)
  n_upper <- sum(is_upper[fast_flagged_idx])
  prop_upper <- n_upper / n_tested
  posterior_alpha <- prior_alpha + n_upper
  posterior_beta <- prior_beta + (n_tested - n_upper)

  hdi <- .compute_beta_hdi(posterior_alpha, posterior_beta, credible_mass = credible_mass)
  bf_01 <- stats::dbeta(guess_prob, posterior_alpha, posterior_beta) /
    stats::dbeta(guess_prob, prior_alpha, prior_beta)
  guess_in_hdi <- guess_prob >= hdi$lower && guess_prob <= hdi$upper

  nlist(
    method = "bayesian", prop_upper, hdi_lower = hdi$lower, hdi_upper = hdi$upper, bf_01, guess_in_hdi,
    bf_evidence = .categorize_bf(bf_01), posterior_alpha, posterior_beta, n_tested,
    rt_threshold = as.numeric(rt_threshold), threshold_type,
    credible_mass, mean_rt_tested = mean(rt_data[fast_flagged_idx])
  )
}

# Jeffreys (1961) scale: BF > 10 (strong), > 3 (moderate), > 1 (anecdotal)
.categorize_bf <- function(bf_01) {
  if (bf_01 > 10) "strong_for_guessing"
  else if (bf_01 > 3) "moderate_for_guessing"
  else if (bf_01 > 1) "anecdotal_for_guessing"
  else if (bf_01 > 1 / 3) "anecdotal_against_guessing"
  else if (bf_01 > 1 / 10) "moderate_against_guessing"
  else "strong_against_guessing"
}

.compute_beta_hdi <- function(alpha, beta, credible_mass = 0.95) {
  # Fall back to equal-tailed interval for extreme (>1000) or U-shaped (both <1) cases
  if (alpha > 1000 || beta > 1000 || (alpha < 1 && beta < 1)) {
    tail_prob <- (1 - credible_mass) / 2
    list(
      lower = stats::qbeta(tail_prob, alpha, beta),
      upper = stats::qbeta(1 - tail_prob, alpha, beta)
    )
  } else {
    lower_percentiles <- seq(0, 1 - credible_mass, length.out = 1000)
    upper_percentiles <- lower_percentiles + credible_mass

    lowers <- stats::qbeta(lower_percentiles, alpha, beta)
    uppers <- stats::qbeta(upper_percentiles, alpha, beta)
    widths <- uppers - lowers

    min_idx <- which.min(widths)
    list(lower = lowers[min_idx], upper = uppers[min_idx])
  }
}


############################################################################# !
# RT MIXTURE (EM) HELPERS                                                 ####
############################################################################# !

# Ex-Gaussian density function
# @param x Numeric vector of values
# @param mu Mean of the Gaussian component
# @param sigma Standard deviation of the Gaussian component
# @param tau Rate parameter of the exponential component
# @param log Logical, return log density if TRUE
dexgauss <- function(x, mu, sigma, tau, log = FALSE) {
  # Ensure positive parameters
  if (sigma <= 0 || tau <= 0) {
    return(rep(if (log) -Inf else 0, length(x)))
  }

  # Ex-Gaussian density: convolution of Gaussian and Exponential
  # Using the standard formula with numerical stability
  z <- (x - mu) / sigma - sigma / tau
  log_dens <- -log(tau) + (sigma^2) / (2 * tau^2) - (x - mu) / tau +
    pnorm(z, log.p = TRUE)

  if (log) {
    return(log_dens)
  }
  exp(log_dens)
}

# Inverse Gaussian (Wald) density function
# @param x Numeric vector of values
# @param mu Mean parameter
# @param lambda Shape parameter
# @param log Logical, return log density if TRUE
dinvgauss <- function(x, mu, lambda, log = FALSE) {
  # Ensure positive parameters and values
  if (mu <= 0 || lambda <= 0) {
    return(rep(if (log) -Inf else 0, length(x)))
  }

  valid <- x > 0
  log_dens <- rep(-Inf, length(x))

  if (any(valid)) {
    xv <- x[valid]
    log_dens[valid] <- 0.5 * (log(lambda) - log(2 * pi) - 3 * log(xv)) -
      (lambda * (xv - mu)^2) / (2 * mu^2 * xv)
  }

  if (log) {
    return(log_dens)
  }
  exp(log_dens)
}

# Compute log-likelihood for a distribution
# @param x Numeric vector of RT values
# @param params Named list of distribution parameters
# @param distribution Character specifying the distribution type
# @param weights Optional numeric vector of observation weights
# @return Log-likelihood value
neg_loglik <- function(x, params, distribution, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(x))
  }

  log_dens <- switch(distribution,
    exgaussian = dexgauss(x, params["mu"], params["sigma"], params["tau"], log = TRUE),
    lognormal = dlnorm(x, params["mu"], params["sigma"], log = TRUE),
    invgaussian = dinvgauss(x, params["mu"], params["lambda"], log = TRUE)
  )

  -sum(weights * log_dens, na.rm = TRUE)
}


# Extract mean and variance from fitted distribution parameters
# @param x Named list of distribution parameters
# @param distribution Character specifying the distribution type
# @return List with mean and var components
.dist_moments <- function(x, distribution = c("exgaussian", "lognormal", "invgaussian")) {
  distribution <- match.arg(distribution)
  switch(distribution,
    exgaussian = list(
      mean = x["mu"] + x["tau"],
      var = x["sigma"]^2 + x["tau"]^2
    ),
    lognormal = list(
      mean = exp(x["mu"] + x["sigma"]^2 / 2),
      var = exp(2 * x["mu"] + x["sigma"]^2) * (exp(x["sigma"]^2) - 1)
    ),
    invgaussian = list(
      mean = x["mu"],
      var = x["mu"]^3 / x["lambda"]
    )
  )
}

# Initialize distribution parameters using method of moments
# @param x Numeric vector of RT values
# @param distribution Character specifying the distribution type
# @return Named list of initial parameter estimates
.init_dist_params <- function(x, distribution) {
  m <- mean(x)
  v <- var(x)
  s <- sd(x)

  switch(distribution,
    exgaussian = {
      # Method of moments for ex-Gaussian
      # Skewness = 2 * tau^3 / (sigma^2 + tau^2)^(3/2)
      # Use simple heuristic: tau captures about 1/3 of the variance
      tau <- max(s / 3, 0.01)
      sigma <- max(sqrt(max(v - tau^2, 0.0001)), 0.01)
      mu <- max(m - tau, 0.01)
      c(mu = mu, sigma = sigma, tau = tau)
    },
    lognormal = {
      # Method of moments for lognormal
      sigma2 <- log(1 + v / m^2)
      sigma <- sqrt(max(sigma2, 0.01))
      mu <- log(m) - sigma2 / 2
      c(mu = mu, sigma = sigma)
    },
    invgaussian = {
      # Method of moments for inverse Gaussian
      mu <- max(m, 0.01)
      lambda <- max(mu^3 / v, 0.01)
      c(mu = mu, lambda = lambda)
    }
  )
}

# Fit distribution parameters using weighted MLE
# @param x Numeric vector of RT values
# @param distribution Character specifying the distribution type
# @param weights Numeric vector of observation weights
# @param init_params Initial parameter estimates
# @return Named list of fitted parameters
.fit_dist_params <- function(x, distribution, weights, init_params) {
  bounds <- .get_param_bounds(distribution)
  result <- tryCatch(
    stats::optim(
      par = init_params,
      fn = \(par) neg_loglik(x, par, distribution, weights),
      method = "L-BFGS-B",
      lower = bounds$lower,
      upper = bounds$upper
    ),
    error = function(e) NULL
  )

  if (is.null(result) || result$convergence != 0) {
    return(init_params)
  }

  result$par
}

# Get parameter bounds for optimization
.get_param_bounds <- function(distribution) {
  switch(distribution,
    exgaussian = list(
      lower = c(-Inf, 1e-6, 1e-6),
      upper = c(Inf, Inf, Inf)
    ),
    lognormal = list(
      lower = c(-Inf, 1e-6),
      upper = c(Inf, Inf)
    ),
    invgaussian = list(
      lower = c(1e-6, 1e-6),
      upper = c(Inf, Inf)
    )
  )
}

# Fit RT mixture model using EM algorithm
# @param x Numeric vector of RT values
# @param distribution Character specifying the parametric distribution
# @param contaminant_bound Numeric vector of length 2 for uniform bounds
# @param init_contaminant Initial contaminant proportion
# @param max_contaminant Maximum allowed contaminant proportion (clipping)
# @param maxit Maximum EM iterations
# @param tol Convergence tolerance
# @return List with fitted params, contaminant proportion, convergence info
.fit_rt_mixture <- function(x, distribution, contaminant_bound,
                            init_contaminant, max_contaminant, maxit, tol) {
  n <- length(x)

  # Filter to valid range for fitting
  x_valid <- x[x >= contaminant_bound[1] & x <= contaminant_bound[2]]
  n_valid <- length(x_valid)

  if (n_valid < 5) {
    return(list(
      params = NULL,
      contaminant_prop = NA,
      converged = FALSE,
      iterations = 0,
      message = "Too few observations in valid range"
    ))
  }

  # Initialize parameters
  pi_c <- init_contaminant # contaminant proportion
  pi_rt <- 1 - pi_c # RT distribution proportion
  dist_params <- .init_dist_params(x_valid, distribution)

  # Uniform density (constant for contaminant component)
  uniform_dens <- 1 / (contaminant_bound[2] - contaminant_bound[1])

  prev_loglik <- -Inf
  converged <- FALSE

  for (iter in seq_len(maxit)) {
    # E-step: compute responsibilities
    dens_rt <- switch(distribution,
      exgaussian = dexgauss(x_valid, dist_params["mu"], dist_params["sigma"],
        dist_params["tau"],
        log = FALSE
      ),
      lognormal = dlnorm(x_valid, dist_params["mu"], dist_params["sigma"]),
      invgaussian = dinvgauss(x_valid, dist_params["mu"], dist_params["lambda"],
        log = FALSE
      )
    )

    # Ensure numerical stability
    dens_rt <- pmax(dens_rt, 1e-300)

    # Posterior probabilities
    numer_rt <- pi_rt * dens_rt
    numer_c <- pi_c * uniform_dens
    denom <- numer_rt + numer_c

    # Responsibilities (prob of being from RT distribution)
    gamma_rt <- numer_rt / denom
    gamma_c <- 1 - gamma_rt

    # Handle numerical issues (NA or NaN in responsibilities)
    if (any(is.na(gamma_rt)) || any(is.nan(gamma_rt))) {
      # Fall back to previous iteration's estimates
      break
    }

    # Compute log-likelihood
    loglik <- sum(log(denom))

    # Handle numerical issues in log-likelihood
    if (is.na(loglik) || is.nan(loglik) || is.infinite(loglik)) {
      break
    }

    # Check convergence
    if (abs(loglik - prev_loglik) < tol) {
      converged <- TRUE
      break
    }
    prev_loglik <- loglik

    # M-step: update parameters
    # Update mixing proportions
    pi_rt <- mean(gamma_rt, na.rm = TRUE)
    pi_c <- 1 - pi_rt

    # Handle edge case where pi_c is NA
    if (is.na(pi_c)) {
      pi_c <- init_contaminant
      pi_rt <- 1 - pi_c
    }

    # Clip contaminant proportion to maximum allowed value
    if (pi_c > max_contaminant) {
      pi_c <- max_contaminant
      pi_rt <- 1 - pi_c
    }

    # Update distribution parameters using weighted MLE
    dist_params <- .fit_dist_params(
      x_valid, distribution, gamma_rt,
      dist_params
    )
  }

  # Warn if contaminant proportion hit the maximum bound
  if (pi_c >= max_contaminant) {
    warning2("Contaminant proportion was clipped to max_contaminant \\
             ({max_contaminant}). This may indicate data quality issues.",
      env.frame = -1
    )
  }

  list(
    params = dist_params,
    contaminant_prop = pi_c,
    converged = converged,
    iterations = iter,
    loglik = if (converged) loglik else NA
  )
}
