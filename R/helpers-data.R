############################################################################# !
# CHECK_DATA METHODS                                                     ####
############################################################################# !

#' @title Generic S3 method for checking data based on model type
#' @description Called by [bmm()] to automatically perform checks on the
#'   data depending on the model type. It will call the appropriate check_data
#'   methods based on the list of classes defined in the .model_* functions. For
#'   models with several classes listed, it will call the functions in the order
#'   they are listed. Thus, any operations that are common to a group of models
#'   should be defined in the appropriate check_data.* function, where \*
#'   corresponds to the shared class. For example, for the .model_imm_abc model,
#'   this corresponds to the following order of check_data.* functions:
#'   check_data() -> check_data.circular(), check_data.non_targets() the output of the
#'   final function is returned to bmm().
#' @param model A model list object returned from check_model()
#' @param data The user supplied data.frame containing the data to be checked
#' @param formula The user supplied formula
#' @return A data.frame with the same number of rows as the input data, but with
#'   additional columns added as necessary, any necessary transformations
#'   applied, and attributes added to the data.frame for later use. If you need
#'   to reuse variables created by the check_data.* functions in subsequent
#'   stages (e.g. in configure_model()), you can store and access them using the
#'   attr() function.
#' @export
#'
#' @keywords internal developer
#'
#' @examples
#' data <- oberauer_lin_2017
#' model <- sdmSimple(resp_error = "dev_rad")
#' formula <- bmf(c ~ 1, kappa ~ 1)
#' checked_data <- check_data(model, data, formula)
check_data <- function(model, data, formula) {
  UseMethod("check_data")
}

#' @export
check_data.default <- function(model, data, formula) {
  data
}

#' @export
check_data.bmmodel <- function(model, data, formula) {
  stopif(missing(data), "Data must be specified using the 'data' argument.")
  data <- try(as.data.frame(data), silent = TRUE)
  stopif(is_try_error(data), "Argument 'data' must be coercible to a data.frame.")
  stopif(!isTRUE(nrow(data) > 0L), "Argument 'data' does not contain observations.")

  attr(data, "data_name") <- substitute_name(data, envir = eval(parent.frame()))
  attr(data, "checked") <- TRUE
  NextMethod("check_data")
}

#' @export
check_data.circular <- function(model, data, formula) {
  resp_name <- model$resp_vars[[1]]
  stopif(
    not_in(resp_name, colnames(data)),
    "The response variable '{resp_name}' is not present in the data."
  )
  warnif(
    max(abs(data[[resp_name]]), na.rm = TRUE) > 2 * pi,
    "It appears your response variable is in degrees.
    The model requires the response variable to be in radians.
    The model will continue to run, but the results may be compromised."
  )

  NextMethod("check_data")
}

#' @export
check_data.non_targets <- function(model, data, formula) {
  nt_features <- model$other_vars$nt_features
  warnif(
    max(abs(data[, nt_features]), na.rm = TRUE) > 2 * pi,
    "It appears at least one of your non_target variables are in degrees.
    The model requires these variable to be in radians.
    The model will continue to run, but the results may be compromised."
  )

  ss <- check_var_set_size(model$other_vars$set_size, data)
  max_set_size <- ss$max_set_size
  ss_numeric <- ss$ss_numeric

  stopif(
    !isTRUE(all.equal(length(nt_features), max_set_size - 1)),
    "The number of columns for non-target values in the argument \\
    'nt_features' should equal max(set_size)-1"
  )

  # create index variables for nt_features and correction variable for theta due to set_size
  lure_idx_vars <- paste0("LureIdx", 1:(max_set_size - 1))
  for (i in 1:(max_set_size - 1)) {
    data[[lure_idx_vars[i]]] <- ifelse(ss_numeric >= (i + 1), 1, 0)
  }
  data$ss_numeric <- ss_numeric
  data$inv_ss <- 1 / (ss_numeric - 1)
  data$inv_ss <- ifelse(is.infinite(data$inv_ss), 1, data$inv_ss)
  data[, nt_features][is.na(data[, nt_features])] <- 0

  # save some variables for later use
  attr(data, "max_set_size") <- max_set_size
  attr(data, "lure_idx_vars") <- lure_idx_vars

  NextMethod("check_data")
}

check_var_set_size <- function(set_size, data) {
  stopif(
    length(set_size) > 1,
    "The set_size variable '{set_size}' must be a single numeric value or \\
    a single variable in your data. You provided a vector of length \\
    {length(set_size)}"
  )

  if (is_data_var(set_size, data)) {
    ss_numeric <- try(as_numeric_vector(data[[set_size]]), silent = TRUE)

    stopif(
      is_try_error(ss_numeric),
      "The set_size variable '{set_size}' must be coercible to a numeric \\
      vector. Did you code your set size as a character vector?"
    )

    max_set_size <- max(ss_numeric, na.rm = TRUE)
  } else {
    max_set_size <- try(as_one_integer(set_size), silent = TRUE)

    stopif(
      is_try_error(max_set_size) | is.logical(set_size),
      "The set_size variable '{set_size}' must be either a variable in your \\
       data or a single numeric value"
    )

    ss_numeric <- rep(max_set_size, nrow(data))
  }

  stopif(
    any(ss_numeric < 1 | ss_numeric %% 1 != 0, na.rm = TRUE),
    "Values of the set_size variable '{set_size}' must be positive whole numbers"
  )

  list(max_set_size = max_set_size, ss_numeric = ss_numeric)
}

############################################################################# !
# HELPER FUNCTIONS                                                       ####
############################################################################# !
#' Calculate response error relative to non-target values
#'
#' @description Given a vector of responses, and the values of non-targets, this
#'   function computes the error relative to each of the non-targets.
#' @param data A `data.frame` object where each row is a single observation
#' @param response Character. The name of the column in `data` which contains
#'   the response
#' @param nt_features Character vector. The names of the columns in `data` which
#'   contain the values of the non-targets
#' @keywords transform
#' @return A `data.frame` with n*m rows, where n is the number of rows of `data`
#'   and m is the number of non-target variables. It preserves all other columns
#'   of `data`, except for the non-target locations, and adds a column `y_nt`,
#'   which contains the transformed response error relative to the non-targets
#'
#' @export
#'
#' @examples
#' data <- oberauer_lin_2017
#' data <- calc_error_relative_to_nontargets(data, "dev_rad", paste0("col_nt", 1:7))
#' hist(data$y_nt, breaks = 100)
#'
calc_error_relative_to_nontargets <- function(data, response, nt_features) {
  y <- y_nt <- non_target_name <- non_target_value <- NULL
  stopif(
    !requireNamespace("tidyr", quietly = TRUE),
    'The "tidyr" package is required for this functionality'
  )
  data <- tidyr::gather(data, non_target_name, non_target_value, eval(nt_features))
  data$y_nt <- wrap(data[[response]] - data[["non_target_value"]])
  data
}

#' @title Wrap angles that extend beyond (-pi;pi)
#' @description On the circular space, angles can be only in the range (-pi;pi
#'   or -180;180). When subtracting angles, this can result in values outside of
#'   this range. For example, when calculating the difference between a value of
#'   10 degrees minus 340 degrees, this results in a difference of 330 degrees.
#'   However, the true difference between these two values is -30 degrees. This
#'   function wraps such values, so that they occur in the circle
#' @param x A numeric vector, matrix or data.frame of angles to be wrapped. In
#'   radians (default) or degrees.
#' @param radians Logical. Is x in radians (default=TRUE) or degrees (FALSE)
#' @return An object of the same type as x
#' @keywords transform
#' @export
#' @examples
#' x <- runif(1000, -pi, pi)
#' y <- runif(1000, -pi, pi)
#' diff <- x - y
#' hist(diff)
#' wrapped_diff <- wrap(x - y)
#' hist(wrapped_diff)
#'
wrap <- function(x, radians = TRUE) {
  stopifnot(is.logical(radians))
  if (radians) {
    return(((x + pi) %% (2 * pi)) - pi)
  }
  ((x + 180) %% (2 * 180)) - 180
}

#' @title Convert degrees to radians or radians to degrees.
#' @description The helper functions `deg2rad` and `rad2deg` should add
#' convenience in transforming data from degrees to radians and from radians to
#' degrees.
#'
#' @name circle_transform
#' @param deg A numeric vector of values in degrees.
#' @param rad A numeric vector of values in radians.
#' @return A numeric vector of the same length as `deg` or `rad`.
#' @keywords transform
#' @export
#' @examples
#' degrees <- runif(100, min = 0, max = 360)
#' radians <- deg2rad(degrees)
#' degrees_again <- rad2deg(radians)
deg2rad <- function(deg) {
  deg * pi / 180
}

#' @rdname circle_transform
#' @export
rad2deg <- function(rad) {
  rad * 180 / pi
}

#' @title Stan data for `bmm` models
#' @description Given the `model`, the `data` and the `formula` for the model,
#'   this function will return the combined stan data generated by `bmm` and
#'   `brms`
#' @inheritParams bmm
#' @aliases standata
#' @param object A `bmmformula` object
#' @param ... Further arguments passed to [brms::standata()]. See the
#'   description of [brms::standata()] for more details
#' @return A named list of objects containing the required data to fit a bmm
#'   model with Stan.
#' @seealso [supported_models()], [brms::standata()]
#' @keywords extract_info
#' @examples
#' sdata1 <- standata(bmf(c ~ 1, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#' str(sdata1)
#' @importFrom brms standata
#' @export
standata.bmmformula <- function(object, data, model, ...) {
  # check model, formula and data, and transform data if necessary
  formula <- object
  configure_options(list(...))
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)

  # generate the model specification to pass to brms later
  config_args <- configure_model(model, data, formula)

  # extract stan data
  dots <- list(...)
  fit_args <- combine_args(nlist(config_args, dots))
  fit_args$object <- fit_args$formula
  fit_args$formula <- NULL
  brms::do_call(brms::standata, fit_args)
}

# check if the data is sorted by the predictors
is_data_ordered <- function(data, formula) {
  dpars <- names(formula)
  predictors <- rhs_vars(formula)
  predictors <- predictors[not_in(predictors, dpars)]
  predictors <- predictors[predictors %in% colnames(data)]
  data <- data[predictors]
  if (length(predictors) > 1) {
    gr_idx <- do.call(paste, c(data, list(sep = "_")))
  } else {
    gr_idx <- unlist(data)
  }
  !has_nonconsecutive_duplicates(gr_idx)
}

# checks if all repetitions of a given value are consecutive in a vector
# by iterating over unique values and checking if all their positions are
# consecutive
has_nonconsecutive_duplicates <- function(vec) {
  unique_vals <- unique(vec)
  cond <- TRUE
  for (val in unique_vals) {
    positions <- which(vec == val)
    cond <- cond & all(diff(positions) == 1)
  }
  !cond
}

############################################################################# !
# EZDM SUMMARY STATISTICS                                                 ####
############################################################################# !

# Validate contaminant_bound argument (shared by ezdm_summary_stats and flag_contaminant_rts)
# @param contaminant_bound Vector of length 2 (numeric or "min"/"max")
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
  invisible(NULL)
}

# Validate init_contaminant and max_contaminant arguments (shared helper)
# @param init_contaminant Numeric. Initial proportion of contaminants
# @param max_contaminant Numeric. Maximum allowed contaminant proportion
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
  invisible(NULL)
}

# Apply a function to each group defined by .by variables (shared helper)
# @param data Data.frame
# @param .by Character vector of grouping variable names
# @param fn Function taking a group data.frame, returning any value
# @param ... Additional arguments passed to fn
# @return List of lists, each with $group_values (1-row data.frame) and $result
.apply_by_groups <- function(data, .by, fn, ...) {
  split_factor <- if (length(.by) == 1) data[[.by]] else interaction(data[.by], drop = TRUE)
  split_data <- split(data, split_factor)
  lapply(names(split_data), function(grp_name) {
    grp_data <- split_data[[grp_name]]
    list(
      group_values = unique(grp_data[.by])[1L, , drop = FALSE],
      result = fn(grp_data, ...)
    )
  })
}

#' Compute Robust Summary Statistics for EZ-Diffusion Model
#'
#' @description Computes robust summary statistics for the EZ-Diffusion Model
#'   by fitting mixture models to raw trial-level RT data, separating
#'   contaminant responses from true responses.
#'
#' @param data A `data.frame` containing trial-level data with RT and accuracy
#'   columns
#' @param rt Character. The name of the column containing reaction times (in
#'   seconds)
#' @param response Character. The name of the column containing response
#'   indicators. Accepts multiple formats:
#'   \itemize{
#'     \item Numeric: 1 = upper/correct, 0 = lower/error
#'     \item Logical: TRUE = upper/correct, FALSE = lower/error
#'     \item Character/Factor: "upper"/"lower", "correct"/"error",
#'       "acc"/"err", "hit"/"miss", "yes"/"no" (case-insensitive)
#'   }
#' @param .by A character vector of column names to group by before computing
#'   summary statistics (e.g., `.by = c("subject", "condition")`). If NULL
#'   (default), computes statistics across all data without grouping.
#' @param version Character. Either "3par" (default) for pooled RTs or "4par"
#'   for separate upper/lower boundary RTs
#' @param distribution Character. The parametric distribution for the RT
#'   component. One of "exgaussian" (default), "lognormal", or "invgaussian"
#' @param method Character. One of "mixture" (default) for robust estimation
#'   via mixture modeling, "robust" for non-parametric robust estimation using
#'   median and IQR/MAD-based variance, or "simple" for standard moment
#'   calculation. The "robust" method is faster and requires no distributional
#'   assumptions, but note that the EZ equations were derived for mean and
#'   variance, so using median may introduce some bias for skewed distributions.
#' @param robust_scale Character. Scale estimator for robust method. Either
#'   "iqr" (default) for IQR-based variance estimation (variance = (IQR/1.349)^2)
#'   or "mad" for MAD-based estimation (variance = MAD^2, where MAD is scaled
#'   to be consistent with SD for normal data). Only used when method = "robust".
#' @param contaminant_bound Vector of length 2 specifying the bounds (in
#'   seconds) for the uniform contaminant distribution. Can be numeric values
#'   or the special strings "min" and "max" to use data-driven bounds:
#'   \itemize{
#'     \item Numeric: Fixed bounds, e.g., c(0.1, 3.0) (default)
#'     \item "min": Use the minimum RT in each group, minus a 50\% buffer
#'     \item "max": Use the maximum RT in each group, plus a 50\% buffer
#'   }
#'   The buffer extends data-driven bounds to ensure conservative estimates.
#'   Examples: c(0.1, 3.0), c("min", "max"), c(0.1, "max"), c("min", 3.0)
#' @param min_trials Integer. Minimum number of trials required for fitting.
#'   Groups with fewer trials will return NA. Default is 10
#' @param init_contaminant Numeric. Initial proportion of contaminants for EM
#'   algorithm. Default is 0.05
#' @param max_contaminant Numeric. Maximum allowed contaminant proportion
#'   (0 < max <= 1). Estimates are clipped to this value to prevent inflated
#'   contaminant proportions. Default is 0.5
#' @param maxit Integer. Maximum number of EM iterations. Default is 100
#' @param tol Numeric. Convergence tolerance for EM algorithm. Default is 1e-6
#' @param adjust_accuracy Logical. If TRUE and method = "mixture", adjust
#'   accuracy counts by removing estimated contaminant guesses using binomial
#'   sampling. Default is FALSE
#' @param guess_rate Numeric. Assumed accuracy rate for contaminant trials
#'   (random guessing). Default is 0.5 (appropriate for 2AFC tasks)
#'
#' @return A `data.frame` with summary statistics. For version = "3par":
#'   grouping variables, mean_rt, var_rt, n_upper, n_trials, contaminant_prop.
#'   When adjust_accuracy = TRUE, also includes n_upper_adj and n_trials_adj.
#'   For version = "4par": grouping variables, mean_rt_upper, mean_rt_lower,
#'   var_rt_upper, var_rt_lower, n_upper, n_trials, contaminant_prop_upper,
#'   contaminant_prop_lower.
#'
#' @details RT outliers and contaminant responses (fast guesses, lapses of
#'   attention) can distort the mean and variance estimates used as input to
#'   the EZ-Diffusion equations. This function addresses this by fitting a
#'   mixture model with two components: a uniform distribution for
#'   contaminants and a parametric RT distribution for true responses.
#'   Robust moments are then extracted from the fitted parametric component.
#'
#' @keywords transform
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' test_data <- data.frame(
#'   subject = rep(1:3, each = 100),
#'   condition = rep(c("A", "B"), 150),
#'   rt = rgamma(300, shape = 5, rate = 10) + 0.3,
#'   correct = rbinom(300, 1, 0.8)
#' )
#'
#' # Compute summary statistics grouped by subject
#' result <- ezdm_summary_stats(test_data,
#'   rt = "rt", response = "correct",
#'   .by = "subject"
#' )
#' print(result)
#'
#' # Group by multiple variables using simple method
#' result_multi <- ezdm_summary_stats(test_data,
#'   rt = "rt",
#'   response = "correct",
#'   .by = c("subject", "condition"),
#'   method = "simple"
#' )
#'
ezdm_summary_stats <- function(
    data,
    rt,
    response,
    .by = NULL,
    version = c("3par", "4par"),
    distribution = c("exgaussian", "lognormal", "invgaussian"),
    method = c("mixture", "simple", "robust"),
    robust_scale = c("iqr", "mad"),
    contaminant_bound = c(0.1, 3.0),
    min_trials = 10,
    init_contaminant = 0.05,
    max_contaminant = 0.5,
    maxit = 100,
    tol = 1e-6,
    adjust_accuracy = FALSE,
    guess_rate = 0.5) {
  stop_missing_args()
  version <- match.arg(version)
  distribution <- match.arg(distribution)
  method <- match.arg(method)
  robust_scale <- match.arg(robust_scale)

  data <- try(as.data.frame(data), silent = TRUE)
  stopif(
    is_try_error(data),
    "Argument 'data' must be coercible to a data.frame"
  )
  stopif(
    !isTRUE(nrow(data) > 0L),
    "Argument 'data' does not contain observations"
  )
  stopif(not_in(rt, colnames(data)), "RT variable '{rt}' not found in data")
  stopif(
    not_in(response, colnames(data)),
    "Response variable '{response}' not found in data"
  )
  .validate_contaminant_bounds(contaminant_bound)
  stopif(
    !is.numeric(min_trials) || min_trials < 1,
    "min_trials must be a positive integer"
  )
  .validate_contaminant_params(init_contaminant, max_contaminant)
  stopif(!is.logical(adjust_accuracy), "adjust_accuracy must be TRUE or FALSE")
  stopif(
    !is.numeric(guess_rate) || guess_rate < 0 || guess_rate > 1,
    "guess_rate must be between 0 and 1"
  )

  # Warn if adjust_accuracy is TRUE but method is "simple"
  warnif(
    adjust_accuracy && method == "simple",
    "adjust_accuracy has no effect with method='simple' (no contaminant estimate)"
  )

  # Warnings for potential data issues
  warnif(
    any(data[[rt]] > 10),
    "Some RT values > 10. Ensure RTs are in seconds, not milliseconds."
  )

  # Filter out non-positive RTs with warning
  warnif(
    any(data[[rt]] <= 0),
    "Non-positive RT founds. These values will be excluded."
  )
  data <- data[data[[rt]] > 0 & !is.na(data[[rt]]), ]

  # Validate grouping variables exist
  missing_by <- setdiff(.by, colnames(data))
  stopif(
    length(missing_by) > 0,
    "Grouping variable(s) not found in data: {paste(missing_by, collapse = ', ')}"
  )

  process_group <- function(grp_data) {
    .process_rt_group(
      rt_data = grp_data[[rt]],
      response_data = grp_data[[response]],
      version = version,
      distribution = distribution,
      method = method,
      robust_scale = robust_scale,
      contaminant_bound = contaminant_bound,
      min_trials = min_trials,
      init_contaminant = init_contaminant,
      max_contaminant = max_contaminant,
      maxit = maxit,
      tol = tol,
      adjust_accuracy = adjust_accuracy,
      guess_rate = guess_rate
    )
  }

  if (length(.by) == 0) {
    result_df <- as.data.frame(process_group(data))
  } else {
    groups <- .apply_by_groups(data, .by, process_group)
    results_list <- lapply(groups, function(g) {
      cbind(g$group_values, as.data.frame(g$result))
    })
    result_df <- do.call(rbind, results_list)
    rownames(result_df) <- NULL
  }

  # Add attributes for diagnostics
  attr(result_df, "version") <- version
  attr(result_df, "distribution") <- distribution
  attr(result_df, "method") <- method

  result_df
}

# Simple aggregation (standard moments)
# @param x Numeric vector
# @return List with mean, var, and n
.simple_aggregation <- function(x) {
  list(
    mean = mean(x, na.rm = TRUE),
    var = var(x, na.rm = TRUE)
  )
}

# Robust aggregation using median and IQR/MAD-based variance
# @param x Numeric vector
# @param scale_method Character, either "iqr" or "mad"
# @return List with mean (median), var, and n
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

# Convert response data to logical indicator for upper boundary
# Handles: numeric (0 or 1), logical (TRUE or FALSE), character/factor
# ("upper" or "lower", "correct" or "error", "acc" or "err", etc.)
# @param x Vector of response values
# @return Logical vector where TRUE = upper boundary response
.convert_response_to_upper <- function(x) {
  # Numeric or logical: treat 1/TRUE as upper
  if (is.numeric(x) || is.logical(x)) {
    return(as.logical(x))
  }

  # Character: match common patterns for upper boundary
  if (is.character(x) || is.factor(x)) {
    x <- tolower(x)
    upper_patterns <- c("upper", "correct", "acc", "1", "true", "yes", "hit")
    lower_patterns <- c(
      "lower", "error", "err", "incorrect", "0", "false",
      "no", "miss", "fa"
    )

    is_upper <- x %in% upper_patterns
    is_lower <- x %in% lower_patterns

    # Check if all responses are recognized
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
        Got class: {class(response_data)[1]}")
}

# Process RT data for a single group
# @param rt_data Numeric vector of RTs
# @param response_data Vector of responses (various formats accepted)
# @param version "3par" or "4par"
# @param distribution Distribution type
# @param method "mixture", "robust", or "simple"
# @param contaminant_bound Bounds for uniform distribution
# @param min_trials Minimum trials required
# @param init_contaminant Initial contaminant proportion
# @param max_contaminant Maximum allowed contaminant proportion
# @param maxit Maximum EM iterations
# @param tol Convergence tolerance
# @param adjust_accuracy Whether to adjust accuracy counts
# @param guess_rate Assumed accuracy rate for contaminants
# @return Named list with summary statistics
.process_rt_group <- function(rt_data, response_data, version, distribution,
                              method, robust_scale, contaminant_bound,
                              min_trials, init_contaminant, max_contaminant,
                              maxit, tol, adjust_accuracy, guess_rate) {
  n_trials <- length(rt_data)
  is_upper <- .convert_response_to_upper(response_data)
  n_upper <- sum(is_upper, na.rm = TRUE)
  resolved_bounds <- .resolve_contaminant_bounds(contaminant_bound, rt_data)

  # Helper to compute moments for a single RT vector
  compute_moments <- function(rt_vec, n_obs) {
    if (n_obs < min_trials) {
      list(mean = NA_real_, var = NA_real_, contaminant_prop = NA_real_)
    } else if (method == "simple") {
      c(.simple_aggregation(rt_vec), contaminant_prop = NA_real_)
    } else if (method == "robust") {
      c(.robust_aggregation(rt_vec, scale_method = robust_scale), contaminant_prop = NA_real_)
    } else if (method == "mixture") {
      fit <- .fit_rt_mixture(
        rt_vec, distribution, resolved_bounds,
        init_contaminant, max_contaminant, maxit, tol
      )

      if (!fit$converged || is.null(fit$params)) {
        warning2("EM did not converge. Using simple moments.", env.frame = -2)
        c(.simple_aggregation(rt_vec), contaminant_prop = NA_real_)
      } else {
        c(.dist_moments(fit$params, distribution), contaminant_prop = fit$contaminant_prop)
      }
    }
  }

  if (version == "3par") {
    moments <- compute_moments(rt_data, n_trials)
    result <- list(
      mean_rt = moments$mean,
      var_rt = moments$var,
      n_upper = n_upper,
      n_trials = n_trials,
      contaminant_prop = moments$contaminant_prop
    )

    if (adjust_accuracy && !is.na(moments$contaminant_prop)) {
      adj <- .adjust_accuracy_counts(
        n_upper, n_trials, moments$contaminant_prop, guess_rate
      )
      result <- c(result, adj)
    }

    return(result)
  }

  # version == "4par": separate by response

  rt_upper <- rt_data[is_upper]
  rt_lower <- rt_data[!is_upper]
  n_lower <- length(rt_lower)

  moments_upper <- compute_moments(rt_upper, n_upper)
  moments_lower <- compute_moments(rt_lower, n_lower)

  result <- list(
    mean_rt_upper = moments_upper$mean,
    mean_rt_lower = moments_lower$mean,
    var_rt_upper = moments_upper$var,
    var_rt_lower = moments_lower$var,
    n_upper = n_upper,
    n_trials = n_trials,
    contaminant_prop_upper = moments_upper$contaminant_prop,
    contaminant_prop_lower = moments_lower$contaminant_prop
  )

  if (adjust_accuracy) {
    contam_upper <- moments_upper$contaminant_prop
    contam_lower <- moments_lower$contaminant_prop
    contam_upper_safe <- if (is.na(contam_upper)) 0 else contam_upper
    contam_lower_safe <- if (is.na(contam_lower)) 0 else contam_lower
    avg_contam <- if (n_upper + n_lower > 0) {
      (n_upper * contam_upper_safe + n_lower * contam_lower_safe) /
        (n_upper + n_lower)
    } else {
      0
    }
    adj <- .adjust_accuracy_counts(n_upper, n_trials, avg_contam, guess_rate)
    result <- c(result, adj)
  }

  result
}

# Resolve contaminant bounds from data
# @param contaminant_bound Vector of length 2 (numeric or "min"/"max")
# @param rt_data Numeric vector of RT values
# @param bound_buffer Proportion to extend data-driven bounds (default 0.5)
# @return Numeric vector of length 2 with resolved bounds
.resolve_contaminant_bounds <- function(contaminant_bound, rt_data,
                                        bound_buffer = 0.5) {
  resolved <- numeric(2)
  data_min <- min(rt_data, na.rm = TRUE)
  data_max <- max(rt_data, na.rm = TRUE)
  data_range <- data_max - data_min

  # Track whether each bound is data-driven
  lower_is_data_driven <- FALSE
  upper_is_data_driven <- FALSE

  # Helper to resolve a single bound
  resolve_single <- function(bound_val, is_lower) {
    if (is.numeric(bound_val)) {
      return(list(value = bound_val, data_driven = FALSE))
    }
    # Character value
    if (tolower(bound_val) == "min") {
      return(list(value = data_min, data_driven = TRUE))
    } else if (tolower(bound_val) == "max") {
      return(list(value = data_max, data_driven = TRUE))
    } else {
      # Try to convert as numeric string (e.g., "0.1" from c(0.1, "max"))
      return(list(value = as.numeric(bound_val), data_driven = FALSE))
    }
  }

  lower_result <- resolve_single(contaminant_bound[1], TRUE)
  upper_result <- resolve_single(contaminant_bound[2], FALSE)

  resolved[1] <- lower_result$value
  resolved[2] <- upper_result$value
  lower_is_data_driven <- lower_result$data_driven
  upper_is_data_driven <- upper_result$data_driven

  # Apply buffer to data-driven bounds to reduce uniform density
  # This improves mixture identifiability by making uniform less competitive
  # Use 50% of range or at least 100ms to ensure conservative estimation
  if (lower_is_data_driven && bound_buffer > 0) {
    buffer_amount <- max(bound_buffer * data_range, 0.1) # At least 100ms
    resolved[1] <- max(0.001, resolved[1] - buffer_amount) # Keep positive
  }
  if (upper_is_data_driven && bound_buffer > 0) {
    buffer_amount <- max(bound_buffer * data_range, 0.1) # At least 100ms
    resolved[2] <- resolved[2] + buffer_amount
  }

  # Ensure lower < upper (swap if necessary due to data-driven bounds)
  if (resolved[1] >= resolved[2]) {
    warning2("Resolved contaminant bounds are invalid (lower >= upper). \\
             Using data range with buffer.", env.frame = -1)
    buffer_amount <- max(bound_buffer * data_range, 0.1)
    resolved[1] <- max(0.001, data_min - buffer_amount)
    resolved[2] <- data_max + buffer_amount
  }

  resolved
}

# Compute adjusted accuracy counts based on contaminant proportion
# Uses binomial sampling to produce integer counts
# @param n_upper Raw count of upper boundary responses
# @param n_trials Total number of trials
# @param contaminant_prop Estimated proportion of contaminants
# @param guess_rate Assumed accuracy rate for contaminants (default 0.5)
# @return List with n_upper_adj and n_trials_adj (integers)
.adjust_accuracy_counts <- function(n_upper, n_trials, contaminant_prop,
                                    guess_rate) {
  if (is.na(contaminant_prop) || contaminant_prop <= 0) {
    return(list(
      n_upper_adj = as.integer(n_upper),
      n_trials_adj = as.integer(n_trials)
    ))
  }

  # Sample number of contaminant trials from binomial distribution
  n_contam <- rbinom(1, size = n_trials, prob = contaminant_prop)

  # Sample number of contaminant upper responses from binomial
  # (contaminants that happened to be "correct" by chance)
  n_contam_upper <- rbinom(1, size = n_contam, prob = guess_rate)

  # Adjusted counts (integers)
  n_trials_adj <- n_trials - n_contam
  n_upper_adj <- n_upper - n_contam_upper

  # Bound to valid range
  n_upper_adj <- max(0L, min(n_upper_adj, n_trials_adj))

  list(
    n_upper_adj = as.integer(n_upper_adj),
    n_trials_adj = as.integer(n_trials_adj)
  )
}

############################################################################# !
# FLAG_CONTAMINANT_RTS                                                   ####
############################################################################# !

#' Flag contaminant reaction times using mixture modeling
#'
#' @description Identifies contaminant RTs (fast guesses, attention lapses) at
#'   the trial level using mixture modeling. Unlike [ezdm_summary_stats()] which
#'   computes aggregated statistics, this function returns per-trial
#'   contamination probabilities or flags to enable filtering before model fitting.
#'
#'   The function fits a mixture of a uniform distribution (representing
#'   contaminants) and an RT distribution (exgaussian, lognormal, or invgaussian)
#'   using an EM algorithm. For each trial, it computes the posterior probability
#'   of being a contaminant.
#'
#' @param data A `data.frame` containing trial-level RT data
#' @param rt Character. Name of the RT column in `data` (in seconds)
#' @param response Character or NULL. Name of the response column. Required when
#'   `version = "separate"`. Can be numeric (0 or 1), logical, character, or factor.
#'   See Details for accepted response formats.
#' @param .by Character vector or NULL. Grouping variables for separate mixture
#'   fits (e.g., `c("subject", "condition")`). Default NULL (single group).
#' @param version Character. Either "aggregate" (default) for single contamination
#'   estimate across all RTs, or "separate" for separate estimates by upper/lower
#'   boundary. Requires `response` when "separate". Mirrors the version argument
#'   in [ezdm_summary_stats()].
#' @param distribution Character. RT distribution for the mixture model:
#'   "exgaussian" (default), "lognormal", or "invgaussian".
#' @param contaminant_bound Numeric vector of length 2. Bounds `[lower, upper]`
#'   for the uniform contaminant distribution. Can be numeric values or
#'   "min"/"max" for data-driven bounds. Default `c(0.1, 3.0)`.
#' @param init_contaminant Numeric. Initial contaminant proportion for EM
#'   algorithm. Must be in (0, 1). Default 0.05.
#' @param max_contaminant Numeric. Maximum allowed contaminant proportion. Values
#'   exceeding this are clipped with a warning. Must be in (0, 1]. Default 0.5.
#' @param maxit Integer. Maximum EM iterations. Default 100.
#' @param tol Numeric. Convergence tolerance for log-likelihood. Default 1e-6.
#' @param output Character. Type of contamination metric to compute:
#'   - "probability" (default): Posterior probability P(contaminant | RT)
#'   - "likelihood_ratio": Ratio of contaminant to RT likelihood
#'   - "flag": Binary indicator based on probabilistic binomial sampling
#' @param validate_fast_guessing Logical. If TRUE and response data is provided,
#'   tests whether fast flagged contaminants show random guessing behavior (~50%
#'   accuracy for 2AFC tasks). By default, tests trials with RT below the 25th
#'   percentile (controlled by `threshold_type = "quantile"` and `rt_threshold = 0.25`
#'   passed via `...`). A warning is issued if responses show evidence against
#'   random guessing. Default FALSE. See [validate_fast_guesses()] for all options.
#' @param ... Additional arguments passed to [validate_fast_guesses()] when
#'   `validate_fast_guessing = TRUE`. See [validate_fast_guesses()] documentation
#'   for supported arguments and their default values
#' @param what_return Character. Type of output to return:
#'   - "data" (default): Returns data.frame with contamination columns only
#'   - "diagnostics": Returns data.frame with mixture fit diagnostics only
#'   - "all": Returns list with both `$data` and `$diagnostics`
#'
#' @details
#'
#' ## Response Formats
#'
#' The `response` variable (when provided) can be in several formats:
#' - **Numeric/Logical**: 1 or TRUE = upper boundary, 0 or FALSE = lower boundary
#' - **Character/Factor**: "upper", "lower", "correct", "error", "acc", "err",
#'   "1", "0", "TRUE", "FALSE", "yes", "no", "hit", "miss", "fa"
#'
#' ## Mixture Model
#'
#' The function fits: `f(RT) = pi_c * Uniform(a,b) + (1-pi_c) * f_RT(RT|theta)`
#'
#' where pi_c is the contaminant proportion, Uniform(a,b) is the contaminant
#' distribution over `contaminant_bound`, and f_RT is the specified RT
#' distribution with parameters theta.
#'
#' ## Probabilistic Flagging
#'
#' When `output = "flag"`, contamination flags are generated probabilistically
#' using binomial sampling based on each trial's contamination probability.
#' This avoids the deterministic threshold issue where all long RTs above a
#' fixed threshold would be flagged, even if they're part of the RT distribution.
#' Each trial has a chance of being flagged proportional to its contamination
#' probability, allowing for more nuanced detection. Use the `seed` argument
#' for reproducible flagging.
#'
#' ## Version Behavior
#'
#' - **version = "3par"**: Fits single mixture to all RTs. Returns
#'   `contam_prob`, `contam_lr`, `contam_flag` columns.
#' - **version = "4par"**: Fits separate mixtures for upper/lower boundary
#'   RTs. Returns `contam_prob_upper`/`_lower` columns. Upper trials have NA
#'   for `_lower` metrics and vice versa.
#'
#' ## Small Groups
#'
#' Groups with fewer than 5 valid observations return NA for contamination
#' metrics with a warning.
#'
#' @return Return type depends on `what_return`:
#'
#' **what_return = "data"**: Data.frame with original columns plus:
#' - `contam_prob`: Contamination probability (version="3par")
#' - `contam_lr`: Likelihood ratio (if output != "flag")
#' - `contam_flag`: Binary contamination flag (if output != "probability")
#' - For version="4par": `_upper` and `_lower` suffixed columns
#'
#' **what_return = "diagnostics"**: Data.frame with one row per group:
#' - Grouping variables (if `.by` specified)
#' - `mixture_params`: List-column with fitted distribution parameters
#' - `contaminant_prop`: Estimated contamination proportion
#' - `converged`: EM convergence status
#' - `iterations`, `loglik`: Fit quality metrics
#' - `n_trials`: Number of trials
#' - For version="4par": `_upper` and `_lower` suffixed columns
#'
#' **what_return = "all"**: List with `$data` and `$diagnostics` components
#'
#' @seealso [ezdm_summary_stats()] for aggregated RT statistics with contamination
#'   handling
#'
#' @keywords transform
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with contaminants
#' library(bmm)
#' set.seed(123)
#' n <- 200
#' rt_clean <- rgamma(150, shape = 5, rate = 10)  # Right-skewed RTs
#' rt_contam <- runif(50, 0.1, 0.2)  # Fast guesses
#'
#' data <- data.frame(
#'   rt = c(rt_clean, rt_contam),
#'   subject = 1,
#'   trial = 1:n
#' )
#'
#' # Basic usage - flag contaminants
#' result <- flag_contaminant_rts(data, rt = "rt")
#' data_clean <- result[!result$contam_flag, ]
#'
#' # Get diagnostics to check mixture fit
#' result_all <- flag_contaminant_rts(data, rt = "rt", what_return = "all")
#' print(result_all$diagnostics)
#'
#' # With grouping
#' data_multi <- rbind(
#'   cbind(data, condition = "A"),
#'   cbind(data, condition = "B")
#' )
#' result_grouped <- flag_contaminant_rts(
#'   data_multi,
#'   rt = "rt",
#'   .by = c("subject", "condition")
#' )
#'
#' # Separate by response boundary (for DDM)
#' data$response <- sample(c(0, 1), n, replace = TRUE)
#' result_sep <- flag_contaminant_rts(
#'   data,
#'   rt = "rt",
#'   response = "response",
#'   version = "4par"
#' )
#' }
flag_contaminant_rts <- function(
    data,
    rt,
    response = NULL,
    .by = NULL,
    version = c("3par", "4par"),
    distribution = c("exgaussian", "lognormal", "invgaussian"),
    contaminant_bound = c(0.1, 3.0),
    init_contaminant = 0.05,
    max_contaminant = 0.5,
    maxit = 100,
    tol = 1e-6,
    output = c("probability", "likelihood_ratio", "flag"),
    validate_fast_guessing = FALSE,
    ...,
    what_return = c("data", "diagnostics", "all")) {

  stop_missing_args()
  version <- match.arg(version)
  distribution <- match.arg(distribution)
  output <- match.arg(output)
  what_return <- match.arg(what_return)

  data <- try(as.data.frame(data), silent = TRUE)
  stopif(is_try_error(data), "Argument 'data' must be coercible to a data.frame")
  stopif(!isTRUE(nrow(data) > 0L), "Argument 'data' does not contain observations")

  stopif(not_in(rt, colnames(data)), "RT variable '{rt}' not found in data")

  stopif(
    version == "4par" && is.null(response),
    "Argument 'response' is required when version = '4par'"
  )
  stopif(
    version == "4par" && not_in(response, colnames(data)),
    "Response variable '{response}' not found in data"
  )

  missing_by <- setdiff(.by, colnames(data))
  stopif(
    length(missing_by) > 0,
    "Grouping variable(s) not found in data: {collapse_comma(missing_by)}"
  )

  .validate_contaminant_bounds(contaminant_bound)
  .validate_contaminant_params(init_contaminant, max_contaminant)

  warnif(any(data[[rt]] > 10), "Some RT values > 10. Ensure RTs are in seconds, not milliseconds.")
  stopif(
    any(data[[rt]] <= 0, na.rm = TRUE),
    "Non-positive RT values found in column '{rt}'. Please remove or correct these before proceeding."
  )
  data <- data[!is.na(data[[rt]]), , drop = FALSE]

  flag_group <- function(grp_data) {
    .flag_rt_group(
      data = grp_data,
      rt = rt,
      response = response,
      version = version,
      distribution = distribution,
      contaminant_bound = contaminant_bound,
      init_contaminant = init_contaminant,
      max_contaminant = max_contaminant,
      maxit = maxit,
      tol = tol,
      output = output,
      validate_fast_guessing = validate_fast_guessing,
      ...
    )
  }

  if (length(.by) == 0) {
    result <- flag_group(data)
  } else {
    groups <- .apply_by_groups(data, .by, flag_group)
    results_list <- lapply(groups, function(g) {
      g$result$diagnostics <- cbind(g$group_values, g$result$diagnostics, row.names = NULL)
      g$result
    })
    result <- list(
      data = do.call(rbind, lapply(results_list, `[[`, "data")),
      diagnostics = do.call(rbind, lapply(results_list, `[[`, "diagnostics"))
    )
    rownames(result$data) <- NULL
    rownames(result$diagnostics) <- NULL
  }

  switch(what_return,
    data = result$data,
    diagnostics = result$diagnostics,
    result
  )
}

# Helper function to process a single group
# @param data Data.frame for a single group
# @param rt Character. Name of RT column
# @param response Character or NULL. Name of response column
# @param version Character. "3par" or "4par"
# @param distribution Character. RT distribution type
# @param contaminant_bound Numeric(2). Bounds for uniform distribution
# @param init_contaminant Numeric. Initial contaminant proportion
# @param max_contaminant Numeric. Maximum contaminant proportion
# @param maxit Integer. Maximum EM iterations
# @param tol Numeric. Convergence tolerance
# @param output Character. Output type
# @param validate_fast_guessing Logical. Whether to validate fast guessing behavior
# @param ... Additional arguments passed to validate_fast_guesses()
# @return List with $data (augmented data.frame) and $diagnostics (fit info)
.flag_rt_group <- function(data, rt, response, version, distribution,
                           contaminant_bound, init_contaminant, max_contaminant,
                           maxit, tol, output, validate_fast_guessing,
                           ...) {

  rt_data <- data[[rt]]
  n_trials <- length(rt_data)

  resolved_bounds <- .resolve_contaminant_bounds(contaminant_bound, rt_data)

  if (version == "3par") {
    fit <- .fit_rt_mixture(
      rt_data, distribution, resolved_bounds,
      init_contaminant, max_contaminant, maxit, tol
    )

    if (!fit$converged || is.null(fit$params)) {
      warning2("EM did not converge for group. Returning NA for contamination metrics.",
               env.frame = -1)
      contam_metrics <- list(
        contam_prob = rep(NA_real_, n_trials),
        contam_lr = rep(NA_real_, n_trials),
        contam_flag = rep(NA, n_trials)
      )
    } else {
      contam_metrics <- .compute_trial_contam_metrics(
        rt_data, fit, distribution, resolved_bounds
      )
    }

    data_out <- data
    data_out$contam_prob <- contam_metrics$contam_prob
    if (output %in% c("likelihood_ratio", "flag")) {
      data_out$contam_lr <- contam_metrics$contam_lr
    }
    if (output %in% c("flag", "probability")) {
      data_out$contam_flag <- contam_metrics$contam_flag
    }

    # Perform fast guess validation if requested and response data available
    # Validation tests the actual flagged trials (from probabilistic flagging)
    fast_guess_result <- NULL
    if (validate_fast_guessing && !is.null(response)) {
      fast_guess_result <- validate_fast_guesses(
        contam_flag = contam_metrics$contam_flag,  # Use actual flags
        rt_data = rt_data,
        response = data[[response]],
        ...  # Pass additional arguments (rt_threshold, prior_alpha, prior_beta, etc.)
      )

      # Issue warning if evidence against guessing is moderate or strong
      # Using Jeffreys scale: BF < 1/3 indicates moderate evidence against H0
      if (!is.na(fast_guess_result$bf_01) &&
          fast_guess_result$bf_01 < 1/3) {
        # Use actual credible mass from result
        cred_mass_pct <- round(fast_guess_result$credible_mass * 100)
        warning2(
          "Fast flagged contaminants (RT < {round(fast_guess_result$rt_threshold, 2)}s) show ",
          "evidence against random guessing. ",
          "Observed upper proportion: {round(fast_guess_result$prop_upper, 2)} ",
          "(guessing would be 0.50). ",
          "{cred_mass_pct}% HDI: ",
          "[{round(fast_guess_result$hdi_lower, 2)}, {round(fast_guess_result$hdi_upper, 2)}]. ",
          "BF_01 = {round(fast_guess_result$bf_01, 2)} ({fast_guess_result$bf_evidence}). ",
          "This suggests flagged contaminants may not be fast guesses, or there is response bias.",
          env.frame = -1
        )
      }
    }

    diagnostics <- data.frame(
      mixture_params = I(list(fit$params)),
      contaminant_prop = fit$contaminant_prop,
      converged = fit$converged,
      iterations = fit$iterations,
      loglik = if (fit$converged) fit$loglik else NA_real_,
      n_trials = n_trials,
      distribution = distribution,
      method = "mixture_em",
      fast_guess_validated = !is.null(fast_guess_result),
      fast_guess_prop_upper = if (!is.null(fast_guess_result)) fast_guess_result$prop_upper else NA_real_,
      fast_guess_hdi_lower = if (!is.null(fast_guess_result)) fast_guess_result$hdi_lower else NA_real_,
      fast_guess_hdi_upper = if (!is.null(fast_guess_result)) fast_guess_result$hdi_upper else NA_real_,
      fast_guess_bf_01 = if (!is.null(fast_guess_result)) fast_guess_result$bf_01 else NA_real_,
      fast_guess_in_hdi = if (!is.null(fast_guess_result)) fast_guess_result$guess_in_hdi else NA,
      fast_guess_evidence = if (!is.null(fast_guess_result)) fast_guess_result$bf_evidence else NA_character_,
      fast_guess_n_tested = if (!is.null(fast_guess_result)) fast_guess_result$n_tested else NA_integer_,
      fast_guess_rt_threshold = if (!is.null(fast_guess_result)) fast_guess_result$rt_threshold else NA_real_,
      fast_guess_mean_rt = if (!is.null(fast_guess_result)) fast_guess_result$mean_rt_tested else NA_real_,
      stringsAsFactors = FALSE
    )

  } else {
    response_data <- data[[response]]
    is_upper <- .convert_response_to_upper(response_data)

    rt_upper <- rt_data[is_upper]
    rt_lower <- rt_data[!is_upper]
    n_upper <- length(rt_upper)
    n_lower <- length(rt_lower)

    fit_upper <- .fit_rt_mixture(
      rt_upper, distribution, resolved_bounds,
      init_contaminant, max_contaminant, maxit, tol
    )
    fit_lower <- .fit_rt_mixture(
      rt_lower, distribution, resolved_bounds,
      init_contaminant, max_contaminant, maxit, tol
    )

    if (!fit_upper$converged || is.null(fit_upper$params)) {
      warning2("EM did not converge for upper boundary. Returning NA.",
               env.frame = -1)
      contam_metrics_upper <- list(
        contam_prob = rep(NA_real_, n_upper),
        contam_lr = rep(NA_real_, n_upper),
        contam_flag = rep(NA, n_upper)
      )
    } else {
      contam_metrics_upper <- .compute_trial_contam_metrics(
        rt_upper, fit_upper, distribution, resolved_bounds
      )
    }

    if (!fit_lower$converged || is.null(fit_lower$params)) {
      warning2("EM did not converge for lower boundary. Returning NA.",
               env.frame = -1)
      contam_metrics_lower <- list(
        contam_prob = rep(NA_real_, n_lower),
        contam_lr = rep(NA_real_, n_lower),
        contam_flag = rep(NA, n_lower)
      )
    } else {
      contam_metrics_lower <- .compute_trial_contam_metrics(
        rt_lower, fit_lower, distribution, resolved_bounds
      )
    }

    data_out <- data
    data_out$contam_prob_upper <- NA_real_
    data_out$contam_prob_lower <- NA_real_
    data_out$contam_prob_upper[is_upper] <- contam_metrics_upper$contam_prob
    data_out$contam_prob_lower[!is_upper] <- contam_metrics_lower$contam_prob

    if (output %in% c("likelihood_ratio", "flag")) {
      data_out$contam_lr_upper <- NA_real_
      data_out$contam_lr_lower <- NA_real_
      data_out$contam_lr_upper[is_upper] <- contam_metrics_upper$contam_lr
      data_out$contam_lr_lower[!is_upper] <- contam_metrics_lower$contam_lr
    }

    if (output %in% c("flag", "probability")) {
      data_out$contam_flag_upper <- NA
      data_out$contam_flag_lower <- NA
      data_out$contam_flag_upper[is_upper] <- contam_metrics_upper$contam_flag
      data_out$contam_flag_lower[!is_upper] <- contam_metrics_lower$contam_flag
    }

    # Perform fast guess validation if requested
    # For version='4par', test upper and lower separately
    # Validation uses the actual flagged trials
    fast_guess_upper <- NULL
    fast_guess_lower <- NULL
    if (validate_fast_guessing) {
      # Test upper boundary flagged contaminants
      fast_guess_upper <- validate_fast_guesses(
        contam_flag = contam_metrics_upper$contam_flag,  # Use actual flags
        rt_data = rt_upper,
        response = rep(TRUE, n_upper),  # All are upper responses
        ...  # Pass additional arguments
      )

      # Test lower boundary flagged contaminants
      fast_guess_lower <- validate_fast_guesses(
        contam_flag = contam_metrics_lower$contam_flag,  # Use actual flags
        rt_data = rt_lower,
        response = rep(FALSE, n_lower),  # All are lower responses
        ...  # Pass additional arguments
      )
    }

    diagnostics <- data.frame(
      mixture_params_upper = I(list(fit_upper$params)),
      mixture_params_lower = I(list(fit_lower$params)),
      contaminant_prop_upper = fit_upper$contaminant_prop,
      contaminant_prop_lower = fit_lower$contaminant_prop,
      converged_upper = fit_upper$converged,
      converged_lower = fit_lower$converged,
      iterations_upper = fit_upper$iterations,
      iterations_lower = fit_lower$iterations,
      loglik_upper = if (fit_upper$converged) fit_upper$loglik else NA_real_,
      loglik_lower = if (fit_lower$converged) fit_lower$loglik else NA_real_,
      n_trials_upper = n_upper,
      n_trials_lower = n_lower,
      distribution = distribution,
      method = "mixture_em",
      fast_guess_validated = !is.null(fast_guess_upper),
      fast_guess_n_upper = if (!is.null(fast_guess_upper)) fast_guess_upper$n_tested else NA_integer_,
      fast_guess_n_lower = if (!is.null(fast_guess_lower)) fast_guess_lower$n_tested else NA_integer_,
      fast_guess_rt_threshold = if (!is.null(fast_guess_upper)) fast_guess_upper$rt_threshold else NA_real_,
      fast_guess_mean_rt_upper = if (!is.null(fast_guess_upper)) fast_guess_upper$mean_rt_tested else NA_real_,
      fast_guess_mean_rt_lower = if (!is.null(fast_guess_lower)) fast_guess_lower$mean_rt_tested else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  list(data = data_out, diagnostics = diagnostics)
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
#' **Note**: This function is called internally by [flag_contaminant_rts()] when
#' `validate_fast_guessing = TRUE`. It can also be used standalone for custom
#' validation workflows.
#'
#' @references
#' Jeffreys, H. (1961). Theory of Probability (3rd ed.). Oxford University Press.
#'
#' @seealso [flag_contaminant_rts()] for contamination detection with optional validation
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
#'   rt_threshold = 0.30  # 30th percentile
#' )
#'
#' # Test using absolute threshold (fixed RT)
#' result2 <- validate_fast_guesses(
#'   contam_flag = contam_flag,
#'   rt_data = rt,
#'   response = response,
#'   threshold_type = "absolute",
#'   rt_threshold = 0.30  # 300ms
#' )
#'
#' print(result1$bf_01)          # Bayes Factor
#' print(result1$bf_evidence)    # Evidence category
#' print(result1$guess_in_hdi)   # Is 0.5 in 95% HDI?
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
  
  # Input validation for maintainability
  stopif(prior_alpha <= 0 || prior_beta <= 0,
         "Beta prior parameters must be positive")
  stopif(guess_prob <= 0 || guess_prob >= 1,
         "guess_prob must be between 0 and 1 (exclusive)")
  
  # Validate rt_threshold based on threshold_type
  if (threshold_type == "quantile") {
    stopif(rt_threshold <= 0 || rt_threshold >= 1,
           "rt_threshold must be a quantile between 0 and 1 (exclusive) when threshold_type='quantile'")
  } else {
    stopif(rt_threshold <= 0,
           "rt_threshold must be positive when threshold_type='absolute'")
  }

  is_upper <- .convert_response_to_upper(response)
  
  # Compute actual RT threshold based on threshold_type
  if (threshold_type == "quantile") {
    rt_threshold_value <- stats::quantile(rt_data, rt_threshold, na.rm = TRUE)
  } else {
    rt_threshold_value <- rt_threshold
  }

  # Select fast flagged trials (vectorized)
  fast_flagged_idx <- contam_flag & (rt_data < rt_threshold_value) &
                      !is.na(contam_flag) & !is.na(rt_data)
  n_tested <- sum(fast_flagged_idx)

  # Posterior: Beta(alpha + n_upper, beta + n_lower)
  # Note: HDI provides useful uncertainty information even with low trial numbers
  n_upper <- sum(is_upper[fast_flagged_idx])
  prop_upper <- n_upper / n_tested
  posterior_alpha <- prior_alpha + n_upper
  posterior_beta <- prior_beta + (n_tested - n_upper)

  # HDI and Bayes Factor
  hdi <- .compute_beta_hdi(posterior_alpha, posterior_beta,
                           credible_mass = credible_mass)
  bf_01 <- stats::dbeta(guess_prob, posterior_alpha, posterior_beta) /
           stats::dbeta(guess_prob, prior_alpha, prior_beta)
  guess_in_hdi <- guess_prob >= hdi$lower && guess_prob <= hdi$upper
  bf_evidence <- .categorize_bf(bf_01)

  method <- "bayesian"
  hdi_lower <- hdi$lower
  hdi_upper <- hdi$upper
  mean_rt_tested <- mean(rt_data[fast_flagged_idx])

  nlist(method, prop_upper, hdi_lower, hdi_upper, bf_01, guess_in_hdi,
        bf_evidence, posterior_alpha, posterior_beta, n_tested,
        rt_threshold = as.numeric(rt_threshold_value), threshold_type,
        credible_mass, mean_rt_tested)
}

# Categorize Bayes Factor using Jeffreys scale
# See: Jeffreys, H. (1961). Theory of Probability (3rd ed.). Oxford University Press.
# Thresholds: BF > 10 (strong), > 3 (moderate), > 1 (anecdotal)
# @param bf_01 Numeric. Bayes factor for H0 vs H1
# @return Character. Evidence category on Jeffreys scale
.categorize_bf <- function(bf_01) {
  if (bf_01 > 10) return("strong_for_guessing")
  if (bf_01 > 3) return("moderate_for_guessing")
  if (bf_01 > 1) return("anecdotal_for_guessing")
  if (bf_01 > 1/3) return("anecdotal_against_guessing")
  if (bf_01 > 1/10) return("moderate_against_guessing")
  "strong_against_guessing"
}

# Compute Highest Density Interval (HDI) for Beta distribution
# Uses grid search to find narrowest interval containing credible_mass
# For extreme parameters, falls back to equal-tailed credible interval
# @param alpha,beta Numeric. Beta distribution parameters
# @param credible_mass Numeric. Probability mass (default from constant)
# @return List with lower and upper bounds
.compute_beta_hdi <- function(alpha, beta,
                              credible_mass = 0.95) {

  # Input validation for numerical stability
  stopif(alpha <= 0 || beta <= 0,
         "Beta parameters alpha and beta must be positive")
  stopif(credible_mass <= 0 || credible_mass >= 1,
         "credible_mass must be between 0 and 1 (exclusive)")

  # Use quantile-based interval for extreme or undefined cases
  # Extreme: alpha or beta > 1000 (numerical precision issues)
  # Undefined: both < 1 (U-shaped, HDI not meaningful)
  if (alpha > 1000 ||
      beta > 1000 ||
      (alpha < 1 && beta < 1)) {
    tail_prob <- (1 - credible_mass) / 2
    return(nlist(
      lower = stats::qbeta(tail_prob, alpha, beta),
      upper = stats::qbeta(1 - tail_prob, alpha, beta)
    ))
  }

  # Grid search: find shortest interval containing credible_mass
  # Higher grid resolution = more accurate HDI but slower computation
  # Grid size = 1000 balances accuracy vs speed
  lower_percentiles <- seq(0, 1 - credible_mass,
                          length.out = 1000)
  upper_percentiles <- lower_percentiles + credible_mass

  # Vectorized quantile computation (efficient)
  lowers <- stats::qbeta(lower_percentiles, alpha, beta)
  uppers <- stats::qbeta(upper_percentiles, alpha, beta)
  widths <- uppers - lowers

  # Return narrowest interval (HDI property)
  min_idx <- which.min(widths)
  nlist(lower = lowers[min_idx], upper = uppers[min_idx])
}

# Compute trial-level contamination metrics from mixture fit
# @param rt_vec Numeric vector of RTs
# @param mixture_fit List returned from .fit_rt_mixture
# @param distribution Character. Distribution type
# @param contaminant_bound Numeric(2). Bounds for uniform
# @return List with contam_prob, contam_lr, contam_flag vectors
.compute_trial_contam_metrics <- function(rt_vec, mixture_fit, distribution,
                                          contaminant_bound) {

  params <- mixture_fit$params
  pi_c <- mixture_fit$contaminant_prop
  pi_rt <- 1 - pi_c

  uniform_dens <- 1 / (contaminant_bound[2] - contaminant_bound[1])

  dens_rt <- switch(distribution,
    exgaussian = dexgauss(rt_vec, params["mu"], params["sigma"],
                          params["tau"], log = FALSE),
    lognormal = dlnorm(rt_vec, params["mu"], params["sigma"]),
    invgaussian = dinvgauss(rt_vec, params["mu"], params["lambda"],
                            log = FALSE)
  )

  dens_rt <- pmax(dens_rt, 1e-300)

  numer_c <- pi_c * uniform_dens
  numer_rt <- pi_rt * dens_rt
  denom <- numer_c + numer_rt

  contam_prob <- numer_c / denom
  contam_lr <- numer_c / numer_rt

  # Probabilistic flagging using binomial sampling
  # Each trial is flagged with probability equal to its contamination probability
  contam_flag <- stats::rbinom(n = length(contam_prob), size = 1, prob = contam_prob) == 1

  list(
    contam_prob = contam_prob,
    contam_lr = contam_lr,
    contam_flag = contam_flag
  )
}
