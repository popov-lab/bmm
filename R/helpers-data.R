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
#' @param seed Integer. Random seed for reproducibility of accuracy adjustment.
#'   If NULL (default), results will vary across runs
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
#' result <- ezdm_summary_stats(test_data, rt = "rt", response = "correct",
#'                              .by = "subject")
#' print(result)
#'
#' # Group by multiple variables using simple method
#' result_multi <- ezdm_summary_stats(test_data, rt = "rt",
#'                                    response = "correct",
#'                                    .by = c("subject", "condition"),
#'                                    method = "simple")
#'
ezdm_summary_stats <- function(
  data,
  rt,
  response,
  .by = NULL,
  version = "3par",
  distribution = "exgaussian",
  method = "mixture",
  robust_scale = "iqr",
  contaminant_bound = c(0.1, 3.0),
  min_trials = 10,
  init_contaminant = 0.05,
  max_contaminant = 0.5,
  maxit = 100,
  tol = 1e-6,
  adjust_accuracy = FALSE,
  guess_rate = 0.5,
  seed = NULL
) {
  # Input validation
  stopif(missing(data), "Argument 'data' is required")
  stopif(missing(rt), "Argument 'rt' is required")
  stopif(missing(response), "Argument 'response' is required")

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

  stopif(
    not_in(version, c("3par", "4par")),
    "version must be '3par' or '4par'"
  )
  stopif(
    not_in(distribution, c("exgaussian", "lognormal", "invgaussian")),
    "distribution must be 'exgaussian', 'lognormal', or 'invgaussian'"
  )
  stopif(
    not_in(method, c("mixture", "simple", "robust")),
    "method must be 'mixture', 'simple', or 'robust'"
  )
  stopif(
    not_in(robust_scale, c("iqr", "mad")),
    "robust_scale must be 'iqr' or 'mad'"
  )
  stopif(
    length(contaminant_bound) != 2,
    "contaminant_bound must be a vector of length 2"
  )
  # Validate each element is either numeric (or coercible) or "min"/"max"
  valid_bound <- function(x) {
    if (is.numeric(x)) return(TRUE)
    if (is.character(x)) {
      # Check if it's "min" or "max"
      if (tolower(x) %in% c("min", "max")) return(TRUE)
      # Check if it's a numeric string
      num_val <- suppressWarnings(as.numeric(x))
      if (!is.na(num_val)) return(TRUE)
    }
    FALSE
  }

  stopif(
    !valid_bound(contaminant_bound[1]) || !valid_bound(contaminant_bound[2]),
    "contaminant_bound elements must be numeric or 'min'/'max'"
  )
  # If both are numeric (or coercible), check order
  is_numeric_val <- function(x) {
    is.numeric(x) ||
      (is.character(x) && !is.na(suppressWarnings(as.numeric(x))) &&
         !(tolower(x) %in% c("min", "max")))
  }
  if (is_numeric_val(contaminant_bound[1]) &&
        is_numeric_val(contaminant_bound[2])) {
    stopif(
      as.numeric(contaminant_bound[1]) >= as.numeric(contaminant_bound[2]),
      "contaminant_bound[1] must be less than contaminant_bound[2]"
    )
  }
  stopif(
    !is.numeric(min_trials) || min_trials < 1,
    "min_trials must be a positive integer"
  )
  stopif(
    !is.numeric(init_contaminant) || init_contaminant <= 0 ||
      init_contaminant >= 1,
    "init_contaminant must be between 0 and 1 (exclusive)"
  )
  stopif(
    !is.numeric(max_contaminant) || max_contaminant <= 0 ||
      max_contaminant > 1,
    "max_contaminant must be between 0 (exclusive) and 1 (inclusive)"
  )
  stopif(
    init_contaminant >= max_contaminant,
    "init_contaminant must be less than max_contaminant"
  )
  stopif(
    !is.logical(adjust_accuracy),
    "adjust_accuracy must be TRUE or FALSE"
  )
  stopif(
    !is.numeric(guess_rate) || guess_rate < 0 || guess_rate > 1,
    "guess_rate must be between 0 and 1"
  )
  stopif(
    !is.null(seed) && (!is.numeric(seed) || length(seed) != 1),
    "seed must be NULL or a single integer value"
  )

  # Warn if adjust_accuracy is TRUE but method is "simple"
  warnif(
    adjust_accuracy && method == "simple",
    "adjust_accuracy has no effect with method='simple' (no contaminant estimate)"
  )

  # Set seed for reproducibility if provided
  if (adjust_accuracy && !is.null(seed)) {
    set.seed(seed)
  }

  # Warnings for potential data issues
  warnif(
    any(data[[rt]] > 10, na.rm = TRUE),
    "Some RT values > 10. Ensure RTs are in seconds, not milliseconds."
  )

  # Filter out non-positive RTs with warning
  non_positive <- sum(data[[rt]] <= 0, na.rm = TRUE)
  warnif(
    non_positive > 0,
    "{non_positive} non-positive RT values will be excluded."
  )
  data <- data[data[[rt]] > 0 & !is.na(data[[rt]]), ]

  # Process grouping variables
  group_vars <- if (is.null(.by)) character(0) else .by

  # Validate grouping variables exist
  for (gv in group_vars) {
    stopif(
      not_in(gv, colnames(data)),
      "Grouping variable '{gv}' not found in data"
    )
  }

  # Process data by groups
  if (length(group_vars) == 0) {
    # No grouping - process all data
    result <- .process_rt_group(
      rt_data = data[[rt]],
      response_data = data[[response]],
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
    result_df <- as.data.frame(result)
  } else {
    # Split by grouping variables
    if (length(group_vars) == 1) {
      split_factor <- data[[group_vars]]
    } else {
      split_factor <- interaction(data[group_vars], drop = TRUE)
    }

    split_data <- split(data, split_factor)

    results_list <- lapply(names(split_data), function(grp_name) {
      grp_data <- split_data[[grp_name]]
      grp_result <- .process_rt_group(
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

      # Extract group values
      grp_values <- unique(grp_data[group_vars])
      cbind(grp_values[1, , drop = FALSE], as.data.frame(grp_result))
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
    var = var(x, na.rm = TRUE),
    n = sum(!is.na(x))
  )
}

# Robust aggregation using median and IQR/MAD-based variance
# @param x Numeric vector
# @param scale_method Character, either "iqr" or "mad"
# @return List with mean (median), var, and n
.robust_aggregation <- function(x, scale_method = "iqr") {
  x <- x[!is.na(x)]
  med <- median(x)

  if (scale_method == "iqr") {
    # IQR-based variance estimation
    # For normal distribution: sigma = IQR / (2 * qnorm(0.75)) ≈ IQR / 1.349
    iqr_val <- IQR(x)
    sd_est <- iqr_val / 1.349
    var_est <- sd_est^2
  } else if (scale_method == "mad") {
    # MAD-based variance estimation
    # mad() already scales by 1.4826 to be consistent with SD for normal data
    mad_val <- mad(x)
    var_est <- mad_val^2
  } else {
    stop2("scale_method must be 'iqr' or 'mad'")
  }

  list(
    mean = med,
    var = var_est,
    n = length(x)
  )
}

# Convert response data to logical indicator for upper boundary
# Handles: numeric (0/1), logical (TRUE/FALSE), character/factor
# ("upper"/"lower", "correct"/"error", "acc"/"err", etc.)
# @param response_data Vector of response values
# @return Logical vector where TRUE = upper boundary response
.convert_response_to_upper <- function(response_data) {
  # Handle factors by converting to character
  if (is.factor(response_data)) {
    response_data <- as.character(response_data)
  }

  # Numeric or logical: treat 1/TRUE as upper
  if (is.numeric(response_data) || is.logical(response_data)) {
    return(as.logical(response_data))
  }

  # Character: match common patterns for upper boundary
  if (is.character(response_data)) {
    response_lower <- tolower(response_data)
    # Define patterns for upper boundary responses
    upper_patterns <- c("upper", "correct", "acc", "1", "true", "yes", "hit")
    lower_patterns <- c("lower", "error", "err", "incorrect", "0", "false",
                        "no", "miss", "fa")

    is_upper <- response_lower %in% upper_patterns
    is_lower <- response_lower %in% lower_patterns

    # Check if all responses are recognized
    unrecognized <- !is_upper & !is_lower & !is.na(response_data)
    if (any(unrecognized)) {
      unique_unrec <- unique(response_data[unrecognized])
      stop2("Unrecognized response values: \\
            {paste(unique_unrec, collapse = ', ')}. Expected values like \\
            'upper'/'lower', 'correct'/'error', 1/0, or TRUE/FALSE.")
    }

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
# @param method "mixture" or "simple"
# @param contaminant_bound Bounds for uniform distribution
# @param min_trials Minimum trials required
# @param init_contaminant Initial contaminant proportion
# @param maxit Maximum EM iterations
# @param tol Convergence tolerance
# @return Named list with summary statistics
.process_rt_group <- function(rt_data, response_data, version, distribution,
                              method, robust_scale, contaminant_bound,
                              min_trials, init_contaminant, max_contaminant,
                              maxit, tol, adjust_accuracy, guess_rate) {
  n_trials <- length(rt_data)

  # Convert response to logical indicator (TRUE = upper boundary)
  is_upper <- .convert_response_to_upper(response_data)
  n_upper <- sum(is_upper, na.rm = TRUE)

  # Resolve contaminant bounds from data if "min" or "max" specified
  resolved_bounds <- .resolve_contaminant_bounds(contaminant_bound, rt_data)

  # Check minimum trials
  if (n_trials < min_trials) {
    if (version == "3par") {
      return(list(
        mean_rt = NA_real_,
        var_rt = NA_real_,
        n_upper = n_upper,
        n_trials = n_trials,
        contaminant_prop = NA_real_
      ))
    } else {
      return(list(
        mean_rt_upper = NA_real_,
        mean_rt_lower = NA_real_,
        var_rt_upper = NA_real_,
        var_rt_lower = NA_real_,
        n_upper = n_upper,
        n_trials = n_trials,
        contaminant_prop_upper = NA_real_,
        contaminant_prop_lower = NA_real_
      ))
    }
  }

  if (version == "3par") {
    # Pool all RTs
    if (method == "simple") {
      agg <- .simple_aggregation(rt_data)
      result <- list(
        mean_rt = agg$mean,
        var_rt = agg$var,
        n_upper = n_upper,
        n_trials = n_trials,
        contaminant_prop = NA_real_
      )
      return(result)
    }

    # Robust method (median + IQR/MAD-based variance)
    if (method == "robust") {
      agg <- .robust_aggregation(rt_data, scale_method = robust_scale)
      result <- list(
        mean_rt = agg$mean,
        var_rt = agg$var,
        n_upper = n_upper,
        n_trials = n_trials,
        contaminant_prop = NA_real_
      )
      return(result)
    }

    # Mixture method
    fit <- .fit_rt_mixture(
      rt_data, distribution, resolved_bounds,
      init_contaminant, max_contaminant, maxit, tol
    )

    if (!fit$converged || is.null(fit$params)) {
      # Fall back to simple moments with warning
      warning2("EM did not converge. Using simple moments.", env.frame = -1)
      agg <- .simple_aggregation(rt_data)
      result <- list(
        mean_rt = agg$mean,
        var_rt = agg$var,
        n_upper = n_upper,
        n_trials = n_trials,
        contaminant_prop = fit$contaminant_prop
      )
      if (adjust_accuracy) {
        adj <- .adjust_accuracy_counts(
          n_upper, n_trials, fit$contaminant_prop, guess_rate
        )
        result$n_upper_adj <- adj$n_upper_adj
        result$n_trials_adj <- adj$n_trials_adj
      }
      return(result)
    }

    moments <- .dist_moments(fit$params, distribution)
    result <- list(
      mean_rt = moments$mean,
      var_rt = moments$var,
      n_upper = n_upper,
      n_trials = n_trials,
      contaminant_prop = fit$contaminant_prop
    )
    if (adjust_accuracy) {
      adj <- .adjust_accuracy_counts(
        n_upper, n_trials, fit$contaminant_prop, guess_rate
      )
      result$n_upper_adj <- adj$n_upper_adj
      result$n_trials_adj <- adj$n_trials_adj
    }
    return(result)

  } else {
    # version == "4par": separate by response
    rt_upper <- rt_data[is_upper]
    rt_lower <- rt_data[!is_upper]
    n_lower <- length(rt_lower)

    if (method == "simple") {
      agg_upper <- if (n_upper >= min_trials) {
        .simple_aggregation(rt_upper)
      } else {
        list(mean = NA_real_, var = NA_real_)
      }
      agg_lower <- if (n_lower >= min_trials) {
        .simple_aggregation(rt_lower)
      } else {
        list(mean = NA_real_, var = NA_real_)
      }
      return(list(
        mean_rt_upper = agg_upper$mean,
        mean_rt_lower = agg_lower$mean,
        var_rt_upper = agg_upper$var,
        var_rt_lower = agg_lower$var,
        n_upper = n_upper,
        n_trials = n_trials,
        contaminant_prop_upper = NA_real_,
        contaminant_prop_lower = NA_real_
      ))
    }

    # Robust method (median + IQR/MAD-based variance)
    if (method == "robust") {
      agg_upper <- if (n_upper >= min_trials) {
        .robust_aggregation(rt_upper, scale_method = robust_scale)
      } else {
        list(mean = NA_real_, var = NA_real_)
      }
      agg_lower <- if (n_lower >= min_trials) {
        .robust_aggregation(rt_lower, scale_method = robust_scale)
      } else {
        list(mean = NA_real_, var = NA_real_)
      }
      return(list(
        mean_rt_upper = agg_upper$mean,
        mean_rt_lower = agg_lower$mean,
        var_rt_upper = agg_upper$var,
        var_rt_lower = agg_lower$var,
        n_upper = n_upper,
        n_trials = n_trials,
        contaminant_prop_upper = NA_real_,
        contaminant_prop_lower = NA_real_
      ))
    }

    # Mixture method for upper boundary
    if (n_upper >= min_trials) {
      fit_upper <- .fit_rt_mixture(
        rt_upper, distribution, resolved_bounds,
        init_contaminant, max_contaminant, maxit, tol
      )
      if (!fit_upper$converged || is.null(fit_upper$params)) {
        warning2("EM for upper boundary did not converge. Using simple moments.",
                 env.frame = -1)
        agg_upper <- .simple_aggregation(rt_upper)
        moments_upper <- list(mean = agg_upper$mean, var = agg_upper$var)
        contam_upper <- fit_upper$contaminant_prop
      } else {
        moments_upper <- .dist_moments(fit_upper$params, distribution)
        contam_upper <- fit_upper$contaminant_prop
      }
    } else {
      moments_upper <- list(mean = NA_real_, var = NA_real_)
      contam_upper <- NA_real_
    }

    # Mixture method for lower boundary
    if (n_lower >= min_trials) {
      fit_lower <- .fit_rt_mixture(
        rt_lower, distribution, resolved_bounds,
        init_contaminant, max_contaminant, maxit, tol
      )
      if (!fit_lower$converged || is.null(fit_lower$params)) {
        warning2("EM for lower boundary did not converge. Using simple moments.",
                 env.frame = -1)
        agg_lower <- .simple_aggregation(rt_lower)
        moments_lower <- list(mean = agg_lower$mean, var = agg_lower$var)
        contam_lower <- fit_lower$contaminant_prop
      } else {
        moments_lower <- .dist_moments(fit_lower$params, distribution)
        contam_lower <- fit_lower$contaminant_prop
      }
    } else {
      moments_lower <- list(mean = NA_real_, var = NA_real_)
      contam_lower <- NA_real_
    }

    result <- list(
      mean_rt_upper = moments_upper$mean,
      mean_rt_lower = moments_lower$mean,
      var_rt_upper = moments_upper$var,
      var_rt_lower = moments_lower$var,
      n_upper = n_upper,
      n_trials = n_trials,
      contaminant_prop_upper = contam_upper,
      contaminant_prop_lower = contam_lower
    )

    if (adjust_accuracy) {
      # Use weighted average of contaminant proportions
      contam_upper_safe <- ifelse(is.na(contam_upper), 0, contam_upper)
      contam_lower_safe <- ifelse(is.na(contam_lower), 0, contam_lower)
      if (n_upper + n_lower > 0) {
        avg_contam <- (n_upper * contam_upper_safe +
                         n_lower * contam_lower_safe) / (n_upper + n_lower)
      } else {
        avg_contam <- 0
      }
      adj <- .adjust_accuracy_counts(n_upper, n_trials, avg_contam, guess_rate)
      result$n_upper_adj <- adj$n_upper_adj
      result$n_trials_adj <- adj$n_trials_adj
    }

    result
  }
}

# Ex-Gaussian density function
# @param x Numeric vector of values
# @param mu Mean of the Gaussian component
# @param sigma Standard deviation of the Gaussian component
# @param tau Rate parameter of the exponential component
# @param log Logical, return log density if TRUE
.dexgauss <- function(x, mu, sigma, tau, log = FALSE) {
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
.dinvgauss <- function(x, mu, lambda, log = FALSE) {
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

# Extract mean and variance from fitted distribution parameters
# @param params Named list of distribution parameters
# @param distribution Character specifying the distribution type
# @return List with mean and var components
.dist_moments <- function(params, distribution) {
  switch(distribution,
    exgaussian = {
      mu <- params$mu
      sigma <- params$sigma
      tau <- params$tau
      list(
        mean = mu + tau,
        var = sigma^2 + tau^2
      )
    },
    lognormal = {
      mu <- params$mu
      sigma <- params$sigma
      list(
        mean = exp(mu + sigma^2 / 2),
        var = exp(2 * mu + sigma^2) * (exp(sigma^2) - 1)
      )
    },
    invgaussian = {
      mu <- params$mu
      lambda <- params$lambda
      list(
        mean = mu,
        var = mu^3 / lambda
      )
    }
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
      list(mu = mu, sigma = sigma, tau = tau)
    },
    lognormal = {
      # Method of moments for lognormal
      sigma2 <- log(1 + v / m^2)
      sigma <- sqrt(max(sigma2, 0.01))
      mu <- log(m) - sigma2 / 2
      list(mu = mu, sigma = sigma)
    },
    invgaussian = {
      # Method of moments for inverse Gaussian
      mu <- max(m, 0.01)
      lambda <- max(mu^3 / v, 0.01)
      list(mu = mu, lambda = lambda)
    }
  )
}

# Compute log-likelihood for a distribution
# @param x Numeric vector of RT values
# @param params Named list of distribution parameters
# @param distribution Character specifying the distribution type
# @param weights Optional numeric vector of observation weights
# @return Log-likelihood value
.dist_loglik <- function(x, params, distribution, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(x))
  }

  log_dens <- switch(distribution,
    exgaussian = .dexgauss(x, params$mu, params$sigma, params$tau, log = TRUE),
    lognormal = dlnorm(x, params$mu, params$sigma, log = TRUE),
    invgaussian = .dinvgauss(x, params$mu, params$lambda, log = TRUE)
  )

  sum(weights * log_dens, na.rm = TRUE)
}

# Fit distribution parameters using weighted MLE
# @param x Numeric vector of RT values
# @param distribution Character specifying the distribution type
# @param weights Numeric vector of observation weights
# @param init_params Initial parameter estimates
# @return Named list of fitted parameters
.fit_dist_params <- function(x, distribution, weights, init_params) {
  # Objective function (negative log-likelihood)
  neg_loglik <- function(par) {
    params <- .par_to_list(par, distribution)
    -1 * .dist_loglik(x, params, distribution, weights)
  }

  # Convert initial params to vector
  init_par <- .list_to_par(init_params, distribution)

  # Set bounds based on distribution
  bounds <- .get_param_bounds(distribution)

  # Optimize
  result <- tryCatch(
    stats::optim(
      par = init_par,
      fn = neg_loglik,
      method = "L-BFGS-B",
      lower = bounds$lower,
      upper = bounds$upper
    ),
    error = function(e) NULL
  )

  if (is.null(result) || result$convergence != 0) {
    # Return initial params if optimization fails
    return(init_params)
  }

  .par_to_list(result$par, distribution)
}

# Convert parameter vector to named list
.par_to_list <- function(par, distribution) {
  switch(distribution,
    exgaussian = list(mu = par[1], sigma = par[2], tau = par[3]),
    lognormal = list(mu = par[1], sigma = par[2]),
    invgaussian = list(mu = par[1], lambda = par[2])
  )
}

# Convert named list to parameter vector
.list_to_par <- function(params, distribution) {
  switch(distribution,
    exgaussian = c(params$mu, params$sigma, params$tau),
    lognormal = c(params$mu, params$sigma),
    invgaussian = c(params$mu, params$lambda)
  )
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
  pi_c <- init_contaminant  # contaminant proportion
  pi_rt <- 1 - pi_c         # RT distribution proportion
  dist_params <- .init_dist_params(x_valid, distribution)

  # Uniform density (constant for contaminant component)
  uniform_dens <- 1 / (contaminant_bound[2] - contaminant_bound[1])

  prev_loglik <- -Inf
  converged <- FALSE

  for (iter in seq_len(maxit)) {
    # E-step: compute responsibilities
    dens_rt <- switch(distribution,
      exgaussian = .dexgauss(x_valid, dist_params$mu, dist_params$sigma,
                            dist_params$tau, log = FALSE),
      lognormal = dlnorm(x_valid, dist_params$mu, dist_params$sigma),
      invgaussian = .dinvgauss(x_valid, dist_params$mu, dist_params$lambda,
                               log = FALSE)
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
    dist_params <- .fit_dist_params(x_valid, distribution, gamma_rt,
                                    dist_params)
  }

  # Warn if contaminant proportion hit the maximum bound
  if (pi_c >= max_contaminant) {
    warning2("Contaminant proportion was clipped to max_contaminant \\
             ({max_contaminant}). This may indicate data quality issues.",
             env.frame = -1)
  }

  list(
    params = dist_params,
    contaminant_prop = pi_c,
    converged = converged,
    iterations = iter,
    loglik = if (converged) loglik else NA
  )
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
    buffer_amount <- max(bound_buffer * data_range, 0.1)  # At least 100ms
    resolved[1] <- max(0.001, resolved[1] - buffer_amount)  # Keep positive
  }
  if (upper_is_data_driven && bound_buffer > 0) {
    buffer_amount <- max(bound_buffer * data_range, 0.1)  # At least 100ms
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
