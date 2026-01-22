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

#' Compute Summary Statistics for EZ-Diffusion Model
#'
#' @description Computes summary statistics for the EZ-Diffusion Model
#'   from raw trial-level RT data.
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
#' @param ... Grouping variables (unquoted column names). Summary statistics
#'   will be computed separately for each combination of these variables.
#' @param version Character. Either "3par" (default) for pooled RTs or "4par"
#'   for separate upper/lower boundary RTs
#' @param min_trials Integer. Minimum number of trials required for fitting.
#'   Groups with fewer trials will return NA. Default is 10
#'
#' @return A `data.frame` with summary statistics. For version = "3par":
#'   grouping variables, mean_rt, var_rt, n_upper, n_trials.
#'   For version = "4par": grouping variables, mean_rt_upper, mean_rt_lower,
#'   var_rt_upper, var_rt_lower, n_upper, n_trials.
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
#'                              subject)
#' print(result)
#'
#' # Group by multiple variables
#' result_multi <- ezdm_summary_stats(test_data, rt = "rt",
#'                                    response = "correct",
#'                                    subject, condition)
#'
ezdm_summary_stats <- function(
  data,
  rt,
  response,
  ...,
  version = "3par",
  min_trials = 10
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
    !is.numeric(min_trials) || min_trials < 1,
    "min_trials must be a positive integer"
  )

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
  group_vars <- as.character(substitute(list(...)))[-1]

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
      min_trials = min_trials
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
        min_trials = min_trials
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
# @param min_trials Minimum trials required
# @return Named list with summary statistics
.process_rt_group <- function(rt_data, response_data, version, min_trials) {
  n_trials <- length(rt_data)

  # Convert response to logical indicator (TRUE = upper boundary)
  is_upper <- .convert_response_to_upper(response_data)
  n_upper <- sum(is_upper, na.rm = TRUE)

  # Check minimum trials
  if (n_trials < min_trials) {
    if (version == "3par") {
      return(list(
        mean_rt = NA_real_,
        var_rt = NA_real_,
        n_upper = n_upper,
        n_trials = n_trials
      ))
    } else {
      return(list(
        mean_rt_upper = NA_real_,
        mean_rt_lower = NA_real_,
        var_rt_upper = NA_real_,
        var_rt_lower = NA_real_,
        n_upper = n_upper,
        n_trials = n_trials
      ))
    }
  }

  if (version == "3par") {
    # Pool all RTs
    agg <- .simple_aggregation(rt_data)
    result <- list(
      mean_rt = agg$mean,
      var_rt = agg$var,
      n_upper = n_upper,
      n_trials = n_trials
    )
    return(result)
  } else {
    # version == "4par": separate by response
    rt_upper <- rt_data[is_upper]
    rt_lower <- rt_data[!is_upper]
    n_lower <- length(rt_lower)

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
      n_trials = n_trials
    ))
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
