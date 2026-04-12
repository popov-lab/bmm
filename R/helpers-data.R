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

.pull_sdt_column <- function(data, var, arg_name) {
  stopif(is.null(var), "Argument `{arg_name}` must not be NULL")
  stopif(length(var) != 1 || !is.character(var),
    "Argument `{arg_name}` must be a single column name"
  )
  stopif(!var %in% names(data),
    "Column '{var}' supplied to `{arg_name}` is missing in the data"
  )
  data[[var]]
}

.as_binary01 <- function(x, arg_name) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  stopif(any(!x %in% c(0, 1), na.rm = TRUE),
    "`{arg_name}` must be coded as 0 and 1"
  )
  as.integer(x)
}

.map_binary_levels <- function(x, arg_name, levels = c(0, 1)) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  stopif(length(levels) != 2,
    "`levels` for `{arg_name}` must have length 2"
  )
  if (all(x %in% c(0, 1), na.rm = TRUE)) {
    return(as.integer(x))
  }

  out <- ifelse(x == levels[1], 0L, ifelse(x == levels[2], 1L, NA_integer_))
  stopif(any(is.na(out)),
    "`{arg_name}` must be coded as 0/1 or using `response_levels = c({levels[1]}, {levels[2]})`"
  )
  out
}

.prepare_binary_response <- function(data, stimulus, response = NULL,
                                     accuracy = NULL,
                                     response_levels = c(0, 1)) {
  stopif(!is.null(response) && !is.null(accuracy),
    "Provide only one of `response` or `accuracy`"
  )

  if (!is.null(response)) {
    resp <- .pull_sdt_column(data, response, "response")
    return(.map_binary_levels(resp, "response", response_levels))
  }

  stopif(is.null(accuracy),
    "Provide `response` or `accuracy` to prepare SDT data"
  )
  acc <- .as_binary01(.pull_sdt_column(data, accuracy, "accuracy"), "accuracy")
  ifelse(acc == 1L, stimulus, 1L - stimulus)
}

.sdt_group_rows <- function(data, group_cols) {
  if (length(group_cols) == 0L) {
    group_id <- factor(rep("all", nrow(data)), levels = "all")
    group_data <- data.frame()
  } else {
    group_id <- interaction(data[group_cols], drop = TRUE, lex.order = TRUE)
    first_rows <- match(levels(group_id), group_id)
    group_data <- data[first_rows, group_cols, drop = FALSE]
    rownames(group_data) <- NULL
  }

  nlist(group_id, group_data)
}

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
#' @param n_levels Integer. Number of confidence levels per response side.
#' @param response Optional integer vector (0/1) indicating the observed
#'   response side: 0 = "new"/"noise", 1 = "old"/"signal". If omitted, the
#'   historical behaviour is used and categories are derived from `stimulus`.
#' @param accuracy Optional integer vector (0/1) indicating whether the
#'   observed response was correct. If supplied, `response` is derived from
#'   `stimulus` and `accuracy`.
#'
#' @return Integer vector of combined response categories (1 to 2*n_levels).
#' @export
#' @examples
#' # 3 confidence levels per stimulus type -> K=6 combined categories
#' stim <- c(0, 0, 0, 1, 1, 1)
#' conf <- c(3, 2, 1, 1, 2, 3)
#' combine_sdt_response(stim, conf, n_levels = 3)
#' # Returns: 1, 2, 3, 4, 5, 6
combine_sdt_response <- function(stimulus, confidence, n_levels,
                                 response = NULL, accuracy = NULL) {
  stopif(any(!stimulus %in% c(0, 1)),
    "stimulus must be 0 (noise) or 1 (signal)"
  )
  stopif(any(confidence < 1 | confidence > n_levels),
    "confidence must be between 1 and {n_levels}"
  )
  stopif(!is.null(response) && !is.null(accuracy),
    "Provide only one of `response` or `accuracy`"
  )

  if (!is.null(response)) {
    stopif(any(!response %in% c(0, 1)),
      "response must be 0 (new/noise) or 1 (old/signal)"
    )
  }

  if (!is.null(accuracy)) {
    stopif(any(!accuracy %in% c(0, 1)),
      "accuracy must be coded as 0 and 1"
    )
    response <- ifelse(accuracy == 1L, stimulus, 1L - stimulus)
  }

  if (is.null(response)) {
    response <- stimulus
  }

  ifelse(response == 0,
    n_levels - confidence + 1L,
    n_levels + confidence
  )
}

#' @title Prepare Trial-Level Data for SDT Models
#' @description Aggregates trial-level SDT data into the wide count formats
#'   expected by [sdt_binary()] and the rating-family SDT models.
#' @param data A data frame in trial-level format.
#' @param stimulus Column name coding the stimulus type. Must be coded as
#'   0 (noise/new) and 1 (signal/old).
#' @param response Optional column name. For binary outcomes, this should code
#'   the observed old/new response. For rating outcomes with no `confidence`
#'   column, it should contain the combined rating category (1 to `n_ratings`).
#' @param confidence Optional column name with confidence ratings. When
#'   supplied, combined SDT rating categories are created via
#'   [combine_sdt_response()].
#' @param accuracy Optional column name with correctness coded as 0/1. This can
#'   be used instead of `response`; the observed response side is then derived
#'   from `stimulus` and `accuracy`.
#' @param id_cols Optional character vector of grouping columns (for example
#'   subject or condition identifiers). The output is aggregated over
#'   `id_cols` and `stimulus`.
#' @param outcome Character. Either `"rating"` (default) or `"binary"`.
#' @param n_ratings Optional integer. For rating outcomes, the total number of
#'   response categories. If omitted, it is inferred from the data.
#' @param response_levels Optional length-2 vector defining the coding of the
#'   observed response variable when `response` is binary but not coded as 0/1.
#' @return A data frame in the count format required by the corresponding SDT
#'   model family.
#' @keywords transform
#' @export
prepare_sdt_data <- function(data, stimulus, response = NULL, confidence = NULL,
                             accuracy = NULL, id_cols = NULL,
                             outcome = c("rating", "binary"),
                             n_ratings = NULL,
                             response_levels = c(0, 1)) {
  data <- as.data.frame(data)
  outcome <- match.arg(outcome)
  stopif(any(!id_cols %in% names(data)),
    "Grouping columns {collapse_comma(setdiff(id_cols, names(data)))} are missing in the data"
  )

  stim <- .as_binary01(.pull_sdt_column(data, stimulus, "stimulus"), "stimulus")
  data[[stimulus]] <- stim
  group_cols <- unique(c(id_cols, stimulus))
  grouping <- .sdt_group_rows(data, group_cols)
  group_rows <- split(seq_len(nrow(data)), grouping$group_id)

  if (outcome == "binary") {
    response_side <- .prepare_binary_response(
      data, stimulus = stim, response = response,
      accuracy = accuracy, response_levels = response_levels
    )

    out <- grouping$group_data
    out$n_old <- vapply(group_rows, function(idx) {
      sum(response_side[idx], na.rm = TRUE)
    }, numeric(1))
    out$n_trials <- vapply(group_rows, length, integer(1))
    out$n_old <- as.integer(out$n_old)
    rownames(out) <- NULL
    return(out)
  }

  if (!is.null(response) && !is.null(accuracy)) {
    stop2("Provide only one of `response` or `accuracy` for rating outcomes.")
  }

  category <- if (!is.null(confidence)) {
    conf <- .pull_sdt_column(data, confidence, "confidence")
    if (is.factor(conf)) {
      conf <- as.integer(as.character(conf))
    }
    stopif(any(conf != round(conf), na.rm = TRUE),
      "`confidence` must contain integer values"
    )

    if (is.null(n_ratings)) {
      n_ratings <- 2L * max(conf, na.rm = TRUE)
    }
    stopif(n_ratings %% 2 != 0,
      "`n_ratings` must be even when `confidence` is supplied"
    )

    response_side <- if (!is.null(response)) {
      .map_binary_levels(
        .pull_sdt_column(data, response, "response"),
        "response", response_levels
      )
    } else {
      NULL
    }

    combine_sdt_response(
      stimulus = stim,
      confidence = conf,
      n_levels = n_ratings / 2L,
      response = response_side,
      accuracy = if (!is.null(accuracy)) {
        .as_binary01(.pull_sdt_column(data, accuracy, "accuracy"), "accuracy")
      } else {
        NULL
      }
    )
  } else {
    stopif(is.null(response),
      "For rating outcomes, provide either `confidence` or a combined-category `response`"
    )
    category <- .pull_sdt_column(data, response, "response")
    if (is.factor(category)) {
      category <- as.integer(as.character(category))
    }
    stopif(any(category != round(category), na.rm = TRUE),
      "Combined rating categories must be integers"
    )
    if (is.null(n_ratings)) {
      n_ratings <- max(category, na.rm = TRUE)
    }
    as.integer(category)
  }

  stopif(any(category < 1 | category > n_ratings, na.rm = TRUE),
    "Rating categories must be between 1 and {n_ratings}"
  )

  count_list <- lapply(group_rows, function(idx) {
    as.integer(tabulate(category[idx], nbins = n_ratings))
  })
  count_mat <- matrix(
    unlist(count_list, use.names = FALSE),
    ncol = n_ratings,
    byrow = TRUE
  )

  out <- grouping$group_data
  count_df <- as.data.frame(count_mat)
  names(count_df) <- paste0("r", seq_len(n_ratings))
  out <- cbind(out, count_df, nTrials = as.integer(rowSums(count_mat)))
  rownames(out) <- NULL
  out
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
  stopif(
    length(contaminant_bound) != 2,
    "contaminant_bound must be a vector of length 2"
  )

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
  stopif(
    !is.numeric(min_trials) || min_trials < 1,
    "min_trials must be a positive integer"
  )
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
  for (gv in .by) {
    stopif(
      not_in(gv, colnames(data)),
      "Grouping variable '{gv}' not found in data"
    )
  }

  if (length(.by) == 0) {
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
    if (length(.by) == 1) {
      split_factor <- data[[.by]]
    } else {
      split_factor <- interaction(data[.by], drop = TRUE)
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
      grp_values <- unique(grp_data[.by])
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
# Handles: numeric (0/1), logical (TRUE/FALSE), character/factor
# ("upper"/"lower", "correct"/"error", "acc"/"err", etc.)
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
      'upper'/'lower', 'correct'/'error', 1/0, or TRUE/FALSE."
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
