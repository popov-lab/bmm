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

# Formulas may use the columns check_data() adds as predictors, so the formula
# checks and bmm_data_check() must not call them unknown although the user's
# data lacks them. Methods chain with NextMethod(), as for revert_check_data()
built_data_columns <- function(model) {
  UseMethod("built_data_columns")
}

#' @exportS3Method
built_data_columns.default <- function(model) {
  character(0)
}

# inv_ss and ss_numeric are built whatever the set sizes; the LureIdx step
# functions are one per non-target feature
#' @exportS3Method
built_data_columns.non_targets <- function(model) {
  c(
    paste0("LureIdx", seq_along(model$other_vars$nt_features)),
    "inv_ss", "ss_numeric",
    NextMethod("built_data_columns")
  )
}

#' @export
check_data.bmmodel <- function(model, data, formula) {
  stopif(missing(data), "Data must be specified using the 'data' argument.")
  data_name <- attr(data, "data_name") %||% substitute_name(data, envir = eval(parent.frame()))
  # plyr registers as.data.frame.function, so with its namespace loaded (a
  # pp_check() loads it through bayesplot) a function coerces to a function
  # instead of failing; only a data frame counts as coerced
  data <- try(as.data.frame(data), silent = TRUE)
  stopif(is_try_error(data) || !is.data.frame(data),
         "Argument 'data' must be coercible to a data.frame.")
  stopif(!isTRUE(nrow(data) > 0L), "Argument 'data' does not contain observations.")

  attr(data, "data_name") <- data_name
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
  # rows at set size 1 hold NA in every non-target column
  warnif(
    any(abs(data[, nt_features]) > 2 * pi, na.rm = TRUE),
    "It appears at least one of your non_target variables are in degrees.
    The model requires these variable to be in radians.
    The model will continue to run, but the results may be compromised."
  )

  ss <- check_var_set_size(model$other_vars$set_size, data)
  # brms drops NA-response rows from a fit's stored frame, so the frame's maximum
  # set size can be below the fit's (#459), and newdata for a prediction can
  # stop at any set size; carry_fit_bounds() hands over the fit's bound. Only
  # revert_check_data() and check_newdata() set it, so data from bmm() and
  # update(newdata = ) is still checked against its own maximum
  max_set_size <- attr(data, "fit_max_set_size") %||% ss$max_set_size
  attr(data, "fit_max_set_size") <- NULL
  ss_numeric <- ss$ss_numeric

  stopif(
    ss$max_set_size > max_set_size,
    "The data hold a set size of {ss$max_set_size}, above the {max_set_size} \\
    the model was fitted with."
  )
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

# The checks on trial-level response times that models with a decision share.
# Rows with a missing RT or response are dropped; fast_rt_hint is the model's
# advice after the 0.100 s warning.
check_rt_data <- function(data, rt, response, fast_rt_hint) {
  stopif(
    not_in(rt, colnames(data)),
    "The RT variable '{rt}' is not present in the data."
  )
  stopif(
    not_in(response, colnames(data)),
    "The response variable '{response}' is not present in the data."
  )

  for (var in c(rt, response)) {
    is_missing <- is.na(data[[var]])
    if (any(is_missing)) {
      data <- data[!is_missing, , drop = FALSE]
      warning2("Some values in {var} were NA. These were removed from the analysis.")
    }
  }
  stopif(
    nrow(data) == 0L,
    "Every row has a missing value in '{rt}' or '{response}', so no \\
    observations are left."
  )

  stopif(
    !is.numeric(data[[rt]]),
    "The RT variable '{rt}' needs to be of type double or integer."
  )
  stopif(
    any(data[[rt]] <= 0),
    "Some reaction times are zero or negative, please check your data. \\
    The likelihood is zero for a response time that is not strictly \\
    greater than the non-decision time."
  )
  warnif(
    any(data[[rt]] > 10),
    "Your data contains reaction times larger than 10 seconds.\n
    Either you have passed reaction times in milliseconds, then please \\
    recode them to seconds and rerun the model.\n
    Or you have very long RTs in your data in which case you might want \\
    to consider outlier filtering."
  )
  warnif(
    any(data[[rt]] < 0.100),
    "Your data contains reaction times smaller than 0.100 seconds.\n
    {fast_rt_hint}"
  )

  data
}

# The decision as integer 0 (lower boundary) and 1 (upper boundary); a factor
# is read by its labels
code_dec_response <- function(x, response) {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.numeric(x)) {
    stopif(
      any(!x %in% c(0, 1)),
      "The response variable '{response}' contains values other than 0 and 1.\n
      Please pass responses coded as 0 (lower boundary) and 1 (upper boundary)."
    )
    as.integer(x)
  } else if (is.logical(x)) {
    warning2(
      "The response variable is boolean and will be internally transformed ",
      "to an integer variable with values 0 for FALSE and 1 for TRUE."
    )
    as.integer(x)
  } else if (is.character(x)) {
    x <- tolower(x)
    stopif(
      any(!x %in% c("upper", "lower")),
      "The response variable '{response}' contains invalid character values.\n
      Please pass only 'upper' or 'lower' as response values, or use \\
      numeric coding (0 = lower, 1 = upper)."
    )
    warning2(
      "The response variable is a character variable and will be internally ",
      "transformed to an integer variable with 0 for 'lower' and 1 for 'upper'."
    )
    as.integer(x == "upper")
  } else {
    stop2(
      "The response variable '{response}' is of type '{typeof(x)}'.\n
      Please provide responses as integer (0/1), logical, character \\
      ('upper'/'lower'), or factor."
    )
  }
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
  stopif(!isTRUE(radians) && !isFALSE(radians), "Argument 'radians' must be TRUE or FALSE.")
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
#' @seealso [bmm_models()], [brms::standata()]
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
  dots <- list(...)
  local_brms_threads(dots)
  configure_options(dots)
  call_brms_extractor(
    brms::standata,
    configure_fit(object, data, model, until = "model"),
    dots
  )
}

#' @title Stan data of a fitted `bmm` model
#' @description The prediction methods of `brms`, such as [brms::log_lik()],
#'   [brms::posterior_predict()], [brms::posterior_epred()],
#'   [brms::predict.brmsfit()], [brms::fitted.brmsfit()] and [brms::loo()],
#'   turn `newdata` into Stan data through `standata()`. For a `bmmfit`,
#'   `newdata` can be in the layout of the data you gave to [bmm()]: the
#'   columns `bmm` built from those data when it fitted the model are built
#'   again. Rows of the fit's model frame,
#'   `fit$data`, are used as they are.
#' @param object A `bmmfit` object
#' @param newdata An optional data frame, in the layout of the data given to
#'   [bmm()] or a subset of the rows of `object$data`. Its rows keep their
#'   order. When `newdata` is in the layout given to [bmm()], the columns
#'   `bmm` builds itself are built again, and a column of `newdata` with one
#'   of those names is replaced.
#' @param ... Further arguments passed to [brms::standata()]
#' @return A named list of objects containing the data to evaluate the model
#'   with Stan.
#' @seealso [brms::standata()], [bmm()]
#' @keywords extract_info
#' @export
standata.bmmfit <- function(object, newdata = NULL, ...) {
  if (lacks_built_columns(object, newdata)) {
    newdata <- check_newdata(object, newdata)
  }
  NextMethod()
}

# The model frame holds the columns check_data() built, and data in the layout
# given to bmm() lack some of them. Slices of the frame carry them all: those
# of kfold() and reloo(), and the grids that conditional_effects() and
# native_parameters() cut from it. So does any newdata for a model that
# builds no columns, which brms then checks as it always did
lacks_built_columns <- function(object, newdata) {
  !is.null(newdata) && any(not_in(
    intersect(built_data_columns(object$bmm$model), names(object$data)),
    names(newdata)
  ))
}

# brms returns predictions in the rows of newdata, so the data are not sorted,
# and a check_data() method that dropped or reordered rows would misalign them
# without an error
check_newdata <- function(object, newdata) {
  withr::local_options(bmm.sort_data = FALSE)
  newdata <- carry_fit_bounds(object$bmm$model, as.data.frame(newdata))
  checked <- check_data(object$bmm$model, newdata, object$bmm$user_formula)
  stopif(
    !identical(row.names(checked), row.names(newdata)),
    "Preparing 'newdata' for this model changed its rows. Please report this \\
    as a bug of bmm."
  )
  checked
}

# Data checked against a fit rather than on its own: the fit's stored frame in
# update() and newdata for a prediction. Both can lack the largest set size
# the fit was made with, which check_data() would then take from the data
carry_fit_bounds <- function(model, data) {
  UseMethod("carry_fit_bounds")
}

#' @exportS3Method
carry_fit_bounds.default <- function(model, data) {
  data
}

# nt_features passed the check when the fit was made, so the fit's largest set
# size is one more than their number. Derived rather than stored on the fit, so
# fits saved before #459 was fixed are covered too
#' @exportS3Method
carry_fit_bounds.non_targets <- function(model, data) {
  attr(data, "fit_max_set_size") <- length(model$other_vars$nt_features) + 1
  data
}

# check if the data is sorted by the predictors
is_data_ordered <- function(data, formula) {
  predictors <- data_predictor_vars(data, formula)
  if (length(predictors) == 0) {
    return(TRUE)
  }
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
