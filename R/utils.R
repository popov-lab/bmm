#' Softmax function and its inverse
#'
#' `softmax` returns the value of the softmax function
#' `softmaxinv` returns the value of the inverse-softmax function
#'
#' The softmax function is a bijective function that maps a real vector with length `m` to a probability vector
#' with length `m` with all non-zero probabilities.  The present functions define the softmax function and its inverse, both with a tuning
#' parameter.
#'
#' The current functions define the softmax as:
#'
#' \deqn{\Large P(\eta_i) = \frac{e^{\lambda \eta_i}}{\sum_{j=1}^m e^{\lambda \eta_j}}}
#'
#' @param eta A numeric vector input
#' @param lambda Tuning parameter (a single positive value)
#' @param p A probability vector (i.e., numeric vector of non-negative values that sum to one)
#' @param eta A numeric vector input
#' @param ref_position The reference position that should be used to calculate the inverse softmax function. The default is the last position.
#' @param ref_value The value the reference position will be set to. The default is 0.
#' @return Value of the softmax function or its inverse
#' @keywords transform
#'
#' @export
#' @examples
#' softmax(5:7)
#' softmaxinv(softmax(5:7), ref_position = 1, ref_value = 5)
softmax <- function(eta, lambda = 1) {
  stopifnot(requireNamespace("matrixStats", quietly = TRUE))
  denom <- matrixStats::logSumExp(lambda * eta)
  exp(c(lambda * eta) - denom)
}

#' @rdname softmax
#' @export
softmaxinv <- function(p, lambda = 1, ref_position = length(p), ref_value = 0) {
  len <- length(p)
  if (len <= 1) {
    return(numeric(0))
  }

  stopif(length(ref_position) > 1, "Please provide a single reference value.")
  stopif(ref_position > len, "The reference value must be less or equal than the length of the probability vector.")

  weights <- (log(p) - log(p[len])) / lambda
  weights + (ref_value - weights[ref_position])
}

#' @title Configure local options during model fitting
#' @description Currently it serves to set local options for parallel processing
#'   and update number of chains if the number of chains is greater than the
#'   number of cores.
#'
#' @param opts A list of options
#' @param env The environment in which to set the options - when set to
#'   parent.frame() the changes would apply to the environment of the function
#'   that called it. In our case, this is the environment of the bmm()
#'   function. Changes will not be propagated to the user environment.
#' @keywords internal developer
#' @noRd
#' @return A list of options to pass to brm()
configure_options <- function(opts, env = parent.frame()) {
  if (isTRUE(opts$parallel)) {
    # for backward compatibility - we depricated the parallel argument
    cores <- parallel::detectCores()
    chains <- opts$chains %||% 4
  } else {
    cores <- opts$cores %||% getOption("mc.cores", 1)
  }
  if (not_in_list("silent", opts)) {
    opts$silent <- getOption("bmm.silent", 1)
  }
  if (is.null(opts$backend) && requireNamespace("cmdstanr", quietly = TRUE)) {
    opts$backend <- "cmdstanr"
  }
  withr::local_options(
    list(
      mc.cores = cores,
      bmm.silent = opts$silent,
      bmm.sort_data = opts$sort_data %||% getOption("bmm.sort_data", "check")
    ),
    .local_envir = env
  )

  # return only options that can be passed to brms/rstan/cmdstanr
  exclude_args <- c("parallel", "sort_data", "cores")
  opts[not_in(names(opts), exclude_args)]
}

# check if a value is not in a vector
not_in <- function(value, vector) {
  !(value %in% vector)
}

# check if a key is not in a list
not_in_list <- function(key, list) {
  !(key %in% names(list))
}

#' wrappers to construct a brms nlf and lf formulas from multiple string
#' arguments
#' @param ... string parts of the formula separated by commas
#' @examples
#' kappa_nts <- paste0("kappa_nt", 1:4)
#' glue_nlf(kappa_nts[i], " ~ kappa") ## same as brms::nlf(kappa_nt1 ~ kappa)
#' @noRd
glue_nlf <- function(..., env.frame = -1) {
  brms::nlf(stats::as.formula(glue(..., .envir = sys.frame(env.frame))))
}

# like glue_nlf but for lf formulas
glue_lf <- function(..., env.frame = -1) {
  brms::lf(stats::as.formula(glue(..., .envir = sys.frame(env.frame))))
}

# function to ensure that if the user wants to overwrite an argument (such as
# init), they can. args$prior is NULL by default or is a user-provided prior
# any argument in args$dots is potentially overwrite a default argument in config_args
combine_args <- function(args) {
  config_args <- args$config_args
  dots <- args$dots
  stopif("family" %in% names(dots), "Unsupported argument 'family'. Use the model argument instead.")
  config_args$prior <- args$prior %||% config_args$prior
  config_args[names(dots)] <- dots
  c(config_args, args$opts)
}

############################
#' @description
#'  stop2, warning2, and message2 are wrappers to the builting functions stop, warning and message.
#'  They automatically suppress the call stack and allow you to use glue syntax for the message.
#' @param ... arguments to pass to glue() - these will form the message to the user
#' @param env.frame the frame in which variables passed to glue will be evaluated. Default is the environment
#'  of the function that calls stop2, message2, warning2
#' @examples
#' name <- "John"
#'
#' stop2("Hi {name}, this is an error")
#'
#' warning2("Hi {name}, this is a long warning \\
#'  that is split on multiple rows \\
#'  but no new lines are added in the output message \\
#'  because of the double slashes at the end of each row.")
#'
#' message2("Hi {name}, this is a long message
#'      that is split on multiple rows
#'      and the printed message shows every new line
#'      because there are not back slashes on any line.
#'
#'      You don't have to worry about identation or adding quotes
#'      or commas to any row. Makes typing long messages a breeze!")
#' @noRd
stop2 <- function(..., env.frame = -1) {
  msg <- glue(..., .envir = sys.frame(env.frame))
  stop(msg, call. = FALSE)
}

warning2 <- function(..., env.frame = -1) {
  msg <- glue(..., .envir = sys.frame(env.frame))
  warning(msg, call. = FALSE)
}

message2 <- function(..., env.frame = -1) {
  silent <- getOption("bmm.silent", 1)
  msg <- glue(..., .envir = sys.frame(env.frame))
  if (silent < 2) {
    message(msg)
  }
  invisible()
}

# stopif and warnif work like the builtin stopifnot except that
# the stoping condition is a single boolean expression,
# and you can provide a custom message
stopif <- function(condition, message) {
  if (condition) {
    stop2(message, env.frame = -2)
  }
}

warnif <- function(condition, message) {
  if (condition) {
    warning2(message, env.frame = -2)
  }
}

# function to ensure proper reading of stan files
read_lines2 <- function(con) {
  lines <- readLines(con, n = -1L, warn = FALSE)
  paste(lines, collapse = "\n")
}

# for testing purposes
install_and_load_bmm_version <- function(version, path) {
  stopif(missing(path))
  if ("package:bmm" %in% search()) {
    detach("package:bmm", unload = TRUE)
  }
  path <- file.path(path, paste0("bmm-", version))
  if (!dir.exists(path) || length(list.files(path)) == 0 ||
    length(list.files(file.path(path, "bmm"))) == 0) {
    dir.create(path)
    remotes::install_github(paste0("venpopov/bmm@", version), lib = path)
  }
  library(bmm, lib.loc = path)
}

#' Extract information from a brmsfit object
#' @param fit A brmsfit object, or a list of brmsfit objects
#' @param what String. What to return:
#'  - "time" for the sampling time per chain
#'  - "time_mean" for the mean sampling time
#' @return Depends on `what` and the class of `fit`. For `brmsfit` objects,
#'   information about the single fit is returned. For `brmsfit_list` objects, a
#'   list or data.frame with the information for each fit is returned.
#'  - "time": A data.frame with the sampling time per chain
#'  - "time_mean": A named numeric vector with the mean sampling time
#' @keywords extract_info
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' fit <- bmm(
#'   formula = bmmformula(c ~ 1, kappa ~ 1),
#'   data = data.frame(y = rsdm(1000)),
#'   model = sdm(resp_error = "y")
#' )
#'
#' fit_info(fit, "time")
fit_info <- function(fit, what) {
  UseMethod("fit_info")
}

#' @export
fit_info.brmsfit <- function(fit, what) {
  fit_attr <- attributes(fit$fit)
  metadata <- fit_attr$metadata
  switch(what,
    time = metadata$time$chains,
    time_mean = colMeans(metadata$time$chains),
    inv_metric = metadata$inv_metric,
  )
}

#' @export
fit_info.brmsfit_list <- function(fit, what) {
  .NotYetImplemented()
}

is_data_var <- function(x, data) {
  is.character(x) && length(x) == 1 && x %in% names(data)
}

is_try_warning <- function(x) {
  inherits(x, "warning")
}

is_bmmodel <- function(x) {
  inherits(x, "bmmodel")
}

is_supported_bmmodel <- function(x) {
  valid_models <- supported_models(print_call = FALSE)
  is_bmmodel(x) && inherits(x, valid_models)
}

is_bmmfit <- function(x) {
  inherits(x, "bmmfit")
}

is_namedlist <- function(x) {
  nms <- names(x)
  is.list(x) && !is.null(nms) && all(nzchar(nms))
}

as_numeric_vector <- function(x) {
  out <- tryCatch(as.numeric(as.character(x)), warning = function(w) w)
  if (is_try_warning(out)) {
    stop2("Cannot coerce '{x}' to a numeric vector")
  }
  out
}

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# for some models it is faster to sample if the normalizing constant is calculated
# only once for all trials that have the same value for the predictors. Currently
# this is only used in the sdm model, and to work it requires that the
# data is ordered by the predictor variables. This function checks if the data is
# ordered by the predictors, and if not, it suggests to the user to sort the data
order_data_query <- function(model, data, formula) {
  sort_data <- getOption("bmm.sort_data", "check")
  dpars <- names(formula)
  predictors <- rhs_vars(formula)
  predictors <- predictors[not_in(predictors, dpars)]
  predictors <- predictors[predictors %in% colnames(data)]

  if (length(predictors) == 0) {
    return(data)
  }

  if (sort_data == "check" && !is_data_ordered(data, formula)) {
    message2(
      "\n\nData is not ordered by predictors.
      You can speed up the model estimation up to several times (!) by ordering \\
      the data by all your predictor columns.\n\n"
    )
    caution_msg <- paste(strwrap("* caution: if you chose Option 2, you need to be careful
      when using brms postprocessing methods that rely on the data order, such as
      generating predictions. Assuming you assigned the result of bmm() to a
      variable called `fit`, you can extract the sorted data from the fitted object
      with:\n\n   data_sorted <- fit$data", width = 80), collapse = "\n")
    caution_msg <- crayon::red(caution_msg)
    disable_msg <- glue("To disable this check and query, use options('bmm.sort_data' \\
      = TRUE) to always sort or options('bmm.sort_data' = FALSE) to never check nor \\
      prompt this question.")

    if (interactive()) {
      var <- utils::menu(
        c(
          "Yes (note: you will receive code to sort your data)",
          paste0(
            "Let bmm sort the data for you and continue with the faster model fitting ",
            crayon::red("(*)")
          ),
          paste0("No, I want to continue with the slower estimation\n\n",
            caution_msg,
            collapse = "\n"
          )
        ),
        title = "Do you want to stop and sort your data? (y/n): "
      )
      if (var == 1) {
        message("Please sort your data by all predictors and then re-run the model.")
        data_name <- attr(data, "data_name")
        if (is.null(data_name)) {
          data_name <- deparse(substitute(data))
        }
        message("To sort your data, use the following code:\n\n")
        message(crayon::green("library(dplyr)"))
        message(crayon::green(data_name, "_sorted <- ", data_name, " %>% arrange(",
          paste(predictors, collapse = ", "),
          ")\n\n",
          sep = ""
        ))
        message("Then re-run the model with the newly sorted data.\n\n", disable_msg)
        stop_quietly()
      } else if (var == 2) {
        message(
          "Your data has been sorted by the following predictors: ",
          paste(predictors, collapse = ", "), "\n"
        )
        data <- data[do.call(order, data[predictors]), ]
      }
    }
    message("\n\n", disable_msg)
  } else if (isTRUE(sort_data)) {
    data <- data[do.call(order, data[predictors]), ]
    message(
      "\nYour data has been sorted by the following predictors: ",
      paste(predictors, collapse = ", "), "\n"
    )
    caution_msg <- paste(strwrap("* caution: you have set `sort_data=TRUE`. You need to be careful
        when using brms postprocessing methods that rely on the data order, such as
        generating predictions. Assuming you assigned the result of bmm() to a
        variable called `fit`, you can extract the sorted data from the fitted object
        with:\n\n   data_sorted <- fit$data", width = 80), collapse = "\n")
    caution_msg <- crayon::red(caution_msg)
    message(caution_msg)
  }
  data
}

# when called from another function, it will return a vector of arguments that are
# missing from the call
missing_args <- function(which = -1) {
  parent_objects <- as.list(sys.frame(which))
  call_head <- sys.call(which)[[1]]
  parent_fun <- if (is.function(call_head)) {
    call_head
  } else if (is.symbol(call_head) || is.character(call_head)) {
    get(as.character(call_head), mode = "function", envir = sys.frame(which))
  } else {
    sys.function(which)
  }
  parent_args <- names(formals(parent_fun))
  parent_args <- parent_args[!parent_args %in% c("...", "")]
  symbols <- names(parent_objects)[sapply(parent_objects, is.symbol)]
  missing <- symbols[symbols %in% parent_args]
  missing
}

# when called from another function, it will stop the execution if any of the
# required arguments are missing
stop_missing_args <- function() {
  missing <- missing_args(-2)
  fun <- deparse(sys.call(-1)[[1]])[1]
  stopif(
    length(missing) > 0,
    "The following required arguments are missing in {fun}(): \\
          {paste(missing, collapse = ', ')}"
  )
}

# custom method form printing nicely formatted character values via cat instead of print
#' @export
print.message <- function(x, ...) {
  cat(x, ...)
}

# returns either x, or all variables that match the regular expression x
# @param x character vector or regular expression
# @param all_variables character vector of all variables within which to search
# @param regex logical. If TRUE, x is treated as a regular expression
get_variables <- function(x, all_variables, regex = FALSE) {
  if (!regex) {
    return(x)
  }

  variables <- all_variables[grep(x, all_variables)]
  stopif(
    length(variables) == 0,
    "No variables found that match the regular expression '{x}'"
  )
  variables
}

# replace the base identical function with a S3 method
identical <- function(x, y, ...) {
  UseMethod("identical")
}

#' @export
identical.default <- function(x, y, ...) {
  base::identical(x, y, ...)
}

#' @export
identical.brmsformula <- function(x, y, ...) {
  res <- waldo::compare(x, y, ignore_formula_env = TRUE)
  length(res) == 0
}

#' @export
identical.bmmformula <- function(x, y, ...) {
  res <- waldo::compare(x, y, ignore_formula_env = TRUE)
  length(res) == 0
}

#' @export
identical.formula <- function(x, y, ...) {
  res <- waldo::compare(x, y, ignore_formula_env = TRUE)
  length(res) == 0
}

#' View or change global bmm options
#' @param sort_data logical. If TRUE, the data will be sorted by the predictors.
#'   If FALSE, the data will not be sorted, but sampling will be slower. If
#'   "check" (the default), `bmm()` will check if the data is sorted, and
#'   ask you via a console prompt if it should be sorted. **Default: "check"**
#' @param parallel logical. If TRUE, chains will be run in parallel. If FALSE,
#'   chains will be run sequentially. You can also set these value for each
#'   model separately via the argument `parallel` in `bmm()`. **Default:
#'   FALSE**
#' @param default_priors logical. If TRUE (default), the default bmm priors will
#'   be used. If FALSE, only the basic `brms` priors will be used. **Default:
#'   TRUE**
#' @param silent numeric. Verbosity level between 0 and 2. If 1 ( the default),
#'   most of the informational messages of compiler and sampler are suppressed.
#'   If 2, even more messages are suppressed. The actual sampling progress is
#'   still printed. **Default: 1**
#' @param color_summary logical. If TRUE, the summary of the model will be
#'  printed in color. **Default: TRUE**
#' @param reset_options logical. If TRUE, the options will be reset to their
#'   default values **Default: FALSE**
#' @param file_refit logical. If TRUE, bmm() will refit the model even if the
#'  file argument is specified. **Default: FALSE**
#' @details The `bmm_options` function is used to view or change the current bmm
#'   options. If no arguments are provided, the function will return the current
#'   options. If arguments are provided, the function will change the options
#'   and return the old options invisibly. If you provide only some of the
#'   arguments, the other options will not be changed. The options are stored in
#'   the global options list and will be used by `bmm()` and other
#'   functions in the `bmm` package. Each of these options can also be set
#'   manually using the built-in `options()` function, by setting the
#'   `bmm.sort_data`,  `bmm.default_priors`, and `bmm.silent` options.
#'
#' @return A message with the current bmm options and their values, and
#'   invisibly returns the old options for use with on.exit() and friends.
#' @examples
#'
#' # view the current options
#' bmm_options()
#'
#' # change the options to always sort the data and to use parallel sampling
#' bmm_options(sort_data = TRUE, parallel = TRUE)
#'
#' # restore the default options
#' bmm_options(reset_options = TRUE)
#'
#' # you can change the options using the options() function as well
#' options(bmm.sort_data = TRUE, bmm.parallel = TRUE)
#' bmm_options()
#'
#' # reset the options to their default values
#' bmm_options(reset_options = TRUE)
#'
#' # bmm_options(sort_data = TRUE, parallel = TRUE) will also return the old options
#' # so you can use it with on.exit()
#' old_op <- bmm_options(sort_data = TRUE, parallel = TRUE)
#' on.exit(bmm_options(old_op))
#'
#' bmm_options(reset_options = TRUE)
#' @export
bmm_options <- function(sort_data, parallel, default_priors, silent,
                        color_summary, file_refit, reset_options = FALSE) {
  opts <- ls()
  stopif(
    !missing(sort_data) && sort_data != "check" && !is.logical(sort_data),
    "sort_data must be one of TRUE, FALSE, or 'check'"
  )
  stopif(
    !missing(parallel) && !is.logical(parallel),
    "parallel must either TRUE or FALSE"
  )
  stopif(
    !missing(default_priors) && !is.logical(default_priors),
    "default_priors must be either TRUE or FALSE"
  )
  stopif(
    !missing(silent) && (!is.numeric(silent) || silent < 0 || silent > 2),
    "silent must be one of 0, 1, or 2"
  )
  stopif(
    !missing(color_summary) && !is.logical(color_summary),
    "color_summary must be either TRUE or FALSE"
  )
  stopif(
    !missing(file_refit) && !is.logical(file_refit),
    "file_refit must be either TRUE or FALSE"
  )

  # set default options if function is called for the first time or if reset_options is TRUE
  if (reset_options) {
    options(
      bmm.sort_data = "check",
      bmm.parallel = FALSE,
      bmm.default_priors = TRUE,
      bmm.silent = 1,
      bmm.color_summary = TRUE,
      bmm.file_refit = FALSE
    )
  }

  # change options if arguments are provided. get argument name and loop over non-missing arguments
  op <- list()
  non_missing_args <- names(match.call())[-1]
  non_missing_args <- non_missing_args[!non_missing_args %in% "reset_options"]
  for (i in non_missing_args) {
    op[[paste0("bmm.", i)]] <- get(i)
  }

  old_op <- options(op)
  message2(
    "\nCurrent bmm options:\n",
    crayon::green(paste0(
      "  sort_data = ", getOption("bmm.sort_data"), "",
      "\n  parallel = ", getOption("bmm.parallel"),
      "\n  default_priors = ", getOption("bmm.default_priors"),
      "\n  silent = ", getOption("bmm.silent"),
      "\n  file_refit = ", getOption("bmm.file_refit"),
      "\n  color_summary = ", getOption("bmm.color_summary"), "\n"
    )),
    "For more information on these options or how to change them, see help(bmm_options).\n"
  )
  invisible(old_op)
}

# an improved version of tryCatch that captures messages as well
# modified version of https://github.com/cran/admisc/blob/master/R/tryCatchWEM.R
tryCatch2 <- function(expr, capture = FALSE) {
  toreturn <- list()
  output <- withVisible(withCallingHandlers(
    tryCatch(expr, error = function(e) {
      toreturn$error <<- e$message
      NULL
    }),
    warning = function(w) {
      toreturn$warning <<- c(toreturn$warning, w$message)
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      toreturn$message <<- paste(toreturn$message, m$message, sep = "")
      invokeRestart("muffleMessage")
    }
  ))
  if (capture && output$visible && !is.null(output$value)) {
    toreturn$output <- utils::capture.output(output$value)
    toreturn$value <- output$value
  }
  if (length(toreturn) > 0) {
    return(toreturn)
  }
}

# resets the environments stored within an objects
reset_env <- function(object, env = NULL, ...) {
  UseMethod("reset_env")
}

#' @export
reset_env.bmmfit <- function(object, env = globalenv(), ...) {
  object$formula <- reset_env(object$formula, env)
  object$family <- reset_env(object$family, env)
  object$bmm$user_formula <- reset_env(object$bmm$user_formula, env)
  object
}

#' @export
reset_env.bmmformula <- function(object, env = globalenv(), ...) {
  # we loop this way instead of lapply or map because the latter
  # will drop the attributes and class of a bmmformula
  for (par in names(object)) {
    object[[par]] <- reset_env(object[[par]], env)
  }
  object
}

#' @export
reset_env.brmsformula <- function(object, env = globalenv(), ...) {
  object$formula <- reset_env(object$formula, env)
  for (par in names(object$pforms)) {
    object$pforms[[par]] <- reset_env(object$pforms[[par]], env)
  }
  if (!is.null(object$family)) {
    object$family <- reset_env(object$family, env)
  }
  object
}

#' @export
reset_env.formula <- function(object, env = globalenv(), ...) {
  environment(object) <- env
  object
}

#' @export
reset_env.brmsfamily <- function(object, env = globalenv(), ...) {
  if (!is.null(object$env)) {
    object$env <- env
  }
  object
}

#' @export
reset_env.default <- function(object, env = NULL, ...) {
  object
}

# Remove all attributes of an object except those specified as protected
# @param x an R object
# @param protect a character vector of attribute names to keep. Default is
#   `c("names", "row.names", "class")`, which are the attributes that a
#   data.frame has by default.
# @return An R object with all attributes removed except those specified in
#   `protect`.
strip_attributes <- function(x, protect = c("names", "row.names", "class"),
                             recursive = FALSE) {
  to_remove <- names(attributes(x))
  to_remove <- to_remove[!to_remove %in% protect]
  attributes(x)[to_remove] <- NULL
  if (recursive && is.list(x)) {
    for (i in seq_along(x)) {
      x[[i]] <- strip_attributes(x[[i]], protect, recursive)
    }
  }
  x
}

deprecated_args <- function(...) {
  dots <- list(...)
  stopif(
    "model_type" %in% names(dots),
    'The "model_type" argument was deprecated on Feb 3, 2024. Either:
         - See ?bmm for the new usage;
         - or install the old version of the package with:
           devtools::install_github("venpopov/bmm@v0.0.1")'
  )
  warnif(
    "parallel" %in% names(dots),
    'The "parallel" argument is deprecated. Please use cores instead.
         See `help("brm")` for more information.'
  )
}

try_read_bmmfit <- function(file, file_refit) {
  if (is.character(file_refit)) {
    stopif(
      !tolower(file_refit) %in% c("never", "always", "on_change"),
      'You have provided an invalid option for the file_refit argument.
      Valid options are: "never" or "always"
      The "on_change" option available in brms is not yet implemented'
    )

    warnif(
      tolower(file_refit) == "on_change",
      'The "on_change" option for the file_refit argument available in brms,
      is currently not implemented for bmm.
      To avoid overwriting an already saved bmmfit object, file_refit was set to "never".'
    )
    file_refit <- ifelse(file_refit == "always", TRUE, FALSE)
  }
  file <- check_rds_file(file)
  if (is.null(file) || file_refit) {
    return(NULL)
  }
  dir <- try(fs::dir_create(dirname(file)))
  stopif(is_try_error(dir), "Cannot create directory for file.")

  out <- suppressWarnings(try(readRDS(file), silent = TRUE))
  if (is_try_error(out)) {
    return(NULL)
  }

  stopif(!is_bmmfit(out), "Object loaded via 'file' is not of class 'bmmfit'.")
  out$file <- file
  out
}

try_save_bmmfit <- function(object, file, compress) {
  file <- check_rds_file(file)
  object$file <- file
  if (!is.null(file)) {
    saveRDS(object, file, compress = compress)
  }
  object
}

check_rds_file <- function(file) {
  if (is.null(file)) {
    return(NULL)
  }
  stopif(!is.character(file) || length(file) > 1, "'file' must be a single string.")
  ext <- fs::path_ext(file)
  if (ext != "rds") {
    file <- paste0(file, ".rds")
  }
  file
}

`%||%` <- function(l, r) {
  if (!is.null(l)) l else r
}

# like unlist, but keeps the final outcome a list of all
# elements of nested lists. Only works 1-level deep
unnest_list <- function(list_of_lists) {
  do.call(c, list_of_lists)
}
