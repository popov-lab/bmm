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
# PRE-FIT DATA CHECK                                                     ####
############################################################################# !

#' @title Human-readable pre-fit data report for bmm models
#' @description Inspect whether your data is coded the way a `bmmodel` expects
#'   *before* committing to Stan compilation and sampling. The function runs
#'   the same validation pipeline that [bmm()] runs internally
#'   ([check_data()] and [check_formula()]), but captures all errors,
#'   warnings and messages instead of throwing them, and combines them with a
#'   readable summary of:
#'
#'   - the response variables, their observed values and the coding the model
#'     expects
#'   - the data columns each parameter formula uses, and how they are coded
#'     (factor levels, character vectors, numeric variables with few unique
#'     values)
#'   - the number of observations per design cell (crossing the grouping
#'     variables from random effects terms with the categorical predictors),
#'     flagging sparse or empty cells
#'   - model-specific diagnostics for common data mistakes (e.g. responses in
#'     degrees rather than radians for circular models, or misplaced `NA`
#'     values in the non-target features of set size varying designs).
#'     These are provided by [data_check_findings()] methods, which new models
#'     can extend
#'
#'   Problems are reported as findings rather than errors, so a single call
#'   shows everything that needs fixing at once.
#' @inheritParams bmm
#' @param min_trials Numeric. Design cells with fewer observations are flagged
#'   in the report. Defaults to 10. Ignored for models fit to aggregated data
#'   (e.g. `m3`, `ezdm`), where one row summarizes many trials.
#' @return An object of class `bmm_data_check` with a print method. The object
#'   is a list containing the summarized `response`, `predictors` and `cells`
#'   information (including the full cell count table in `$cells$counts`), all
#'   `findings`, and the captured results of the validation `pipeline`.
#' @seealso [bmm()], [check_data()], [data_check_findings()]
#' @keywords extract_info
#' @examples
#' # three-parameter mixture model: radians check, NA padding of the
#' # non-target features across set sizes, trials per participant x set size
#' bmm_data_check(
#'   bmf(kappa ~ 0 + set_size + (0 + set_size | ID), thetat ~ 0 + set_size),
#'   data = oberauer_lin_2017,
#'   model = mixture3p(
#'     resp_error = "dev_rad",
#'     nt_features = paste0("col_nt", 1:7),
#'     set_size = "set_size"
#'   )
#' )
#'
#' # interference measurement model: additionally checks the NA padding
#' # of the non-target distances
#' bmm_data_check(
#'   bmf(c ~ 1, a ~ 1, kappa ~ 1, s ~ 1),
#'   data = oberauer_lin_2017,
#'   model = imm(
#'     resp_error = "dev_rad",
#'     nt_features = paste0("col_nt", 1:7),
#'     nt_distances = paste0("dist_nt", 1:7),
#'     set_size = "set_size",
#'     version = "full"
#'   )
#' )
#'
#' # signal discrimination model
#' bmm_data_check(
#'   bmf(c ~ 0 + setsize, kappa ~ 0 + setsize),
#'   data = zhang_luck_2008,
#'   model = sdm(resp_error = "response_error")
#' )
#'
#' # memory measurement model: response category counts as response variables
#' bmm_data_check(
#'   bmf(c ~ 1 + (1 | ID), a ~ 1 + (1 | ID)),
#'   data = oberauer_lewandowsky_2019_e1,
#'   model = m3(
#'     resp_cats = c("corr", "other", "npl"),
#'     num_options = c("n_corr", "n_other", "n_npl"),
#'     choice_rule = "simple",
#'     version = "ss"
#'   )
#' )
#'
#' # diffusion decision model: reaction time plausibility and response coding
#' bmm_data_check(
#'   bmf(drift ~ condition + (condition | ID)),
#'   data = data_color_judgement_task,
#'   model = ddm(rt = "rt", response = "response_correct")
#' )
#' @export
bmm_data_check <- function(formula, data, model, min_trials = 10) {
  stop_missing_args()
  stopif(
    !is_bmmformula(formula),
    "The formula argument must be a bmmformula created with bmmformula() or bmf()"
  )
  stopif(
    !is.numeric(min_trials) || length(min_trials) != 1 || min_trials < 1,
    "Argument 'min_trials' must be a single positive number"
  )
  withr::local_options(bmm.sort_data = FALSE)

  data_name <- substitute_name(data)
  data <- try(as.data.frame(data), silent = TRUE)
  stopif(is_try_error(data), "Argument 'data' must be coercible to a data.frame.")
  stopif(!isTRUE(nrow(data) > 0L), "Argument 'data' does not contain observations.")

  model <- check_model(model, data, formula)
  pipeline <- capture_check_conditions({
    checked_data <- check_data(model, data, formula)
    check_formula(model, checked_data, formula)
  })

  response <- summarise_response_vars(model, data)
  predictors <- summarise_predictor_vars(model, data, formula)
  cells <- summarise_design_cells(data, formula)
  findings <- c(
    response_var_findings(response),
    generic_data_findings(
      predictors, cells,
      min_trials = if (uses_aggregate_data(model)) NULL else min_trials
    ),
    data_check_findings(model, data, formula)
  )

  structure(
    nlist(
      model, data_name,
      n_obs = nrow(data),
      response, predictors, cells, findings, pipeline
    ),
    class = "bmm_data_check"
  )
}

capture_check_conditions <- function(expr) {
  warnings <- character(0)
  messages <- character(0)
  error <- NULL
  withCallingHandlers(
    tryCatch(expr, error = function(e) error <<- conditionMessage(e)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )
  nlist(error, warnings, messages)
}

summarise_response_vars <- function(model, data) {
  annotations <- response_annotations(model)
  info <- lapply(names(model$resp_vars), function(arg) {
    vars <- model$resp_vars[[arg]]
    data.frame(
      variable = vars,
      summary = vapply(vars, function(v) describe_data_column(data[[v]]), character(1)),
      expected = annotations[[arg]] %||% NA_character_,
      row.names = NULL
    )
  })
  do.call(rbind, info)
}

describe_data_column <- function(x, max_values = 6) {
  if (is.null(x)) {
    return("not found in the data")
  }
  n_na <- sum(is.na(x))
  na_info <- if (n_na > 0) glue(", {n_na} NA") else ""
  if (all(is.na(x))) {
    return(glue("{class(x)[1]}, all values NA"))
  }
  if (is.numeric(x)) {
    values <- sort(unique(x[!is.na(x)]))
    if (length(values) <= max_values) {
      return(glue("numeric with {length(values)} unique values: {paste(values, collapse = ', ')}{na_info}"))
    }
    rng <- signif(range(values), 3)
    return(glue("numeric, range [{rng[1]}, {rng[2]}]{na_info}"))
  }
  if (is.logical(x)) {
    return(glue("logical, {sum(x, na.rm = TRUE)} TRUE / {sum(!x, na.rm = TRUE)} FALSE{na_info}"))
  }
  values <- if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
  shown <- collapse_comma(as.character(utils::head(values, max_values)))
  more <- if (length(values) > max_values) ", ..." else ""
  glue("{class(x)[1]} with {length(values)} levels: {shown}{more}{na_info}")
}

summarise_predictor_vars <- function(model, data, formula) {
  par_names <- unique(c(names(formula), names(model$parameters)))
  pred_map <- lapply(rhs_vars(formula, collapse = FALSE), function(vars) {
    setdiff(vars, par_names)
  })
  pred_map <- pred_map[lengths(pred_map) > 0]

  used_vars <- unique(unlist(pred_map))
  group_vars <- intersect(re_group_vars(formula), colnames(data))
  data_vars <- intersect(used_vars, colnames(data))

  coding <- data.frame(
    variable = data_vars,
    role = ifelse(data_vars %in% group_vars, "grouping", "predictor"),
    class = vapply(data_vars, function(v) class(data[[v]])[1], character(1)),
    n_unique = vapply(data_vars, function(v) {
      length(unique(data[[v]][!is.na(data[[v]])]))
    }, integer(1)),
    summary = vapply(data_vars, function(v) describe_data_column(data[[v]]), character(1)),
    row.names = NULL
  )
  nlist(pred_map, coding, unknown_vars = setdiff(used_vars, colnames(data)), group_vars)
}

re_group_vars <- function(formula) {
  vars <- lapply(formula, function(f) {
    if (!is_formula(f) || length(f) == 0) {
      return(character(0))
    }
    extract_bar_group_vars(f[[length(f)]])
  })
  as.character(unique(unlist(vars)))
}

extract_bar_group_vars <- function(expr) {
  if (!is.call(expr)) {
    return(character(0))
  }
  if (identical(expr[[1]], quote(`|`)) || identical(expr[[1]], quote(`||`))) {
    return(all.vars(expr[[3]]))
  }
  unique(unlist(lapply(as.list(expr)[-1], extract_bar_group_vars)))
}

summarise_design_cells <- function(data, formula) {
  group_vars <- intersect(re_group_vars(formula), colnames(data))
  fixed_vars <- setdiff(intersect(rhs_vars(formula), colnames(data)), group_vars)
  cell_vars <- fixed_vars[vapply(fixed_vars, function(v) {
    is_discrete_var(data[[v]])
  }, logical(1))]

  vars <- c(group_vars, cell_vars)
  if (length(vars) == 0) {
    return(nlist(group_vars, cell_vars, counts = NULL))
  }

  n_combos <- prod(vapply(vars, function(v) {
    length(unique(data[[v]]))
  }, numeric(1)))
  if (n_combos > 1e5) {
    return(nlist(group_vars, cell_vars, counts = NULL))
  }

  counts <- as.data.frame(
    table(data[, vars, drop = FALSE]),
    responseName = "n", stringsAsFactors = FALSE
  )
  nlist(group_vars, cell_vars, counts)
}

is_discrete_var <- function(x, max_unique = 10) {
  if (is.numeric(x)) {
    return(length(unique(x[!is.na(x)])) <= max_unique)
  }
  is.factor(x) || is.character(x) || is.logical(x)
}

data_check_finding <- function(severity, message) {
  nlist(severity, message = as.character(message))
}

response_var_findings <- function(response) {
  missing_vars <- response$variable[response$summary == "not found in the data"]
  if (length(missing_vars) == 0) {
    return(list())
  }
  list(data_check_finding("warning", glue(
    "The response variable(s) {collapse_comma(missing_vars)} are not present \\
    in the data."
  )))
}

generic_data_findings <- function(predictors, cells, min_trials) {
  findings <- list()
  if (length(predictors$unknown_vars) > 0) {
    findings <- c(findings, list(data_check_finding("warning", glue(
      "The formula uses {collapse_comma(predictors$unknown_vars)}, which are \\
      neither columns in the data nor model parameters."
    ))))
  }

  coding <- predictors$coding
  pred_rows <- coding[coding$role == "predictor", ]
  chr_vars <- pred_rows$variable[pred_rows$class == "character"]
  if (length(chr_vars) > 0) {
    findings <- c(findings, list(data_check_finding("note", glue(
      "Character predictor(s) {collapse_comma(chr_vars)} will be coerced to \\
      factors with alphabetical level order. Use factor() to control the \\
      reference level."
    ))))
  }
  num_vars <- pred_rows$variable[pred_rows$class %in% c("numeric", "integer") &
    pred_rows$n_unique <= 6]
  if (length(num_vars) > 0) {
    findings <- c(findings, list(data_check_finding("note", glue(
      "Numeric predictor(s) {collapse_comma(num_vars)} have few unique values \\
      and will be treated as continuous (a single slope). If they code \\
      discrete conditions, convert them to factors for one estimate per level."
    ))))
  }

  c(findings, cell_count_findings(cells$counts, min_trials))
}

# aggregate-response models code many trials per row, so row counts per
# design cell say nothing about trial numbers and min_trials does not apply
uses_aggregate_data <- function(model) {
  inherits(model, "m3") || inherits(model, "ezdm")
}

cell_count_findings <- function(counts, min_trials) {
  if (is.null(counts)) {
    return(list())
  }
  findings <- list()
  low <- if (is.null(min_trials)) {
    counts[0, , drop = FALSE]
  } else {
    counts[counts$n > 0 & counts$n < min_trials, , drop = FALSE]
  }
  if (nrow(low) > 0) {
    findings <- c(findings, list(data_check_finding("note", glue(
      "{nrow(low)} design cell(s) have fewer than {min_trials} observations: \\
      {format_cell_labels(utils::head(low, 3))}\\
      {if (nrow(low) > 3) ', ...' else ''}"
    ))))
  }
  empty <- counts[counts$n == 0, , drop = FALSE]
  if (nrow(empty) > 0) {
    findings <- c(findings, list(data_check_finding("note", glue(
      "{nrow(empty)} design cell combination(s) contain no observations - \\
      fine if the design is intentionally incomplete."
    ))))
  }
  findings
}

format_cell_labels <- function(counts) {
  vars <- setdiff(names(counts), "n")
  labels <- apply(counts[vars], 1, function(row) {
    paste(paste0(vars, "=", row), collapse = ", ")
  })
  paste0(labels, " (n=", counts$n, ")", collapse = "; ")
}

#' @title Generic S3 method for model-specific pre-fit data findings
#' @description Called by [bmm_data_check()] to collect model-specific
#'   diagnostics of common data coding mistakes. Like [check_data()], methods
#'   dispatch over the model class chain from general to specific, and every
#'   method should concatenate its own findings with `NextMethod()`. Each
#'   finding is a list with elements `severity` ("warning" or "note") and
#'   `message`. Unlike [check_data()], these diagnostics never throw - they
#'   only describe likely problems, so a report can always be produced.
#' @param model A `bmmodel` object
#' @param data The user supplied data.frame, before any transformations by
#'   [check_data()]
#' @param formula The user supplied `bmmformula`
#' @return A list of findings (possibly empty)
#' @keywords internal developer
#' @export
data_check_findings <- function(model, data, formula) {
  UseMethod("data_check_findings")
}

#' @export
data_check_findings.default <- function(model, data, formula) {
  list()
}

#' @export
data_check_findings.circular <- function(model, data, formula) {
  c(
    circular_range_findings(data, model$resp_vars[[1]], "response"),
    NextMethod()
  )
}

#' @export
data_check_findings.non_targets <- function(model, data, formula) {
  nt_features <- model$other_vars$nt_features
  c(
    circular_range_findings(data, nt_features, "non-target feature"),
    set_size_padding_findings(data, nt_features, model$other_vars$set_size, "nt_features"),
    NextMethod()
  )
}

#' @export
data_check_findings.imm_bsc <- function(model, data, formula) {
  c(
    set_size_padding_findings(
      data, model$other_vars$nt_distances, model$other_vars$set_size, "nt_distances"
    ),
    NextMethod()
  )
}

#' @export
data_check_findings.imm_full <- function(model, data, formula) {
  c(
    set_size_padding_findings(
      data, model$other_vars$nt_distances, model$other_vars$set_size, "nt_distances"
    ),
    NextMethod()
  )
}

circular_range_findings <- function(data, vars, label) {
  vars <- intersect(vars, colnames(data))
  if (length(vars) == 0 || all(is.na(data[, vars]))) {
    return(list())
  }
  max_abs <- max(abs(data[, vars]), na.rm = TRUE)
  if (max_abs > 2 * pi) {
    return(list(data_check_finding("warning", glue(
      "The {label} variable(s) {collapse_comma(vars)} contain values up to \\
      {round(max_abs, 1)} - almost certainly degrees. The model requires \\
      radians in [-pi, pi]; convert with deg2rad() and wrap()."
    ))))
  }
  if (max_abs > pi) {
    return(list(data_check_finding("warning", glue(
      "The {label} variable(s) {collapse_comma(vars)} contain values in \\
      (pi, 2*pi], suggesting angles coded on [0, 2*pi) instead of the \\
      required [-pi, pi]. Recode with wrap()."
    ))))
  }
  list()
}

# in set size varying designs the nt columns must contain values at positions
# up to set_size - 1 and NA padding after that; check_data silently replaces
# NAs, so misplaced ones are only visible here
set_size_padding_findings <- function(data, cols, set_size, arg) {
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    return(list(data_check_finding("warning", glue(
      "The '{arg}' variable(s) {collapse_comma(missing_cols)} are not present in the data."
    ))))
  }
  ss <- try(check_var_set_size(set_size, data), silent = TRUE)
  if (is_try_error(ss) || length(cols) != ss$max_set_size - 1) {
    return(list())
  }

  na_mat <- as.matrix(is.na(data[, cols, drop = FALSE]))
  position <- matrix(seq_along(cols), nrow(data), length(cols), byrow = TRUE)
  active <- position <= (ss$ss_numeric - 1)
  padding <- position > (ss$ss_numeric - 1)
  active[is.na(active)] <- FALSE
  padding[is.na(padding)] <- FALSE

  findings <- list()
  n_na_active <- sum(na_mat & active)
  if (n_na_active > 0) {
    rows <- which(rowSums(na_mat & active) > 0)
    findings <- c(findings, list(data_check_finding("warning", glue(
      "{n_na_active} NA value(s) in '{arg}' at positions within the set size \\
      of the trial ({length(rows)} row(s), first at row {rows[1]}). These \\
      positions must contain values - during fitting NAs are silently \\
      replaced by placeholder values, which biases the estimates."
    ))))
  }
  n_value_padding <- sum(!na_mat & padding)
  if (n_value_padding > 0) {
    rows <- which(rowSums(!na_mat & padding) > 0)
    findings <- c(findings, list(data_check_finding("note", glue(
      "{n_value_padding} value(s) in '{arg}' at positions beyond set_size - 1 \\
      ({length(rows)} row(s)). Positions beyond the set size of a trial are \\
      ignored during fitting - make sure the set_size variable is correct."
    ))))
  }
  findings
}

#' @export
print.bmm_data_check <- function(x, color = getOption("bmm.color_summary", TRUE), ...) {
  withr::local_options(bmm.color_summary = color)
  cat("Pre-fit data check\n")
  cat(style("purple1")("Model:"), gsub("\\s+", " ", construct_model_call(x$model)), "\n")
  cat(style("purple1")("Data: "), glue("{x$n_obs} observations ('{x$data_name}')"), "\n\n")
  print_response_section(x$response)
  print_predictor_section(x$predictors)
  print_cells_section(x$cells)
  print_findings_section(x$findings)
  print_pipeline_section(x$pipeline)
  invisible(x)
}

print_response_section <- function(response) {
  cat(style("green")("Response variables:\n"))
  labels <- format(response$variable)
  for (i in seq_len(nrow(response))) {
    expected <- if (is.na(response$expected[i])) {
      ""
    } else {
      glue(" (expected: {response$expected[i]})")
    }
    cat("  ", labels[i], "  ", response$summary[i], expected, "\n", sep = "")
  }
  cat("\n")
}

print_predictor_section <- function(predictors) {
  if (length(predictors$pred_map) == 0) {
    return(invisible())
  }
  cat(style("green")("Data columns used by each parameter formula:\n"))
  pars <- format(names(predictors$pred_map))
  for (i in seq_along(predictors$pred_map)) {
    cat("  ", pars[i], "  ", paste(predictors$pred_map[[i]], collapse = ", "), "\n", sep = "")
  }
  coding <- predictors$coding
  if (nrow(coding) > 0) {
    cat(style("green")("\nPredictor coding:\n"))
    labels <- format(coding$variable)
    roles <- format(paste0("[", coding$role, "]"))
    for (i in seq_len(nrow(coding))) {
      cat("  ", labels[i], "  ", roles[i], " ", coding$summary[i], "\n", sep = "")
    }
  }
  cat("\n")
}

print_cells_section <- function(cells) {
  if (is.null(cells$counts)) {
    return(invisible())
  }
  observed <- cells$counts$n[cells$counts$n > 0]
  cat(style("green")("Observations per design cell:\n"))
  cat("  ", paste(c(cells$group_vars, cells$cell_vars), collapse = " x "),
    ": ", length(observed), " non-empty cell(s)\n",
    sep = ""
  )
  cat("  observations per cell: min ", min(observed),
    " / median ", stats::median(observed), " / max ", max(observed), "\n\n",
    sep = ""
  )
}

print_findings_section <- function(findings) {
  if (length(findings) == 0) {
    cat("Findings: ", style("green")("No issues detected."), "\n\n", sep = "")
    return(invisible())
  }
  cat("Findings:\n")
  severities <- vapply(findings, `[[`, character(1), "severity")
  for (f in findings[order(severities != "warning")]) {
    if (f$severity == "warning") {
      cat_wrapped(f$message, "  ! ", color = "red")
    } else {
      cat_wrapped(f$message, "  - ")
    }
  }
  cat("\n")
}

print_pipeline_section <- function(pipeline) {
  cat("Hard checks: ")
  if (!is.null(pipeline$error)) {
    cat("\n")
    cat_wrapped(glue("FAILED: {pipeline$error}"), "  ! ", color = "red")
  } else if (length(pipeline$warnings) == 0) {
    cat(style("green")("passed"), "\n", sep = "")
  } else {
    cat("passed with ", length(pipeline$warnings), " warning(s):\n", sep = "")
  }
  for (w in pipeline$warnings) {
    cat_wrapped(w, "  ! ", color = "red")
  }
  for (m in pipeline$messages) {
    cat_wrapped(m, "  - ")
  }
  invisible()
}

cat_wrapped <- function(text, prefix, color = NULL) {
  text <- gsub("\\s+", " ", text)
  lines <- strwrap(text, width = 80, initial = prefix, exdent = nchar(prefix))
  if (!is.null(color)) {
    lines <- style(color)(lines)
  }
  cat(lines, sep = "\n")
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
  dots <- list(...)
  local_brms_threads(dots)

  # check model, formula and data, and transform data if necessary
  formula <- object
  configure_options(dots)
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)

  # generate the model specification to pass to brms later
  config_args <- configure_model(model, data, formula)

  # extract stan data
  fit_args <- combine_args(nlist(config_args, dots))
  fit_args$object <- fit_args$formula
  fit_args$formula <- NULL
  brms::do_call(brms::standata, fit_args)
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
#'   assumptions, but note that the EZ equations were derived for mean and
#'   variance, so using median may introduce some bias for skewed distributions.
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
#' @param min_trials Integer. Minimum number of trials required for fitting.
#'   Returns NA if fewer trials are available. Default is 10
#' @param init_contaminant Numeric. Initial proportion of contaminants for EM
#'   algorithm. Default is 0.05
#' @param max_contaminant Numeric. Maximum allowed contaminant proportion
#'   (0 < max <= 1). Estimates are clipped to this value to prevent inflated
#'   contaminant proportions. Default is 0.5
#' @param maxit Integer. Maximum number of EM iterations. Default is 100
#' @param tol Numeric. Convergence tolerance for EM algorithm. Default is 1e-6
#'
#' @return A 1-row `data.frame`. For version = "3par": `mean_rt`, `var_rt`,
#'   `n_upper`, `n_trials`, `contaminant_prop`. For version = "4par":
#'   `mean_rt_upper`, `mean_rt_lower`, `var_rt_upper`, `var_rt_lower`,
#'   `n_upper`, `n_trials`, `contaminant_prop_upper`, `contaminant_prop_lower`.
#'
#' @details RT outliers and contaminant responses (fast guesses, lapses of
#'   attention) can distort the mean and variance estimates used as input to
#'   the EZ-Diffusion equations. This function addresses this by fitting a
#'   mixture model with two components: a uniform distribution for
#'   contaminants and a parametric RT distribution for true responses.
#'   Robust moments are then extracted from the fitted parametric component.
#'
#'   This function is designed to work with [dplyr::group_by()] and
#'   [dplyr::reframe()] for grouped operations. Use [adjust_ezdm_accuracy()]
#'   as a separate step if you need to adjust accuracy counts for
#'   contamination.
#'
#' @seealso [adjust_ezdm_accuracy()] for adjusting accuracy counts,
#'   [flag_contaminant_rts()] for trial-level contamination probabilities,
#'   [ezdm()] for fitting the EZ-Diffusion Model
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
    maxit = 100,
    tol = 1e-6) {
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
        warning2("EM did not converge. Using robust moments.")
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
    data.frame(
      mean_rt = moments$mean,
      var_rt = moments$var,
      n_upper = n_upper,
      n_trials = n_trials,
      contaminant_prop = moments$contaminant_prop
    )
  } else {
    rt_upper <- rt[is_upper]
    rt_lower <- rt[!is_upper]
    moments_upper <- compute_moments(rt_upper, length(rt_upper))
    moments_lower <- compute_moments(rt_lower, length(rt_lower))
    data.frame(
      mean_rt_upper = moments_upper$mean,
      mean_rt_lower = moments_lower$mean,
      var_rt_upper = moments_upper$var,
      var_rt_lower = moments_lower$var,
      n_upper = n_upper,
      n_trials = n_trials,
      contaminant_prop_upper = moments_upper$contaminant_prop,
      contaminant_prop_lower = moments_lower$contaminant_prop
    )
  }
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

#' Adjust Accuracy Counts for Contamination
#'
#' @description Adjusts accuracy counts (`n_upper`, `n_trials`) by removing
#'   estimated contaminant trials using binomial sampling. Contaminant trials
#'   are assumed to produce correct responses at a fixed guess rate (e.g., 0.5
#'   for 2AFC tasks).
#'
#' @param n_upper Numeric. Count of upper boundary (correct) responses.
#' @param n_trials Numeric. Total number of trials.
#' @param contaminant_prop Numeric. Estimated proportion of contaminant trials
#'   (e.g., from the `contaminant_prop` column of [ezdm_summary_stats()]).
#' @param guess_rate Numeric. Assumed accuracy rate for contaminant trials
#'   (random guessing). Default is 0.5 (appropriate for 2AFC tasks).
#'
#' @return A 1-row `data.frame` with columns `n_upper_adj` and `n_trials_adj`
#'   (integers). When `contaminant_prop` is `NA` or <= 0, returns the original
#'   counts unchanged.
#'
#' @details Uses binomial sampling to estimate the number of contaminant trials
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
#' # Adjust accuracy for estimated 10% contamination
#' set.seed(42)
#' adjust_ezdm_accuracy(n_upper = 80, n_trials = 100, contaminant_prop = 0.1)
#'
#' # In a pipeline with ezdm_summary_stats
#' # library(dplyr)
#' # mydata |>
#' #   group_by(subject) |>
#' #   reframe(ezdm_summary_stats(rt, response)) |>
#' #   mutate(adjust_ezdm_accuracy(n_upper, n_trials, contaminant_prop))
#'
adjust_ezdm_accuracy <- function(n_upper, n_trials, contaminant_prop,
                                 guess_rate = 0.5) {
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
#' @param maxit Integer. Maximum EM iterations. Default 100.
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
#'   | RT), with a `"diagnostics"` attribute containing a one-row data.frame
#'   with columns: `mixture_params` (list), `contaminant_prop`, `converged`,
#'   `iterations`, `loglik`, `n_trials`, `distribution`, `method`.
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
    maxit = 100,
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
    warning2("EM did not converge. Returning NA.")
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
