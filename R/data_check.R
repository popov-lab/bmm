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
#'
#'   Defining a method is entirely optional: the default method returns an
#'   empty list, and [bmm_data_check()] produces its full generic report for
#'   models without any method. Add one only when a model has a common data
#'   mistake that the hard checks in [check_data()] deliberately tolerate.
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
  if (is.null(response) || nrow(response) == 0) {
    return(invisible())
  }
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
