############################################################################# !
# TRIAL-LEVEL DATA FOR THE COUNT MODELS                                  ####
############################################################################# !

# The count models take either cells of counts or one row per trial. The
# constructor decides which, from the arguments the user names, and records it
# in other_vars$trial_level; data are never sniffed, because 0/1 counts with
# one trial per cell are legal aggregated data. A model without the marker
# takes cells only (m3, ezdm) or is a trial model anyway (circular, RT models)
is_trial_level <- function(model) {
  isTRUE(model$other_vars$trial_level)
}

# The response of each trial as an integer one-hot matrix with one column per
# category. NA rows stay NA, as brms drops them. The binomial models need only
# the validation and the column of the second category; the multinomial models
# use the matrix as their response
.trial_rows_to_counts <- function(values, categories, what) {
  idx <- match(values, categories)
  bad <- unique(values[is.na(idx) & !is.na(values)])
  stopif(
    length(bad) > 0,
    "With one row per trial, {what} must hold one of {collapse_comma(categories)} \\
    in every row. Found: {collapse_comma(utils::head(bad, 5))}\\
    {if (length(bad) > 5) ', ...' else ''}. For counts aggregated over trials, \\
    construct the model with n_trials."
  )
  counts <- matrix(0L, nrow = length(values), ncol = length(categories),
                   dimnames = list(NULL, as.character(categories)))
  counts[cbind(which(!is.na(idx)), idx[!is.na(idx)])] <- 1L
  counts[is.na(idx), ] <- NA_integer_
  counts
}

# the binomial models count trials in the column the user named, or in the
# constant column check_data() adds for trial rows
.sdt_trials_var <- function(model) {
  model$other_vars$n_trials %||% "nTrials"
}

# The binomial likelihoods take all rows at once (loop = FALSE), and brms
# slices only the response under reduce_sum, so their per-row arrays (vint1,
# trials) are sliced here exactly when brms slices (see
# brms_slices_likelihood). vint2 is dist_type, one code for every row
.sdt_binomial_row_vars <- function() {
  if (brms_slices_likelihood()) {
    return(c("vint1[start:end]", "vint2[1]", "trials[start:end]"))
  }
  c("vint1", "vint2[1]", "trials")
}

# trial rows of the binomial models become cells of one trial each. nTrials is
# the name the multinomial models already reserve for their trial counts; a
# column the model reads under that name would be overwritten with 1s
.add_binomial_trials <- function(data, model, formula) {
  resp_var <- model$resp_vars$response
  stopif(not_in(resp_var, colnames(data)), "Variables '{resp_var}' missing in the data")
  stopif(
    "nTrials" %in% c(unlist(model$resp_vars), unlist(model$other_vars), rhs_vars(formula)),
    "With one row per trial, {model$name} adds a column 'nTrials' to the data. \\
    Rename your column 'nTrials' and the argument or formula term that names it."
  )
  counts <- .trial_rows_to_counts(data[[resp_var]], c(0, 1),
                                  glue("the response column '{resp_var}'"))
  data[[resp_var]] <- unname(counts[, "1"])
  warnif("nTrials" %in% colnames(data),
         "Column 'nTrials' in your data is reserved by {model$name} and will \\
         be overwritten with 1 (one trial per row)")
  data$nTrials <- 1L
  data
}

# the "Trial-level data" part of the help page of every model that carries the
# marker; model_docs() appends it under the trial-level requirements, with the
# stimulus bullet only for models with a stimulus column
trial_level_docs <- paste0(
  "Without trial-level predictors or item terms, trial rows give the same ",
  "posterior as the counts aggregated from them, which are sufficient ",
  "statistics; the fit only takes longer, because its cost grows with the ",
  "number of rows. Crossed random effects such as `(1 | id) + (1 | item)` ",
  "need trial rows, unless the counts are kept per participant and item.\n\n"
)

trial_level_stimulus_docs <- paste0(
  "   - With `criterion` (the default), an item's effect on `d` is informed ",
  "only by the trials on which the item is the signal: on noise trials `d` ",
  "does not enter the likelihood, so the item effects on `d` of items shown ",
  "only as noise stay at their prior. With `bias`, `d` also moves the ",
  "criterion and enters every trial.\n"
)

trial_level_notes_docs <- paste0(
  "   - `loo()` leaves out single trials on trial rows and whole cells on ",
  "aggregated rows. These are different predictions, and the log-likelihoods ",
  "of the trial rows sum to the cell's log-likelihood less the log binomial ",
  "coefficient, so ELPD values compare only fits to the same layout. ",
  "`kfold(fit, group = \"id\")` leaves out whole participants with either ",
  "layout.\n",
  "   - `bmm(..., threads = n)` with `backend = \"cmdstanr\"` evaluates the ",
  "likelihood on `n` threads per chain, which shortens trial-level fits on a ",
  "machine with cores to spare.\n",
  "   - `ppp(fit, participant = \"id\")` names the participant when the formula ",
  "has more than one grouping variable, such as item terms.\n"
)
