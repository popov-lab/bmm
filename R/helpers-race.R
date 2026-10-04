############################################################################# !
# RACE MODELS: SHARED INFRASTRUCTURE                                     ####
############################################################################# !
# The race models (lba, lnr, rdm) code a trial the same way: the response time
# is the brms response, and the winning category and the number of
# accumulators racing for each category travel with it as vint() columns. The
# helpers below build and read that coding, so a race model supplies only its
# finishing-time distribution. `prefix` is the model's name, which names its
# coding columns (.<prefix>_cat, .<prefix>_n1, ...)

race_cat_column <- function(prefix) {
  paste0(".", prefix, "_cat")
}

race_count_columns <- function(prefix, n_cats) {
  paste0(".", prefix, "_n", seq_len(n_cats))
}

# the names a response category cannot take, besides those of the model's own
# parameters: brms declares `real Intercept` for the response's intercept,
# which collides with the `vector[N] Intercept` a category of that name would
# declare ("Identifier "Intercept" is already in use", stanc 2.40)
race_reserved_names <- function(shared_pars) {
  c("mu", shared_pars, "intercept")
}

# The categories of a custom version are the formula's parameters other than
# the model's shared ones, and each becomes a brms dpar and a Stan identifier
race_category_names <- function(formula, shared_pars) {
  cat_pars <- setdiff(names(formula), shared_pars)
  stopif(
    length(cat_pars) == 0,
    "Custom version requires at least one accumulator parameter in the formula."
  )

  bad_internal_names <- cat_pars[tolower(cat_pars) %in% race_reserved_names(shared_pars)]
  stopif(
    length(bad_internal_names) > 0,
    "Category names cannot use reserved internal parameter names: \\
    {collapse_comma(bad_internal_names)}."
  )

  bad_names <- intersect(cat_pars, .stan_reserved)
  stopif(
    length(bad_names) > 0,
    "Category names cannot be Stan reserved words: {collapse_comma(bad_names)}. \\
    Please rename the affected response categories."
  )

  bad_dpar_names <- cat_pars[grepl("[0-9]$", cat_pars)]
  stopif(
    length(bad_dpar_names) > 0,
    "Category names cannot end in a number because brms uses them as \\
    distributional parameters: {collapse_comma(bad_dpar_names)}. \\
    Please rename the affected response categories."
  )

  bad_underscore_names <- cat_pars[grepl("_", cat_pars, fixed = TRUE)]
  stopif(
    length(bad_underscore_names) > 0,
    "Category names cannot contain underscores because brms rejects them as \\
    distributional parameters: {collapse_comma(bad_underscore_names)}. \\
    Please rename the affected response categories."
  )
  cat_pars
}

race_check_data <- function(model, data) {
  rt_var <- model$resp_vars$rt
  response_var <- model$resp_vars$response

  stopif(
    not_in(rt_var, colnames(data)),
    "The RT variable '{rt_var}' is not present in the data."
  )

  n_na_rt <- sum(is.na(data[, rt_var]))
  stopif(
    n_na_rt > 0,
    "The RT variable '{rt_var}' contains {n_na_rt} NA values. \\
    Please remove or impute missing values before fitting the model."
  )

  stopif(
    !typeof(data[, rt_var]) %in% c("double", "integer"),
    "The RT variable '{rt_var}' needs to be of type double or integer."
  )
  stopif(
    any(data[, rt_var] <= 0),
    "Some reaction times are zero or negative, please check your data. \\
    The likelihood is zero for a response time that is not strictly \\
    greater than the non-decision time."
  )
  warnif(
    any(data[, rt_var] > 10),
    "Your data contains reaction times larger than 10 seconds.\n
    Either you have passed reaction times in milliseconds, then please \\
    recode them to seconds and rerun the model.\n
    Or you have very long RTs in your data in which case you might want \\
    to consider outlier filtering."
  )
  warnif(
    any(data[, rt_var] < 0.100),
    "Your data contains reaction times smaller than 0.100 seconds.\n
    It is likely that the model will not be able to sample with the \\
    current settings of the initial values.\n
    Either pass your own initial value function or consider filtering \\
    reaction times below 0.100 seconds."
  )

  stopif(
    not_in(response_var, colnames(data)),
    "The response variable '{response_var}' is not present in the data."
  )

  n_na_resp <- sum(is.na(data[, response_var]))
  stopif(
    n_na_resp > 0,
    "The response variable '{response_var}' contains {n_na_resp} NA values. \\
    Please remove or impute missing values before fitting the model."
  )
}

# The simple version races one correct accumulator against n_choices - 1 error
# accumulators, so the category code is 1 for a correct and 2 for any error
race_code_simple_response <- function(model, data, prefix) {
  response_var <- model$resp_vars$response
  n_alt <- model$other_vars$n_choices

  if (is.factor(data[, response_var]) || is.character(data[, response_var])) {
    labels <- as.character(data[, response_var])
    coded <- suppressWarnings(as.integer(labels))
    # as.integer() turns a label that is not a number into NA, and every check
    # below then compares against NA rather than refusing the label
    stopif(
      anyNA(coded),
      "The response variable '{response_var}' must be integer-coded 1:{n_alt} \\
      (1 = correct, 2:{n_alt} = errors) for version 'simple', but contains the \\
      non-numeric label(s) {collapse_comma(unique(labels[is.na(coded)]))}. Use \\
      version 'custom' for named response categories."
    )
    data[, response_var] <- coded
  }

  stopif(
    !is.numeric(data[, response_var]),
    "The response variable '{response_var}' must be numeric (integer-coded 1:{n_alt})."
  )

  resp_vals <- data[, response_var]
  stopif(
    any(resp_vals < 1) || any(resp_vals > n_alt) ||
      any(resp_vals != round(resp_vals)),
    "The response variable '{response_var}' must contain integers in 1:{n_alt}."
  )
  data[, response_var] <- as.integer(data[, response_var])

  counts <- race_count_columns(prefix, 2)
  data[[race_cat_column(prefix)]] <- ifelse(data[, response_var] == 1L, 1L, 2L)
  data[[counts[1]]] <- 1L
  data[[counts[2]]] <- n_alt - 1L
  data
}

race_code_custom_response <- function(model, data, prefix, shared_pars) {
  response_var <- model$resp_vars$response
  cat_names <- model$other_vars$resp_cats

  if (is.factor(data[, response_var])) {
    data[, response_var] <- as.character(data[, response_var])
  }

  stopif(
    !is.character(data[, response_var]),
    "For version 'custom', the response variable '{response_var}' must \\
    contain character labels matching the formula category names."
  )

  data_levels <- unique(data[, response_var])
  # a level named after one of the model's own parameters cannot have reached the
  # formula as a category, so it would otherwise be reported as unspecified
  bad_levels <- data_levels[tolower(data_levels) %in% race_reserved_names(shared_pars)]
  stopif(
    length(bad_levels) > 0,
    "Response levels cannot use reserved internal parameter names: \\
    {collapse_comma(bad_levels)}."
  )
  missing_in_formula <- setdiff(data_levels, cat_names)
  missing_in_data <- setdiff(cat_names, data_levels)
  stopif(
    length(missing_in_formula) > 0,
    "Response levels {collapse_comma(missing_in_formula)} in the data are \\
    not specified in the formula."
  )
  warnif(
    length(missing_in_data) > 0,
    "Formula categories {collapse_comma(missing_in_data)} are not present \\
    in the data."
  )

  data[[race_cat_column(prefix)]] <- unname(
    stats::setNames(seq_along(cat_names), cat_names)[data[, response_var]]
  )
  race_code_accumulators(model, data, prefix)
}

# `accumulators` is NULL (one per category), a named number per category, or a
# named column per category for counts that vary over trials. A column may give
# a category zero accumulators on a trial it sits out -- the likelihood and the
# simulators skip it -- but a category needs at least one on every trial it
# wins, because the likelihood opens with log(n_win)
race_code_accumulators <- function(model, data, prefix) {
  cat_names <- model$other_vars$resp_cats
  accumulators <- model$other_vars$accumulators
  cols <- race_count_columns(prefix, length(cat_names))

  if (is.null(accumulators)) {
    data[cols] <- 1L
  } else {
    stopif(
      !is.numeric(accumulators) && !is.character(accumulators),
      "accumulators must be NULL, a named numeric vector of positive integers, \\
      or a named character vector of column names."
    )
    stopif(
      is.null(names(accumulators)) || any(names(accumulators) == "") ||
        anyDuplicated(names(accumulators)) > 0,
      "accumulators must be a uniquely named vector with one entry for each \\
      formula category: {collapse_comma(cat_names)}"
    )
    stopif(
      !setequal(names(accumulators), cat_names),
      "accumulators must have exactly the formula categories: \\
      {collapse_comma(cat_names)}"
    )
    counts <- if (is.numeric(accumulators)) {
      race_constant_counts(accumulators[cat_names], nrow(data))
    } else {
      race_column_counts(data, accumulators[cat_names])
    }
    data[cols] <- counts
  }

  winner <- data[[race_cat_column(prefix)]]
  won <- as.matrix(data[cols])[cbind(seq_len(nrow(data)), winner)]
  empty_winners <- unique(cat_names[winner[won < 1L]])
  stopif(
    length(empty_winners) > 0,
    "Category {collapse_comma(empty_winners)} wins on trials where accumulators \\
    gives it no accumulator. A category that can be chosen needs at least one \\
    accumulator on every trial where it wins."
  )
  data
}

race_constant_counts <- function(counts, n_rows) {
  invalid <- counts[!is.finite(counts) | counts < 1 | counts != round(counts)]
  stopif(
    length(invalid) > 0,
    "accumulators must contain a positive integer for each formula category. \\
    Invalid value(s): {collapse_comma(glue('{names(invalid)} = {invalid}'))}"
  )
  lapply(counts, function(n) rep(as.integer(n), n_rows))
}

race_column_counts <- function(data, columns) {
  missing_cols <- setdiff(columns, colnames(data))
  stopif(
    length(missing_cols) > 0,
    "accumulators columns {collapse_comma(missing_cols)} not found in the data."
  )
  lapply(columns, function(col_name) {
    col_vals <- data[, col_name]
    stopif(
      !is.numeric(col_vals),
      "accumulators column '{col_name}' must be numeric."
    )
    stopif(
      anyNA(col_vals) || any(!is.finite(col_vals)),
      "accumulators column '{col_name}' contains NA or non-finite values."
    )
    # fractional or negative counts reach Stan as a NaN target
    stopif(
      any(col_vals < 0 | col_vals != round(col_vals)),
      "accumulators column '{col_name}' must contain integers >= 0."
    )
    as.integer(col_vals)
  })
}

race_bf <- function(model, prefix, n_cats) {
  vint_args <- paste(c(race_cat_column(prefix), race_count_columns(prefix, n_cats)), collapse = ", ")
  brms::bf(glue("{model$resp_vars$rt} | vint({vint_args}) ~ 1"))
}

# The vint() columns hold the winning category and the per-category accumulator
# counts, one value per observation, in the order the Stan signature expects.
# reduce_sum slices Y inside partial_log_lik but passes a custom family's `vars`
# through whole, so a bare "vint1" would pair each slice's response times with
# the top of the data. start/end only exist inside partial_log_lik, so the
# columns are sliced only where brms really threads (see brms_slices_likelihood)
race_family_vars <- function(n_cats) {
  vars <- paste0("vint", seq_len(n_cats + 1))
  if (brms_slices_likelihood()) paste0(vars, "[start:end]") else vars
}

# the accumulator counts of observation i, one per category
race_counts <- function(prep, i, n_cats) {
  vapply(
    seq_len(n_cats),
    function(j) prep$data[[paste0("vint", j + 1)]][i],
    integer(1)
  )
}

# brms keeps only the vint() columns, so the response comes back from the
# winning category and a column form of `accumulators` from the counts. The
# simple version's code is 1 for correct and 2 for any error: which error was
# given is gone, and the likelihood never used it
race_revert_check_data <- function(model, data, prefix) {
  response <- model$resp_vars$response
  cat_names <- model$other_vars$resp_cats
  accumulators <- model$other_vars$accumulators
  codes <- data[[race_cat_column(prefix)]]
  if (not_in(response, colnames(data))) {
    data[[response]] <- if (is.null(cat_names)) codes else cat_names[codes]
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), response)
  }
  if (is.character(accumulators)) {
    for (i in seq_along(cat_names)) {
      col <- accumulators[[cat_names[i]]]
      if (not_in(col, colnames(data))) {
        data[[col]] <- data[[race_count_columns(prefix, i)[i]]]
        attr(data, "rebuilt") <- c(attr(data, "rebuilt"), col)
      }
    }
  }
  data[grep(glue("^\\.{prefix}_(cat|n[0-9]+)$"), colnames(data))] <- NULL
  data
}

# brms::posterior_predict() returns one matrix, so the winning category rides
# along as a second observable. vint1 is the category code that check_data()
# wrote, in the order of the family's accumulator dpars.
race_pp_observables <- function(model) {
  cats <- model$other_vars$resp_cats %||% c("correct", "error")
  coding <- collapse_comma(glue("{seq_along(cats)} = {cats}"))
  list(
    observed = c(rt = "Y", response = "vint1"),
    checks = list(
      rt = .pp_observable(function(d) d$rt, label = "Response time"),
      response = .pp_observable(
        function(d) d$response,
        label = glue("Response category ({coding})"),
        type = "bars"
      )
    )
  )
}

# E[T] = int_0^inf S(t) dt for a non-negative T, evaluated deterministically so
# that two calls on the same draws return the same number. The substitution
# u = log(t) turns a race survivor, which decays over several orders of
# magnitude in t, into an integrand S(e^u) e^u that a uniform grid in u
# resolves. `log_surv` takes the grid times and returns a draws x grid matrix of
# log survivors. Below `t_lo` the survivor is 1 to the caller's chosen
# tolerance, so that stretch contributes t_lo exactly. On lognormal and Wald
# races the error is below 4e-7 relative to stats::integrate(), against up to
# 2e-5 for a trapezoid of 20 points per decade
race_expected_time <- function(log_surv, t_lo, t_hi, n_grid = 1024L) {
  u <- seq(log(t_lo), log(t_hi), length.out = n_grid)
  integrand <- exp(sweep(log_surv(exp(u)), 2, u, `+`))
  du <- u[2] - u[1]
  edge <- (integrand[, 1] + integrand[, n_grid]) / 2
  t_lo + du * (rowSums(integrand) - edge)
}

# The integration range of race_expected_time() for a race whose survivor has
# no closed-form quantile, searched decade by decade from 1 s: t_lo is the
# largest decade where every draw's survivor is still within `tol` of 1, t_hi
# the smallest where every draw's S(t) t is below `tol`. S(t) t, not S(t),
# because it bounds the neglected tail of a survivor that decays like a power
# of t (the normal-drift LBA race decays like 1/t^2). A survivor that decays
# more slowly than 1/t has no finite mean; t_max stops the search for it. A
# draw whose survivor is not finite at a decade counts as not yet inside `tol`
race_time_range <- function(log_surv, tol = 1e-12, t_max = 1e7) {
  t_hi <- 1
  while (t_hi < t_max && !isTRUE(all(exp(log_surv(t_hi)) * t_hi <= tol))) {
    t_hi <- t_hi * 10
  }
  t_lo <- 1
  while (t_lo > t_hi * 1e-12 && !isTRUE(all(-expm1(log_surv(t_lo)) <= tol))) {
    t_lo <- t_lo / 10
  }
  c(t_lo, t_hi)
}
