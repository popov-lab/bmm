############################################################################# !
# SDT SHARED PIPELINE HELPERS                                            ####
# Used by model_sdt_*.R files for check_data, bmf2bf, configure_model    #
############################################################################# !


############################################################################# !
# DISTRIBUTION ID & NAME MAPPING                                          ####
############################################################################# !

# registry position defines the Stan dist_type integer (see .sdt_dists)
.sdt_dist_id <- function(dist) {
  match(dist, names(.sdt_dists))
}

.sdt_dist_names <- names(.sdt_dists)


############################################################################# !
# DATA VALIDATION HELPERS                                                 ####
############################################################################# !

.validate_sdt_counts <- function(data, resp_var, n_trials_var) {
  required <- c(resp_var, n_trials_var)
  missing <- setdiff(required, colnames(data))
  stopif(length(missing) > 0,
    "Variables {collapse_comma(missing)} missing in the data")

  resp_vals <- data[[resp_var]]
  stopif(any(resp_vals < 0, na.rm = TRUE),
    "Response variable '{resp_var}' must contain non-negative counts")
  warnif(any(resp_vals != round(resp_vals), na.rm = TRUE),
    "Response variable '{resp_var}' should contain integer counts")

  trial_vals <- data[[n_trials_var]]
  stopif(any(trial_vals <= 0, na.rm = TRUE),
    "Variable '{n_trials_var}' must contain positive values")
  warnif(any(trial_vals != round(trial_vals), na.rm = TRUE),
    "Variable '{n_trials_var}' should contain integer counts")

  stopif(any(resp_vals > trial_vals, na.rm = TRUE),
    "Response counts in '{resp_var}' must not exceed '{n_trials_var}'")
}

# Multi-column variant for the multinomial models (ranking, rating), whose
# response is one count column per category.
.validate_sdt_count_cols <- function(data, resp_cols) {
  missing <- setdiff(resp_cols, colnames(data))
  stopif(length(missing) > 0,
    "Response columns {collapse_comma(missing)} missing in the data")

  for (col in resp_cols) {
    vals <- data[[col]]
    stopif(any(vals < 0, na.rm = TRUE),
      "Response column '{col}' must contain non-negative counts")
    warnif(any(vals != round(vals), na.rm = TRUE),
      "Response column '{col}' should contain integer counts")
  }
}

# Accept the set size either as a constant or as a column name, so models can
# fit trials with different set sizes jointly. Returns one set size per row.
.sdt_resolve_set_size <- function(m, data) {
  if (is.character(m)) {
    stopif(!m %in% colnames(data),
      "Set-size column '{m}' missing in the data")
    stopif(!is.numeric(data[[m]]),
      "Set-size column '{m}' must be numeric")
    warnif(any(data[[m]] != round(data[[m]]), na.rm = TRUE),
      "Set-size column '{m}' should contain integer values")
    out <- as.integer(data[[m]])
  } else {
    out <- rep.int(as.integer(m), nrow(data))
  }
  stopif(any(out < 2, na.rm = TRUE), "Set size must be an integer >= 2")
  out
}
