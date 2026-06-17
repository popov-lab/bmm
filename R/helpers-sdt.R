############################################################################# !
# SDT SHARED PIPELINE HELPERS                                            ####
# Shared across model_sdt_*.R files (binary, mafc, ranking, rating) for  #
# check_data and configure_model. Rating-only helpers live in            #
# model_sdt_rating.R.                                                     #
############################################################################# !


############################################################################# !
# DISTRIBUTION ID & NAME MAPPING                                          ####
############################################################################# !

.sdt_dist_id <- function(dist) {
  .SDT_DISTS[[dist]]$id
}

# ordered by id so vint2[i] in log_lik/posterior_predict indexes the right name
.sdt_dist_names <- names(.SDT_DISTS)[order(vapply(.SDT_DISTS, `[[`, 0L, "id"))]


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

  stopif(any(resp_vals > trial_vals, na.rm = TRUE),
    "Response counts in '{resp_var}' must not exceed '{n_trials_var}'")
}

# Accept the set size either as a constant or as a column name, so models can
# fit trials with different set sizes jointly. Returns one set size per row.
.sdt_resolve_set_size <- function(m, data) {
  if (is.character(m)) {
    stopif(!m %in% colnames(data),
      "Set-size column '{m}' missing in the data")
    out <- as.integer(data[[m]])
  } else {
    out <- rep.int(as.integer(m), nrow(data))
  }
  stopif(any(out < 2, na.rm = TRUE), "Set size must be an integer >= 2")
  out
}
