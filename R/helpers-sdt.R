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
# FORMULA CHECKS                                                          ####
############################################################################# !

# The likelihood consumes the stimulus variable to set each row's role, so a
# stimulus term in a parameter formula is confounded with that parameter's own
# intercept: moving b_d_stimulus by +0.5 moves d's intercept by -0.25 and
# criterion by +0.125 at identical likelihood. Profiled flat to 6e-11 on a
# 5-condition design that is otherwise richly identified, for d, criterion and
# sdratio alike. Formats without a stimulus other-var (sdt_mafc, sdt_ranking)
# skip this.
#' @export
check_formula.sdt <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  if (is.null(stim_var)) {
    return(NextMethod("check_formula"))
  }
  # rhs_vars() returns random-effect grouping variables too, so (1 | stimulus)
  # is caught here as well -- a random intercept per stimulus level is the same
  # confound
  uses_stim <- vapply(rhs_vars(formula, collapse = FALSE),
                      function(x) stim_var %in% x, logical(1))
  stopif(
    any(uses_stim),
    "The formula for parameter(s) {collapse_comma(names(uses_stim)[uses_stim])} \\
    uses the stimulus variable '{stim_var}'. {model$name} already uses \\
    '{stim_var}' in the likelihood to tell noise rows from signal rows, so a \\
    term in '{stim_var}' -- as a predictor or as a grouping factor -- is \\
    confounded with that parameter's own intercept and the model is not \\
    identified. Drop '{stim_var}' from the formula."
  )
  NextMethod("check_formula")
}


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
