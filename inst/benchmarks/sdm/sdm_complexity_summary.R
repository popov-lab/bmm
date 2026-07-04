sdm_complexity_summary <- function(data,
                                   subject = "subject",
                                   condition = "condition",
                                   task = "task",
                                   run_vars = c(subject, condition, task),
                                   n_subject_level_effects = NA_integer_,
                                   n_pairwise_correlations = NA_integer_) {
  stopifnot(is.data.frame(data))
  required <- unique(c(subject, condition, task, run_vars))
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  runs <- unique(data[run_vars])
  trials_per_run <- as.integer(table(interaction(data[run_vars], drop = TRUE)))

  data.frame(
    n_trials = nrow(data),
    n_subjects = n_unique(data[[subject]]),
    n_conditions = n_unique(data[[condition]]),
    n_tasks = n_unique(data[[task]]),
    run_vars = paste(run_vars, collapse = ","),
    n_sdm_predictor_runs = nrow(runs),
    compression_ratio = nrow(data) / nrow(runs),
    trials_per_run_min = min(trials_per_run),
    trials_per_run_median = stats::median(trials_per_run),
    trials_per_run_max = max(trials_per_run),
    n_subject_level_effects = n_subject_level_effects,
    n_pairwise_correlations = n_pairwise_correlations
  )
}

n_unique <- function(x) {
  length(unique(x))
}

if (identical(environment(), globalenv()) && !interactive()) {
  source("inst/benchmarks/sdm/simulate_sdm_data.R")
  sim <- simulate_sdm_benchmark_data(config = sdm_benchmark_configs()$smoke)
  out <- sdm_complexity_summary(sim$data)
  stopifnot(out$n_trials == 100L)
  stopifnot(out$n_sdm_predictor_runs == 10L)
  stopifnot(out$compression_ratio == 10)
}
