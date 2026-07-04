sdm_benchmark_configs <- function() {
  list(
    smoke = list(
      n_subjects = 5L,
      n_conditions = 2L,
      n_tasks = 1L,
      trials_per_subject_condition_task = 10L
    ),
    subjects30_conditions3 = list(
      n_subjects = 30L,
      n_conditions = 3L,
      n_tasks = 1L,
      trials_per_subject_condition_task = 100L
    ),
    subjects30_conditions5 = list(
      n_subjects = 30L,
      n_conditions = 5L,
      n_tasks = 1L,
      trials_per_subject_condition_task = 100L
    ),
    subjects60_conditions5 = list(
      n_subjects = 60L,
      n_conditions = 5L,
      n_tasks = 1L,
      trials_per_subject_condition_task = 100L
    ),
    recovery = list(
      n_subjects = 30L,
      n_conditions = 3L,
      n_tasks = 1L,
      trials_per_subject_condition_task = 100L,
      recovery = TRUE
    )
  )
}

simulate_sdm_benchmark_data <- function(config = NULL,
                                        n_subjects = 5L,
                                        n_conditions = 2L,
                                        n_tasks = 1L,
                                        trials_per_subject_condition_task = 10L,
                                        seed = 123,
                                        condition_correlation = 0.5,
                                        c_means = NULL,
                                        kappa_means = NULL,
                                        subject_log_c_sd = 0.25,
                                        subject_log_kappa_sd = 0.25,
                                        parametrization = "sqrtexp",
                                        recovery = FALSE) {
  if (!is.null(config)) {
    return(do.call(simulate_sdm_benchmark_data, config))
  }

  check_positive_int(n_subjects, "n_subjects")
  check_positive_int(n_conditions, "n_conditions")
  check_positive_int(n_tasks, "n_tasks")
  check_positive_int(trials_per_subject_condition_task, "trials_per_subject_condition_task")
  if (!is.numeric(condition_correlation) || length(condition_correlation) != 1L ||
      condition_correlation <= -1 || condition_correlation >= 1) {
    stop("condition_correlation must be a single number between -1 and 1.", call. = FALSE)
  }

  with_sdm_seed(seed, {
    condition <- factor(paste0("condition", seq_len(n_conditions)))
    subject <- factor(paste0("subject", seq_len(n_subjects)))
    task <- factor(paste0("task", seq_len(n_tasks)))

    log_c_means <- default_log_means(c_means, n_conditions, c(2.5, 5.5), "c_means")
    log_kappa_means <- default_log_means(kappa_means, n_conditions, c(2, 6), "kappa_means")

    log_c <- draw_subject_condition_values(
      n_subjects, log_c_means, subject_log_c_sd, condition_correlation
    )
    log_kappa <- draw_subject_condition_values(
      n_subjects, log_kappa_means, subject_log_kappa_sd, condition_correlation
    )

    true_parameters <- expand.grid(
      subject = subject,
      condition = condition,
      task = task,
      KEEP.OUT.ATTRS = FALSE
    )
    cond_idx <- as.integer(true_parameters$condition)
    subj_idx <- as.integer(true_parameters$subject)
    true_parameters$mu <- 0
    true_parameters$c <- exp(log_c[cbind(subj_idx, cond_idx)])
    true_parameters$kappa <- exp(log_kappa[cbind(subj_idx, cond_idx)])

    data <- true_parameters[rep(seq_len(nrow(true_parameters)),
      each = trials_per_subject_condition_task
    ), c("subject", "condition", "task")]
    data$trial <- ave(
      seq_len(nrow(data)),
      data$subject, data$condition, data$task,
      FUN = seq_along
    )

    param_idx <- match(
      interaction(data$subject, data$condition, data$task),
      interaction(true_parameters$subject, true_parameters$condition, true_parameters$task)
    )
    data$mu <- true_parameters$mu[param_idx]
    data$c <- true_parameters$c[param_idx]
    data$kappa <- true_parameters$kappa[param_idx]
    data$y <- simulate_grouped_sdm(data, parametrization)
    data <- data[c("subject", "condition", "task", "trial", "y")]
    rownames(data) <- NULL

    list(
      data = data,
      true_parameters = true_parameters,
      metadata = list(
        seed = seed,
        n_subjects = n_subjects,
        n_conditions = n_conditions,
        n_tasks = n_tasks,
        trials_per_subject_condition_task = trials_per_subject_condition_task,
        n_trials = nrow(data),
        condition_correlation = condition_correlation,
        subject_log_c_sd = subject_log_c_sd,
        subject_log_kappa_sd = subject_log_kappa_sd,
        parametrization = parametrization,
        recovery = recovery,
        response_error = "y",
        response_error_sampler = "bmm::rsdm rejection sampling"
      )
    )
  })
}

check_positive_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x %% 1 != 0) {
    stop(name, " must be a positive whole number.", call. = FALSE)
  }
}

default_log_means <- function(x, n, range, name) {
  if (is.null(x)) {
    return(log(seq(range[1], range[2], length.out = n)))
  }
  if (!is.numeric(x) || length(x) != n || any(x <= 0)) {
    stop(name, " must contain one positive native-scale value per condition.", call. = FALSE)
  }
  log(x)
}

draw_subject_condition_values <- function(n_subjects, means, subject_sd, condition_correlation) {
  n_conditions <- length(means)
  cor_mat <- matrix(condition_correlation, n_conditions, n_conditions)
  diag(cor_mat) <- 1
  cov_mat <- cor_mat * subject_sd^2
  z <- matrix(stats::rnorm(n_subjects * n_conditions), n_subjects, n_conditions)
  sweep(z %*% chol(cov_mat), 2, means, "+")
}

simulate_grouped_sdm <- function(data, parametrization) {
  split_idx <- split(seq_len(nrow(data)), interaction(data$subject, data$condition, data$task))
  y <- numeric(nrow(data))
  rsdm_fun <- get_rsdm()
  for (idx in split_idx) {
    y[idx] <- rsdm_fun(
      n = length(idx),
      mu = data$mu[idx[1]],
      c = data$c[idx[1]],
      kappa = data$kappa[idx[1]],
      parametrization = parametrization
    )
  }
  y
}

get_rsdm <- function() {
  if (exists("rsdm", mode = "function")) {
    return(get("rsdm", mode = "function"))
  }
  if (requireNamespace("bmm", quietly = TRUE)) {
    return(getExportedValue("bmm", "rsdm"))
  }
  stop("Load bmm first, for example with devtools::load_all().", call. = FALSE)
}

with_sdm_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

if (identical(environment(), globalenv()) && !interactive()) {
  out <- simulate_sdm_benchmark_data(config = sdm_benchmark_configs()$smoke)
  stopifnot(nrow(out$data) == 5L * 2L * 10L)
  stopifnot(all(out$data$y >= -pi & out$data$y <= pi))
  stopifnot(nrow(out$true_parameters) == 5L * 2L)
}
