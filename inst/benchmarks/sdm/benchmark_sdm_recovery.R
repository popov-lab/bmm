source("inst/benchmarks/sdm/simulate_sdm_data.R")

run_sdm_recovery_check <- function(iter = 3000L,
                                   warmup = 1000L,
                                   chains = 4L,
                                   cores = chains,
                                   threads_per_chain = 1L,
                                   seed = 123,
                                   output_dir = "inst/benchmarks/sdm/results") {
  check_positive_int(iter, "iter")
  check_positive_int(warmup, "warmup")
  check_positive_int(chains, "chains")
  check_positive_int(cores, "cores")
  check_positive_int(threads_per_chain, "threads_per_chain")
  if (cores > chains) {
    stop("cores must be <= chains for local recovery runs.", call. = FALSE)
  }
  estimated_cpu_threads <- cores * threads_per_chain
  if (estimated_cpu_threads > 8L) {
    stop("cores * threads_per_chain must be <= 8 for local recovery runs.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  sim <- simulate_sdm_benchmark_data(config = c(sdm_benchmark_configs()$recovery, seed = seed))
  formula <- bmf(
    c ~ 0 + condition + (0 + condition | subject),
    kappa ~ 0 + condition + (0 + condition | subject)
  )
  threads <- if (threads_per_chain > 1L) {
    brms::threading(threads_per_chain)
  }

  fit_time <- system.time({
    fit <- bmm(
      formula = formula,
      data = sim$data,
      model = sdm(resp_error = "y"),
      backend = "cmdstanr",
      sort_data = TRUE,
      silent = 2,
      chains = chains,
      cores = cores,
      iter = iter,
      warmup = warmup,
      seed = seed,
      refresh = 0,
      threads = threads
    )
  })

  recovery <- fixed_effect_recovery(fit, sim$metadata$n_conditions)
  diagnostics <- collect_recovery_diagnostics(fit, unname(fit_time[["elapsed"]]))
  git <- git_metadata()
  recovery$benchmark_name <- "sdm_recovery"
  recovery$git_branch <- git$branch
  recovery$git_commit <- git$commit
  recovery$total_model_time_sec <- unname(fit_time[["elapsed"]])
  recovery$divergent_transitions <- diagnostics$divergent_transitions
  recovery$total_leapfrog_steps <- diagnostics$total_leapfrog_steps
  recovery$min_bulk_ess <- diagnostics$min_bulk_ess
  recovery$min_bulk_ess_per_total_sec <- diagnostics$min_bulk_ess_per_total_sec
  recovery$chains <- chains
  recovery$cores <- cores
  recovery$threads_per_chain <- threads_per_chain
  recovery$estimated_cpu_threads <- estimated_cpu_threads
  recovery$iter <- iter
  recovery$warmup <- warmup
  recovery$seed <- seed

  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  base <- file.path(output_dir, paste0("sdm_recovery_", git$commit, "_", stamp))
  utils::write.csv(recovery, paste0(base, ".csv"), row.names = FALSE)
  saveRDS(list(recovery = recovery, metadata = sim$metadata), paste0(base, ".rds"))
  recovery
}

fixed_effect_recovery <- function(fit, n_conditions) {
  fe <- as.data.frame(brms::fixef(fit))
  fe$term <- rownames(fe)
  true <- data.frame(
    parameter = rep(c("c", "kappa"), each = n_conditions),
    condition = rep(paste0("condition", seq_len(n_conditions)), 2L),
    true = c(
      seq(2.5, 5.5, length.out = n_conditions),
      seq(2, 6, length.out = n_conditions)
    )
  )
  true$term <- paste0(true$parameter, "_condition", true$condition)
  out <- merge(true, fe, by = "term", all.x = TRUE, sort = FALSE)
  out$estimate_native <- exp(out$Estimate)
  out$q2.5_native <- exp(out$Q2.5)
  out$q97.5_native <- exp(out$Q97.5)
  out$abs_error <- abs(out$estimate_native - out$true)
  out$true_in_95_ci <- out$true >= out$q2.5_native & out$true <= out$q97.5_native
  out
}

collect_recovery_diagnostics <- function(fit, total_time_sec) {
  nuts <- try(brms::nuts_params(fit), silent = TRUE)
  divergent <- NA_integer_
  leapfrog <- NA_real_
  if (!inherits(nuts, "try-error")) {
    divergent <- sum(nuts$Parameter == "divergent__" & nuts$Value == 1, na.rm = TRUE)
    leapfrog <- sum(nuts$Value[nuts$Parameter == "n_leapfrog__"], na.rm = TRUE)
  }

  draws <- try(posterior::summarise_draws(fit, "ess_bulk"), silent = TRUE)
  min_ess_bulk <- NA_real_
  if (!inherits(draws, "try-error") && "ess_bulk" %in% names(draws)) {
    min_ess_bulk <- min(draws$ess_bulk, na.rm = TRUE)
  }

  data.frame(
    divergent_transitions = divergent,
    total_leapfrog_steps = leapfrog,
    min_bulk_ess = min_ess_bulk,
    min_bulk_ess_per_total_sec = min_ess_bulk / max(total_time_sec, .Machine$double.eps)
  )
}

git_metadata <- function() {
  list(
    branch = system2("git", c("branch", "--show-current"), stdout = TRUE),
    commit = system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
  )
}

parse_sdm_recovery_args <- function(args) {
  get_arg <- function(name, default) {
    hit <- grep(paste0("^", name, "="), args, value = TRUE)
    if (length(hit) == 0L) {
      return(default)
    }
    sub(paste0("^", name, "="), "", hit[[1]])
  }
  list(
    iter = as.integer(get_arg("--iter", "3000")),
    warmup = as.integer(get_arg("--warmup", "1000")),
    chains = as.integer(get_arg("--chains", "4")),
    cores = as.integer(get_arg("--cores", "4")),
    threads_per_chain = as.integer(get_arg("--threads-per-chain", "1")),
    seed = as.integer(get_arg("--seed", "123"))
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- parse_sdm_recovery_args(commandArgs(trailingOnly = TRUE))
  print(do.call(run_sdm_recovery_check, args))
}
