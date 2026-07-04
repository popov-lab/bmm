source("inst/benchmarks/sdm/simulate_sdm_data.R")
source("inst/benchmarks/sdm/sdm_complexity_summary.R")

run_sdm_baseline_benchmark <- function(config_name = "smoke",
                                       iter = 100L,
                                       warmup = 50L,
                                       chains = 2L,
                                       seed = 123,
                                       output_dir = "inst/benchmarks/sdm/results") {
  configs <- sdm_benchmark_configs()
  if (!config_name %in% names(configs)) {
    stop("Unknown config: ", config_name, call. = FALSE)
  }
  check_positive_int(iter, "iter")
  check_positive_int(warmup, "warmup")
  check_positive_int(chains, "chains")
  if (chains > 6L) {
    stop("chains must be <= 6 for local benchmark runs.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sim <- simulate_sdm_benchmark_data(config = c(configs[[config_name]], seed = seed))
  formula <- bmf(
    c ~ 0 + condition + (0 + condition | subject),
    kappa ~ 0 + condition + (0 + condition | subject)
  )
  model <- sdm(resp_error = "y")

  translation_time <- system.time({
    stan_code <- stancode(formula, data = sim$data, model = model)
  })
  stan_data_time <- system.time({
    stan_data <- standata(formula, data = sim$data, model = model)
  })
  stan_code_hash <- hash_text(stan_code)

  started_at <- Sys.time()
  fit_time <- system.time({
    fit <- bmm(
      formula = formula,
      data = sim$data,
      model = model,
      backend = "cmdstanr",
      sort_data = TRUE,
      silent = 2,
      chains = chains,
      cores = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      refresh = 0
    )
  })
  finished_at <- Sys.time()

  diagnostics <- collect_sdm_diagnostics(fit, unname(fit_time[["elapsed"]]))
  complexity <- sdm_complexity_summary(
    sim$data,
    n_subject_level_effects = 2L * sim$metadata$n_conditions,
    n_pairwise_correlations = choose(2L * sim$metadata$n_conditions, 2L)
  )
  git <- git_metadata()

  result <- cbind(
    data.frame(
      benchmark_name = paste0("sdm_baseline_", config_name),
      config_name = config_name,
      git_branch = git$branch,
      git_commit = git$commit,
      git_tag = git$tag,
      package_version = as.character(utils::packageVersion("bmm")),
      brms_version = as.character(utils::packageVersion("brms")),
      cmdstanr_version = as.character(utils::packageVersion("cmdstanr")),
      cmdstan_version = cmdstanr::cmdstan_version(error_on_NA = FALSE),
      hardware = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
      detected_cores = parallel::detectCores(logical = TRUE),
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
      finished_at = format(finished_at, "%Y-%m-%d %H:%M:%S %Z"),
      stan_translation_time_sec = unname(translation_time[["elapsed"]]),
      standata_time_sec = unname(stan_data_time[["elapsed"]]),
      cpp_compilation_time_sec = NA_real_,
      total_model_time_sec = unname(fit_time[["elapsed"]]),
      stan_code_hash = stan_code_hash
    ),
    complexity,
    diagnostics
  )

  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  base <- file.path(output_dir, paste0(result$benchmark_name, "_", git$commit, "_", stamp))
  utils::write.csv(result, paste0(base, ".csv"), row.names = FALSE)
  saveRDS(list(result = result, metadata = sim$metadata), paste0(base, ".rds"))

  result
}

collect_sdm_diagnostics <- function(fit, total_time_sec) {
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
    commit = system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
    tag = paste(system2("git", c("tag", "--points-at", "HEAD"), stdout = TRUE), collapse = ";")
  )
}

hash_text <- function(x) {
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(x, tf)
  tools::md5sum(tf)[[1]]
}

parse_sdm_benchmark_args <- function(args) {
  get_arg <- function(name, default) {
    hit <- grep(paste0("^", name, "="), args, value = TRUE)
    if (length(hit) == 0L) {
      return(default)
    }
    sub(paste0("^", name, "="), "", hit[[1]])
  }
  list(
    config_name = get_arg("--config", "smoke"),
    iter = as.integer(get_arg("--iter", "100")),
    warmup = as.integer(get_arg("--warmup", "50")),
    chains = as.integer(get_arg("--chains", "2")),
    seed = as.integer(get_arg("--seed", "123"))
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- parse_sdm_benchmark_args(commandArgs(trailingOnly = TRUE))
  print(do.call(run_sdm_baseline_benchmark, args))
}
