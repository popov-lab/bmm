devtools::load_all()

cmdstan_dir <- file.path("tests", "internal", "cmdstanr_output")
if (!dir.exists(cmdstan_dir)) {
  dir.create(cmdstan_dir, recursive = TRUE)
}
options(cmdstanr_write_stan_file_dir = normalizePath(cmdstan_dir))

rdm_benchmark_metrics <- function(fit, elapsed_sec) {
  draws <- posterior::as_draws_matrix(fit)
  draw_names <- colnames(draws)
  is_sampler <- grepl("__$", draw_names)
  is_constant <- apply(draws, 2, function(x) stats::sd(x) == 0)
  keep <- !is_sampler & !is_constant

  bulk_ess <- posterior::ess_bulk(draws[, keep, drop = FALSE])
  tail_ess <- posterior::ess_tail(draws[, keep, drop = FALSE])
  nuts <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)

  divergences <- NA_integer_
  max_treedepth <- NA_integer_
  if (!is.null(nuts)) {
    divergences <- sum(nuts$Parameter == "divergent__" & nuts$Value > 0)
    treedepth_vals <- nuts$Value[nuts$Parameter == "treedepth__"]
    if (length(treedepth_vals) > 0) {
      max_treedepth <- max(treedepth_vals)
    }
  }

  data.frame(
    parameter = names(bulk_ess),
    bulk_ess = unname(bulk_ess),
    tail_ess = unname(tail_ess),
    ess_per_sec = unname(bulk_ess) / elapsed_sec,
    elapsed_sec = elapsed_sec,
    divergences = divergences,
    max_treedepth = max_treedepth,
    row.names = NULL
  )
}

simulate_rdm_custom <- function(n, drift, counts, gap, ndt, s = 1, sp = 0) {
  resp_names <- names(counts)
  full_drift <- unlist(Map(rep, drift, counts), use.names = FALSE)
  raw <- rrdm(n = n, drift = full_drift, gap = gap, ndt = ndt, s = s, sp = sp)
  idx <- rep(seq_along(counts), counts)
  raw$response <- resp_names[idx[raw$response]]
  raw
}

run_rdm_benchmark <- function(name, data, formula, model,
                              threads = NULL, sort_data = TRUE,
                              iter = 1000, warmup = 500,
                              chains = 2, cores = 2) {
  cat("\n", strrep("-", 80), "\n", sep = "")
  cat("Benchmark:", name, "\n")
  cat("sort_data =", sort_data, "\n")
  cat("threads =", if (is.null(threads)) "none" else deparse(threads), "\n")
  cat(strrep("-", 80), "\n", sep = "")

  fit <- NULL
  elapsed <- system.time({
    fit <- bmm(
      formula = formula,
      data = data,
      model = model,
      backend = "cmdstanr",
      chains = chains,
      cores = cores,
      iter = iter,
      warmup = warmup,
      refresh = 0,
      silent = 2,
      sort_data = sort_data,
      threads = threads
    )
  })[["elapsed"]]

  metrics <- rdm_benchmark_metrics(fit, elapsed)
  print(metrics)

  invisible(list(fit = fit, metrics = metrics))
}

simple_fixed_data <- rrdm(
  n = 500,
  drift = c(3.0, 1.2),
  gap = 1.0,
  ndt = 0.2
)

simple_free_data <- rrdm(
  n = 500,
  drift = c(3.1, 1.4, 1.4, 1.4),
  gap = 1.0,
  ndt = 0.22,
  s = 0.9,
  sp = 0.04
)

custom_free_data <- simulate_rdm_custom(
  n = 500,
  drift = c(corr = 2.8, lure = 1.9, npl = 1.2),
  counts = c(corr = 1, lure = 2, npl = 2),
  gap = 1.05,
  ndt = 0.24,
  s = 0.95,
  sp = 0.03
)

benchmark_cases <- list(
  simple_fixed = list(
    data = simple_fixed_data,
    formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1),
    model = rdm(rt = "rt", response = "response", n_alternatives = 2)
  ),
  simple_free_sp = list(
    data = simple_free_data,
    formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, s ~ 1, sp ~ 1),
    model = rdm(rt = "rt", response = "response", n_alternatives = 4)
  ),
  custom_free_sp = list(
    data = custom_free_data,
    formula = bmf(corr ~ 1, lure ~ 1, npl ~ 1, gap ~ 1, ndt ~ 1, s ~ 1, sp ~ 1),
    model = rdm(
      rt = "rt",
      response = "response",
      version = "custom",
      num_alternatives = c(corr = 1, lure = 2, npl = 2)
    )
  )
)

benchmark_results <- list()
for (name in names(benchmark_cases)) {
  case <- benchmark_cases[[name]]
  benchmark_results[[paste0(name, "_serial")]] <- run_rdm_benchmark(
    name = paste(name, "(serial)"),
    data = case$data,
    formula = case$formula,
    model = case$model,
    threads = NULL,
    sort_data = TRUE
  )

  benchmark_results[[paste0(name, "_threaded")]] <- run_rdm_benchmark(
    name = paste(name, "(threaded, sort_data = FALSE)"),
    data = case$data,
    formula = case$formula,
    model = case$model,
    threads = brms::threading(2),
    sort_data = FALSE
  )
}

benchmark_summary <- do.call(
  rbind,
  Map(function(name, result) {
    metrics <- result$metrics
    data.frame(
      run = name,
      elapsed_sec = unique(metrics$elapsed_sec),
      mean_ess_per_sec = mean(metrics$ess_per_sec, na.rm = TRUE),
      divergences = unique(metrics$divergences),
      max_treedepth = unique(metrics$max_treedepth),
      row.names = NULL
    )
  }, names(benchmark_results), benchmark_results)
)

cat("\nBenchmark summary\n")
print(benchmark_summary)

invisible(list(details = benchmark_results, summary = benchmark_summary))
