# Manual SDT benchmark script
# Run interactively; not part of automated CI.

devtools::load_all(".")

backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

benchmark_fit <- function(name, model, formula, data,
                          iter = 1000, warmup = 500, chains = 2, cores = 2,
                          seed = 1234) {
  start <- Sys.time()
  fit <- bmm(
    formula = formula,
    data = data,
    model = model,
    backend = backend,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    refresh = 0
  )
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  summ <- summary(fit)$fixed
  nuts <- try(brms::nuts_params(fit), silent = TRUE)
  sampler <- if (inherits(nuts, "try-error")) {
    data.frame(Parameter = character())
  } else {
    nuts
  }

  data.frame(
    case = name,
    backend = backend,
    elapsed_sec = elapsed,
    divergences = sum(sampler$Parameter == "divergent__" & sampler$Value > 0),
    max_treedepth = sum(sampler$Parameter == "treedepth__" & sampler$Value >= 10),
    min_bulk_ess = if ("Bulk_ESS" %in% colnames(summ)) min(summ[, "Bulk_ESS"]) else NA_real_,
    min_tail_ess = if ("Tail_ESS" %in% colnames(summ)) min(summ[, "Tail_ESS"]) else NA_real_,
    row.names = NULL
  )
}

cases <- list(
  list(
    name = "rating_uv_hard",
    data = rsdt_rating(
      n_per_cell = 150, n_subjects = 20,
      dprime = 2.2, criterion = 0.1, sdratio = 1.4,
      n_ratings = 6, spacing = 0.35, threshold_type = "parsimonious"
    ),
    model = sdt(
      version = "rating",
      response = paste0("r", 1:6),
      stimulus = "stimulus",
      threshold_type = "parsimonious"
    ),
    formula = bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  ),
  list(
    name = "dp_uv",
    data = rsdt_dp(
      n_per_cell = 150, n_subjects = 20,
      dprime = 1.4, criterion = 0, Ro = 0.35, Rn = 0.1,
      sdratio = 1.25, n_ratings = 6, spacing = 0.4
    ),
    model = sdt(
      version = "dp",
      response = paste0("r", 1:6),
      stimulus = "stimulus"
    ),
    formula = bmf(
      dprime ~ 1, criterion ~ 1, spacing ~ 1,
      Ro ~ 1, Rn ~ 1, sdratio ~ 1
    )
  ),
  list(
    name = "ranking_normal_uv",
    data = rsdt_ranking(
      n_per_cell = 120, n_subjects = 20,
      dprime = 1.3, m = 4, dist = "normal", sdratio = 1.25
    ),
    model = sdt(version = "ranking", response = "observed", rank = "rank", m = 4,
                dist = "normal"),
    formula = bmf(dprime ~ 1, sdratio ~ 1)
  ),
  list(
    name = "cdp_2way",
    data = rsdt_cdp(
      n_per_cell = 120, n_subjects = 20,
      dprimef = 0.8, dprimer = 1.0,
      criterion = 0, spacing = 0.45, rcrit = 0.25,
      n_ratings = 6
    ),
    model = sdt_cdp(
      new_response = c("n1", "n2", "n3"),
      old_know = c("k4", "k5", "k6"),
      old_remember = c("r4", "r5", "r6"),
      stimulus = "stimulus"
    ),
    formula = bmf(dprimef ~ 1, dprimer ~ 1, criterion ~ 1, spacing ~ 1, rcrit ~ 1)
  )
)

results <- do.call(
  rbind,
  lapply(cases, function(case) {
    benchmark_fit(case$name, case$model, case$formula, case$data)
  })
)

print(results)
