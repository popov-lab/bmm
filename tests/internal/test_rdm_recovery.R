devtools::load_all()

cmdstan_dir <- file.path("tests", "internal", "cmdstanr_output")
if (!dir.exists(cmdstan_dir)) {
  dir.create(cmdstan_dir, recursive = TRUE)
}
options(cmdstanr_write_stan_file_dir = normalizePath(cmdstan_dir))

rdm_recovery_summary <- function(fit, truth) {
  draws <- as.data.frame(fit)

  transformed <- lapply(names(truth), function(par) {
    draw_name <- switch(
      par,
      driftc = "b_driftc_Intercept",
      drifte = "b_drifte_Intercept",
      gap = "b_gap_Intercept",
      ndt = "b_ndt_Intercept",
      s = "b_s_Intercept",
      sp = "b_sp_Intercept",
      corr = "b_corr_Intercept",
      lure = "b_lure_Intercept",
      npl = "b_npl_Intercept",
      stop2("Unsupported recovery parameter '{par}'")
    )

    x <- draws[[draw_name]]
    if (par %in% c("driftc", "drifte", "gap", "s", "sp", "corr", "lure", "npl")) {
      x <- exp(x)
    }
    if (par == "ndt") {
      x <- exp(x)
    }
    x
  })
  names(transformed) <- names(truth)

  out <- do.call(rbind, lapply(names(truth), function(par) {
    x <- transformed[[par]]
    lower <- stats::quantile(x, 0.025)
    upper <- stats::quantile(x, 0.975)
    data.frame(
      parameter = par,
      truth = truth[[par]],
      mean = mean(x),
      q2.5 = lower,
      q97.5 = upper,
      covered = lower <= truth[[par]] && truth[[par]] <= upper,
      row.names = NULL
    )
  }))

  out
}

fit_rdm_recovery_case <- function(name, data, formula, model, truth,
                                  iter = 1000, warmup = 500,
                                  chains = 2, cores = 2) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("Recovery case:", name, "\n")
  cat(strrep("=", 80), "\n", sep = "")

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
    silent = 2
  )

  out <- rdm_recovery_summary(fit, truth)
  print(out)
  cat("Coverage rate:", mean(out$covered), "\n")

  invisible(list(fit = fit, summary = out))
}

simulate_rdm_custom <- function(n, drift, counts, gap, ndt, s = 1, sp = 0) {
  resp_names <- names(counts)
  full_drift <- unlist(Map(rep, drift, counts), use.names = FALSE)
  raw <- rrdm(n = n, drift = full_drift, gap = gap, ndt = ndt, s = s, sp = sp)
  idx <- rep(seq_along(counts), counts)
  raw$response <- resp_names[idx[raw$response]]
  raw
}

simple_2choice_data <- rrdm(
  n = 400,
  drift = c(3.2, 1.2),
  gap = 1.0,
  ndt = 0.2
)

simple_4choice_data <- rrdm(
  n = 500,
  drift = c(3.0, 1.3, 1.3, 1.3),
  gap = 1.0,
  ndt = 0.22,
  s = 0.9,
  sp = 0.04
)

custom_fixed_data <- simulate_rdm_custom(
  n = 500,
  drift = c(corr = 3.1, lure = 1.8, npl = 1.0),
  counts = c(corr = 1, lure = 2, npl = 3),
  gap = 0.95,
  ndt = 0.23
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

recovery_cases <- list(
  simple_2choice = list(
    data = simple_2choice_data,
    formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1),
    model = rdm(rt = "rt", response = "response", n_choices = 2),
    truth = c(driftc = 3.2, drifte = 1.2, gap = 1.0, ndt = 0.2)
  ),
  simple_4choice_free_sp = list(
    data = simple_4choice_data,
    formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, s ~ 1, sp ~ 1),
    model = rdm(rt = "rt", response = "response", n_choices = 4),
    truth = c(driftc = 3.0, drifte = 1.3, gap = 1.0, ndt = 0.22, s = 0.9, sp = 0.04)
  ),
  custom_fixed = list(
    data = custom_fixed_data,
    formula = bmf(corr ~ 1, lure ~ 1, npl ~ 1, gap ~ 1, ndt ~ 1),
    model = rdm(
      rt = "rt",
      response = "response",
      version = "custom",
      accumulators = c(corr = 1, lure = 2, npl = 3)
    ),
    truth = c(corr = 3.1, lure = 1.8, npl = 1.0, gap = 0.95, ndt = 0.23)
  ),
  custom_free_sp = list(
    data = custom_free_data,
    formula = bmf(corr ~ 1, lure ~ 1, npl ~ 1, gap ~ 1, ndt ~ 1, s ~ 1, sp ~ 1),
    model = rdm(
      rt = "rt",
      response = "response",
      version = "custom",
      accumulators = c(corr = 1, lure = 2, npl = 2)
    ),
    truth = c(corr = 2.8, lure = 1.9, npl = 1.2, gap = 1.05, ndt = 0.24, s = 0.95, sp = 0.03)
  )
)

recovery_results <- lapply(names(recovery_cases), function(name) {
  case <- recovery_cases[[name]]
  fit_rdm_recovery_case(
    name = name,
    data = case$data,
    formula = case$formula,
    model = case$model,
    truth = case$truth
  )
})
names(recovery_results) <- names(recovery_cases)

invisible(recovery_results)
