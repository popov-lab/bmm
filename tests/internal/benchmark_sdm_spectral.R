devtools::load_all(quiet = TRUE)

select_benchmark_rows <- function(dat, n_per_set_size = 60L) {
  split_dat <- split(dat, dat$set_size)
  rows <- lapply(split_dat, function(x) utils::head(x, n_per_set_size))
  do.call(rbind, rows)
}

get_sdm_benchmark_params <- function(dat, version) {
  set_size <- as.integer(as.character(dat$set_size))

  out <- list(
    mu = rep(0, nrow(dat)),
    c = 1.8 + 0.18 * set_size,
    kappa = 2.2 + 0.55 * set_size,
    a = rep(0, nrow(dat)),
    s = rep(0, nrow(dat))
  )

  if (version %in% c("abc", "full")) {
    out$a <- 0.18 + 0.08 * set_size
  }
  if (version %in% c("bsc", "full")) {
    out$s <- 0.5 + 0.2 * set_size
  }

  out
}

build_density_arg_rows <- function(dat, version) {
  params <- get_sdm_benchmark_params(dat, version)
  nt_features <- as.matrix(dat[, paste0("col_nt", 1:7), drop = FALSE])
  nt_distances <- as.matrix(dat[, paste0("dist_nt", 1:7), drop = FALSE])
  lure_idx <- as.matrix(dat[, grep("^LureIdx", names(dat), value = TRUE), drop = FALSE])

  lapply(seq_len(nrow(dat)), function(i) {
    list(
      x = dat$dev_rad[i],
      mu = params$mu[i],
      c = params$c[i],
      kappa = params$kappa[i],
      a = params$a[i],
      s = params$s[i],
      nt_features = nt_features[i, ],
      nt_distances = nt_distances[i, ],
      lure_idx = lure_idx[i, ],
      version = version
    )
  })
}

evaluate_density_mode <- function(arg_rows, mode = c("adaptive", "fixed_256", "fixed_512")) {
  mode <- match.arg(mode)
  L <- switch(
    mode,
    adaptive = NULL,
    fixed_256 = 256L,
    fixed_512 = 512L
  )

  t0 <- proc.time()[["elapsed"]]
  log_density <- vapply(
    arg_rows,
    function(args) {
      do.call(
        .dsdm_spectral_single,
        c(
          args,
          list(
            log = TRUE,
            parametrization = "sqrtexp",
            L = L
          )
        )
      )
    },
    numeric(1)
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  data.frame(
    mode = mode,
    elapsed_sec = elapsed,
    mean_log_density = mean(log_density),
    sum_log_density = sum(log_density)
  )
}

summarize_trial_levels <- function(arg_rows) {
  rows <- lapply(arg_rows, function(args) {
    spec <- .sdm_item_spec(
      mu = args$mu,
      c = args$c,
      a = args$a,
      s = args$s,
      nt_features = args$nt_features,
      nt_distances = args$nt_distances,
      lure_idx = args$lure_idx,
      version = args$version
    )
    diag <- .sdm_trial_diagnostics(spec$weights, spec$rel_angles, args$kappa)
    data.frame(
      L = .sdm_get_L_general(
        kappa = diag$kappa,
        J = diag$J,
        w_sum = diag$w_sum,
        w_max = diag$w_max,
        R = diag$R,
        delta_min = diag$delta_min
      ),
      J = diag$J,
      kappa = diag$kappa
    )
  })

  out <- do.call(rbind, rows)
  attr(out, "table") <- table(out$L)
  out
}

run_short_fit <- function(version, dat, iter = 250, warmup = 125) {
  formula <- switch(
    version,
    abc = bmf(
      c ~ 0 + set_size + (0 + set_size || ID),
      a ~ 0 + set_size + (0 + set_size || ID),
      kappa ~ 0 + set_size + (0 + set_size || ID)
    ),
    bsc = bmf(
      c ~ 0 + set_size + (0 + set_size || ID),
      s ~ 0 + set_size + (0 + set_size || ID),
      kappa ~ 0 + set_size + (0 + set_size || ID)
    ),
    full = bmf(
      c ~ 0 + set_size + (0 + set_size || ID),
      a ~ 0 + set_size + (0 + set_size || ID),
      s ~ 0 + set_size + (0 + set_size || ID),
      kappa ~ 0 + set_size + (0 + set_size || ID)
    )
  )
  model <- switch(
    version,
    abc = sdm(
      resp_error = "dev_rad",
      nt_features = paste0("col_nt", 1:7),
      set_size = "set_size",
      version = "abc"
    ),
    bsc = sdm(
      resp_error = "dev_rad",
      nt_features = paste0("col_nt", 1:7),
      nt_distances = paste0("dist_nt", 1:7),
      set_size = "set_size",
      version = "bsc"
    ),
    full = sdm(
      resp_error = "dev_rad",
      nt_features = paste0("col_nt", 1:7),
      nt_distances = paste0("dist_nt", 1:7),
      set_size = "set_size",
      version = "full"
    )
  )

  t0 <- proc.time()[["elapsed"]]
  fit <- bmm(
    formula = formula,
    data = dat,
    model = model,
    backend = "cmdstanr",
    chains = 1,
    cores = 1,
    iter = iter,
    warmup = warmup,
    refresh = 0,
    seed = 123,
    silent = 2
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  diag <- fit$fit$sampler_diagnostics()
  divergences <- sum(diag[, , "divergent__"])
  max_treedepth <- max(diag[, , "treedepth__"])

  data.frame(
    version = version,
    elapsed_sec = elapsed,
    divergences = divergences,
    max_treedepth = max_treedepth
  )
}

run_benchmarks <- function(n_density_per_set_size = 60L,
                           n_fit_per_set_size = 35L,
                           fit_iter = 250,
                           fit_warmup = 125) {
  dat_density <- select_benchmark_rows(oberauer_lin_2017, n_density_per_set_size)
  dat_fit <- select_benchmark_rows(oberauer_lin_2017, n_fit_per_set_size)
  versions <- c("abc", "bsc", "full")

  density_results <- lapply(versions, function(version) {
    args <- build_density_arg_rows(dat_density, version)
    density_bench <- do.call(rbind, lapply(
      c("adaptive", "fixed_256", "fixed_512"),
      function(mode) evaluate_density_mode(args, mode = mode)
    ))
    level_info <- summarize_trial_levels(args)
    density_bench$version <- version
    density_bench$mean_L <- mean(level_info$L)
    density_bench$median_L <- stats::median(level_info$L)
    density_bench$max_L <- max(level_info$L)
    attr(density_bench, "L_table") <- attr(level_info, "table")
    density_bench
  })

  fit_results <- do.call(
    rbind,
    lapply(versions, run_short_fit, dat = dat_fit, iter = fit_iter, warmup = fit_warmup)
  )

  list(
    density = density_results,
    fit = fit_results
  )
}

results <- run_benchmarks()

cat("Density benchmark results\n")
print(do.call(rbind, results$density))
cat("\nAdaptive L distributions\n")
for (res in results$density) {
  cat("\nVersion:", unique(res$version), "\n")
  print(attr(res, "L_table"))
}
cat("\nShort fit results\n")
print(results$fit)
