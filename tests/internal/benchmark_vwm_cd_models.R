# Sampling speed of the change detection models on comparable data.
#
# All four share the same decision rule and differ only in the retrieval
# distribution they integrate over the decision arc, so the differences here are
# the cost of that integral: 32 Gauss-Legendre nodes per von Mises component,
# against an adaptive split rule for sdm_cd whose density is far more peaked.
#
# Run with:  Rscript tests/internal/benchmark_vwm_cd_models.R

devtools::load_all()

# mean per-chain warmup and sampling time, for either backend
chain_times <- function(fit) {
  if (inherits(fit$fit, "CmdStanMCMC")) {
    tm <- fit$fit$time()$chains
    return(c(warmup = mean(tm$warmup), sampling = mean(tm$sampling)))
  }
  tm <- rstan::get_elapsed_time(fit$fit)
  c(warmup = mean(tm[, "warmup"]), sampling = mean(tm[, "sample"]))
}

benchmark_vwm_cd_models <- function(n = 800, iter = 1000, warmup = 500,
                                    chains = 2, backend = "cmdstanr",
                                    seed = 20260805) {
  set.seed(seed)

  base <- data.frame(
    target = 0,
    probe = runif(n, -pi, pi),
    nt1 = runif(n, -pi, pi),
    nt2 = runif(n, -pi, pi),
    d1 = runif(n, 0.2, 2.5),
    d2 = runif(n, 0.2, 2.5),
    ss = 3
  )
  nt <- cbind(base$nt1, base$nt2)
  ntd <- cbind(base$d1, base$d2)

  cases <- list(
    mixture2p = list(
      resp = rmixture2p_cd(n, base$probe, kappa = 6, p_mem = 0.75),
      formula = bmf(kappa ~ 1, thetat ~ 1),
      model = mixture2p_cd(response = "resp", probe = "probe", target = "target")
    ),
    mixture3p = list(
      resp = rmixture3p_cd(n, base$probe, nt_features = nt, kappa = 6,
                           thetat = 1.2, thetant = 0.2),
      formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1),
      model = mixture3p_cd(response = "resp", probe = "probe", target = "target",
                           nt_features = c("nt1", "nt2"), set_size = "ss")
    ),
    imm = list(
      resp = rimm_cd(n, base$probe, nt_features = nt, nt_distances = ntd,
                     kappa = 6, c = 5, a = 0.6, s = 1.35),
      formula = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
      model = imm_cd(response = "resp", probe = "probe", target = "target",
                     nt_features = c("nt1", "nt2"),
                     nt_distances = c("d1", "d2"), set_size = "ss")
    ),
    sdm = list(
      resp = rsdm_cd(n, base$probe, c = 5, kappa = 4),
      formula = bmf(c ~ 1, kappa ~ 1),
      model = sdm_cd(response = "resp", probe = "probe", target = "target")
    )
  )

  out <- lapply(names(cases), function(nm) {
    case <- cases[[nm]]
    dat <- base
    dat$resp <- case$resp

    fit <- bmm(case$formula, dat, case$model, backend = backend,
               chains = chains, cores = chains, iter = iter, warmup = warmup,
               refresh = 0, silent = 2)

    tm <- chain_times(fit)
    ess <- suppressWarnings(min(brms::neff_ratio(fit), na.rm = TRUE)) *
      (iter - warmup) * chains

    data.frame(
      model = nm,
      sampling_s = round(tm[["sampling"]], 1),
      warmup_s = round(tm[["warmup"]], 1),
      min_ess = round(ess),
      ess_per_s = round(ess / tm[["sampling"]], 1),
      max_rhat = round(suppressWarnings(max(brms::rhat(fit), na.rm = TRUE)), 4)
    )
  })

  res <- do.call(rbind, out)
  res$relative <- round(res$sampling_s / res$sampling_s[res$model == "mixture2p"], 2)
  res
}

if (sys.nframe() == 0L) {
  print(benchmark_vwm_cd_models(), row.names = FALSE)
}
