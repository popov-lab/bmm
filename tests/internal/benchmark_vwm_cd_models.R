devtools::load_all()

benchmark_vwm_cd_models <- function(iter = 400, warmup = 200, backend = "rstan") {
  set.seed(1)

  simulate_cases <- list(
    sdm = list(
      data = {
        dat <- data.frame(target = 0, probe = runif(80, -pi, pi))
        dat$resp <- rsdm_cd(80, probe = dat$probe, c = 4, kappa = 3)
        dat
      },
      formula = bmf(c ~ 1, kappa ~ 1),
      model = sdm_cd(response = "resp", probe = "probe", target = "target")
    ),
    mixture2p = list(
      data = {
        dat <- data.frame(target = 0, probe = runif(80, -pi, pi))
        dat$resp <- rmixture2p_cd(80, probe = dat$probe, kappa = 5, p_target = 0.75)
        dat
      },
      formula = bmf(kappa ~ 1, thetat ~ 1),
      model = mixture2p_cd(response = "resp", probe = "probe", target = "target")
    ),
    mixture3p = list(
      data = {
        dat <- data.frame(
          target = 0,
          probe = runif(80, -pi, pi),
          nt1 = runif(80, -pi, pi),
          nt2 = runif(80, -pi, pi),
          ss = 3
        )
        dat$resp <- vapply(seq_len(nrow(dat)), function(i) {
          rmixture3p_cd(
            1,
            probe = dat$probe[i],
            nt_features = c(dat$nt1[i], dat$nt2[i]),
            lure_idx = c(1, 1),
            p_target = 0.6,
            p_nontarget = 0.25,
            kappa = 5
          )
        }, numeric(1))
        dat
      },
      formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1),
      model = mixture3p_cd(
        response = "resp",
        probe = "probe",
        target = "target",
        nt_features = c("nt1", "nt2"),
        set_size = "ss"
      )
    ),
    imm = list(
      data = {
        dat <- data.frame(
          target = 0,
          probe = runif(80, -pi, pi),
          nt1 = runif(80, -pi, pi),
          nt2 = runif(80, -pi, pi),
          d1 = runif(80, 0.3, 1.3),
          d2 = runif(80, 0.4, 1.6),
          ss = 3
        )
        dat$resp <- vapply(seq_len(nrow(dat)), function(i) {
          rimm_cd(
            1,
            probe = dat$probe[i],
            nt_features = c(dat$nt1[i], dat$nt2[i]),
            nt_distances = c(dat$d1[i], dat$d2[i]),
            lure_idx = c(1, 1),
            kappa = 5,
            c = 2,
            a = 0.5,
            s = 2
          )
        }, numeric(1))
        dat
      },
      formula = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
      model = imm_cd(
        response = "resp",
        probe = "probe",
        target = "target",
        nt_features = c("nt1", "nt2"),
        nt_distances = c("d1", "d2"),
        set_size = "ss"
      )
    )
  )

  results <- lapply(names(simulate_cases), function(name) {
    case <- simulate_cases[[name]]
    dat <- case$data
    t0 <- proc.time()[["elapsed"]]
    fit <- bmm(
      formula = case$formula,
      data = dat,
      model = case$model,
      backend = backend,
      chains = 2,
      cores = 2,
      iter = iter,
      warmup = warmup,
      seed = 123,
      refresh = 0
    )
    diag <- cd_fit_diagnostics(fit)
    elapsed <- proc.time()[["elapsed"]] - t0
    min_ess_per_sec <- diag$min_neff_ratio * brms::ndraws(fit) / elapsed
    cbind(model = name, elapsed = elapsed, min_ess_per_sec = min_ess_per_sec, diag)
  })

  do.call(rbind, results)
}

benchmark_vwm_cd_models()
