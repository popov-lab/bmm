devtools::load_all(quiet = TRUE)

softmax_logz_reference <- function(mu, c, kappa, a = 0, s = 1,
                                   nt_features = numeric(),
                                   nt_distances = numeric(),
                                   lure_idx = rep(1, length(nt_features)),
                                   version = "full",
                                   parametrization = "sqrtexp",
                                   L = 4096L) {
  case <- .sdm_spectral_case(
    mu = mu,
    c = c,
    kappa = kappa,
    a = a,
    s = s,
    nt_features = nt_features,
    nt_distances = nt_distances,
    lure_idx = lure_idx,
    version = version,
    parametrization = parametrization,
    L = L
  )
  case$logz
}

generate_sdm_scenario <- function() {
  K <- sample(0:7, 1)
  version <- sample(c("abc", "bsc", "full"), 1)
  out <- list(
    version = version,
    mu = 0,
    c = runif(1, 1.5, 6.0),
    kappa = runif(1, 1.5, 15),
    a = runif(1, 0.05, 1.5),
    s = runif(1, 0.3, 3.0),
    nt_features = if (K > 0) runif(K, -pi, pi) else numeric(),
    nt_distances = if (K > 0) runif(K, 0.2, 2.0) else numeric(),
    lure_idx = if (K > 0) rep(1, K) else numeric()
  )

  hard_case <- runif(1)
  if (K >= 2 && hard_case < 0.33) {
    center <- runif(1, -0.35, 0.35)
    offsets <- c(0, runif(1, -0.12, 0.12), runif(1, -0.25, 0.25),
                 runif(1, -0.6, 0.6), runif(1, -1, 1), runif(1, -1.4, 1.4),
                 runif(1, -1.8, 1.8))
    out$nt_features[seq_len(K)] <- center + offsets[seq_len(K)]
  } else if (K >= 2 && hard_case < 0.66) {
    out$nt_features <- seq(-pi, pi - 2 * pi / K, length.out = K)
  }

  if (version == "abc") {
    out$nt_distances <- numeric()
  }
  if (version == "bsc") {
    out$a <- 0
  }

  out
}

find_minimal_L <- function(scenario, target_tol = 1e-6,
                           levels = c(32L, 64L, 128L, 256L, 512L)) {
  ref <- do.call(softmax_logz_reference, c(scenario, list(L = 4096L)))
  err <- vapply(
    levels,
    function(L) {
      est <- do.call(softmax_logz_reference, c(scenario, list(L = L)))
      abs(est - ref)
    },
    numeric(1)
  )
  idx <- which(err <= target_tol)[1]

  list(
    minimal_L = if (is.na(idx)) NA_integer_ else levels[idx],
    err = err,
    ref = ref
  )
}

assess_selector <- function(n_scenarios = 5000, seed = 123) {
  set.seed(seed)

  rows <- lapply(seq_len(n_scenarios), function(i) {
    scenario <- generate_sdm_scenario()
    diag <- .sdm_trial_diagnostics(
      weights = .sdm_item_spec(
        mu = scenario$mu,
        c = scenario$c,
        a = scenario$a,
        s = scenario$s,
        nt_features = scenario$nt_features,
        nt_distances = scenario$nt_distances,
        lure_idx = scenario$lure_idx,
        version = scenario$version
      )$weights,
      rel_angles = .sdm_item_spec(
        mu = scenario$mu,
        c = scenario$c,
        a = scenario$a,
        s = scenario$s,
        nt_features = scenario$nt_features,
        nt_distances = scenario$nt_distances,
        lure_idx = scenario$lure_idx,
        version = scenario$version
      )$rel_angles,
      kappa = scenario$kappa
    )
    minimal <- find_minimal_L(scenario)
    adaptive_L <- .sdm_get_L_general(
      kappa = diag$kappa,
      J = diag$J,
      w_sum = diag$w_sum,
      w_max = diag$w_max,
      R = diag$R,
      delta_min = diag$delta_min
    )
    adaptive_err <- do.call(softmax_logz_reference, c(scenario, list(L = adaptive_L)))

    data.frame(
      version = scenario$version,
      minimal_L = minimal$minimal_L,
      adaptive_L = adaptive_L,
      adaptive_abs_err = abs(adaptive_err - minimal$ref),
      kappa = diag$kappa,
      J = diag$J,
      w_sum = diag$w_sum,
      w_max = diag$w_max,
      R = diag$R,
      delta_min = diag$delta_min
    )
  })

  do.call(rbind, rows)
}

results <- assess_selector()

cat("Adaptive L counts\n")
print(table(results$adaptive_L, useNA = "ifany"))
cat("\nMinimal L counts\n")
print(table(results$minimal_L, useNA = "ifany"))
cat("\nAdaptive selector error summary\n")
print(summary(results$adaptive_abs_err))
cat("\nProportion within 1e-6:", mean(results$adaptive_abs_err <= 1e-6), "\n")
cat("Proportion within 1e-5:", mean(results$adaptive_abs_err <= 1e-5), "\n")
