  // cache normalization constant across sorted observations
  {
    real log_Z_cache;
    for (n in 1:N) {
      if (n == 1 || c[n] != c[n-1] || kappa[n] != kappa[n-1]) {
        log_Z_cache = sdm_cd_log_Z(c[n], kappa[n], SDM_CD_GRID);
      }
      target += sdm_cd_log_prob(Y[n], mu[n], c[n], kappa[n], beta[n],
                                 vreal1[n], log_Z_cache, SDM_CD_GRID);
    }
  }
