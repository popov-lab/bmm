  {
    for (n in 1:N) {
      target += mixture2p_cd_log_prob(
        Y[n], probe_cd[n], mu[n], kappa[n], thetat[n], beta[n], CD_GRID, CD_DX
      );
    }
  }
