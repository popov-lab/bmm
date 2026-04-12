  {
    for (n in 1:N) {
      target += mixture3p_cd_log_prob(
        Y[n], probe_cd[n], mu[n], kappa[n], thetat[n], thetant[n], beta[n],
        cd_nt_features[n], cd_lure_idx[n], CD_GRID, CD_DX
      );
    }
  }
