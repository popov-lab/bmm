  {
    for (n in 1:N) {
      target += imm_full_cd_log_prob(
        Y[n], probe_cd[n], mu[n], kappa[n], c[n], a[n], s[n], beta[n],
        cd_nt_features[n], cd_nt_distances[n], cd_lure_idx[n], CD_GRID, CD_DX
      );
    }
  }
