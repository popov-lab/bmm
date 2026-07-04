  // adaptive calculation of the normalization constant
  for (g in 1:G_sdm_runs) {
    if (sdm_run_start[g] >= start && sdm_run_start[g] <= end) {
      int n = sdm_run_start[g] - start + 1;
      real z = sdm_simple_ldenom_chquad_adaptive(c[n], kappa[n], COSN);
      target += sdm_run_count[g] * z;
    }
  }
  target += -(log2()+log(pi()))*N;
