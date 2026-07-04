  // adaptive calculation of the normalization constant
  for (g in 1:G_sdm_runs) {
    int n = sdm_run_start[g];
    real z = sdm_simple_ldenom_chquad_adaptive(c[n], kappa[n], COSN);
    target += sdm_run_count[g] * z;
  }
  target += -(log2()+log(pi()))*N;
