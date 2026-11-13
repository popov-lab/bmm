  target += sdm_cont_run_ldenom(c, kappa, sdm_activation, sdm_scaling, sdm_gl_x, sdm_gl_w,
                                sdm_wide_cosm1, sdm_wide_w, G_sdm_runs, sdm_run_start,
                                sdm_run_count);
  target += -(log2()+log(pi()))*N;
