  target += sdm_simple_run_ldenom_slice(c, kappa, sdm_gl_x, sdm_gl_w,
                                        sdm_wide_cosm1, sdm_wide_w, start,
                                        end, G_sdm_runs, sdm_run_start,
                                        sdm_run_count);
  target += -(log2()+log(pi()))*N;
