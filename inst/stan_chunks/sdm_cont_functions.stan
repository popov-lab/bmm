  // log of the numerator of the likelihood of the sdm() cells other than the
  // default; sdm_kernel_functions.stan defines the cell
  real sdm_cont_lpdf(vector y, vector mu, vector c, vector kappa, data int sdm_activation,
                     data int sdm_scaling) {
    return sum(exp(sdm_log_cstar(c, kappa, sdm_activation, sdm_scaling)
                   + sdm_log_s(y - mu, kappa, sdm_activation)));
  }

  real sdm_cont_ldenom(real c, real kappa, int activation, int scaling, data vector gl_x,
                       data vector gl_w, data vector wide_cosm1, data vector wide_w) {
    real log_peak = sdm_log_cstar(c, kappa, activation, scaling);
    if (activation == 1) {
      return sdm_vm_ldenom_peak(log_peak, kappa, gl_x, gl_w, wide_cosm1, wide_w);
    }
    return sdm_laplace_ldenom(log_peak, kappa, gl_x, gl_w);
  }

  real sdm_cont_run_ldenom(vector c, vector kappa, int activation, int scaling,
                           data vector gl_x, data vector gl_w, data vector wide_cosm1,
                           data vector wide_w, int G_runs, array[] int run_start,
                           array[] int run_count) {
    real out = 0;
    for (g in 1:G_runs) {
      int n = run_start[g];
      out += run_count[g] * sdm_cont_ldenom(c[n], kappa[n], activation, scaling, gl_x, gl_w,
                                            wide_cosm1, wide_w);
    }
    return out;
  }

  real sdm_cont_run_ldenom_slice(vector c, vector kappa, int activation, int scaling,
                                 data vector gl_x, data vector gl_w, data vector wide_cosm1,
                                 data vector wide_w, int start, int end, int G_runs,
                                 array[] int run_start, array[] int run_count) {
    real out = 0;
    for (g in 1:G_runs) {
      if (run_start[g] >= start && run_start[g] <= end) {
        int n = run_start[g] - start + 1;
        out += run_count[g] * sdm_cont_ldenom(c[n], kappa[n], activation, scaling, gl_x, gl_w,
                                              wide_cosm1, wide_w);
      }
    }
    return out;
  }
