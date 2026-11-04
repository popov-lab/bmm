  // IMM (Oberauer & Lin, 2017); R twin .imm_log_weights() in R/distributions.R

  real imm_abc_core(real y, real mu, real kappa, real tau, real c, real a,
                    real b, int ss, vector nt, int nodes, data vector logk,
                    data vector dlogk, data real logJ_min, data real dlogJ) {
    int n_nt = ss - 1;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    real log_c = log(c);
    real log_a = log(a);

    cosd[1] = cos(y - mu);
    logw[1] = log_sum_exp(log_c, log_a);
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = log_a;
    }

    real total = log_sum_exp(append_row(logw, log(b)));
    return circmix_vp_ld(cosd, logw - total, log(b) - total, kappa, tau, nodes,
                         logk, dlogk, logJ_min, dlogJ);
  }

  real imm_bsc_core(real y, real mu, real kappa, real tau, real c, real s,
                    real b, int ss, vector nt, vector dist, int nodes,
                    data vector logk, data vector dlogk, data real logJ_min,
                    data real dlogJ) {
    int n_nt = ss - 1;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    real log_c = log(c);

    cosd[1] = cos(y - mu);
    logw[1] = log_c;
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = log_c - s * dist[j];
    }

    real total = log_sum_exp(append_row(logw, log(b)));
    return circmix_vp_ld(cosd, logw - total, log(b) - total, kappa, tau, nodes,
                         logk, dlogk, logJ_min, dlogJ);
  }

  real imm_full_core(real y, real mu, real kappa, real tau, real c, real a,
                     real s, real b, int ss, vector nt, vector dist, int nodes,
                     data vector logk, data vector dlogk, data real logJ_min,
                     data real dlogJ) {
    int n_nt = ss - 1;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    real log_c = log(c);
    real log_a = log(a);

    cosd[1] = cos(y - mu);
    logw[1] = log_sum_exp(log_c, log_a);
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = log_sum_exp(log_c - s * dist[j], log_a);
    }

    real total = log_sum_exp(append_row(logw, log(b)));
    return circmix_vp_ld(cosd, logw - total, log(b) - total, kappa, tau, nodes,
                         logk, dlogk, logJ_min, dlogJ);
  }
