  /* Interference measurement model (Oberauer & Lin, 2017). The weights are a
   * Luce choice over the activations of the target, the ss - 1 active
   * non-targets, and a background b, which is the theory rather than a link, so
   * these versions keep it and gain only variable precision.
   *
   * c, a, s and b arrive on their natural scale because brms applies the
   * inverse link before calling the likelihood. The activations are assembled
   * on the log scale with log_sum_exp rather than as log(c * exp(-s * d) + a),
   * which is the same quantity without the intermediate exponentials.
   */

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
