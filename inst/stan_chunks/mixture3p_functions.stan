  // mixture3p cores; R twins are .dmixture3p_*() in R/distributions.R
  real mixture3p_simple_core(real y, real mu, real kappa, real tau, real thetat,
                             real thetant, int ss, vector nt, int nodes,
                             data vector logk, data vector dlogk,
                             data real logJ_min, data real dlogJ) {
    int n_nt = ss - 1;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    cosd[1] = cos(y - mu);
    logw[1] = thetat;
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = thetant - log(n_nt);
    }
    // guessing is the softmax reference category (log weight 0)
    real total = log_sum_exp(append_row(logw, 0.0));
    return circmix_vp_ld(cosd, logw - total, -total, kappa, tau, nodes,
                         logk, dlogk, logJ_min, dlogJ);
  }

  // capacity rule sets p_mem, ptarget sets which stored item is reported;
  // at set size 1 the target is the only stored item
  real mixture3p_slot_core(real y, real mu, real kappa, real tau, real K,
                           real ptarget, int ss, vector nt, int nodes,
                           data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    int n_nt = ss - 1;
    real p_mem = fmin(1.0, K / ss);
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    cosd[1] = cos(y - mu);
    logw[1] = log(p_mem) + (n_nt > 0 ? log(ptarget) : 0.0);
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = log(p_mem) + log1m(ptarget) - log(n_nt);
    }
    return circmix_vp_ld(cosd, logw, log1m(p_mem), kappa, tau, nodes,
                         logk, dlogk, logJ_min, dlogJ);
  }

  // slot averaging crossed with swapping; swap weights shared by both slot counts
  real mixture3p_slot_averaging_core(real y, real mu, real kappa, real tau,
                                     real K, real ptarget, int ss, vector nt,
                                     int nodes, data vector logk,
                                     data vector dlogk, data real logJ_min,
                                     data real dlogJ) {
    int n_nt = ss - 1;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    cosd[1] = cos(y - mu);
    logw[1] = n_nt > 0 ? log(ptarget) : 0.0;
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = log1m(ptarget) - log(n_nt);
    }
    return circmix_slot_averaging_ld(cosd, logw, K, ss, kappa, tau, nodes,
                                     logk, dlogk, logJ_min, dlogJ);
  }
