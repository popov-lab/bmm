  // mixture2p (Zhang & Luck, 2008); parameters arrive on the natural scale, tau = 0 means constant precision
  real mixture2p_simple_core(real y, real mu, real kappa, real tau, real thetat,
                             int nodes, data vector logk, data vector dlogk,
                             data real logJ_min, data real dlogJ) {
    return circmix_vp_ld([cos(y - mu)]', [log(thetat)]', log1m(thetat),
                         kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // fixed-resolution slots: p_mem = min(1, K / ss)
  real mixture2p_slot_core(real y, real mu, real kappa, real tau, real K, int ss,
                           int nodes, data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    real p_mem = fmin(1.0, K / ss);
    return circmix_vp_ld([cos(y - mu)]', [log(p_mem)]', log1m(p_mem),
                         kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // slots plus averaging: the single item's weight given that it is held is one
  real mixture2p_slot_averaging_core(real y, real mu, real kappa, real tau,
                                     real K, int ss, int nodes, data vector logk,
                                     data vector dlogk, data real logJ_min,
                                     data real dlogJ) {
    return circmix_slot_averaging_ld([cos(y - mu)]', [0.0]', K, ss, kappa, tau,
                                     nodes, logk, dlogk, logJ_min, dlogJ);
  }
