  // mixture2p (Zhang & Luck, 2008): log weights of memory and guessing responses
  vector mixture2p_simple_logw(real thetat) {
    return [log(thetat), log1m(thetat)]';
  }

  // log weights with K memory slots: P(memory) = min(1, K / set size)
  vector mixture2p_slot_logw(real K, int ss) {
    real p_mem = fmin(1.0, K / ss);
    return [log(p_mem), log1m(p_mem)]';
  }

  // log density of a response under the two-parameter mixture model
  real mixture2p_simple_core(real y, real mu, real kappa, real tau, real thetat,
                             int nodes, data vector logk, data vector dlogk,
                             data real logJ_min, data real dlogJ) {
    vector[2] w = mixture2p_simple_logw(thetat);
    return circmix_vp_ld([cos(y - mu)]', head(w, 1), w[2],
                         kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // log density of a response under the slot model
  real mixture2p_slot_core(real y, real mu, real kappa, real tau, real K, int ss,
                           int nodes, data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    vector[2] w = mixture2p_slot_logw(K, ss);
    return circmix_vp_ld([cos(y - mu)]', head(w, 1), w[2],
                         kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // log density of a response under the slots-plus-averaging model
  real mixture2p_slot_averaging_core(real y, real mu, real kappa, real tau,
                                     real K, int ss, int nodes, data vector logk,
                                     data vector dlogk, data real logJ_min,
                                     data real dlogJ) {
    return circmix_slot_averaging_ld([cos(y - mu)]', [0.0]', K, ss, kappa, tau,
                                     nodes, logk, dlogk, logJ_min, dlogJ);
  }
