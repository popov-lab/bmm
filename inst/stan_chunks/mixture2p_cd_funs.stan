  // mixture2p_cd: the mixture2p retrieval weights judged by the change-detection rule; y = 1 is "change"
  real mixture2p_cd_simple_core(int y, real mu, real kappa, real tau, real thetat,
                                real criterion, real probe, int nodes, int rich,
                                data vector logk, data vector dlogk, data real logJ_min,
                                data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[2] w = mixture2p_simple_logw(thetat);
    real p = circmix_cd_vp_psame([probe - mu]', head(w, 1), w[2], kappa, tau, criterion,
                                 exp(w[1]), rich, nodes, logk, dlogk, logJ_min, dlogJ,
                                 gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }

  real mixture2p_cd_slot_core(int y, real mu, real kappa, real tau, real K,
                              real criterion, int ss, real probe, int nodes, int rich,
                              data vector logk, data vector dlogk, data real logJ_min,
                              data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[2] w = mixture2p_slot_logw(K, ss);
    real p = circmix_cd_vp_psame([probe - mu]', head(w, 1), w[2], kappa, tau, criterion,
                                 exp(w[1]), rich, nodes, logk, dlogk, logJ_min, dlogJ,
                                 gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }

  real mixture2p_cd_slot_averaging_core(int y, real mu, real kappa, real tau, real K,
                                        real criterion, int ss, real probe, int nodes,
                                        int rich, data vector logk, data vector dlogk,
                                        data real logJ_min, data real dlogJ,
                                        data vector gl_x, data vector gl_w) {
    vector[5] branch = circmix_slot_averaging_branches(K, ss, kappa, logk, dlogk,
                                                       logJ_min, dlogJ);
    real p = circmix_cd_slot_averaging_psame([probe - mu]', [0.0]', branch, tau,
                                             criterion, rich, nodes, logk, dlogk,
                                             logJ_min, dlogJ, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }
