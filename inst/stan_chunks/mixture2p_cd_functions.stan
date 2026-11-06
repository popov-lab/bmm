  // mixture2p_cd: the mixture2p retrieval weights (mixture2p_functions.stan) judged by the change-detection rule
  // (circmix_cd_functions.stan); y = 1 is "change". The observer's prior of storage is the target's retrieval weight,
  // except under slot averaging, where circmix_cd_slot_averaging_psame() sets it
  real mixture2p_cd_simple_core(int y, real mu, real kappa, real thetat, real criterion,
                                real probe, int rich, data vector logk, data vector dlogk,
                                data real logJ_min, data real dlogJ, data vector gl_x,
                                data vector gl_w) {
    vector[2] w = mixture2p_simple_logw(thetat);
    return cd_bernoulli_lpmf(y | circmix_cd_psame([probe - mu]', head(w, 1), w[2], kappa,
                                                  criterion, thetat, gl_x, gl_w));
  }

  real mixture2p_cd_slot_core(int y, real mu, real kappa, real K, real criterion, int ss,
                              real probe, int rich, data vector logk, data vector dlogk,
                              data real logJ_min, data real dlogJ, data vector gl_x,
                              data vector gl_w) {
    vector[2] w = mixture2p_slot_logw(K, ss);
    return cd_bernoulli_lpmf(y | circmix_cd_psame([probe - mu]', head(w, 1), w[2], kappa,
                                                  criterion, exp(w[1]), gl_x, gl_w));
  }

  real mixture2p_cd_slot_averaging_core(int y, real mu, real kappa, real K, real criterion,
                                        int ss, real probe, int rich, data vector logk,
                                        data vector dlogk, data real logJ_min,
                                        data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[4] branch = circmix_slot_averaging_branches(K, ss, kappa, logk, dlogk, logJ_min,
                                                       dlogJ);
    return cd_bernoulli_lpmf(y | circmix_cd_slot_averaging_psame([probe - mu]', [0.0]', branch,
                                                                 criterion, rich, gl_x, gl_w));
  }
