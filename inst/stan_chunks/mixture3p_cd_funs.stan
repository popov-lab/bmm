  // mixture3p change-detection cores; R twins are .mixture3p_cd_psame_*() in R/model_mixture3p_cd.R
  real mixture3p_cd_simple_core(int y, real mu, real kappa, real tau, real thetat,
                                real thetant, real criterion, int ss, real probe,
                                vector nt, int nodes, int rich, data vector logk,
                                data vector dlogk, data real logJ_min, data real dlogJ,
                                data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = mixture3p_simple_logw(thetat, thetant, ss);
    real p = circmix_cd_vp_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                 w[ss + 1], kappa, tau, criterion, exp(w[1]), rich,
                                 nodes, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }

  real mixture3p_cd_slot_core(int y, real mu, real kappa, real tau, real K, real pnt,
                              real criterion, int ss, real probe, vector nt,
                              int nodes, int rich, data vector logk,
                              data vector dlogk, data real logJ_min, data real dlogJ,
                              data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = mixture3p_slot_logw(K, pnt, ss);
    real p = circmix_cd_vp_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                 w[ss + 1], kappa, tau, criterion, exp(w[1]), rich,
                                 nodes, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }

  real mixture3p_cd_slot_averaging_core(int y, real mu, real kappa, real tau, real K,
                                        real pnt, real criterion, int ss, real probe,
                                        vector nt, int nodes, int rich,
                                        data vector logk, data vector dlogk,
                                        data real logJ_min, data real dlogJ,
                                        data vector gl_x, data vector gl_w) {
    vector[5] branch = circmix_slot_averaging_branches(K, ss, kappa, logk, dlogk,
                                                       logJ_min, dlogJ);
    real p = circmix_cd_slot_averaging_psame(probe - circmix_locations(mu, nt, ss),
                                             mixture3p_slot_averaging_logw(pnt, ss),
                                             branch, tau, criterion, rich, nodes,
                                             logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }
