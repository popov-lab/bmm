  // mixture3p_cd: the mixture3p retrieval weights (mixture3p_functions.stan) judged by the change-detection rule
  // (circmix_cd_functions.stan); y = 1 is "change". R twin: .mixture3p_cd_psame() in R/distributions.R. The observer's
  // prior p_s is the target's retrieval weight, except under slot averaging, where circmix_cd_slot_averaging_psame()
  // sets it
  real mixture3p_cd_simple_core(int y, real mu, real kappa, real thetat, real thetant,
                                real criterion, int ss, real probe, vector nt, int rich,
                                data vector logk, data vector dlogk, data real logJ_min,
                                data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = mixture3p_simple_logw(thetat, thetant, ss);
    return cd_bernoulli_lpmf(y | circmix_cd_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                                  w[ss + 1], kappa, criterion, exp(w[1]), gl_x, gl_w));
  }

  real mixture3p_cd_slot_core(int y, real mu, real kappa, real K, real ptarget,
                              real criterion, int ss, real probe, vector nt, int rich,
                              data vector logk, data vector dlogk, data real logJ_min,
                              data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = mixture3p_slot_logw(K, ptarget, ss);
    return cd_bernoulli_lpmf(y | circmix_cd_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                                  w[ss + 1], kappa, criterion, exp(w[1]), gl_x, gl_w));
  }

  // the swap weights are given that the item is held; how many slots it holds is independent of which item is
  // reported, so they serve both slot counts
  real mixture3p_cd_slot_averaging_core(int y, real mu, real kappa, real K, real ptarget,
                                        real criterion, int ss, real probe, vector nt, int rich,
                                        data vector logk, data vector dlogk, data real logJ_min,
                                        data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[4] branch = circmix_slot_averaging_branches(K, ss, kappa, logk, dlogk, logJ_min, dlogJ);
    return cd_bernoulli_lpmf(y | circmix_cd_slot_averaging_psame(probe - circmix_locations(mu, nt, ss),
                                                                 mixture3p_slot_averaging_logw(ptarget, ss),
                                                                 branch, criterion, rich, gl_x, gl_w));
  }
