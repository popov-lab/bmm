  // imm_cd: the imm retrieval weights (imm_functions.stan) judged by the change-detection rule
  // (circmix_cd_functions.stan); y = 1 is "change". R twin: .imm_cd_psame() in R/distributions.R. The observer's
  // prior p_s is the target's retrieval weight; precision is one state, so rich does not enter
  real imm_cd_full_core(int y, real mu, real kappa, real c, real a, real s, real b,
                        real criterion, int ss, real probe, vector nt, vector dist, int rich,
                        data vector logk, data vector dlogk, data real logJ_min,
                        data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = imm_full_logw(c, a, s, b, ss, dist);
    return cd_bernoulli_lpmf(y | circmix_cd_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                                  w[ss + 1], kappa, criterion, exp(w[1]), gl_x, gl_w));
  }

  real imm_cd_bsc_core(int y, real mu, real kappa, real c, real s, real b,
                       real criterion, int ss, real probe, vector nt, vector dist, int rich,
                       data vector logk, data vector dlogk, data real logJ_min,
                       data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = imm_bsc_logw(c, s, b, ss, dist);
    return cd_bernoulli_lpmf(y | circmix_cd_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                                  w[ss + 1], kappa, criterion, exp(w[1]), gl_x, gl_w));
  }

  real imm_cd_abc_core(int y, real mu, real kappa, real c, real a, real b,
                       real criterion, int ss, real probe, vector nt, int rich,
                       data vector logk, data vector dlogk, data real logJ_min,
                       data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = imm_abc_logw(c, a, b, ss);
    return cd_bernoulli_lpmf(y | circmix_cd_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                                  w[ss + 1], kappa, criterion, exp(w[1]), gl_x, gl_w));
  }
