  // IMM change detection (Lin & Oberauer, 2022): imm weights through the generic evaluator; R twin .dimm_cd_version()
  real imm_cd_full_core(int y, real mu, real kappa, real tau, real c, real a,
                        real s, real b, real criterion, int ss, real probe,
                        vector nt, vector dist, int nodes, int rich,
                        data vector logk, data vector dlogk, data real logJ_min,
                        data real dlogJ, data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = imm_full_logw(c, a, s, b, ss, dist);
    real p = circmix_cd_vp_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                 w[ss + 1], kappa, tau, criterion, exp(w[1]), rich,
                                 nodes, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }

  real imm_cd_bsc_core(int y, real mu, real kappa, real tau, real c, real s,
                       real b, real criterion, int ss, real probe, vector nt,
                       vector dist, int nodes, int rich, data vector logk,
                       data vector dlogk, data real logJ_min, data real dlogJ,
                       data vector gl_x, data vector gl_w) {
    vector[ss + 1] w = imm_bsc_logw(c, s, b, ss, dist);
    real p = circmix_cd_vp_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                 w[ss + 1], kappa, tau, criterion, exp(w[1]), rich,
                                 nodes, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }

  real imm_cd_abc_core(int y, real mu, real kappa, real tau, real c, real a,
                       real b, real criterion, int ss, real probe, vector nt,
                       int nodes, int rich, data vector logk, data vector dlogk,
                       data real logJ_min, data real dlogJ, data vector gl_x,
                       data vector gl_w) {
    vector[ss + 1] w = imm_abc_logw(c, a, b, ss);
    real p = circmix_cd_vp_psame(probe - circmix_locations(mu, nt, ss), head(w, ss),
                                 w[ss + 1], kappa, tau, criterion, exp(w[1]), rich,
                                 nodes, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
    return cd_bernoulli_lpmf(y | p);
  }
