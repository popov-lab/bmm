  // IMM-abc (Oberauer & Lin, 2017): normalised log weights of target, non-targets and guessing
  vector imm_abc_logw(real c, real a, real b, int ss) {
    int n_nt = ss - 1;
    vector[ss + 1] logw;
    real log_a = log(a);
    logw[1] = log_sum_exp(log(c), log_a);
    for (j in 1:n_nt) {
      logw[j + 1] = log_a;
    }
    logw[ss + 1] = log(b);
    return logw - log_sum_exp(logw);
  }

  // IMM-bsc: the same weights, with non-target weights that fall with the distance to the target
  vector imm_bsc_logw(real c, real s, real b, int ss, vector dist) {
    int n_nt = ss - 1;
    vector[ss + 1] logw;
    real log_c = log(c);
    logw[1] = log_c;
    for (j in 1:n_nt) {
      logw[j + 1] = log_c - s * dist[j];
    }
    logw[ss + 1] = log(b);
    return logw - log_sum_exp(logw);
  }

  // full IMM: the same weights, combining the abc and the bsc version
  vector imm_full_logw(real c, real a, real s, real b, int ss, vector dist) {
    int n_nt = ss - 1;
    vector[ss + 1] logw;
    real log_c = log(c);
    real log_a = log(a);
    logw[1] = log_sum_exp(log_c, log_a);
    for (j in 1:n_nt) {
      logw[j + 1] = log_sum_exp(log_c - s * dist[j], log_a);
    }
    logw[ss + 1] = log(b);
    return logw - log_sum_exp(logw);
  }

  // log density of a response under the IMM-abc
  real imm_abc_core(real y, real mu, real kappa, real tau, real c, real a,
                    real b, int ss, vector nt, int nodes, data vector logk,
                    data vector dlogk, data real logJ_min, data real dlogJ) {
    vector[ss + 1] w = imm_abc_logw(c, a, b, ss);
    return circmix_vp_ld(cos(y - circmix_locations(mu, nt, ss)), head(w, ss),
                         w[ss + 1], kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // log density of a response under the IMM-bsc
  real imm_bsc_core(real y, real mu, real kappa, real tau, real c, real s,
                    real b, int ss, vector nt, vector dist, int nodes,
                    data vector logk, data vector dlogk, data real logJ_min,
                    data real dlogJ) {
    vector[ss + 1] w = imm_bsc_logw(c, s, b, ss, dist);
    return circmix_vp_ld(cos(y - circmix_locations(mu, nt, ss)), head(w, ss),
                         w[ss + 1], kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // log density of a response under the full IMM
  real imm_full_core(real y, real mu, real kappa, real tau, real c, real a,
                     real s, real b, int ss, vector nt, vector dist, int nodes,
                     data vector logk, data vector dlogk, data real logJ_min,
                     data real dlogJ) {
    vector[ss + 1] w = imm_full_logw(c, a, s, b, ss, dist);
    return circmix_vp_ld(cos(y - circmix_locations(mu, nt, ss)), head(w, ss),
                         w[ss + 1], kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }
