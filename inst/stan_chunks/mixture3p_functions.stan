  // mixture3p: normalised log weights of target, ss - 1 non-targets and guessing
  vector mixture3p_simple_logw(real thetat, real thetant, int ss) {
    int n_nt = ss - 1;
    vector[ss + 1] logw;
    logw[1] = thetat;
    for (j in 1:n_nt) {
      logw[j + 1] = thetant - log(n_nt);
    }
    logw[ss + 1] = 0;
    return logw - log_sum_exp(logw);
  }

  // log weights with K memory slots; ptarget is the probability that a memory response reports the target
  vector mixture3p_slot_logw(real K, real ptarget, int ss) {
    int n_nt = ss - 1;
    real p_mem = fmin(1.0, K / ss);
    vector[ss + 1] logw;
    logw[1] = log(p_mem) + (n_nt > 0 ? log(ptarget) : 0.0);
    for (j in 1:n_nt) {
      logw[j + 1] = log(p_mem) + log1m(ptarget) - log(n_nt);
    }
    logw[ss + 1] = log1m(p_mem);
    return logw;
  }

  // log weights of target and non-targets among memory responses, without guessing
  vector mixture3p_slot_averaging_logw(real ptarget, int ss) {
    int n_nt = ss - 1;
    vector[ss] logw;
    logw[1] = n_nt > 0 ? log(ptarget) : 0.0;
    for (j in 1:n_nt) {
      logw[j + 1] = log1m(ptarget) - log(n_nt);
    }
    return logw;
  }

  // log density of a response under the three-parameter mixture model
  real mixture3p_simple_core(real y, real mu, real kappa, real tau, real thetat,
                             real thetant, int ss, vector nt, int nodes,
                             data vector logk, data vector dlogk,
                             data real logJ_min, data real dlogJ) {
    vector[ss + 1] w = mixture3p_simple_logw(thetat, thetant, ss);
    return circmix_vp_ld(cos(y - circmix_locations(mu, nt, ss)), head(w, ss),
                         w[ss + 1], kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // log density of a response under the slot model with non-target responses
  real mixture3p_slot_core(real y, real mu, real kappa, real tau, real K,
                           real ptarget, int ss, vector nt, int nodes,
                           data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    vector[ss + 1] w = mixture3p_slot_logw(K, ptarget, ss);
    return circmix_vp_ld(cos(y - circmix_locations(mu, nt, ss)), head(w, ss),
                         w[ss + 1], kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // log density of a response under the slots-plus-averaging model with non-target responses
  real mixture3p_slot_averaging_core(real y, real mu, real kappa, real tau,
                                     real K, real ptarget, int ss, vector nt,
                                     int nodes, data vector logk,
                                     data vector dlogk, data real logJ_min,
                                     data real dlogJ) {
    return circmix_slot_averaging_ld(cos(y - circmix_locations(mu, nt, ss)),
                                     mixture3p_slot_averaging_logw(ptarget, ss),
                                     K, ss, kappa, tau, nodes,
                                     logk, dlogk, logJ_min, dlogJ);
  }
