  /* Three-parameter mixture model (Bays, Catalao & Husain, 2009) and its
   * capacity-limited versions. Each core builds the weights of the target and
   * the ss - 1 active non-targets and hands them to the shared circmix
   * functions. nt holds the non-target locations padded to max_set_size - 1;
   * only the first ss - 1 entries are read, so a trial never pays for
   * components its set size does not have.
   *
   * At set size 1 the loop is empty and no swap can occur, which is why these
   * versions need no constant prior pinning the non-target weight at set size 1.
   */

  real mixture3p_simple_core(real y, real mu, real kappa, real tau, real thetat,
                             real thetant, int ss, vector nt, int nodes,
                             data vector logk, data vector dlogk,
                             data real logJ_min, data real dlogJ) {
    int n_nt = ss - 1;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    cosd[1] = cos(y - mu);
    logw[1] = thetat;
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = thetant - log(n_nt);
    }
    // guessing is the reference category of the softmax, i.e. weight zero on
    // the log scale before normalising
    real total = log_sum_exp(append_row(logw, 0.0));
    return circmix_vp_ld(cosd, logw - total, -total, kappa, tau, nodes,
                         logk, dlogk, logJ_min, dlogJ);
  }

  /* Storage and swapping are orthogonal here: the capacity rule says how often
   * a response comes from memory, and pnt says which stored item is reported.
   * A softmax over the mixture weights cannot express this, because p_mem and
   * the swap probability are then jointly constrained and neither is separately
   * the quantity a storage model predicts.
   */
  real mixture3p_slot_core(real y, real mu, real kappa, real tau, real K,
                           real pnt, int ss, vector nt, int nodes,
                           data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    int n_nt = ss - 1;
    real p_mem = fmin(1.0, K / ss);
    real swap = n_nt > 0 ? pnt : 0.0;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    cosd[1] = cos(y - mu);
    logw[1] = log(p_mem) + log1m(swap);
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = log(p_mem) + log(swap) - log(n_nt);
    }
    return circmix_vp_ld(cosd, logw, log1m(p_mem), kappa, tau, nodes,
                         logk, dlogk, logJ_min, dlogJ);
  }

  /* Slot averaging crossed with swapping. Which item is reported and how many
   * slots it holds are independent, so the swap weights serve both slot counts
   * and only the concentration differs between them.
   */
  real mixture3p_slot_averaging_core(real y, real mu, real kappa, real tau,
                                     real K, real pnt, int ss, vector nt,
                                     int nodes, data vector logk,
                                     data vector dlogk, data real logJ_min,
                                     data real dlogJ) {
    int n_nt = ss - 1;
    real swap = n_nt > 0 ? pnt : 0.0;
    vector[n_nt + 1] cosd;
    vector[n_nt + 1] logw;
    cosd[1] = cos(y - mu);
    logw[1] = log1m(swap);
    for (j in 1:n_nt) {
      cosd[j + 1] = cos(y - nt[j]);
      logw[j + 1] = log(swap) - log(n_nt);
    }
    return circmix_slot_averaging_ld(cosd, logw, K, ss, kappa, tau, nodes,
                                     logk, dlogk, logJ_min, dlogJ);
  }
