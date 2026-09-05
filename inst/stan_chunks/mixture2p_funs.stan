  /* Two-parameter mixture model (Zhang & Luck, 2008) and its capacity-limited
   * versions. Each core turns its version's storage rule into mixture weights
   * and hands them to the shared circmix functions.
   *
   * kappa, thetat and K arrive on their natural scale, because brms applies the
   * inverse link before calling the likelihood. tau = 0 is what the versions
   * without variable precision pass, and collapses the quadrature in
   * circmix_vp_ld() to the constant-precision density.
   */

  real mixture2p_simple_core(real y, real mu, real kappa, real tau, real thetat,
                             int nodes, data vector logk, data vector dlogk,
                             data real logJ_min, data real dlogJ) {
    return circmix_vp_ld([cos(y - mu)]', [log(thetat)]', log1m(thetat),
                         kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  /* Fixed-resolution slots (Zhang & Luck, 2008): an item is either held at full
   * precision or not held at all, so p_mem = min(1, K / ss) and precision does
   * not vary with set size. Below capacity nothing is guessed, and that boundary
   * value of exactly one is what a softmax over mixture weights cannot reach.
   */
  real mixture2p_slot_core(real y, real mu, real kappa, real tau, real K, int ss,
                           int nodes, data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    real p_mem = fmin(1.0, K / ss);
    return circmix_vp_ld([cos(y - mu)]', [log(p_mem)]', log1m(p_mem),
                         kappa, tau, nodes, logk, dlogk, logJ_min, dlogJ);
  }

  // Slots and averaging (Zhang & Luck, 2008): the single item is held with
  // whatever slots it receives, so the item weight given holding is one
  real mixture2p_slot_averaging_core(real y, real mu, real kappa, real tau,
                                     real K, int ss, int nodes, data vector logk,
                                     data vector dlogk, data real logJ_min,
                                     data real dlogJ) {
    return circmix_slot_averaging_ld([cos(y - mu)]', [0.0]', K, ss, kappa, tau,
                                     nodes, logk, dlogk, logJ_min, dlogJ);
  }
