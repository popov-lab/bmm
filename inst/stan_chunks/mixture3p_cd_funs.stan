  #include 'fun_tan_half.stan'

  // Change detection likelihood for mixture3p model
  // Based on Lin & Oberauer (2022, Cognitive Psychology)
  //
  // Extends mixture2p_cd by adding non-target von Mises components
  // to the retrieval distribution.
  //
  // mu is included because brms always passes it as the first dpar.
  // nt_features and lure_idx are passed via vreal/vint.
  //
  // Args layout (brms order: dpars then vars):
  //   mu       - unused (void_mu = TRUE)
  //   kappa    - concentration parameter
  //   thetat   - target mixture weight (log-odds scale via identity link)
  //   thetant  - non-target mixture weight (log-odds scale)
  //   beta     - decision criterion
  //   vreal1   - probe_centered
  //   vreal2..vreal(K+1) - nt_features[1..K] (K = max_set_size - 1)
  //   vint1..vintK       - lure_idx[1..K]

  // Note: This function handles a FIXED max_set_size. A separate Stan file
  // is generated per max_set_size value, or we use the maximum and mask
  // inactive non-targets via lure_idx.
