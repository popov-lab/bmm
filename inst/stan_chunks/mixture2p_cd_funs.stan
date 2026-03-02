  #include 'fun_tan_half.stan'

  // Change detection likelihood for mixture2p model
  // Based on Lin & Oberauer (2022, Cognitive Psychology)
  //
  // P(change) = integral sigma(LLR(x) - beta) * P_retrieve(x) dx
  // P_retrieve(x) = thetat * vM(x | 0, kappa) + (1 - thetat) / (2*pi)
  // LLR(x) = log P_retrieve(x | target) - log P_retrieve(x | probe)
  //
  // The indicator I(LLR > beta) is smoothed with a logistic sigmoid to ensure
  // differentiability for HMC sampling.

  real mixture2p_cd_lpmf(int y, real mu, real kappa, real thetat, real beta, real probe) {
    int n_quad = 101;
    real dx = 2 * pi() / (n_quad - 1);
    real p_change = 0;
    real log_uniform = -log(2 * pi());
    real sharpness = 5;

    // precompute von Mises normalization -- depends only on kappa
    real log_vm_norm = log(2 * pi()) + log_modified_bessel_first_kind(0, kappa);

    for (i in 1:n_quad) {
      real x = -pi() + (i - 1) * dx;

      real vm_ret = kappa * cos(x - mu) - log_vm_norm;
      real vm_same = kappa * cos(x - probe - mu) - log_vm_norm;

      real log_p_retrieve = log_mix(thetat, vm_ret, log_uniform);

      real log_p_x_given_same = log_mix(thetat, vm_same, log_uniform);
      real llr = log_p_retrieve - log_p_x_given_same;
      real w = inv_logit(sharpness * (llr - beta));

      p_change += w * exp(log_p_retrieve) * dx;
    }

    p_change = fmin(fmax(p_change, 1e-10), 1 - 1e-10);
    if (y == 1) return log(p_change);
    return log1m(p_change);
  }
