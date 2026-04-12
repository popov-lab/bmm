  real mixture2p_cd_lpmf(array[] int y, vector mu, vector kappa,
                         vector thetat, vector beta) {
    return 0;
  }

  real mixture2p_cd_log_prob(int y, real probe, real mu, real kappa,
                             real thetat, real beta, vector grid, real dx) {
    int n_quad = size(grid);
    real log_uniform = -log(2 * pi());
    vector[n_quad] log_p_ret;
    vector[n_quad] log_p_same;
    vector[n_quad] log_integrand;
    real log_vm_norm = log(2 * pi()) + log_modified_bessel_first_kind(0, kappa);

    for (i in 1:n_quad) {
      real x = grid[i];

      real vm_ret = kappa * cos(x - mu) - log_vm_norm;
      real vm_same = kappa * cos(x - probe - mu) - log_vm_norm;

      log_p_ret[i] = log_mix(thetat, vm_ret, log_uniform);
      log_p_same[i] = log_mix(thetat, vm_same, log_uniform);
      log_integrand[i] = log_inv_logit(5 * (log_p_ret[i] - log_p_same[i] - beta)) +
        log_p_ret[i];
    }

    {
      real log_p_change = log_sum_exp(log_integrand) + log(dx);
      log_p_change = fmin(fmax(log_p_change, log(1e-10)), log1m(1e-10));
      if (y == 1) return log_p_change;
      return log1m_exp(log_p_change);
    }
  }
