  // dummy _lpmf — brms requires this but actual likelihood is in the likelihood block
  real sdm_simple_cd_lpmf(array[] int y, vector mu, vector c_par,
                           vector kappa, vector beta) {
    return 0;
  }

  // compute log-normalization constant via quadrature (cached in likelihood block)
  real sdm_cd_log_Z(real c_par, real kappa, vector grid) {
    int n_quad = size(grid);
    real base = c_par * sqrt(kappa) * inv(sqrt2()) * inv_sqrt(pi());
    vector[n_quad] eta;
    for (j in 1:n_quad) {
      eta[j] = base * exp(kappa * (cos(grid[j]) - 1));
    }
    return log_sum_exp(eta);
  }

  // compute log P(y | params) for a single observation, given cached log_Z
  real sdm_cd_log_prob(int y, real mu, real c_par, real kappa, real beta,
                       real probe, real log_Z, vector grid) {
    int n_quad = size(grid);
    real base = c_par * sqrt(kappa) * inv(sqrt2()) * inv_sqrt(pi());
    vector[n_quad] log_integrand;

    for (i in 1:n_quad) {
      real x = grid[i];
      real eta_ret = base * exp(kappa * (cos(x - mu) - 1));
      real eta_same = base * exp(kappa * (cos(x - probe - mu) - 1));
      log_integrand[i] = log_inv_logit(5 * (eta_ret - eta_same - beta)) + eta_ret;
    }

    {
      real log_p_change = log_sum_exp(log_integrand) - log_Z;
      log_p_change = fmin(fmax(log_p_change, log(1e-10)), log1m(1e-10));
      if (y == 1) return log_p_change;
      return log1m_exp(log_p_change);
    }
  }
