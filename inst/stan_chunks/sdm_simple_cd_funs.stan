  #include 'fun_tan_half.stan'

  // dummy _lpmf — brms requires this but actual likelihood is in the likelihood block
  real sdm_simple_cd_lpmf(array[] int y, vector mu, vector c_par,
                           vector kappa, vector beta) {
    return 0;
  }

  // compute log-normalization constant via quadrature (cached in likelihood block)
  real sdm_cd_log_Z(real c_par, real kappa, vector grid) {
    int n_quad = size(grid);
    real base = c_par * sqrt(kappa) * inv(sqrt2()) * inv_sqrt(pi());
    real norm_sum = 0;
    for (j in 1:n_quad) {
      norm_sum += exp(base * exp(kappa * (cos(grid[j]) - 1)));
    }
    return log(norm_sum);
  }

  // compute log P(y | params) for a single observation, given cached log_Z
  real sdm_cd_log_prob(int y, real mu, real c_par, real kappa, real beta,
                       real probe, real log_Z, vector grid) {
    int n_quad = size(grid);
    real base = c_par * sqrt(kappa) * inv(sqrt2()) * inv_sqrt(pi());
    real sharpness = 5;
    real weighted_sum = 0;

    for (i in 1:n_quad) {
      real x = grid[i];
      real ln_ret = base * exp(kappa * (cos(x - mu) - 1));
      real ln_same = base * exp(kappa * (cos(x - probe - mu) - 1));

      real exp_ln_ret = exp(ln_ret);
      real llr = ln_ret - ln_same;
      real w = inv_logit(sharpness * (llr - beta));
      weighted_sum += w * exp_ln_ret;
    }

    real p_change = weighted_sum / exp(log_Z);
    p_change = fmin(fmax(p_change, 1e-10), 1 - 1e-10);
    if (y == 1) return log(p_change);
    return log1m(p_change);
  }
