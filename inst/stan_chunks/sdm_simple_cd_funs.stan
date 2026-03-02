  #include 'fun_tan_half.stan'

  real sdm_log_numer(real x, real mu, real c_par, real kappa) {
    return c_par * exp(kappa * (cos(x - mu) - 1)) * sqrt(kappa) * inv(sqrt2()) * inv_sqrt(pi());
  }

  real sdm_log_density(real x, real mu, real c_par, real kappa, int n_norm, real norm_dx) {
    real log_numer = sdm_log_numer(x, mu, c_par, kappa);
    real norm_sum = 0;
    for (j in 1:n_norm) {
      real t = -pi() + (j - 1) * norm_dx;
      norm_sum += exp(sdm_log_numer(t, mu, c_par, kappa));
    }
    real log_denom = log(norm_sum * norm_dx);
    return log_numer - log_denom;
  }

  real sdm_simple_cd_lpmf(int y, real mu, real c_par, real kappa, real beta, real probe) {
    int n_quad = 101;
    real dx = 2 * pi() / (n_quad - 1);
    real p_change = 0;
    real sharpness = 5;
    int n_norm = 101;
    real norm_dx = 2 * pi() / (n_norm - 1);

    for (i in 1:n_quad) {
      real x = -pi() + (i - 1) * dx;

      real log_p_retrieve = sdm_log_density(x, mu, c_par, kappa, n_norm, norm_dx);
      real log_p_x_given_same = sdm_log_density(x, probe + mu, c_par, kappa, n_norm, norm_dx);

      real llr = log_p_retrieve - log_p_x_given_same;
      real w = inv_logit(sharpness * (llr - beta));

      p_change += w * exp(log_p_retrieve) * dx;
    }

    p_change = fmin(fmax(p_change, 1e-10), 1 - 1e-10);
    if (y == 1) return log(p_change);
    return log1m(p_change);
  }
