  real sdm_abc_kernel(real y, real mu, real kappa) {
    return exp(kappa * (cos(y - mu) - 1)) * sqrt(kappa) * inv_sqrt(2 * pi());
  }

  real sdm_abc_activation_obs(real y, int n, real mu, real c, real a, real kappa,
                              matrix nt_features, matrix lure_idx) {
    int K = cols(nt_features);
    real exp_c = exp(c);
    real exp_a = exp(a);
    real out = (exp_c + exp_a) * sdm_abc_kernel(y, mu, kappa);

    for (k in 1:K) {
      if (lure_idx[n, k] > 0.5) {
        out += exp_a * sdm_abc_kernel(y, mu + nt_features[n, k], kappa);
      }
    }

    return out;
  }

  real sdm_abc_logz_obs(int n, real mu, real c, real a, real kappa, vector grid, real grid_step,
                        matrix nt_features, matrix lure_idx) {
    int M = num_elements(grid);
    vector[M] act;

    for (m in 1:M) {
      act[m] = log(sdm_abc_activation_obs(grid[m], n, mu, c, a, kappa, nt_features, lure_idx));
    }

    return log_sum_exp(act) + log(grid_step);
  }
