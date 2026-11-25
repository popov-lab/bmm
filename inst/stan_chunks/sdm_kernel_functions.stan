  #include 'fun_tan_half.stan'

  // log peak activation; activation: 1 = von Mises, 2 = Laplace; scaling: 1 = density, 2 = peak, 3 = range
  real sdm_log_cstar(real c, real kappa, int activation, int scaling) {
    if (scaling == 2) {
      return c;
    }
    if (scaling == 1) {
      return activation == 1 ? c + 0.5 * log(kappa / (2 * pi()))
                             : c + log(kappa) - log2() - log1m_exp(-kappa * pi());
    }
    return c - log1m_exp(activation == 1 ? -2 * kappa : -kappa * pi());
  }

  // vectorized version of sdm_log_cstar
  vector sdm_log_cstar(vector c, vector kappa, int activation, int scaling) {
    if (scaling == 2) {
      return c;
    }
    if (scaling == 1) {
      return activation == 1 ? c + 0.5 * log(kappa / (2 * pi()))
                             : c + log(kappa) - log2() - log1m_exp(-kappa * pi());
    }
    return c - log1m_exp(activation == 1 ? -2 * kappa : -kappa * pi());
  }

  // log of the activation at angular distance delta from the target, relative to its peak
  vector sdm_log_s(vector delta, vector kappa, int activation) {
    if (activation == 1) {
      return -2 * kappa .* square(sin(0.5 * delta));
    }
    vector[rows(delta)] d = abs(delta);
    return -kappa .* fmin(d, 2 * pi() - d);
  }

  // normalising term of the SDM log likelihood with von Mises activation, given the log peak
  real sdm_vm_ldenom_peak(real log_peak, real kappa, data vector gl_x, data vector gl_w,
                          data vector wide_cosm1, data vector wide_w) {
    int n_nodes = num_elements(gl_x);
    real peak = exp(log_peak);
    if (kappa * (peak + 1) <= 4000) {
      return log(pi()) - peak - log(dot_product(wide_w, exp(peak * expm1(kappa * wide_cosm1))));
    }
    real w = fmin(pi(), sqrt(2 * fmax(log_peak, 1) / kappa));
    real s = exp(-0.5 * (log_peak + log(kappa)));
    vector[9] cuts = [0, pi(), 0.5 * w, w, 2 * w, 5 * w, s, 3 * s, 9 * s]';
    real log_int = negative_infinity();
    for (i in 1:9) {
      cuts[i] = fmin(fmax(cuts[i], 0), pi());
    }
    cuts = sort_asc(cuts);
    for (i in 1:8) {
      if (cuts[i + 1] > cuts[i]) {
        real half = (cuts[i + 1] - cuts[i]) / 2;
        vector[n_nodes] u = (cuts[i] + cuts[i + 1]) / 2 + half * gl_x;
        vector[n_nodes] d = peak * expm1(-2 * kappa * square(sin(u / 2)));
        real d_max = max(d);
        log_int = log_sum_exp(log_int, log(half) + d_max + log(dot_product(gl_w, exp(d - d_max))));
      }
    }
    return log(pi()) - peak - log_int;
  }

  // normalising term of the SDM log likelihood with Laplace activation, given the log peak
  real sdm_laplace_ldenom(real log_peak, real kappa, data vector gl_x, data vector gl_w) {
    int n_nodes = num_elements(gl_x);
    real peak = exp(log_peak);
    real v0 = fmax(log_peak, 1);
    vector[7] cuts = append_row([0, pi()]', [3 / peak, 30 / peak, v0, v0 + 6, v0 + 30]' / kappa);
    real log_int = negative_infinity();
    for (i in 1:7) {
      cuts[i] = fmin(fmax(cuts[i], 0), pi());
    }
    cuts = sort_asc(cuts);
    for (i in 1:6) {
      if (cuts[i + 1] > cuts[i]) {
        real half = (cuts[i + 1] - cuts[i]) / 2;
        vector[n_nodes] u = (cuts[i] + cuts[i + 1]) / 2 + half * gl_x;
        vector[n_nodes] d = peak * expm1(-kappa * u);
        real d_max = max(d);
        log_int = log_sum_exp(log_int, log(half) + d_max + log(dot_product(gl_w, exp(d - d_max))));
      }
    }
    return log(pi()) - peak - log_int;
  }
