  #include 'fun_tan_half.stan'

  // The activation and normaliser of the sdm() cells other than the default
  // (#558), which keeps sdm_simple_funs.stan. The density is
  // exp(c* s(y - mu)) / Z with s(0) = 1, so c* is the peak of the activation,
  // and c is on the log scale. activation: 1 von Mises, 2 Laplace; scaling:
  // 1 density, 2 peak, 3 range. "density" divides c by the integral of s over
  // the circle (sqrt(2 pi / kappa), the normal approximation, for von Mises;
  // exactly 2 (1 - exp(-kappa pi)) / kappa for Laplace), "peak" takes c as the
  // activation of the target, and "range" as the contrast c* (1 - s(pi)),
  // whose additive constant cancels against Z and is dropped
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

  // log s(delta); cos(delta) - 1 as -2 sin(delta / 2)^2, which keeps its
  // digits near the target. y and mu lie in [-pi, pi], so |delta| < 2 pi
  vector sdm_log_s(vector delta, vector kappa, int activation) {
    if (activation == 1) {
      return -2 * kappa .* square(sin(0.5 * delta));
    }
    vector[rows(delta)] d = abs(delta);
    return -kappa .* fmin(d, 2 * pi() - d);
  }

  // The body of sdm_simple_ldenom() (sdm_simple_funs.stan, where the reasons
  // for each branch are written out) taking the log peak instead of c;
  // test-model_sdm-stan.R pins the two to agree
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

  // Minus the log of the mean of exp(peak exp(-kappa u)) over the circle,
  // from Gauss-Legendre panels on [0, pi]; the kink at u = 0 is an end, so
  // the integrand is smooth on each panel. In v = kappa u, exp(eta) falls by
  // e per 1 / peak near the top, which the cuts at 3 and 30 / peak resolve,
  // and past v0 = log peak, where the activation is 1, the integrand over its
  // floor is exp(exp(-(v - v0))) whatever the peak, which the cuts at v0, v0 +
  // 6 and v0 + 30 resolve. Within 1.7e-12 of the exponential-integral series
  // to log peak 15 and kappa 3000 (#558)
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
