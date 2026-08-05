  // Half-width of the spike in the sdm density. Near its peak the activation
  // behaves like eta ~ A * exp(-kappa u^2 / 2) with A the peak height, so the
  // density has settled onto its floor once A * exp(-kappa u^2 / 2) drops
  // below 1. The integrand spans up to 27 orders of magnitude inside that
  // spike, which is why a fixed quadrature rule over the whole arc fails.
  real sdm_cd_spike_width(real c, real kappa) {
    real peak = c * sqrt(kappa * inv(2 * pi()));
    return fmin(pi(), sqrt(2 * log(fmax(peak, e())) * inv(kappa)));
  }

  // log of the integral of the unnormalised sdm density over [lo, hi], with
  // the density centred on mu. The interval is split at every periodic image
  // of the peak (the arc may span the whole circle) and at multiples of the
  // spike width, so the spike always sits at a sub-interval endpoint where
  // Gauss-Legendre clusters its nodes. Terms are accumulated by log-sum-exp
  // because the integrand overflows well inside the sampled parameter range.
  real sdm_cd_log_int(real lo, real hi, real mu, real c, real kappa,
                      data vector gl_x, data vector gl_w) {
    int n_nodes = num_elements(gl_x);
    real w = sdm_cd_spike_width(c, kappa);
    real base = c * sqrt(kappa * inv(2 * pi()));
    array[9] real mult = {-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3};
    vector[29] cuts;
    vector[n_nodes] log_gl_w = log(gl_w);
    real out = negative_infinity();
    int k = 3;

    cuts[1] = lo;
    cuts[2] = hi;
    for (p in -1:1) {
      for (j in 1:9) {
        cuts[k] = fmin(fmax(mu + p * 2 * pi() + mult[j] * w, lo), hi);
        k += 1;
      }
    }
    cuts = sort_asc(cuts);

    for (i in 1:28) {
      real a = cuts[i];
      real b = cuts[i + 1];
      if (b > a) {
        real half = (b - a) / 2;
        vector[n_nodes] u = (a + b) / 2 + half * gl_x - mu;
        out = log_sum_exp(
          out,
          log_sum_exp(log(half) + log_gl_w + base * exp(kappa * (cos(u) - 1)))
        );
      }
    }

    return out;
  }

  // Half-width of the arc around the probe within which a "same" response is
  // given. The rule is the one of Lin & Oberauer (2022): respond "change" when
  // the memory density at the probe falls below the uniform, scaled by the
  // criterion. The sdm density decreases monotonically in the distance from
  // its centre, so that region is again an arc. Unlike the von Mises models
  // the boundary needs the normalising constant, and no parameter drops out of
  // it at an unbiased criterion.
  real sdm_cd_crit_angle(real c, real kappa, real criterion, real log_z) {
    real peak = c * sqrt(kappa * inv(2 * pi()));
    real thresh = log_z - criterion - log(2 * pi());
    if (thresh >= peak) return 0;
    if (thresh <= peak * exp(-2 * kappa)) return pi();
    return acos(fmin(fmax(1 + log(thresh / peak) * inv(kappa), -1), 1));
  }

  real sdm_cd_lpmf(array[] int y, vector mu, vector c, vector kappa,
                   vector criterion, data array[] real probe, data vector gl_x,
                   data vector gl_w) {
    int N = size(y);
    real out = 0;
    real log_z = 0;

    for (n in 1:N) {
      real hw;
      real p_same;

      // the normalising constant is rotation invariant, so it only has to be
      // recomputed when c or kappa change; the data are sorted by predictors
      if (n == 1 || c[n] != c[n - 1] || kappa[n] != kappa[n - 1]) {
        log_z = sdm_cd_log_int(-pi(), pi(), 0, c[n], kappa[n], gl_x, gl_w);
      }

      hw = sdm_cd_crit_angle(c[n], kappa[n], criterion[n], log_z);

      if (hw <= 0) {
        p_same = 0;
      } else if (hw >= pi()) {
        p_same = 1;
      } else {
        p_same = exp(
          sdm_cd_log_int(probe[n] - hw, probe[n] + hw, mu[n], c[n], kappa[n],
                         gl_x, gl_w) - log_z
        );
      }

      out += cd_bernoulli_lpmf(y[n] | p_same);
    }

    return out;
  }
