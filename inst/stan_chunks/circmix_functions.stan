  // 1 - I1(kappa) / I0(kappa) for large kappa
  real circmix_one_minus_A(real kappa) {
    real u = inv(kappa);
    return u * (0.5 + u * (0.125 + u * (0.125 + u * (25.0 / 128 + u * 13.0 / 32))));
  }

  // Fisher information J of a von Mises distribution with concentration kappa
  real circmix_J(real kappa) {
    if (kappa <= 0) {
      return 0;
    }
    if (kappa > 1e3) {
      return kappa * (1 - circmix_one_minus_A(kappa));
    }
    return kappa * exp(log_modified_bessel_first_kind(1, kappa)
                       - log_modified_bessel_first_kind(0, kappa));
  }

  // index of the grid cell that contains x on an evenly spaced grid
  int circmix_cell(real x, data real x_min, data real dx, int n) {
    int lo = 1;
    int hi = n;
    while (hi - lo > 1) {
      int mid = (lo + hi) %/% 2;
      if (x_min + (mid - 1) * dx <= x) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return lo;
  }

  // concentration kappa for a Fisher information J (inverse of circmix_J), from a lookup table
  real circmix_kappa(real J, data vector logk, data vector dlogk,
                     data real logJ_min, data real dlogJ) {
    int n = num_elements(logk);
    real t = log(J);
    if (t <= logJ_min) {
      return sqrt(2 * J) * (1 + J / 8);
    }
    if (t >= logJ_min + (n - 1) * dlogJ) {
      return J + 0.5;
    }
    int i = circmix_cell(t, logJ_min, dlogJ, n);
    real s = (t - logJ_min) / dlogJ - (i - 1);
    real s2 = s * s;
    real s3 = s2 * s;
    return exp((2 * s3 - 3 * s2 + 1) * logk[i]
               + (-2 * s3 + 3 * s2) * logk[i + 1]
               + dlogJ * ((s3 - 2 * s2 + s) * dlogk[i]
                          + (s3 - s2) * dlogk[i + 1]));
  }

  // log(I0(kappa) * exp(-kappa)), the scaled modified Bessel function of order 0
  real circmix_log_besselI0_scaled(real kappa) {
    if (kappa < 2000) {
      return log_modified_bessel_first_kind(0, kappa) - kappa;
    }
    return -0.5 * (log(2 * pi()) + log(kappa))
           + log1p(inv(8 * kappa) + 9 / (128 * square(kappa)) + 225 / (3072 * kappa^3));
  }

  // log density of a von Mises mixture with one kappa plus uniform guessing; cosd = cos(y - mu)
  real circmix_ld(vector cosd, vector logw, real logw_guess, real kappa) {
    return log_sum_exp(log_sum_exp(logw + kappa * (cosd - 1))
                         - circmix_log_besselI0_scaled(kappa),
                       logw_guess)
           - log(2 * pi());
  }

  // the same mixture with a separate kappa for each component
  real circmix_het_ld(vector cosd, vector logw, real logw_guess, vector kappa) {
    int k = num_elements(cosd);
    vector[k] lp;
    for (j in 1:k) {
      lp[j] = logw[j] + kappa[j] * (cosd[j] - 1)
              - circmix_log_besselI0_scaled(kappa[j]);
    }
    return log_sum_exp(log_sum_exp(lp), logw_guess) - log(2 * pi());
  }

  // log of the Simpson's rule weight of node i out of n
  real circmix_simpson_lw(int i, int n) {
    if (i == 1 || i == n) {
      return 0;
    }
    return i % 2 == 0 ? log(4) : log(2);
  }

  // 1 if memory precision is so low that responses from memory are uniform
  int circmix_vp_uniform(real J, real tau) {
    return J < 1e-17 * sqrt(tau);
  }

  // shape of the gamma distribution of precision; infinite when precision need not be integrated out
  real circmix_vp_shape(real J, real tau) {
    if (tau <= 0 || circmix_vp_uniform(J, tau)) {
      return positive_infinity();
    }
    real shape = J / tau;
    if (8 * sqrt(trigamma(shape)) < 1e-6) {
      return positive_infinity();
    }
    return shape;
  }

  // log normalising term of the gamma density in circmix_vp_log_gamma
  real circmix_gamma_log_norm(real shape) {
    if (shape < 10) {
      return shape * log(shape) - shape - lgamma(shape);
    }
    real u = inv(shape);
    real u2 = square(u);
    return 0.5 * (log(shape) - log(2 * pi()))
           - u * (1.0 / 12 - u2 * (1.0 / 360 - u2 * (1.0 / 1260 - u2 * (1.0 / 1680 - u2 / 1188))));
  }

  // log density of t = log(J) when J is gamma distributed with the given shape and scale tau
  real circmix_vp_log_gamma(real t, real shape, real tau) {
    real delta = t - log(shape) - log(tau);
    if (delta > 700) {
      return negative_infinity();
    }
    return circmix_gamma_log_norm(shape) - shape * (expm1(delta) - delta);
  }

  // first and second derivative, in log(J), of the log integrand in circmix_vp_component
  vector circmix_vp_slopes(real t, real c, real shape, real tau,
                           data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    real J = exp(t);
    real kappa = circmix_kappa(J, logk, dlogk, logJ_min, dlogJ);
    real A;
    real one_minus_A;
    real curvature;
    if (kappa > 1e3) {
      real u = inv(kappa);
      one_minus_A = circmix_one_minus_A(kappa);
      A = 1 - one_minus_A;
      curvature = (0.5 + u * (0.25 + u * (0.375 + u * (25.0 / 32 + u * 65.0 / 32))))
                  / square(kappa * one_minus_A);
    } else {
      A = exp(log_modified_bessel_first_kind(1, kappa)
              - log_modified_bessel_first_kind(0, kappa));
      one_minus_A = 1 - A;
      curvature = (1 - A / kappa - square(A)) / square(one_minus_A);
    }
    real drift = ((c - 1) + one_minus_A) * A / (one_minus_A * (1 + A));
    return [shape - J / tau + drift,
            -J / tau - curvature * square(A / (1 + A))
            + drift * curvature * (1 + square(A)) / square(1 + A)]';
  }

  // circmix_vp_component for gamma shapes below 2
  real circmix_vp_component_tail(real c, real shape, real tau, int nodes,
                                 data vector logk, data vector dlogk,
                                 data real logJ_min, data real dlogJ) {
    real J0 = 1e-4;
    real rate = inv(tau);
    real tail = exp(gamma_lcdf(J0 | shape, rate))
                + c * sqrt(2) * sqrt(tau) * exp(lgamma(shape + 0.5) - lgamma(shape)) * exp(gamma_lcdf(J0 | shape + 0.5, rate))
                + (square(c) - 0.5) * shape * tau * exp(gamma_lcdf(J0 | shape + 1, rate));
    real lo = log(J0);
    real hi = fmax(fmin(log(tau) + log(fmax(shape * (1 - inv(9 * shape) + 8.5 / (3 * sqrt(shape)))^3, 40)), 709),
                   lo + 1);
    real step = (hi - lo) / (nodes - 1);
    vector[nodes] lp;
    for (i in 1:nodes) {
      real ti = lo + (i - 1) * step;
      real J = exp(ti);
      real kappa = circmix_kappa(J, logk, dlogk, logJ_min, dlogJ);
      lp[i] = circmix_simpson_lw(i, nodes) + circmix_vp_log_gamma(ti, shape, tau)
              + kappa * (c - 1) - circmix_log_besselI0_scaled(kappa);
    }
    return log_sum_exp(log_sum_exp(lp) + log(step / 3), log(tail)) - log(2 * pi());
  }

  // log density of one von Mises component with gamma-distributed precision; c = cos(y - mu)
  real circmix_vp_component(real c, real shape, real tau, int nodes,
                            data vector logk, data vector dlogk,
                            data real logJ_min, data real dlogJ) {
    if (shape < 2) {
      return circmix_vp_component_tail(c, shape, tau, nodes, logk, dlogk, logJ_min, dlogJ);
    }
    real rate = inv(tau);
    real t = log((shape + 0.5) / (rate + (1 - c)));
    for (i in 1:4) {
      vector[2] d = circmix_vp_slopes(t, c, shape, tau, logk, dlogk, logJ_min, dlogJ);
      real newton;
      if (d[2] < 0) {
        newton = -d[1] / d[2];
      } else {
        newton = d[1] > 0 ? 1 : -1;
      }
      t += fmin(fmax(newton, -1), 1);
    }
    real d2 = circmix_vp_slopes(t, c, shape, tau, logk, dlogk, logJ_min, dlogJ)[2];
    real scale = sqrt(trigamma(shape));
    if (d2 < 0) {
      scale = fmax(scale, inv_sqrt(-d2));
    }
    real step = 18 * scale / (nodes - 1);
    real lo = t - 12 * scale;
    vector[nodes] lp;
    for (i in 1:nodes) {
      real ti = lo + (i - 1) * step;
      real J = exp(ti);
      real kappa = circmix_kappa(J, logk, dlogk, logJ_min, dlogJ);
      lp[i] = circmix_simpson_lw(i, nodes) + circmix_vp_log_gamma(ti, shape, tau)
              + kappa * (c - 1) - circmix_log_besselI0_scaled(kappa);
    }
    return log_sum_exp(lp) + log(step / 3) - log(2 * pi());
  }

  // mixture log density with variable precision; tau = 0 means constant precision
  real circmix_vp_ld(vector cosd, vector logw, real logw_guess,
                     real kappa, real tau, int nodes,
                     data vector logk, data vector dlogk,
                     data real logJ_min, data real dlogJ) {
    if (tau <= 0) {
      return circmix_ld(cosd, logw, logw_guess, kappa);
    }
    real J = circmix_J(kappa);
    if (tau == positive_infinity() || circmix_vp_uniform(J, tau)) {
      return circmix_ld(cosd, logw, logw_guess, 0);
    }
    real shape = circmix_vp_shape(J, tau);
    if (is_inf(shape)) {
      return circmix_ld(cosd, logw, logw_guess, kappa);
    }
    int k = num_elements(cosd);
    vector[k] lp = rep_vector(negative_infinity(), k);
    for (j in 1:k) {
      if (logw[j] > negative_infinity()) {
        lp[j] = logw[j] + circmix_vp_component(cosd[j], shape, tau, nodes,
                                               logk, dlogk, logJ_min, dlogJ);
      }
    }
    return log_sum_exp(log_sum_exp(lp), logw_guess - log(2 * pi()));
  }

  // slots per item for capacity K and set size ss: [whole slots, probability of one more slot]
  vector circmix_slots(real K, int ss) {
    real q = K * inv(ss);
    real f = floor(q);
    return [f, q - f]';
  }

  // slot averaging: [any slot held, kappa with fewer slots, kappa with one more slot, its probability]
  vector circmix_slot_averaging_branches(real K, int ss, real kappa,
                                         data vector logk, data vector dlogk,
                                         data real logJ_min, data real dlogJ) {
    vector[2] allocation = circmix_slots(K, ss);
    real slots = allocation[1];
    real J1 = circmix_J(kappa);
    int held = slots >= 0.5;
    real kappa_lo = 0;
    if (held) {
      kappa_lo = circmix_kappa(slots * J1, logk, dlogk, logJ_min, dlogJ);
    }
    real kappa_hi = circmix_kappa((slots + 1) * J1, logk, dlogk, logJ_min, dlogJ);
    return [held, kappa_lo, kappa_hi, allocation[2]]';
  }

  // mixture log density under slot averaging (Zhang & Luck, 2008): precision grows with the slots held
  real circmix_slot_averaging_ld(vector cosd, vector logw, real K, int ss,
                                 real kappa, real tau, int nodes,
                                 data vector logk, data vector dlogk,
                                 data real logJ_min, data real dlogJ) {
    vector[4] branch = circmix_slot_averaging_branches(K, ss, kappa, logk, dlogk,
                                                       logJ_min, dlogJ);
    real extra = branch[4];
    real ld_hi = circmix_vp_ld(cosd, logw, negative_infinity(), branch[3], tau,
                               nodes, logk, dlogk, logJ_min, dlogJ);
    real ld_lo = -log(2 * pi());
    if (branch[1] == 1) {
      ld_lo = circmix_vp_ld(cosd, logw, negative_infinity(), branch[2], tau,
                            nodes, logk, dlogk, logJ_min, dlogJ);
    }
    if (extra > 0) {
      return log_sum_exp(log1m(extra) + ld_lo, log(extra) + ld_hi);
    }
    return ld_lo + extra * expm1(ld_hi - ld_lo);
  }

  // locations of the target followed by the ss - 1 non-targets of a trial
  vector circmix_locations(real mu, vector nt, int ss) {
    return append_row(mu, head(nt, ss - 1));
  }
