  // Shared numerics for mixture2p, mixture3p and imm; derivations sit by the R twins in R/helpers-circmix.R

  // Fisher information of a von Mises about its location
  real circmix_J(real kappa) {
    if (kappa <= 0) {
      return 0;
    }
    return kappa * exp(log_modified_bessel_first_kind(1, kappa)
                       - log_modified_bessel_first_kind(0, kappa));
  }

  // index of the grid cell containing x, for a grid uniform in x
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

  // inverse of circmix_J: cubic Hermite on the log J table from .circmix_kappa_table(), exact asymptotics outside it
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

  // log(I0(kappa) e^-kappa), as .circmix_log_besselI0_scaled(): log I0(kappa) - kappa loses up to 2.2e-13 below kappa 2000 and 1.6e-12 below 1e4, the series 7.1e-15 from 2000 on; log(2 pi kappa) would overflow above kappa 2.9e307
  real circmix_log_besselI0_scaled(real kappa) {
    if (kappa < 2000) {
      return log_modified_bessel_first_kind(0, kappa) - kappa;
    }
    return -0.5 * (log(2 * pi()) + log(kappa))
           + log1p(inv(8 * kappa) + 9 / (128 * square(kappa)) + 225 / (3072 * kappa^3));
  }

  // von Mises mixture sharing one kappa plus uniform guessing; normalised log weights, cosd = cos(y - mu_k)
  real circmix_ld(vector cosd, vector logw, real logw_guess, real kappa) {
    return log_sum_exp(log_sum_exp(logw + kappa * cosd)
                         - log_modified_bessel_first_kind(0, kappa),
                       logw_guess)
           - log(2 * pi());
  }

  // as circmix_ld, but every component carries its own concentration
  real circmix_het_ld(vector cosd, vector logw, real logw_guess, vector kappa) {
    int k = num_elements(cosd);
    vector[k] lp;
    for (j in 1:k) {
      lp[j] = logw[j] + kappa[j] * cosd[j]
              - log_modified_bessel_first_kind(0, kappa[j]);
    }
    return log_sum_exp(log_sum_exp(lp), logw_guess) - log(2 * pi());
  }

  // composite Simpson log weight for node i of n, n odd
  real circmix_simpson_lw(int i, int n) {
    if (i == 1 || i == n) {
      return 0;
    }
    return i % 2 == 0 ? log(4) : log(2);
  }

  // 1 when memory is uniform to double precision, J(kappa) / sqrt(tau) < 1e-17; see .circmix_vp_uniform()
  int circmix_vp_uniform(real J, real tau) {
    return J < 1e-17 * sqrt(tau);
  }

  // gamma shape J(kappa) / tau of the precision distribution; +Inf when no grid is needed: constant precision (tau = 0 or a grid narrower than 1e-6) or uniform memory (circmix_vp_uniform()), as .circmix_vp_rows(); any other positive shape is integrable, see circmix_vp_component()
  real circmix_vp_shape(real kappa, real tau) {
    if (tau <= 0) {
      return positive_infinity();
    }
    real J = circmix_J(kappa);
    if (circmix_vp_uniform(J, tau)) {
      return positive_infinity();
    }
    real shape = J / tau;
    if (8 * sqrt(trigamma(shape)) < 1e-6) {
      return positive_infinity();
    }
    return shape;
  }

  // first and second derivative in t = log J of the log integrand in circmix_vp_component(), as .circmix_vp_slopes()
  vector circmix_vp_slopes(real t, real c, real shape, real tau,
                           data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ) {
    real J = exp(t);
    real kappa = circmix_kappa(J, logk, dlogk, logJ_min, dlogJ);
    real A;
    real one_minus_A;
    real dA;
    if (kappa > 1e3) {
      one_minus_A = inv(2 * kappa) + inv(8 * square(kappa)) + inv(8 * kappa^3);
      A = 1 - one_minus_A;
      dA = inv(2 * square(kappa)) + inv(4 * kappa^3) + 3 * inv(8 * kappa^4);
    } else {
      A = exp(log_modified_bessel_first_kind(1, kappa)
              - log_modified_bessel_first_kind(0, kappa));
      one_minus_A = 1 - A;
      dA = 1 - A / kappa - square(A);
    }
    real one_minus_A2 = one_minus_A * (1 + A);
    real dkappa = A / one_minus_A2;
    return [shape - J / tau + (c - A) * dkappa,
            -J / tau - dA * square(dkappa)
            + (c - A) * dA * (1 + square(A)) * dkappa / square(one_minus_A2)]';
  }

  // circmix_vp_component() below a gamma shape of 2, where sd(log J) ~ 1 / shape outgrows any fixed grid; see .circmix_vp_component_tail()
  // Below J0 = 1e-4, kappa = sqrt(2 J) (1 + J / 8) and 2 pi times the von Mises density is
  // 1 + sqrt(2 J) c + J (c^2 - 1/2) + O(J^(3/2)), so that part of the integral is three incomplete
  // gamma moments, accurate to J0^(3/2) = 1e-6 of its mass, and positive since the c terms move it
  // by at most sqrt(2 J0) + J0 / 2. The Simpson grid covers [J0, J_hi] in log J, J_hi the
  // Wilson-Hilferty quantile at z = 8.5, at least 40 tau; a J_hi below J0 (tau < 2.5e-6, under
  // 1e-15 of the mass above J0) spans [J0, e J0] instead. J_hi is capped at e^709 and tau enters as
  // log(tau) and sqrt(tau), because 40 tau and 2 tau overflow past tau 4.5e306 and 9e307.
  real circmix_vp_component_tail(real c, real shape, real tau, int nodes,
                                 data vector logk, data vector dlogk,
                                 data real logJ_min, data real dlogJ) {
    real J0 = 1e-4;
    real rate = inv(tau);
    // not gamma_p(): its reverse mode drops the derivative once J0 / tau > 10 shape
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
      // kappa c - log I0(kappa) would cancel the gamma weight away once kappa nears 1e17, see .circmix_vp_log_integrand()
      lp[i] = circmix_simpson_lw(i, nodes) + ti + gamma_lpdf(J | shape, rate)
              + kappa * (c - 1) - circmix_log_besselI0_scaled(kappa);
    }
    return log_sum_exp(log_sum_exp(lp) + log(step / 3), log(tail)) - log(2 * pi());
  }

  // log integral over J ~ gamma(shape, scale = tau) of one von Mises with c = cos(y - mu); grid centred on the integrand's mode, see .circmix_vp_component_mode(); below shape 2 the analytic tail of circmix_vp_component_tail()
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
      lp[i] = circmix_simpson_lw(i, nodes) + ti + gamma_lpdf(J | shape, rate)
              + kappa * (c - 1) - circmix_log_besselI0_scaled(kappa);
    }
    return log_sum_exp(lp) + log(step / 3) - log(2 * pi());
  }

  // circmix_ld with J marginalised over gamma(J(kappa) / tau, scale = tau), one grid per component; tau = 0 is constant precision
  real circmix_vp_ld(vector cosd, vector logw, real logw_guess,
                     real kappa, real tau, int nodes,
                     data vector logk, data vector dlogk,
                     data real logJ_min, data real dlogJ) {
    // tau -> Inf sends the shape to 0 with the mean fixed, leaving J at 0: memory is uniform, as in .circmix_vp_ld(), and close to that limit integrating would hit a subnormal shape
    if (tau == positive_infinity() || circmix_vp_uniform(circmix_J(kappa), tau)) {
      return circmix_ld(cosd, logw, logw_guess, 0);
    }
    real shape = circmix_vp_shape(kappa, tau);
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

  // slot allocation of capacity K over set size ss: [floor(K / ss), probability of one more slot]
  vector circmix_slots(real K, int ss) {
    real q = K * inv(ss);
    real f = floor(q);
    return [f, q - f]';
  }

  // slot-averaging branches: [held, kappa_lo, kappa_hi, extra]; kappa_lo is 0 when no slot is held, extra is the probability of one more slot
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
    // extra stays on the probability scale: log(extra) at a whole K / ss is log(0), whose NaN adjoint reaches K even unused
    return [held, kappa_lo, kappa_hi, allocation[2]]';
  }

  // slot averaging (Zhang & Luck, 2008): a j-slot item has J = j * J(kappa); composes with variable precision via tau
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
    // whole K / ss: the same value without log(0), whose NaN adjoint would poison K; the gradient is the right-hand slope
    return ld_lo + extra * expm1(ld_hi - ld_lo);
  }

  // target first, then the ss - 1 active non-targets of a trial
  vector circmix_locations(real mu, vector nt, int ss) {
    return append_row(mu, head(nt, ss - 1));
  }
