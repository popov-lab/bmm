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

  // log integral over J ~ gamma(shape, scale = tau) of one von Mises with c = cos(y - mu); grid centred on the integrand's mode, see .circmix_vp_component()
  real circmix_vp_component(real c, real shape, real tau, int nodes,
                            data vector logk, data vector dlogk,
                            data real logJ_min, data real dlogJ) {
    real rate = inv(tau);
    real t = log((shape + 0.5) / (rate + 1 - c));
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
              + kappa * c - log_modified_bessel_first_kind(0, kappa);
    }
    return log_sum_exp(lp) + log(step / 3) - log(2 * pi());
  }

  // circmix_ld with J marginalised over gamma(J(kappa) / tau, scale = tau), one grid per component; tau = 0 is constant precision
  real circmix_vp_ld(vector cosd, vector logw, real logw_guess,
                     real kappa, real tau, int nodes,
                     data vector logk, data vector dlogk,
                     data real logJ_min, data real dlogJ) {
    // kappa = 0 has J = 0: the gamma on J is a point mass there, as in R
    if (tau <= 0 || kappa <= 0) {
      return circmix_ld(cosd, logw, logw_guess, kappa);
    }
    real shape = circmix_J(kappa) / tau;
    if (shape < 60.0 / nodes) {
      reject("bmm: the variable-precision quadrature needs a gamma shape ",
             "J(kappa)/tau of at least 60 / vp_nodes = ", 60.0 / nodes,
             ", but got ", shape,
             ". Raise vp_nodes above ", 60.0 / shape,
             ", or use a more informative prior on tau.");
    }
    if (8 * sqrt(trigamma(shape)) < 1e-6) {
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

  // slot averaging (Zhang & Luck, 2008): a j-slot item has J = j * J(kappa); composes with variable precision via tau
  real circmix_slot_averaging_ld(vector cosd, vector logw, real K, int ss,
                                 real kappa, real tau, int nodes,
                                 data vector logk, data vector dlogk,
                                 data real logJ_min, data real dlogJ) {
    vector[2] allocation = circmix_slots(K, ss);
    real slots = allocation[1];
    real extra = allocation[2];
    real J1 = circmix_J(kappa);

    real kappa_hi = circmix_kappa((slots + 1) * J1, logk, dlogk, logJ_min, dlogJ);
    real ld_hi = circmix_vp_ld(cosd, logw, negative_infinity(), kappa_hi, tau,
                               nodes, logk, dlogk, logJ_min, dlogJ);
    real ld_lo = -log(2 * pi());
    if (slots >= 0.5) {
      real kappa_lo = circmix_kappa(slots * J1, logk, dlogk, logJ_min, dlogJ);
      ld_lo = circmix_vp_ld(cosd, logw, negative_infinity(), kappa_lo, tau,
                            nodes, logk, dlogk, logJ_min, dlogJ);
    }
    if (extra > 0) {
      return log_sum_exp(log1m(extra) + ld_lo, log(extra) + ld_hi);
    }
    // whole K / ss: the same value without log(0), whose NaN adjoint would poison K; the gradient is the right-hand slope
    return ld_lo + extra * expm1(ld_hi - ld_lo);
  }
