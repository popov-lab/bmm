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
      return sqrt(2 * J);
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

  // circmix_ld with J marginalised over gamma(J(kappa) / tau, scale = tau) on a Simpson grid in log J; tau = 0 is constant precision
  real circmix_vp_ld(vector cosd, vector logw, real logw_guess,
                     real kappa, real tau, int nodes,
                     data vector logk, data vector dlogk,
                     data real logJ_min, data real dlogJ) {
    if (tau <= 0) {
      return circmix_ld(cosd, logw, logw_guess, kappa);
    }
    real shape = circmix_J(kappa) / tau;
    if (shape < 40.0 / nodes) {
      reject("bmm: the variable-precision quadrature holds the log density to ",
             "about 1e-5 while the gamma shape J(kappa)/tau stays above ",
             40.0 / nodes, ", but got ", shape,
             ". Raise vp_nodes, or use a more informative prior on tau.");
    }
    real half_width = 8 * sqrt(trigamma(shape));
    if (half_width < 1e-6) {
      return circmix_ld(cosd, logw, logw_guess, kappa);
    }
    real step = 2 * half_width / (nodes - 1);
    real centre = digamma(shape) + log(tau);
    real rate = inv(tau);
    vector[nodes] lp;
    for (i in 1:nodes) {
      real t = centre - half_width + (i - 1) * step;
      real J = exp(t);
      lp[i] = circmix_simpson_lw(i, nodes) + t + gamma_lpdf(J | shape, rate)
              + circmix_ld(cosd, logw, logw_guess,
                           circmix_kappa(J, logk, dlogk, logJ_min, dlogJ));
    }
    return log_sum_exp(lp) + log(step / 3);
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
    real lp_hi = log(extra)
                 + circmix_vp_ld(cosd, logw, negative_infinity(), kappa_hi, tau,
                                 nodes, logk, dlogk, logJ_min, dlogJ);

    if (slots < 0.5) {
      return log_sum_exp(log1m(extra) - log(2 * pi()), lp_hi);
    }

    real kappa_lo = circmix_kappa(slots * J1, logk, dlogk, logJ_min, dlogJ);
    return log_sum_exp(
      log1m(extra) + circmix_vp_ld(cosd, logw, negative_infinity(), kappa_lo,
                                   tau, nodes, logk, dlogk, logJ_min, dlogJ),
      lp_hi
    );
  }
