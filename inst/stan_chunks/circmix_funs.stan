  /* Shared numerics for the circular mixture models (mixture2p, mixture3p, imm).
   *
   * Component weights always arrive normalised on the log scale: exp(logw)
   * together with exp(logw_guess) sums to one. The cosines cos(y - mu_k) are
   * passed in rather than computed here because they do not depend on kappa,
   * which is what makes the variable-precision quadrature affordable -- the
   * whole grid re-uses one set of cosines.
   */

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

  /* Inverse of circmix_J, by cubic Hermite interpolation on a grid uniform in
   * log J. Outside the tabulated range the exact asymptotics are used:
   * J -> kappa^2 / 2 as kappa -> 0, and J -> kappa - 1/2 as kappa -> infinity.
   * logk and dlogk hold log kappa and d log kappa / d log J at the nodes, and
   * are built by .circmix_kappa_table() in R and passed in as data, so that the
   * R and Stan implementations of the inverse cannot drift apart.
   */
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

  // log density of a mixture of von Mises sharing one concentration, plus a
  // uniform guessing component
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

  /* Variable precision: the same mixture, with the Fisher information of the
   * memory components marginalised over J ~ gamma(shape = Jbar / tau, scale =
   * tau), where Jbar = circmix_J(kappa). tau = 0 leaves a point mass at Jbar,
   * i.e. the constant-precision model, which is why tau is the parameter the
   * variable_precision argument frees.
   *
   * Nodes sit on a composite Simpson grid in log J, centred on
   * E[log J] = digamma(shape) + log(tau) with half width
   * 8 sd(log J) = 8 sqrt(trigamma(shape)). The offsets are constants, so only
   * the centre and the width depend on parameters, and they do so smoothly.
   */
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

  /* Slot allocation of a capacity K over set size ss. An item receives
   * floor(K / ss) slots with probability 1 - r and one more with probability r,
   * which is continuous in K across the integer crossings. Returns
   * [floor(K / ss), r].
   */
  vector circmix_slots(real K, int ss) {
    real q = K * inv(ss);
    real f = floor(q);
    return [f, q - f]';
  }
