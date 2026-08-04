  // Half-width of the arc around the probe within which the Bayesian decision
  // rule of Lin & Oberauer (2022, Eq. 8) yields a "same" response. Their
  // Appendix B shows that thetat cancels out of the boundary when the criterion
  // is 0, so free_crit gates the extra term and its contribution to the gradient.
  real mixture2p_cd_crit_angle(real kappa, real criterion, real thetat,
                               int free_crit) {
    real cos_d = log_modified_bessel_first_kind(0, kappa) / kappa;
    if (free_crit) {
      real scaled = (exp(-criterion) - 1 + thetat) / thetat;
      if (scaled <= 0) return pi();
      cos_d += log(scaled) / kappa;
    }
    if (cos_d >= 1) return 0;
    if (cos_d <= -1) return pi();
    return acos(cos_d);
  }

  real mixture2p_cd_lpmf(array[] int y, vector mu, vector kappa, vector thetat,
                         vector criterion, data array[] real probe, data vector gl_x,
                         data vector gl_w, data int free_crit) {
    int N = size(y);
    int n_nodes = num_elements(gl_x);
    real out = 0;

    for (n in 1:N) {
      real half_width = mixture2p_cd_crit_angle(kappa[n], criterion[n],
                                                thetat[n], free_crit);
      real p_same;

      if (half_width <= 0) {
        p_same = 0;
      } else if (half_width >= pi()) {
        p_same = 1;
      } else {
        // mass of the retrieval distribution inside the arc, by Gauss-Legendre
        vector[n_nodes] x = half_width * gl_x + probe[n] - mu[n];
        vector[n_nodes] dens = exp(kappa[n] * cos(x) -
          log_modified_bessel_first_kind(0, kappa[n])) * inv(2 * pi());
        p_same = thetat[n] * half_width * dot_product(dens, gl_w) +
          (1 - thetat[n]) * half_width * inv(pi());
      }

      p_same = fmin(fmax(p_same, machine_precision()), 1 - machine_precision());
      out += y[n] == 1 ? log1m(p_same) : log(p_same);
    }

    return out;
  }
