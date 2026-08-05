  // Shared pieces of the Bayesian decision rule for single-probe change
  // detection (Lin & Oberauer, 2022). The observer retrieves a feature x and
  // responds "change" when the log-likelihood ratio of their Equation 8
  // exceeds the criterion. That ratio increases monotonically with the
  // distance between x and the probe, so the region yielding a "same"
  // response is always an arc centred on the probe. The probability of a
  // "change" response is then the retrieval mass outside that arc, and each
  // model only has to supply its own retrieval distribution.

  // Half-width of the "same" arc: the angular distance at which the von Mises
  // density used by the decision process crosses the uniform, scaled by the
  // criterion. Their Appendix B shows that p_s cancels when the criterion is
  // 0, so free_crit keeps it out of the gradient in the unbiased model.
  real cd_crit_angle(real kappa, real criterion, real p_s, int free_crit) {
    real cos_d = log_modified_bessel_first_kind(0, kappa) / kappa;
    if (free_crit) {
      real scaled = (exp(-criterion) - 1 + p_s) / p_s;
      if (scaled <= 0) return pi();
      cos_d += log(scaled) / kappa;
    }
    if (cos_d >= 1) return 0;
    if (cos_d <= -1) return pi();
    return acos(cos_d);
  }

  // Mass of a von Mises(centre, kappa) inside the arc of half-width hw around
  // the probe, by Gauss-Legendre quadrature.
  real cd_vm_arc_mass(real centre, real probe, real hw, real kappa,
                      data vector gl_x, data vector gl_w) {
    int n_nodes = num_elements(gl_x);
    vector[n_nodes] x = hw * gl_x + probe - centre;
    vector[n_nodes] dens = exp(kappa * cos(x) -
      log_modified_bessel_first_kind(0, kappa)) * inv(2 * pi());
    return hw * dot_product(dens, gl_w);
  }

  // Binary log-likelihood from the probability of a "same" response. Kept on
  // the "same" scale so that rare "same" responses do not lose precision.
  real cd_bernoulli_lpmf(int y, real p_same) {
    real p = fmin(fmax(p_same, machine_precision()), 1 - machine_precision());
    return y == 1 ? log1m(p) : log(p);
  }
