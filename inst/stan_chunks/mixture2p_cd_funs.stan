  #include 'fun_tan_half.stan'

  // Change detection likelihood for mixture2p model
  // Based on Lin & Oberauer (2022, Cognitive Psychology)
  //
  // P(change) = integral P[LLR(x, probe) > beta] * P_retrieve(x) dx
  // P_retrieve(x) = thetat * vM(x | 0, kappa) + (1 - thetat) / (2*pi)
  // LLR(x, probe) = log P(x, probe | change) - log P(x, probe | same)

  real mixture2p_cd_lpmf(int y, real mu, real kappa, real thetat, real beta, real probe) {
    int n_quad = 101;
    real dx = 2 * pi() / (n_quad - 1);
    real p_change = 0;
    real log_uniform = -log(2 * pi());

    for (i in 1:n_quad) {
      real x = -pi() + (i - 1) * dx;

      real log_p_retrieve = log_mix(thetat,
                                     von_mises_lpdf(x | 0, kappa),
                                     log_uniform);
      real p_retrieve = exp(log_p_retrieve);

      real log_p_x_given_same = log_mix(thetat,
                                         von_mises_lpdf(x | probe, kappa),
                                         log_uniform);
      real log_p_same = log_p_x_given_same + log_uniform;
      real log_p_change_hyp = log_p_retrieve + log_uniform;
      real llr = log_p_change_hyp - log_p_same;

      if (llr > beta) {
        p_change += p_retrieve * dx;
      }
    }

    p_change = fmin(fmax(p_change, 1e-10), 1 - 1e-10);
    if (y == 1) return log(p_change);
    return log1m(p_change);
  }
