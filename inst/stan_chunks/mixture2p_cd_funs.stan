  #include 'fun_tan_half.stan'

  // Change detection likelihood for mixture2p model
  // Based on Lin & Oberauer (2022, Cognitive Psychology)
  //
  // The model computes P(response = "change" | probe, kappa, thetat, beta)
  // via numerical integration over possible retrieved features x:
  //   P(change) = integral P[LLR(x, probe) > beta] * P_retrieve(x) dx
  //
  // P_retrieve(x) = thetat * vM(x | 0, kappa) + (1 - thetat) / (2*pi)
  //
  // LLR(x, probe) = log P(x, probe | change) - log P(x, probe | same)
  //   where P(x, probe | same) assumes probe = target, so x ~ retrieval dist centered on probe
  //   and P(x, probe | change) assumes probe != target, so x and probe are independent

  // mu is included in the signature because brms always passes it as the
  // first dpar, but it is unused in CD models (void_mu = TRUE, fixed to 0).
  // vreal1 (probe_centered) is passed last after all dpars.
  real mixture2p_cd_lpmf(int y, real mu, real kappa, real thetat, real beta, real probe) {
    int n_quad = 101;
    real dx = 2 * pi() / (n_quad - 1);
    real p_change = 0;
    real p_guess = 1 - thetat;
    real log_uniform = -log(2 * pi());

    for (i in 1:n_quad) {
      real x = -pi() + (i - 1) * dx;

      // Retrieval density: P(x | memory array) centered at 0 (target)
      real log_p_retrieve = log_mix(thetat,
                                     von_mises_lpdf(x | 0, kappa),
                                     log_uniform);
      real p_retrieve = exp(log_p_retrieve);

      // LLR computation:
      // H_same: probe is the target, so x ~ retrieval dist centered at probe
      //   P(x | same, probe) = thetat * vM(x | probe, kappa) + (1-thetat)/(2pi)
      //   P(probe | same) = 1/(2pi)  [uniform over circle]
      //   joint: P(x, probe | same) = P(x|same,probe) * P(probe|same)
      real log_p_x_given_same = log_mix(thetat,
                                         von_mises_lpdf(x | probe, kappa),
                                         log_uniform);
      real log_p_same = log_p_x_given_same + log_uniform;

      // H_change: probe is random, independent of memory
      //   P(x | change) = P(x | memory array) = retrieval dist centered at 0
      //   P(probe | change) = 1/(2pi)
      //   joint: P(x, probe | change) = P(x|change) * P(probe|change)
      real log_p_change_hyp = log_p_retrieve + log_uniform;

      // LLR = log P(change) - log P(same)
      real llr = log_p_change_hyp - log_p_same;

      // Decision: respond "change" if LLR > beta
      if (llr > beta) {
        p_change += p_retrieve * dx;
      }
    }

    // Clamp to valid probability range
    p_change = fmin(fmax(p_change, 1e-10), 1 - 1e-10);

    if (y == 1) {
      return log(p_change);
    } else {
      return log1m(p_change);
    }
  }
