  // EZ-Diffusion Model likelihood for aggregated data, with the start point
  // midway between the boundaries.
  //
  // The response counts follow a binomial; the two RT summaries follow the
  // sampling distribution of the mean and variance of `trials` decision times,
  // matched to its exact first four cumulants: the mean, the variance, and the
  // third and fourth, which carry the skew and the tail weight. Decision times
  // are right-skewed, so the independent normal and scaled chi-square terms this
  // replaced understated the sampling variance of `vrt` several-fold and
  // ignored its correlation with `mrt`.
  //
  // With a symmetric start point the decision time is independent of which
  // boundary is hit, so all `trials` responses inform one set of cumulants.
  //
  // ezdm_symmetric_lpdf is defined in ezdm_cumulants.stan, which is assembled
  // before this file and whose header maps the code onto the EZ equations.
  //
  // mu is a dummy dpar required by brms and is not used.
  real ezdm_3par_lpdf(real mrt, real mu, real drift, real bound, real ndt,
                      real s, real vrt, int hits, int trials) {
    real s_sq = square(s);
    real k = drift / s_sq;
    // k * bound is the logit of the EZ proportion correct,
    // pC = 1 / (1 + exp(-drift bound / s^2)); on that scale the binomial stays
    // finite where pC itself rounds to 1
    return binomial_logit_lpmf(hits | trials, k * bound)
           + ezdm_symmetric_lpdf(mrt | vrt, trials, ndt, bound, square(k), s_sq);
  }
