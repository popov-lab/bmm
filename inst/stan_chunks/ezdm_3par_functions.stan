  // EZ-Diffusion Model likelihood for aggregated data, symmetric start point.
  //
  // The response counts follow a binomial; the two RT summaries follow the
  // sampling distribution of the mean and variance of `trials` decision times,
  // matched to the exact first four cumulants (issue #407). Decision times are
  // right-skewed, so the independent normal and scaled chi-square terms this
  // replaced understated the sampling variance of `vrt` several-fold and
  // ignored its correlation with `mrt`.
  //
  // With a symmetric start point the decision time is independent of which
  // boundary is hit, so all `trials` responses inform one set of cumulants.
  //
  // mu is a dummy dpar required by brms and is not used.
  real ezdm_3par_lpdf(real mrt, real mu, real drift, real bound, real ndt,
                      real s, real vrt, int hits, int trials) {
    real s_sq = square(s);
    real k = drift / s_sq;
    // logit of expm1(-k bound) / expm1(-2 k bound) = 1 / (1 + exp(-k bound))
    return binomial_logit_lpmf(hits | trials, k * bound)
           + ezdm_symmetric_lpdf(mrt | vrt, trials, ndt, bound, square(k), s_sq);
  }
