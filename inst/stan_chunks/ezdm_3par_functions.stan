  // EZ-diffusion likelihood for aggregated data, with the start point midway between the boundaries
  // mu is a dummy dpar required by brms and is not used
  real ezdm_3par_lpdf(real mrt, real mu, real drift, real bound, real ndt,
                      real s, real vrt, int hits, int trials) {
    real s_sq = square(s);
    real k = drift / s_sq;
    // k * bound is the logit of pC, which stays finite where pC itself rounds to 1
    return binomial_logit_lpmf(hits | trials, k * bound)
           + ezdm_symmetric_lpdf(mrt | vrt, trials, ndt, bound, square(k), s_sq);
  }
