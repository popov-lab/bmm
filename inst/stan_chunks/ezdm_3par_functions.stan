  // EZ-diffusion log likelihood of aggregated data, start point midway between the boundaries; mu is unused
  real ezdm_3par_lpdf(real mrt, real mu, real drift, real bound, real ndt,
                      real s, real vrt, int hits, int trials) {
    real s_sq = square(s);
    real k = drift / s_sq;
    return binomial_logit_lpmf(hits | trials, k * bound)
           + ezdm_symmetric_lpdf(mrt | vrt, trials, ndt, bound, square(k), s_sq);
  }
