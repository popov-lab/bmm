  // EZ-Diffusion Model likelihood for aggregated data with a free start point.
  //
  // A biased start point makes the decision time depend on which boundary is
  // hit, so each boundary carries its own cumulants and its own response count,
  // and the summaries must be supplied separately for the two boundaries. Both
  // sets come from one function evaluated at the distance from the start point
  // to the far boundary: zr * bound above, (1 - zr) * bound below.
  //
  // A boundary with fewer than two responses has no sample variance and
  // contributes only through the binomial term. bmm() rarely gets here: brms
  // drops rows whose summaries are NA, which is how rezdm() and
  // ezdm_summary_stats() code such a boundary.
  //
  // mu is a dummy dpar required by brms and is not used.
  real ezdm_4par_lpdf(real mrt_upper, real mu, real drift, real bound, real ndt,
                      real zr, real s, real mrt_lower, real vrt_upper,
                      real vrt_lower, int hits, int trials) {
    int misses = trials - hits;
    real b_upper = zr * bound;
    real b_lower = bound - b_upper;
    real w = square(drift) / pow(s, 4);

    real lp = binomial_lpmf(hits | trials,
                            ezdm_pc(b_upper, bound, drift / square(s)));

    if (hits >= 2) {
      vector[4] cumulants = ezdm_cumulants(b_upper, bound, w, s);
      vector[4] rt = ezdm_rt_terms(cumulants, hits);
      lp += gamma_lpdf(vrt_upper | rt[1], rt[2])
            + normal_lpdf(mrt_upper | ndt + cumulants[1]
                                      + rt[3] * (vrt_upper - cumulants[2]),
                          rt[4]);
    }

    if (misses >= 2) {
      vector[4] cumulants = ezdm_cumulants(b_lower, bound, w, s);
      vector[4] rt = ezdm_rt_terms(cumulants, misses);
      lp += gamma_lpdf(vrt_lower | rt[1], rt[2])
            + normal_lpdf(mrt_lower | ndt + cumulants[1]
                                      + rt[3] * (vrt_lower - cumulants[2]),
                          rt[4]);
    }

    return lp;
  }
