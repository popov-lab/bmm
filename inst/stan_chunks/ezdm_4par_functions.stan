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
    real s_sq = square(s);
    real k = drift / s_sq;
    real b_upper = zr * bound;
    real b_lower = bound - b_upper;

    // the logit of the EZ proportion correct for a free start point; at
    // zr = 0.5 it is drift bound / s^2, the 3par one
    real lp = binomial_logit_lpmf(hits | trials, ezdm_logit_pc(b_upper, b_lower, k));
    if (hits < 2 && misses < 2) {
      return lp;
    }

    real w = square(k);
    real x0 = square(bound) * w;
    int series = x0 < 0.49;
    real h1;
    real h2;
    real h3;
    real h4;
    if (series) {
      real scale = 2 * square(bound) / s_sq;
      real scale_sq = square(scale);
      h1 = scale * ezdm_log_sinhc_d1(x0);
      h2 = scale_sq * ezdm_log_sinhc_d2(x0);
      h3 = scale_sq * scale * ezdm_log_sinhc_d3(x0);
      h4 = square(scale_sq) * ezdm_log_sinhc_d4(x0);
    } else {
      real t0 = sqrt(x0);
      real e0 = exp(-2 * t0);
      real p0 = ezdm_coth_term(t0, e0);
      real q0 = ezdm_csch_term(t0, e0);
      h1 = p0;
      h2 = -p0 - q0;
      h3 = ezdm_cgf_d3(p0, q0);
      h4 = ezdm_cgf_d4(p0, q0, t0);
    }

    if (hits >= 2) {
      lp += ezdm_boundary_lpdf(mrt_upper | vrt_upper, hits, ndt, b_upper, w,
                               s_sq, series, h1, h2, h3, h4);
    }
    if (misses >= 2) {
      lp += ezdm_boundary_lpdf(mrt_lower | vrt_lower, misses, ndt, b_lower, w,
                               s_sq, series, h1, h2, h3, h4);
    }
    return lp;
  }
