  // EZ-Diffusion Model likelihood for aggregated data with a free start point,
  // zr, given as a proportion of the boundary separation.
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
  // Every ezdm_ function called here is defined in ezdm_cumulants.stan or
  // ezdm_series.stan, which are assembled before this file; the header of the
  // former maps the code onto the EZ equations.
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

    // Each cumulant is a term in b minus the same term in b0 = bound. The b0
    // terms are computed here because both boundaries share them. x0 < 0.49 is
    // t = bound sqrt(w) < 0.7, where ezdm_cumulants.stan switches from its
    // closed forms to the series.
    real w = square(k);
    real x0 = square(bound) * w;
    int series = x0 < 0.49;
    real b0_d1;
    real b0_d2;
    real b0_d3;
    real b0_d4;
    if (series) {
      real scale = 2 * square(bound) / s_sq;
      real scale_sq = square(scale);
      b0_d1 = scale * ezdm_log_sinhc_d1(x0);
      b0_d2 = scale_sq * ezdm_log_sinhc_d2(x0);
      b0_d3 = scale_sq * scale * ezdm_log_sinhc_d3(x0);
      b0_d4 = square(scale_sq) * ezdm_log_sinhc_d4(x0);
    } else {
      real t0 = sqrt(x0);
      real exp_m2t0 = exp(-2 * t0);
      real p0 = ezdm_coth_term(t0, exp_m2t0);
      real q0 = ezdm_csch_term(t0, exp_m2t0);
      b0_d1 = p0;
      b0_d2 = -p0 - q0;
      b0_d3 = ezdm_cgf_d3(p0, q0);
      b0_d4 = ezdm_cgf_d4(p0, q0, t0);
    }

    if (hits >= 2) {
      lp += ezdm_boundary_lpdf(mrt_upper | vrt_upper, hits, ndt, b_upper, w,
                               s_sq, series, b0_d1, b0_d2, b0_d3, b0_d4);
    }
    if (misses >= 2) {
      lp += ezdm_boundary_lpdf(mrt_lower | vrt_lower, misses, ndt, b_lower, w,
                               s_sq, series, b0_d1, b0_d2, b0_d3, b0_d4);
    }
    return lp;
  }
