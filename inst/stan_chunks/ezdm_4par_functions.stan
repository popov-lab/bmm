  // EZ-diffusion likelihood for aggregated data, with a free start point zr as a proportion of bound
  // mu is a dummy dpar required by brms and is not used
  real ezdm_4par_lpdf(real mrt_upper, real mu, real drift, real bound, real ndt,
                      real zr, real s, real mrt_lower, real vrt_upper,
                      real vrt_lower, int hits, int trials, int rt_upper,
                      int rt_lower) {
    int misses = trials - hits;
    // check_data() fills the summaries of an unused boundary with placeholders, so these gate every read
    int use_upper = rt_upper == 1 && hits >= 2;
    int use_lower = rt_lower == 1 && misses >= 2;
    real s_sq = square(s);
    real k = drift / s_sq;
    real b_upper = zr * bound;
    real b_lower = bound - b_upper;

    // the logit of pC for a free start point; at zr = 0.5 it is the 3par one
    real lp = binomial_logit_lpmf(hits | trials, ezdm_logit_pc(b_upper, b_lower, k));
    if (!use_upper && !use_lower) {
      return lp;
    }

    // both boundaries share the b0 = bound terms; x0 < 0.49 is the seam to the series
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

    if (use_upper) {
      lp += ezdm_boundary_lpdf(mrt_upper | vrt_upper, hits, ndt, b_upper, w,
                               s_sq, series, b0_d1, b0_d2, b0_d3, b0_d4);
    }
    if (use_lower) {
      lp += ezdm_boundary_lpdf(mrt_lower | vrt_lower, misses, ndt, b_lower, w,
                               s_sq, series, b0_d1, b0_d2, b0_d3, b0_d4);
    }
    return lp;
  }
