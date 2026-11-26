  // t * coth(t), computed from exp(-2t)
  real ezdm_coth_term(real t, real exp_m2t) {
    return t * (1 + exp_m2t) / (1 - exp_m2t);
  }

  // t^2 / sinh(t)^2, computed from exp(-2t)
  real ezdm_csch_term(real t, real exp_m2t) {
    return 4 * square(t) * exp_m2t * inv_square(1 - exp_m2t);
  }

  // third-derivative term of the cumulant generating function of the decision time
  real ezdm_cgf_d3(real p, real q) {
    return 3 * (p + q) + 2 * p * q;
  }

  // fourth-derivative term of the cumulant generating function of the decision time
  real ezdm_cgf_d4(real p, real q, real t) {
    return -(15 * (p + q) + 12 * p * q + 2 * q * (2 * square(t) + 3 * q));
  }

  // log density of the observed mean and variance of n response times, given their expected moments
  real ezdm_summaries_lpdf(real mrt, real vrt, int n, real mean, real VRT,
                           real k3, real k4) {
    real inv_W = inv(k4 / n + 2 * square(VRT) / (n - 1));
    real cov_mean_var = k3 / n;
    return gamma_lpdf(vrt | square(VRT) * inv_W, VRT * inv_W)
           + normal_lpdf(mrt | mean + cov_mean_var * inv_W * (vrt - VRT),
                         sqrt(VRT / n - square(cov_mean_var) * inv_W));
  }

  // log density of the response time summaries at one boundary, a distance b from the start point
  real ezdm_boundary_lpdf(real mrt, real vrt, int n, real ndt, real b, real w,
                          real s_sq, int series, real b0_d1, real b0_d2,
                          real b0_d3, real b0_d4) {
    if (series) {
      real x = square(b) * w;
      real scale = 2 * square(b) / s_sq;
      real scale_sq = square(scale);
      return ezdm_summaries_lpdf(mrt | vrt, n,
        ndt + b0_d1 - scale * ezdm_log_sinhc_d1(x),
        scale_sq * ezdm_log_sinhc_d2(x) - b0_d2,
        b0_d3 - scale_sq * scale * ezdm_log_sinhc_d3(x),
        square(scale_sq) * ezdm_log_sinhc_d4(x) - b0_d4);
    }
    real t = b * sqrt(w);
    real exp_m2t = exp(-2 * t);
    real p = ezdm_coth_term(t, exp_m2t);
    real q = ezdm_csch_term(t, exp_m2t);
    real inv_s2w = inv(s_sq * w);
    real inv_s2w_sq = square(inv_s2w);
    return ezdm_summaries_lpdf(mrt | vrt, n,
      ndt + (b0_d1 - p) * inv_s2w,
      (-p - q - b0_d2) * inv_s2w_sq,
      (b0_d3 - ezdm_cgf_d3(p, q)) * inv_s2w_sq * inv_s2w,
      (ezdm_cgf_d4(p, q, t) - b0_d4) * square(inv_s2w_sq));
  }

  // log density of the response time summaries with the start point midway between the boundaries
  real ezdm_symmetric_lpdf(real mrt, real vrt, int n, real ndt, real bound,
                           real w, real s_sq) {
    real y = 0.25 * square(bound) * w;
    if (y < 0.1225) {
      real scale = 0.5 * square(bound) / s_sq;
      real scale_sq = square(scale);
      return ezdm_summaries_lpdf(mrt | vrt, n,
        ndt + scale * ezdm_log_cosh_d1(y),
        -scale_sq * ezdm_log_cosh_d2(y),
        scale_sq * scale * ezdm_log_cosh_d3(y),
        -square(scale_sq) * ezdm_log_cosh_d4(y));
    }
    real u = sqrt(y);
    real exp_m2u = exp(-2 * u);
    real p = u * (1 - exp_m2u) / (1 + exp_m2u);
    real q = -4 * y * exp_m2u * inv_square(1 + exp_m2u);
    real inv_s2w = inv(s_sq * w);
    real inv_s2w_sq = square(inv_s2w);
    return ezdm_summaries_lpdf(mrt | vrt, n,
      ndt + p * inv_s2w,
      (p + q) * inv_s2w_sq,
      ezdm_cgf_d3(p, q) * inv_s2w_sq * inv_s2w,
      -ezdm_cgf_d4(p, q, u) * square(inv_s2w_sq));
  }

  // log(expm1(x) / x) for x close to zero
  real ezdm_log_expm1_ratio(real x) {
    return log1p(x * (0.5 + x * (1.0 / 6 + x / 24)));
  }

  // logit of the predicted proportion of correct responses
  real ezdm_logit_pc(real b_upper, real b_lower, real k) {
    real up = 2 * k * b_upper;
    real lo = 2 * k * b_lower;
    if (abs(up + lo) < 1e-4) {
      return log(b_upper / b_lower) + up + ezdm_log_expm1_ratio(-up)
             - ezdm_log_expm1_ratio(-lo);
    }
    if (k > 0) {
      return up + log1m_exp(-up) - log1m_exp(-lo);
    }
    return lo + log1m_exp(up) - log1m_exp(lo);
  }
