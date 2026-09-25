  // Sampling distribution of the EZ-diffusion summaries, shared by the 3par and
  // 4par likelihoods. Both have the shape of the EZ model they extend:
  //   hits ~ binomial(trials, pC)
  //   vrt  ~ gamma(...)          centred on VRT
  //   mrt  ~ normal(...)         centred on ndt + MDT
  // pC, MDT and VRT are the EZ equations. What is added are the third and fourth
  // cumulants k3 and k4 of the decision time, which set the spread of vrt and
  // its correlation with mrt (ezdm_summaries_lpdf). Everything else in this file
  // computes MDT, VRT, k3 and k4 without losing digits at small drift.
  //
  // For X_t = x0 + v t + s W_t absorbing at +-z, the boundary-conditional
  // cumulant generating function is log sinh(q b) - log sinh(q b0) up to a
  // constant, with lambda the argument of the transform,
  // q = sqrt(v^2 + 2 lambda s^2) / s^2, b0 = 2z = bound and b the
  // distance from the starting point to the far boundary. The lower boundary is
  // the upper one at x0 -> -x0, so this covers both, and 3par is b = bound / 2.
  //
  // q^2 = w + 2 lambda / s^2 is linear in lambda, so with w = drift^2 / s^4 the
  // cumulants are (-1)^n (2 / s^2)^n f^(n)(w) for
  // f(w) = log sinh(b sqrt(w)) - log sinh(b0 sqrt(w)). Drift enters only through
  // w, so nothing here needs an absolute value or a zero-drift branch.
  //
  // Two regimes. Below t = b0 sqrt(w) = 0.7 the closed forms lose k4 outright
  // (relative error 1e-4 at t = 0.1), because each derivative cancels to order
  // t^(2n), and the series of ezdm_series.stan take over, which are exact at
  // w = 0. There, with G(x) = log(sinh(sqrt x) / sqrt x), the log terms of b
  // and b0 differ by a constant, so f^(n)(w) = b^(2n) G^(n)(b^2 w)
  // - b0^(2n) G^(n)(b0^2 w). The code calls t = 0.7 the seam; it appears as
  // x < 0.49 for x = b0^2 w and as y < 0.1225 for y = (b0 / 2)^2 w.
  //
  // Everything is scalar on purpose. Each function is called once per
  // observation, and a vector or array built inside a user-defined function is
  // a heap allocation per call that stanc cannot optimize away; returning the
  // cumulants as vectors cost about one and a half times as much per gradient.
  // This is also why the four b0 terms travel as four arguments.
  //
  // The R code in R/distributions.R is the plainer statement of the same math,
  // and the parity tests in tests/testthat/test-model_ezdm.R hold the two
  // together:
  //   ezdm_summaries_lpdf                        .ez_rt_terms()
  //   ezdm_boundary_lpdf, ezdm_symmetric_lpdf    .ezdm_cumulants()
  //   ezdm_logit_pc                              .ezdm_logit_pc()
  // R has no counterpart of ezdm_symmetric_lpdf: it takes 3par through the
  // general path at b = bound / 2.

  // Closed-form regime. The w-derivatives of log sinh(t), t = b sqrt(w), with
  // their powers of w stripped out are polynomials in p = t coth(t) and
  // q = t^2 csch(t)^2:
  //   D1 = p, D2 = -p - q, D3 = 3p + 3q + 2pq,
  //   D4 = -(15p + 15q + 12pq + 2q(2t^2 + 3q)),
  // and k_n = (-1)^n (D_n(b) - D_n(b0)) / (s^2 w)^n. Both come from one
  // exp(-2t), which stays finite where cosh(t) / sinh(t) is Inf / Inf.
  real ezdm_coth_term(real t, real exp_m2t) {
    return t * (1 + exp_m2t) / (1 - exp_m2t);
  }

  real ezdm_csch_term(real t, real exp_m2t) {
    return 4 * square(t) * exp_m2t * inv_square(1 - exp_m2t);
  }

  real ezdm_cgf_d3(real p, real q) {
    return 3 * (p + q) + 2 * p * q;
  }

  real ezdm_cgf_d4(real p, real q, real t) {
    return -(15 * (p + q) + 12 * p * q + 2 * q * (2 * square(t) + 3 * q));
  }

  // Log density of the two RT summaries of n responses, given the mean
  // (ndt + MDT) and the cumulants VRT, k3, k4 of one response time. var_rt
  // follows a Gamma matched to its mean VRT and its variance
  // W = k4 / n + 2 VRT^2 / (n - 1), and mean_rt is normal conditional on it with
  // Cov(mean_rt, var_rt) = k3 / n. With k3 = k4 = 0 these are the terms of the
  // previous likelihood, gamma((n - 1) / 2, (n - 1) / (2 VRT)) and
  // normal(mean, sqrt(VRT / n)).
  real ezdm_summaries_lpdf(real mrt, real vrt, int n, real mean, real VRT,
                           real k3, real k4) {
    real inv_W = inv(k4 / n + 2 * square(VRT) / (n - 1));
    real cov_mean_var = k3 / n;
    return gamma_lpdf(vrt | square(VRT) * inv_W, VRT * inv_W)
           + normal_lpdf(mrt | mean + cov_mean_var * inv_W * (vrt - VRT),
                         sqrt(VRT / n - square(cov_mean_var) * inv_W));
  }

  // RT summaries of the n responses at the boundary a distance b from the start
  // point. b0_d1..b0_d4 are the b0 terms of the four cumulants, which the two
  // boundaries of 4par share: (2 b0^2 / s^2)^n G^(n)(b0^2 w) in the series
  // regime, D_n(b0 sqrt(w)) in the closed one. In both regimes the four
  // arguments handed on are ndt + MDT, VRT, k3 and k4 of this boundary.
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

  // RT summaries of 3par, where the start point is at b = b0 / 2. There
  // log sinh(t / 2) - log sinh(t) = -log(2 cosh(t / 2)), so the difference of two
  // terms becomes one: f(w) = -C(y) - log 2 with C(y) = log cosh(sqrt y) and
  // y = b^2 w. This is the only reason 3par does not call ezdm_boundary_lpdf
  // with b = bound / 2, which took about 1.6 times as long per gradient.
  //
  // The closed forms are the ones above with p = u tanh(u), q = u^2 sech(u)^2
  // and the sign of q reversed, because tanh and -sech^2 obey the same
  // derivative rules as coth and csch^2. With u = bound drift / (2 s^2) the
  // first two are the EZ equations: MDT = bound / (2 drift) tanh(u) and
  // VRT = bound s^2 / (2 drift^3) (tanh(u) - u sech(u)^2).
  real ezdm_symmetric_lpdf(real mrt, real vrt, int n, real ndt, real bound,
                           real w, real s_sq) {
    real y = 0.25 * square(bound) * w;
    // the seam of the general path, t = 2 sqrt(y) = 0.7
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

  // Logit of pC, the probability of the upper boundary, for a start point
  // b_upper below it and b_lower above the lower one. With k = drift / s^2
  // signed, pC = expm1(-2 k b_upper) / expm1(-2 k bound) is the lambda = 0 value
  // of the same transform, and 1 - pC is pC at the mirrored start point and
  // drift. In their ratio the bound terms cancel:
  //   logit(pC) = up + log(1 - exp(-up)) - log(1 - exp(-lo)),
  // up = 2 k b_upper, lo = 2 k b_lower. No log is taken of a probability that
  // has rounded to 0 or 1, mirroring drift and start point flips the sign
  // exactly, and at b_upper = b_lower it is drift bound / s^2, the 3par logit.
  //
  // For k < 0 the same ratio is written in exp(up) and exp(lo), which stay
  // finite there. Near k = 0 it is 0 / 0. Its limit b_upper / b_lower would be
  // a constant, with a zero gradient in drift where the true one is not, and
  // drift = 0 is where a sampler initialized at zero starts. The series
  // expm1(x) / x = 1 + x / 2 + x^2 / 6 + x^3 / 24 is exact to 1e-18 below
  // |x| = 1e-4 and carries the gradient.
  real ezdm_log_expm1_ratio(real x) {
    return log1p(x * (0.5 + x * (1.0 / 6 + x / 24)));
  }

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
