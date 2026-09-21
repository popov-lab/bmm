  // Cumulants of the decision time conditional on hitting a boundary, shared by
  // the 3par and 4par EZ-diffusion likelihoods.
  //
  // For X_t = x0 + v t + s W_t absorbing at +-z, the boundary-conditional
  // cumulant generating function is log sinh(q b) - log sinh(q b0) up to a
  // constant, with q = sqrt(v^2 + 2 lambda s^2) / s^2, b0 = 2z = bound and b the
  // distance from the starting point to the far boundary. The lower boundary is
  // the upper one at x0 -> -x0, so this covers both, and 3par is b = bound / 2.
  //
  // q^2 = w + 2 lambda / s^2 is linear in lambda, so with w = drift^2 / s^4 the
  // cumulants are (-1)^n (2 / s^2)^n f^(n)(w) for
  // f(w) = log sinh(b sqrt(w)) - log sinh(b0 sqrt(w)). Drift enters only through
  // w, so nothing here needs an absolute value or a zero-drift branch.
  //
  // Two regimes, which must match .ezdm_cumulants() in R/distributions.R. Below
  // t = b0 sqrt(w) = 0.7 the closed forms lose k4 outright (relative error 1e-4
  // at t = 0.1), because each derivative cancels to order t^(2n), and a series
  // takes over that is exact at w = 0.
  //
  // Everything is scalar on purpose. Each function is called once per
  // observation, and a vector or array built inside a user-defined function is
  // a heap allocation per call that stanc cannot optimize away; returning the
  // cumulants as vectors cost about one and a half times as much per gradient.

  // Series regime. log sinh(sqrt x) = G(x) + log(x) / 2 with
  // G(x) = log(sinh(sqrt x) / sqrt x) = sum_j a_j x^j, and the log terms of b
  // and b0 differ by a constant, so f^(n)(w) = b^(2n) G^(n)(b^2 w)
  // - b0^(2n) G^(n)(b0^2 w). The G^(n) are generated into ezdm_series.stan.

  // Closed-form regime. The w-derivatives of log sinh(t), t = b sqrt(w), with
  // their powers of w stripped out are polynomials in p = t coth(t) and
  // q = t^2 csch(t)^2:
  //   D1 = p, D2 = -p - q, D3 = 3p + 3q + 2pq,
  //   D4 = -(15p + 15q + 12pq + 2q(2t^2 + 3q)),
  // and k_n = (-1)^n (D_n(b) - D_n(b0)) / (s^2 w)^n. Both come from one
  // e = exp(-2t), which stays finite where cosh(t) / sinh(t) is Inf / Inf.
  real ezdm_coth_term(real t, real e) {
    return t * (1 + e) / (1 - e);
  }

  real ezdm_csch_term(real t, real e) {
    return 4 * square(t) * e * inv_square(1 - e);
  }

  real ezdm_cgf_d3(real p, real q) {
    return 3 * (p + q) + 2 * p * q;
  }

  real ezdm_cgf_d4(real p, real q, real t) {
    return -(15 * (p + q) + 12 * p * q + 2 * q * (2 * square(t) + 3 * q));
  }

  // Log density of the two RT summaries of n responses given the cumulants of
  // one decision time: var_rt follows a Gamma matched to the mean VRT and the
  // variance W = k4 / n + 2 VRT^2 / (n - 1) of a sample variance, and mean_rt is
  // normal conditional on it with Cov(mean_rt, var_rt) = k3 / n. Matches
  // .ez_rt_terms() in R/distributions.R.
  real ezdm_summaries_lpdf(real mrt, real vrt, int n, real mean, real VRT,
                           real k3, real k4) {
    real inv_W = inv(k4 / n + 2 * square(VRT) / (n - 1));
    real cov_mean_var = k3 / n;
    return gamma_lpdf(vrt | square(VRT) * inv_W, VRT * inv_W)
           + normal_lpdf(mrt | mean + cov_mean_var * inv_W * (vrt - VRT),
                         sqrt(VRT / n - square(cov_mean_var) * inv_W));
  }

  // RT summaries of the n responses at the boundary a distance b from the start
  // point. h1..h4 are the b0 terms, which the two boundaries of 4par share:
  // (2 b0^2 / s^2)^n G^(n)(b0^2 w) in the series regime, D_n(b0 sqrt(w)) in the
  // closed one.
  real ezdm_boundary_lpdf(real mrt, real vrt, int n, real ndt, real b, real w,
                          real s_sq, int series, real h1, real h2, real h3,
                          real h4) {
    if (series) {
      real x = square(b) * w;
      real scale = 2 * square(b) / s_sq;
      real scale_sq = square(scale);
      return ezdm_summaries_lpdf(mrt | vrt, n,
        ndt + h1 - scale * ezdm_log_sinhc_d1(x),
        scale_sq * ezdm_log_sinhc_d2(x) - h2,
        h3 - scale_sq * scale * ezdm_log_sinhc_d3(x),
        square(scale_sq) * ezdm_log_sinhc_d4(x) - h4);
    }
    real t = b * sqrt(w);
    real e = exp(-2 * t);
    real p = ezdm_coth_term(t, e);
    real q = ezdm_csch_term(t, e);
    real c = inv(s_sq * w);
    real c_sq = square(c);
    return ezdm_summaries_lpdf(mrt | vrt, n,
      ndt + (h1 - p) * c,
      (-p - q - h2) * c_sq,
      (h3 - ezdm_cgf_d3(p, q)) * c_sq * c,
      (ezdm_cgf_d4(p, q, t) - h4) * square(c_sq));
  }

  // RT summaries of 3par, where the start point is at b = b0 / 2. There
  // log sinh(t / 2) - log sinh(t) = -log(2 cosh(t / 2)), so the difference of two
  // terms becomes one: f(w) = -C(y) - log 2 with C(y) = log cosh(sqrt y) and
  // y = b^2 w; its derivatives C^(n) are generated into ezdm_series.stan. The
  // closed forms are the ones above with p = u tanh(u), q = u^2 sech(u)^2 and
  // the sign of q reversed, because tanh and -sech^2 obey the same derivative
  // rules as coth and csch^2.
  real ezdm_symmetric_lpdf(real mrt, real vrt, int n, real ndt, real bound,
                           real w, real s_sq) {
    real y = 0.25 * square(bound) * w;
    // the seam of the generic path, t = 2 sqrt(y) = 0.7
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
    real e = exp(-2 * u);
    real p = u * (1 - e) / (1 + e);
    real q = -4 * y * e * inv_square(1 + e);
    real c = inv(s_sq * w);
    real c_sq = square(c);
    return ezdm_summaries_lpdf(mrt | vrt, n,
      ndt + p * c,
      (p + q) * c_sq,
      ezdm_cgf_d3(p, q) * c_sq * c,
      -ezdm_cgf_d4(p, q, u) * square(c_sq));
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
