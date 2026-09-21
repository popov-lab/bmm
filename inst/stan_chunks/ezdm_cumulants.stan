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
  // cumulants as vectors cost 1.5 to 1.8 times as much per gradient.

  // Series regime. log sinh(sqrt x) = G(x) + log(x) / 2 with
  // G(x) = log(sinh(sqrt x) / sqrt x) = sum_j a_j x^j, and the log terms of b
  // and b0 differ by a constant, so f^(n)(w) = b^(2n) G^(n)(b^2 w)
  // - b0^(2n) G^(n)(b0^2 w). The functions below are G^(n) by Horner, 16 terms,
  // with a_j j! / (j - n)! formed in exact rational arithmetic and rounded
  // once (generated, not typed: a wrong coefficient makes the truncated series
  // diverge instead of failing loudly). At x <= 0.49 the truncation error is
  // below 1e-13.
  real ezdm_log_sinhc_d1(real x) {
    real acc = -1.2336844022586037e-16;
    acc = fma(acc, x, 1.2175977014591684e-15);
    acc = fma(acc, x, -1.2017207666653852e-14);
    acc = fma(acc, x, 1.1860508700116827e-13);
    acc = fma(acc, x, -1.1705853409912441e-12);
    acc = fma(acc, x, 1.1553216299501312e-11);
    acc = fma(acc, x, -1.1402575602296092e-10);
    acc = fma(acc, x, 1.1253923258404497e-09);
    acc = fma(acc, x, -1.1107304394989839e-08);
    acc = fma(acc, x, 1.0962973925936889e-07);
    acc = fma(acc, x, -1.0822021404031986e-06);
    acc = fma(acc, x, 1.0688899577788467e-05);
    acc = fma(acc, x, -0.00010582010582010582);
    acc = fma(acc, x, 0.0010582010582010583);
    acc = fma(acc, x, -0.011111111111111112);
    acc = fma(acc, x, 0.16666666666666666);
    return acc;
  }

  real ezdm_log_sinhc_d2(real x) {
    real acc = -1.8505266033879055e-15;
    acc = fma(acc, x, 1.704636782042836e-14);
    acc = fma(acc, x, -1.5622369966650008e-13);
    acc = fma(acc, x, 1.4232610440140192e-12);
    acc = fma(acc, x, -1.2876438750903686e-11);
    acc = fma(acc, x, 1.1553216299501312e-10);
    acc = fma(acc, x, -1.0262318042066481e-09);
    acc = fma(acc, x, 9.003138606723598e-09);
    acc = fma(acc, x, -7.775113076492888e-08);
    acc = fma(acc, x, 6.577784355562133e-07);
    acc = fma(acc, x, -5.411010702015993e-06);
    acc = fma(acc, x, 4.275559831115387e-05);
    acc = fma(acc, x, -0.00031746031746031746);
    acc = fma(acc, x, 0.0021164021164021165);
    acc = fma(acc, x, -0.011111111111111112);
    return acc;
  }

  real ezdm_log_sinhc_d3(real x) {
    real acc = -2.5907372447430678e-14;
    acc = fma(acc, x, 2.2160278166556864e-13);
    acc = fma(acc, x, -1.874684395998001e-12);
    acc = fma(acc, x, 1.565587148415421e-11);
    acc = fma(acc, x, -1.2876438750903685e-10);
    acc = fma(acc, x, 1.039789466955118e-09);
    acc = fma(acc, x, -8.209854433653185e-09);
    acc = fma(acc, x, 6.302197024706518e-08);
    acc = fma(acc, x, -4.6650678458957325e-07);
    acc = fma(acc, x, 3.2888921777810667e-06);
    acc = fma(acc, x, -2.1644042808063973e-05);
    acc = fma(acc, x, 0.0001282667949334616);
    acc = fma(acc, x, -0.0006349206349206349);
    acc = fma(acc, x, 0.0021164021164021165);
    return acc;
  }

  real ezdm_log_sinhc_d4(real x) {
    real acc = -3.367958418165988e-13;
    acc = fma(acc, x, 2.659233379986824e-12);
    acc = fma(acc, x, -2.0621528355978012e-11);
    acc = fma(acc, x, 1.5655871484154212e-10);
    acc = fma(acc, x, -1.1588794875813318e-09);
    acc = fma(acc, x, 8.318315735640945e-09);
    acc = fma(acc, x, -5.74689810355723e-08);
    acc = fma(acc, x, 3.781318214823911e-07);
    acc = fma(acc, x, -2.3325339229478664e-06);
    acc = fma(acc, x, 1.3155568711124267e-05);
    acc = fma(acc, x, -6.493212842419192e-05);
    acc = fma(acc, x, 0.0002565335898669232);
    acc = fma(acc, x, -0.0006349206349206349);
    return acc;
  }

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

  // P(hit the boundary at distance b) = expm1(-2 k b) / expm1(-2 k b0), the
  // lambda = 0 value of the same transform with k = drift / s^2 signed. Both
  // expm1 calls overflow when k < 0, so the identity
  // expm1(u) = -exp(u) expm1(-u) moves the evaluation to the finite side.
  real ezdm_pc(real b, real b0, real k) {
    real u = -2 * k * b;
    real u0 = -2 * k * b0;
    if (u0 == 0) {
      return b / b0;
    }
    if (u0 > 0) {
      return exp(u - u0) * expm1(-u) / expm1(-u0);
    }
    return expm1(u) / expm1(u0);
  }
