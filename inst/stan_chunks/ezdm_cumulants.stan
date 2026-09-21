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
  // Regimes and coefficients must match .ezdm_cumulants() in R/distributions.R.

  // The w-derivatives of log sinh(b sqrt(w)) with their powers of w stripped
  // out, as [P, Q, R, S]. Each is a sum of same-sign terms, so nothing cancels
  // here; only the difference the caller takes does. Above t = 30 the csch^2
  // terms are below 1e-26 and are dropped rather than evaluated, because
  // t^4 * csch^2(t) becomes Inf * 0 = NaN once t^4 overflows.
  vector ezdm_cgf_derivatives(real t) {
    if (t > 30) {
      return [t, -t, 3 * t, -15 * t]';
    }
    real coth = 1 / tanh(t);
    real csch_sq = 1 / square(sinh(t));
    real t_sq = square(t);
    real t_cubed = t_sq * t;
    real p = t * coth;
    return [
      p,
      -p - t_sq * csch_sq,
      3 * p + 3 * t_sq * csch_sq + 2 * t_cubed * coth * csch_sq,
      -(15 * p + 15 * t_sq * csch_sq + 12 * t_cubed * coth * csch_sq
        + 2 * t_sq * t_sq * csch_sq * (2 + 3 * csch_sq))
    ]';
  }

  real ezdm_falling_factorial(int j, int n) {
    real out = 1;
    for (i in 0 : n - 1) {
      out *= j - i;
    }
    return out;
  }

  // k_n = (-1)^n (2/s^2)^n sum_{j >= n} a_j (b^2j - b0^2j) j!/(j-n)! w^(j-n),
  // by Horner, with a_j the coefficients of log(sinh(x) / x). The coefficient
  // differences are formed analytically, so the cancellation that defeats the
  // closed forms at small drift never happens; exact at w = 0.
  vector ezdm_cumulants_series(real b, real b0, real w, real s) {
    int j_max = 16;
    // denominators carry a .0 because Stan integer literals stop at 2^31. R
    // without long double (e.g. aarch64 macOS) misparses integers above 2^53,
    // so a_12, a_13, a_15 and a_16 can differ from R's by 1-2 ulp; their weight
    // below the t = 0.7 seam is under 1e-20
    array[16] real a = {
      1.0 / 6.0, -1.0 / 180.0, 1.0 / 2835.0, -1.0 / 37800.0,
      1.0 / 467775.0, -691.0 / 3831077250.0, 2.0 / 127702575.0,
      -3617.0 / 2605132530000.0, 43867.0 / 350813659321125.0,
      -174611.0 / 15313294652906250.0, 155366.0 / 147926426347074375.0,
      -236364091.0 / 2423034863565078262500.0,
      1315862.0 / 144228265688397515625.0,
      -3392780147.0 / 3952575621190533915703125.0,
      6892673020804.0 / 84913182070036240111050234375.0,
      -7709321041217.0 / 999843529136357459316262500000.0
    };

    real b_sq = square(b);
    real b0_sq = square(b0);
    real power_b = b_sq;
    real power_b0 = b0_sq;
    array[16] real coef_diff;
    for (j in 1 : j_max) {
      coef_diff[j] = a[j] * (power_b - power_b0);
      power_b *= b_sq;
      power_b0 *= b0_sq;
    }

    vector[4] out;
    real scale = 2 / square(s);
    for (n in 1 : 4) {
      real acc = coef_diff[j_max] * ezdm_falling_factorial(j_max, n);
      for (step in 1 : j_max - n) {
        int j = j_max - step;
        acc = acc * w + coef_diff[j] * ezdm_falling_factorial(j, n);
      }
      out[n] = (n % 2 == 0 ? 1 : -1) * pow(scale, n) * acc;
    }
    return out;
  }

  // [MDT, VRT, k3, k4] for the boundary at distance b from the start point.
  vector ezdm_cumulants(real b, real b0, real w, real s) {
    real root_w = sqrt(w);
    // below t = 0.7 the closed forms lose k4 outright (relative error 1e-4 at
    // t = 0.1) because the difference cancels to order t^(2n)
    if (b0 * root_w < 0.7) {
      return ezdm_cumulants_series(b, b0, w, s);
    }
    vector[4] at_b = ezdm_cgf_derivatives(b * root_w);
    vector[4] at_b0 = ezdm_cgf_derivatives(b0 * root_w);
    real s_sq = square(s);
    return [
      (at_b0[1] - at_b[1]) / (s_sq * w),
      (at_b[2] - at_b0[2]) / (square(s_sq) * square(w)),
      (at_b0[3] - at_b[3]) / (pow(s_sq, 3) * pow(w, 3)),
      (at_b[4] - at_b0[4]) / (pow(s_sq, 4) * pow(w, 4))
    ]';
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

  // Sampling distribution of the two RT summaries given the cumulants, as
  // [gamma shape, gamma rate, conditional slope, conditional sd]. Matches
  // .ez_rt_terms() in R/distributions.R.
  vector ezdm_rt_terms(vector cumulants, int n) {
    real VRT = cumulants[2];
    real W = cumulants[4] / n + 2 * square(VRT) / (n - 1);
    real shape = square(VRT) / W;
    real cov_mean_var = cumulants[3] / n;
    return [shape, shape / VRT, cov_mean_var / W,
            sqrt(VRT / n - square(cov_mean_var) / W)]';
  }
