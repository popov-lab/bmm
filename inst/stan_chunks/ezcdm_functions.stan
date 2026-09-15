// EZ circular diffusion model: exact von Mises sufficient-statistic term for the
// angle summaries; for RT variance a Gamma matched to the exact mean VRT and
// variance W = k4 / n + 2 VRT^2 / (n - 1) of the sample variance; for mean RT a
// normal conditional on the RT variance with the exact covariance k3 / n, where
// k3 and k4 are the third and fourth cumulants of the decision time.
real ezcdm_lpdf(real mean_angle, real mu, real driftrate, real driftangle,
                real bound, real ndt,
                real var_angle, real mean_rt, real var_rt, int trials) {
  real kappa = bound * driftrate;
  real log_I0 = log_modified_bessel_first_kind(0, kappa);
  real k2 = square(kappa);
  // MRT = ndt + L r, VRT = L^2 v, k3 = L^3 c3, k4 = L^4 c4 with L = bound / driftrate,
  // except below kappa = 0.01, where L = bound^2 and r, v, c3, c4 are the series
  // divided by powers of kappa: there R^2 - 1 + 2R/kappa cancels, and L would
  // overflow when driftrate underflows. Same regimes as .ezcdm_moments().
  real L;
  real r;
  real v;
  real c3;
  real c4;
  if (kappa < 1e-2) {
    L = square(bound);
    r = 0.5 - k2 / 16 + square(k2) / 96;
    v = 1.0 / 8 - k2 / 24;
  } else {
    L = bound / driftrate;
    if (kappa > 1000) {
      // exp(log I1 - log I0) loses ~1e-12, which the cancellation in v amplifies;
      // the asymptotic expansion is exact to double precision from kappa = 1000
      real u = inv(kappa);
      r = (1 - u * (3.0 / 8 + u * (15.0 / 128 + u * (105.0 / 1024 + u * 14175.0 / 98304))))
          / (1 + u * (1.0 / 8 + u * (9.0 / 128 + u * (75.0 / 1024 + u * 11025.0 / 98304))));
    } else {
      r = exp(log_modified_bessel_first_kind(1, kappa) - log_I0);
    }
    v = square(r) - 1 + 2 * r / kappa;
  }
  if (kappa < 0.25) {
    c3 = 1.0 / 12 - k2 * (11.0 / 256 - k2 * (19.0 / 1280
         - k2 * (473.0 / 110592 - k2 * (1145.0 / 1032192
         - k2 * (101369.0 / 377487360 - k2 * 946523.0 / 15288238080.0)))));
    c4 = 11.0 / 128 - k2 * (19.0 / 320 - k2 * (473.0 / 18432
         - k2 * (1145.0 / 129024 - k2 * (101369.0 / 37748736
         - k2 * 946523.0 / 1274019840))));
    if (kappa >= 1e-2) {
      c3 *= kappa * k2;
      c4 *= square(k2);
    }
  } else {
    real u = inv(kappa);
    if (kappa > 100) {
      c3 = square(u) * (3 - u * (4 + u * (15.0 / 8 + u * (3 + u * (875.0 / 128
           + u * (39.0 / 2 + u * (67599.0 / 1024 + u * 515.0 / 2)))))));
      c4 = u * square(u) * (15 - u * (24 + u * (105.0 / 8 + u * (24 + u * (7875.0 / 128
           + u * (195 + u * 743589.0 / 1024))))));
    } else {
      real r2 = square(r);
      real u2 = square(u);
      c3 = r2 * (2 * r + 6 * u) + r * (8 * u2 - 2) - 4 * u;
      c4 = r2 * (6 * r2 + 24 * r * u + 44 * u2 - 8) + r * u * (48 * u2 - 20) + 2 - 24 * u2;
    }
  }
  real VRT = square(L) * v;
  // n W / L^4; the Gamma shape is VRT^2 / W and the conditional normal of
  // mean_rt has slope (k3 / n) / W and variance VRT / n - (k3 / n)^2 / W
  real nW = c4 + 2 * trials * square(v) / (trials - 1);
  real shape = trials * square(v) / nW;
  return trials * (kappa * (1 - var_angle) * cos(mean_angle - driftangle) - log_I0)
         + gamma_lpdf(var_rt | shape, shape / VRT)
         + normal_lpdf(mean_rt | ndt + L * (r + c3 * (var_rt - VRT) / (square(L) * nW)),
                       L * sqrt((v - square(c3) / nW) / trials));
}
