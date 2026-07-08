// log-PDF of the shifted Wald distribution
// Optimized to compute entirely in log-space for numerical stability
real swald_lpdf(real rt, real drift, real bound, real ndt, real sigma) {
  // compute shifted response time
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return negative_infinity();

  // pre-compute common terms
  real sigma_sq = square(sigma);
  real log_t = log(t_shifted);

  // log-normalization: log(bound / sqrt(2*pi*sigma^2*t^3))
  // = log(bound) - 0.5*(log(2*pi) + 2*log(sigma) + 3*log(t))
  real log_norm = log(bound) - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log_t);

  // log-kernel: -((bound - drift*t)^2) / (2*sigma^2*t)
  real residual = bound - drift * t_shifted;
  real log_kernel = -0.5 * square(residual) / (sigma_sq * t_shifted);

  return log_norm + log_kernel;
}

// log shifted Wald survivor function (complementary CDF)
// Optimized for numerical stability using Stan's log-space CDF functions
// and log_diff_exp to avoid overflow/underflow issues
real swald_lccdf(real rt, real drift, real bound, real ndt, real sigma) {
  // compute shifted response time
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return 0;  // process hasn't started, survival = 1

  // pre-compute common terms
  real sqrt_t = sqrt(t_shifted);
  real sigma_sqrt_t = sigma * sqrt_t;
  real sigma_sq = square(sigma);

  // standardized arguments for the normal CDF
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;

  // log of the scaling constant (kept in log-space to avoid overflow)
  real log_c = 2 * bound * drift / sigma_sq;

  // Survival function: S(t) = 1 - Phi(z1) - exp(c)*Phi(z2)
  //                         = Phi(-z1) - exp(c)*Phi(z2)
  //
  // Compute in log-space for numerical stability:
  // log(S(t)) = log(exp(log_Phi(-z1)) - exp(log_c + log_Phi(z2)))
  //           = log_diff_exp(std_normal_lcdf(-z1), log_c + std_normal_lcdf(z2))
  // std_normal_lcdf(-z1) is used instead of the equivalent
  // std_normal_lccdf(z1), which underflows to -Inf for z1 > ~8.3
  real log_term1 = std_normal_lcdf(-z1 | );  // log(1 - Phi(z1)) = log(Phi(-z1))
  real log_term2 = log_c + std_normal_lcdf(z2 | );  // log(exp(c) * Phi(z2))

  // log_diff_exp(a, b) = log(exp(a) - exp(b)), stable when a > b
  return log_diff_exp(log_term1, log_term2);
}

// log-PDF of the shifted Wald with uniform trial-to-trial variability in the
// non-decision time, onset convention NDT ~ uniform(ndt, ndt + sndt):
// f(t) = [S_W(t - ndt - sndt) - S_W(t - ndt)] / sndt (Miller et al. 2018, Eq. 6)
real swald_sndt_lpdf(real rt, real drift, real bound, real ndt, real sndt, real sigma) {
  if (sndt < 0) return negative_infinity();
  // the convolution is continuous at sndt = 0, so tiny sndt can use the plain
  // density (also the code path for the default fixed sndt = 0)
  if (sndt < 1e-8) return swald_lpdf(rt | drift, bound, ndt, sigma);
  if (rt - ndt <= 0) return negative_infinity();

  real ls_lo = swald_lccdf(rt | drift, bound, ndt + sndt, sigma);
  real ls_hi = swald_lccdf(rt | drift, bound, ndt, sigma);

  // for defective (negative-drift) accumulators both survivors converge to
  // the same positive constant in the deep tail, so their log-difference
  // drops below fp precision; the midpoint rule for the density is second
  // order in sndt and stable there. The strip rt - ndt <= sndt is excluded
  // because ls_lo = 0 exactly, making log_diff_exp accurate at any gap.
  if (ls_lo - ls_hi < 1e-8 && rt - ndt - sndt > 0) {
    return swald_lpdf(rt | drift, bound, ndt + sndt / 2, sigma);
  }
  return log_diff_exp(ls_lo, ls_hi) - log(sndt);
}

// log survivor of the shifted Wald + uniform NDT, for censored observations:
// S_conv(t) = [G(x1) - G(x2)] / sndt with x1 = t - ndt, x2 = x1 - sndt, and
// G(x) = int_0^x S_W(u) du = x * S_W(x) + M1(x), extended by G(x) = x for
// x <= 0 where S_W = 1 (this covers the strip ndt < t < ndt + sndt without a
// special case). M1(x) = int_0^x u f_W(u) du is the partial expectation of
// the (possibly defective) Wald. Computed in natural space: for negative
// drift the signs of mu = bound/drift and the Phi-bracket flip together and
// cancel, which a log-space form would need explicit sign tracking for.
real swald_sndt_lccdf(real rt, real drift, real bound, real ndt, real sndt, real sigma) {
  if (sndt < 0) return negative_infinity();
  if (sndt < 1e-8) return swald_lccdf(rt | drift, bound, ndt, sigma);

  real x1 = rt - ndt;
  if (x1 <= 0) return 0;

  real sigma_sq = square(sigma);
  vector[2] xs = [x1, x1 - sndt]';
  vector[2] g;
  for (k in 1:2) {
    real x = xs[k];
    if (x <= 0) {
      g[k] = x;
    } else {
      real sqrt_x = sqrt(x);
      real m1;
      if (abs(drift) < 1e-6 * sigma_sq / bound) {
        // mu = bound/drift diverges as drift -> 0 while the Phi-bracket
        // vanishes; use the exact drift = 0 limit of M1
        real w = bound / (sigma * sqrt_x);
        m1 = 2 * bound * (sqrt_x * exp(std_normal_lpdf(w)) / sigma
                          - (bound / sigma_sq) * Phi(-w));
      } else {
        real z1 = (drift * x - bound) / (sigma * sqrt_x);
        real z2 = -(drift * x + bound) / (sigma * sqrt_x);
        m1 = (bound / drift)
             * (Phi(z1) - exp(2 * drift * bound / sigma_sq + std_normal_lcdf(z2)));
      }
      g[k] = x * exp(swald_lccdf(x | drift, bound, 0, sigma)) + m1;
    }
  }

  real s_conv = (g[1] - g[2]) / sndt;
  // deep-tail cancellation of the G-difference: fall back to the midpoint
  // rule, second order in sndt and irrelevant at ~ -35 nats
  if (s_conv <= 1e-300) {
    return swald_lccdf(rt | drift, bound, ndt + sndt / 2, sigma);
  }
  return log(s_conv);
}
