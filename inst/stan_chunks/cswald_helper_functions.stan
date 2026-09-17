// log(exp(a) - exp(b)) that returns -inf instead of NaN when rounding makes
// b >= a, so one underflowing survivor cannot poison the whole gradient
real swald_log_diff_exp(real a, real b) {
  if (b >= a) return negative_infinity();
  return log_diff_exp(a, b);
}

// log(Phi(z)) with an exact derivative: Phi's derivative is the normal density,
// whereas std_normal_lcdf's gradient is an approximation (relative error up to
// ~1e-4). The lcdf is used only where Phi rounds to 0 (z < ~-37.5)
real swald_log_Phi(real z) {
  real p = Phi(z);
  return p > 0 ? log(p) : std_normal_lcdf(z | );
}

// log shifted Wald survivor S = Phi(-z1) - exp(log_c) * Phi(z2), with
// z1 = (drift*t - bound) / (sigma*sqrt(t)), z2 = -(drift*t + bound) / (sigma*sqrt(t)),
// log_c = 2*bound*drift/sigma^2. Computed in probability space for the exact
// derivative (see swald_log_Phi); exp(log_c)*Phi(z2) <= Phi(-z1) <= 1, so the
// exp cannot overflow. The log-space form takes over when the difference
// underflows or rounds to <= 0
real swald_log_surv(real z1, real z2, real log_c) {
  real surv = Phi(-z1) - exp(log_c + swald_log_Phi(z2));
  if (surv <= 0) {
    return swald_log_diff_exp(std_normal_lcdf(-z1 | ),
                              log_c + std_normal_lcdf(z2 | ));
  }
  return log(surv);
}

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
real swald_lccdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return 0;  // process hasn't started, survival = 1

  real sigma_sqrt_t = sigma * sqrt(t_shifted);
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;
  real log_c = 2 * bound * drift / square(sigma);

  return swald_log_surv(z1, z2, log_c);
}

// Vectorized forms of swald_lpdf / swald_lccdf for the loop = FALSE families.
// Both take the already-shifted time t = rt - ndt. Keep them in step with the
// scalar forms above and with .dwald() / .pwald() in R/distributions.R

// t <= 0 anywhere makes the summed target -inf, exactly as the scalar sum would
vector swald_log_dens_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
  if (min(t) <= 0) return rep_vector(negative_infinity(), n);
  vector[n] residual = bound - drift .* t;
  return log(bound)
         - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log(t))
         - 0.5 * square(residual) ./ (square(sigma) .* t);
}

// requires t > 0 elementwise. Elements whose probability-space survivor
// underflows are recomputed by the scalar swald_log_surv below. The fmax/fmin
// clamps are required: without them a log(0) or exp(inf) in a discarded
// element injects NaN adjoints into the shared parameters even though its
// value is overwritten. Every clamped element is one the loop recomputes
vector swald_log_surv_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
  vector[n] denom = sigma .* sqrt(t);
  vector[n] dxt = drift .* t;
  vector[n] z1 = (dxt - bound) ./ denom;
  vector[n] z2 = -(dxt + bound) ./ denom;
  vector[n] log_c = 2 * (bound .* drift) ./ square(sigma);
  vector[n] p2 = Phi(z2);
  vector[n] surv = Phi(-z1) - exp(fmin(log_c + log(fmax(p2, 1e-300)), 0));
  vector[n] out = log(fmax(surv, 1e-300));
  for (k in 1:n) {
    if (p2[k] < 1e-300 || surv[k] < 1e-300) {
      out[k] = swald_log_surv(z1[k], z2[k], log_c[k]);
    }
  }
  return out;
}

// log-CDF of the shifted Wald, F_W(t) = Phi(z1) + exp(c) * Phi(z2), as a
// log-space sum so that it stays finite deep in the left tail
real swald_lcdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return negative_infinity();

  real sigma_sqrt_t = sigma * sqrt(t_shifted);
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;
  real log_c = 2 * bound * drift / square(sigma);

  return log_sum_exp(swald_log_Phi(z1), log_c + swald_log_Phi(z2));
}

// Integrated survivor G(x) = int_0^x S_W(u) du = x * S_W(x) + M1(x), where
// M1(x) = int_0^x u f_W(u) du is the partial expectation of the (possibly
// defective) Wald; G(x) = x for x <= 0. Mirrors .gwald() in R/distributions.R
real swald_gint(real x, real drift, real bound, real sigma) {
  if (x <= 0) return x;

  real sigma_sq = square(sigma);
  real sqrt_x = sqrt(x);
  real dx = drift * x;
  real z1 = (dx - bound) / (sigma * sqrt_x);
  real z2 = -(dx + bound) / (sigma * sqrt_x);
  // q = exp(c) * Phi(z2) <= F_W(x) <= 1 mathematically, so the exp cannot
  // overflow (same argument as in swald_log_surv)
  real q = exp(2 * bound * drift / sigma_sq + swald_log_Phi(z2));
  real m1;

  if (abs(drift) < 1e-6 * sigma_sq / bound) {
    // mu = bound / drift diverges as drift -> 0 while the Phi-bracket vanishes;
    // use the exact drift = 0 limit of M1 instead of the 0 * inf cancellation
    real w = bound / (sigma * sqrt_x);
    m1 = 2 * bound * (sqrt_x * exp(std_normal_lpdf(w | )) / sigma
                      - (bound / sigma_sq) * Phi(-w));
  } else {
    m1 = (bound / drift) * (Phi(z1) - q);
  }

  return x * exp(swald_lccdf(x | drift, bound, 0, sigma)) + m1;
}

// log of (1/sndt) * int_{x-sndt}^{x} S_W(u) du by composite Simpson on four
// intervals in log space. Mirrors .simpson_log_mean() in R/distributions.R
real swald_log_surv_mean(real x, real drift, real bound, real sigma, real sndt) {
  vector[5] weights = log(to_vector({1, 4, 2, 4, 1}) / 12);
  vector[5] terms;
  for (k in 1:5) {
    terms[k] = swald_lccdf(x - sndt + (k - 1) * sndt / 4 | drift, bound, 0, sigma)
               + weights[k];
  }
  return log_sum_exp(terms);
}

// log-PDF of the shifted Wald with uniform trial-to-trial variability in the
// non-decision time, NDT ~ uniform(ndt, ndt + sndt):
// f(t) = [S_W(t - ndt - sndt) - S_W(t - ndt)] / sndt (Miller et al. 2017, Eq. 6).
// The cutoffs 1e-8 must match .sndt_min and .cancellation_tol in
// R/distributions.R
real swald_sndt_lpdf(real rt, real drift, real bound, real ndt, real sndt, real sigma) {
  if (sndt < 0) return negative_infinity();
  // the convolution is continuous at sndt = 0, so tiny sndt takes the plain
  // density (this is also the path for the default fixed sndt = 0)
  if (sndt < 1e-8) return swald_lpdf(rt | drift, bound, ndt, sigma);

  real t1 = rt - ndt;
  if (t1 <= 0) return negative_infinity();

  // strip ndt < rt <= ndt + sndt: the earlier survivor is exactly 1, so the
  // density is F_W(t1) / sndt; the log-CDF stays finite where log(1 - S)
  // would underflow
  if (t1 <= sndt) return swald_lcdf(rt | drift, bound, ndt, sigma) - log(sndt);

  real surv_early = swald_lccdf(rt | drift, bound, ndt + sndt, sigma);
  real surv_late = swald_lccdf(rt | drift, bound, ndt, sigma);

  // for defective (negative-drift) accumulators both survivors converge to the
  // same constant in the deep tail and their difference cancels; the midpoint
  // rule (second order in sndt) is stable there
  if (surv_early - surv_late < 1e-8) {
    return swald_lpdf(rt | drift, bound, ndt + sndt / 2, sigma);
  }
  return swald_log_diff_exp(surv_early, surv_late) - log(sndt);
}

// log survivor of the shifted Wald + uniform NDT, for censored observations:
// S_conv(t) = [G(x1) - G(x1 - sndt)] / sndt with x1 = t - ndt. Once the
// difference cancels (relative guard) the same integral is recovered by
// Simpson quadrature in log space, whose terms are all positive
real swald_sndt_lccdf(real rt, real drift, real bound, real ndt, real sndt, real sigma) {
  if (sndt < 0) return negative_infinity();
  if (sndt < 1e-8) return swald_lccdf(rt | drift, bound, ndt, sigma);

  real x1 = rt - ndt;
  if (x1 <= 0) return 0;

  real g_hi = swald_gint(x1, drift, bound, sigma);
  real g_lo = swald_gint(x1 - sndt, drift, bound, sigma);
  real delta = g_hi - g_lo;

  if (delta <= 1e-8 * fmax(abs(g_hi), abs(g_lo))) {
    return swald_log_surv_mean(x1, drift, bound, sigma, sndt);
  }
  return log(delta / sndt);
}
