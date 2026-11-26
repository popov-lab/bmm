// log(exp(a) - exp(b)), returning -inf when b >= a
real swald_log_diff_exp(real a, real b) {
  if (b >= a) return negative_infinity();
  return log_diff_exp(a, b);
}

// log of the standard normal CDF, stable in the far lower tail
real swald_log_Phi(real z) {
  real p = Phi(z);
  return p >= 1e-300 ? log(p) : std_normal_lcdf(z | );
}

// log survivor function of the shifted Wald distribution from its standardised arguments
real swald_log_surv(real z1, real z2, real log_c) {
  real surv = Phi(-z1) - exp(log_c + swald_log_Phi(z2));
  if (surv < 1e-300) {
    return swald_log_diff_exp(std_normal_lcdf(-z1 | ),
                              log_c + std_normal_lcdf(z2 | ));
  }
  return log(surv);
}

// log density of the shifted Wald distribution
real swald_lpdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0 || bound <= 0) return negative_infinity();

  real sigma_sq = square(sigma);
  real log_t = log(t_shifted);

  real log_norm = log(bound) - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log_t);

  real residual = bound - drift * t_shifted;
  real log_kernel = -0.5 * square(residual) / (sigma_sq * t_shifted);

  return log_norm + log_kernel;
}

// log survivor function (1 - CDF) of the shifted Wald distribution
real swald_lccdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return 0;
  if (bound <= 0) return negative_infinity();

  real sigma_sqrt_t = sigma * sqrt(t_shifted);
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;
  real log_c = 2 * bound * drift / square(sigma);

  return swald_log_surv(z1, z2, log_c);
}

// shifted Wald log density for a vector of decision times t = rt - ndt
vector swald_log_dens_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
  if (min(t) <= 0 || min(bound) <= 0) return rep_vector(negative_infinity(), n);
  vector[n] residual = bound - drift .* t;
  return log(bound)
         - 0.5 * (log(2 * pi()) + 2 * log(sigma) + 3 * log(t))
         - 0.5 * square(residual) ./ (square(sigma) .* t);
}

// shifted Wald log survivor function for a vector of decision times t = rt - ndt
vector swald_log_surv_vec(vector t, vector drift, vector bound, vector sigma) {
  int n = rows(t);
  if (min(bound) <= 0) return rep_vector(negative_infinity(), n);
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

// log CDF of the shifted Wald distribution
real swald_lcdf(real rt, real drift, real bound, real ndt, real sigma) {
  real t_shifted = rt - ndt;
  if (t_shifted <= 0) return negative_infinity();

  real sigma_sqrt_t = sigma * sqrt(t_shifted);
  real z1 = (drift * t_shifted - bound) / sigma_sqrt_t;
  real z2 = -(drift * t_shifted + bound) / sigma_sqrt_t;
  real log_c = 2 * bound * drift / square(sigma);

  return log_sum_exp(swald_log_Phi(z1), log_c + swald_log_Phi(z2));
}

// integral of the survivor function of the Wald distribution from 0 to x
real swald_gint(real x, real drift, real bound, real sigma) {
  if (x <= 0) return x;

  real sigma_sq = square(sigma);
  real sqrt_x = sqrt(x);
  real dx = drift * x;
  real z1 = (dx - bound) / (sigma * sqrt_x);
  real z2 = -(dx + bound) / (sigma * sqrt_x);
  real q = exp(2 * bound * drift / sigma_sq + swald_log_Phi(z2));
  real m1;

  if (abs(drift) < 3e-8 * sigma_sq / bound) {
    real w = bound / (sigma * sqrt_x);
    m1 = 2 * bound * (sqrt_x * exp(std_normal_lpdf(w | )) / sigma
                      - (bound / sigma_sq) * Phi(-w));
  } else {
    m1 = (bound / drift) * (Phi(z1) - q);
  }

  return x * exp(swald_lccdf(x | drift, bound, 0, sigma)) + m1;
}

// log mean of the Wald survivor function over [x - sndt, x] by Simpson's rule
real swald_log_surv_mean(real x, real drift, real bound, real sigma, real sndt) {
  vector[5] weights = log(to_vector({1, 4, 2, 4, 1}) / 12);
  vector[5] terms;
  for (k in 1:5) {
    terms[k] = swald_lccdf(x - sndt + (k - 1) * sndt / 4 | drift, bound, 0, sigma)
               + weights[k];
  }
  return log_sum_exp(terms);
}

// log density of the shifted Wald distribution with non-decision time uniform on [ndt, ndt + sndt]
real swald_sndt_lpdf(real rt, real drift, real bound, real ndt, real sndt, real sigma) {
  if (sndt < 0) return negative_infinity();
  if (sndt < 1e-8) return swald_lpdf(rt | drift, bound, ndt, sigma);

  real t1 = rt - ndt;
  if (t1 <= 0) return negative_infinity();

  if (t1 <= sndt) return swald_lcdf(rt | drift, bound, ndt, sigma) - log(sndt);

  real surv_early = swald_lccdf(rt | drift, bound, ndt + sndt, sigma);
  real surv_late = swald_lccdf(rt | drift, bound, ndt, sigma);
  if (surv_early - surv_late >= 1e-8) {
    return swald_log_diff_exp(surv_early, surv_late) - log(sndt);
  }

  real cdf_late = swald_lcdf(rt | drift, bound, ndt, sigma);
  real cdf_early = swald_lcdf(rt | drift, bound, ndt + sndt, sigma);
  if (cdf_late - cdf_early >= 1e-8) {
    return swald_log_diff_exp(cdf_late, cdf_early) - log(sndt);
  }

  vector[5] weights = log(to_vector({1, 4, 2, 4, 1}) / 12);
  vector[5] terms;
  for (k in 1:5) {
    terms[k] = swald_lpdf(t1 - sndt + (k - 1) * sndt / 4 | drift, bound, 0, sigma)
               + weights[k];
  }
  return log_sum_exp(terms);
}

// log survivor function of the shifted Wald distribution with uniform non-decision time
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
