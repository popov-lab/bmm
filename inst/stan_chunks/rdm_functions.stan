real rdm_log_normal_cdf_diff(real z_lo, real z_hi) {
  if (z_hi <= z_lo) return negative_infinity();
  if (z_lo >= 0) {
    return log_diff_exp(
      std_normal_lccdf(z_lo | ),
      std_normal_lccdf(z_hi | )
    );
  }
  if (z_hi <= 0) {
    return log_diff_exp(
      std_normal_lcdf(z_hi | ),
      std_normal_lcdf(z_lo | )
    );
  }
  return log(Phi(z_hi) - Phi(z_lo));
}

real rdm_full_pdf_raw(real t, real drift, real bound, real A, real s) {
  real sqrt_t = sqrt(t);
  real s_sqrt_t = s * sqrt_t;
  real alpha = (bound - A - t * drift) / s_sqrt_t;
  real beta = (bound - t * drift) / s_sqrt_t;

  return (1.0 / A) * (
    -drift * Phi(alpha) +
      (s / sqrt_t) * exp(std_normal_lpdf(alpha)) +
      drift * Phi(beta) -
      (s / sqrt_t) * exp(std_normal_lpdf(beta))
  );
}

real rdm_log_full_pdf(real t, real drift, real bound, real A, real s) {
  real sqrt_t = sqrt(t);
  real s_sqrt_t = s * sqrt_t;
  real alpha = (bound - A - t * drift) / s_sqrt_t;
  real beta = (bound - t * drift) / s_sqrt_t;

  // Near ndt the raw-scale expression underflows because both antiderivative
  // terms are extremely close to the same upper limit. In the positive-tail
  // region we can rewrite the difference as a sum of two tiny positive terms.
  if (alpha > 0) {
    real log_term1 = log(drift) + rdm_log_normal_cdf_diff(alpha, beta);
    real log_term2 = log(s / sqrt_t) +
      log_diff_exp(std_normal_lpdf(alpha), std_normal_lpdf(beta));
    return log_sum_exp(log_term1, log_term2) - log(A);
  }

  {
    real pdf = rdm_full_pdf_raw(t, drift, bound, A, s);
    if (pdf <= 0) return negative_infinity();
    return log(pdf);
  }
}

real rdm_full_cdf_raw(real t, real drift, real bound, real A, real s) {
  real s2 = square(s);
  real sqrt_t = sqrt(t);
  real s_sqrt_t = s * sqrt_t;
  real alpha1 = (drift * t - bound) / s_sqrt_t;
  real alpha2 = (drift * t - (bound - A)) / s_sqrt_t;
  real beta1 = -(drift * t + bound) / s_sqrt_t;
  real beta2 = -(drift * t + (bound - A)) / s_sqrt_t;

  return (1.0 / (2 * drift * A)) * (Phi(alpha2) - Phi(alpha1)) +
    (s * sqrt_t / A) * (
      alpha2 * Phi(alpha2) - alpha1 * Phi(alpha1) +
        exp(std_normal_lpdf(alpha2)) - exp(std_normal_lpdf(alpha1))
    ) -
    (1.0 / (2 * drift * A)) * (
      exp(2 * drift * (bound - A) / s2) * Phi(beta2) -
        exp(2 * drift * bound / s2) * Phi(beta1)
    );
}

real rdm_full_surv_antiderivative(real u, real t, real drift, real s) {
  real s2 = square(s);
  real s_sqrt_t = s * sqrt(t);
  real drift_t = drift * t;
  real q = s2 / (2 * drift);
  real y = (u - drift_t) / s_sqrt_t;
  real z = -(drift_t + u) / s_sqrt_t;

  return (u - drift_t - q) * Phi(y) +
    s_sqrt_t * exp(std_normal_lpdf(y)) -
    q * exp(2 * drift * u / s2) * Phi(z);
}

real rdm_full_surv_raw(real t, real drift, real bound, real A, real s) {
  return (
    rdm_full_surv_antiderivative(bound, t, drift, s) -
      rdm_full_surv_antiderivative(bound - A, t, drift, s)
  ) / A;
}

// Log-survival policy: the CDF form is well-conditioned when the CDF is small
// (log1m(cdf)), the direct survival antiderivative when it is large. A CDF
// dipping below 0 is cancellation noise near the support boundary where the
// true CDF ~ 0, so survival ~ 1 (log = 0); a non-positive direct survival in
// the deep tail means the true survival underflowed (log = -inf).
real rdm_log_full_surv(real t, real drift, real bound, real A, real s) {
  real cdf = rdm_full_cdf_raw(t, drift, bound, A, s);

  if (cdf <= 0) {
    return 0;
  }
  if (cdf < 0.5) {
    return log1m(cdf);
  }

  {
    real surv = rdm_full_surv_raw(t, drift, bound, A, s);
    if (surv > 0 && surv < 1) {
      return log(surv);
    }
  }
  return negative_infinity();
}

real rdm_log_lik_one(real rt, array[] real drift, real gap, real ndt,
                     real s, real sp, int response, array[] int n,
                     int use_start_var) {
  int n_cats = size(drift);
  real b = gap + sp;
  real t = rt - ndt;
  real lp;

  if (t <= 0) return negative_infinity();

  if (use_start_var == 0) {
    lp = log(n[response]) + swald_lpdf(rt | drift[response], b, ndt, s);
  } else {
    lp = log(n[response]) + rdm_log_full_pdf(t, drift[response], b, sp, s);
  }

  // single-pass race: the winner contributes n_win - 1 survival copies, each
  // loser n_j copies. The strict reps > 0 guard skips the winner's survival
  // when its category has one accumulator (all correct trials; every trial
  // for K = 2) and keeps zero-count or underflowed terms out of the sum
  // (0 * -inf would otherwise poison the likelihood with NaN).
  for (j in 1:n_cats) {
    int reps = (j == response) ? n[j] - 1 : n[j];
    if (reps > 0) {
      real log_surv;
      if (use_start_var == 0) {
        log_surv = swald_lccdf(rt | drift[j], b, ndt, s);
      } else {
        log_surv = rdm_log_full_surv(t, drift[j], b, sp, s);
      }
      lp += reps * log_surv;
    }
  }

  return lp;
}
