// Racing diffusion model (Tillman et al., 2020); R twins .dwald_full() / .pwald_full() in R/distributions.R

// log Phi(z) via erfc; upper tails as lcdf(-z), never the lccdf of z (-inf from z = 8.26 on Stan Math 5.3)
real rdm_log_Phi(real z) {
  if (z < -37.5) return std_normal_lcdf(z | );
  return log(0.5 * erfc(-0.7071067811865475 * z));
}

// log(Phi(b) - Phi(a)) for b > a; log-differences go through swald_log_diff_exp() to keep the gradient finite
real rdm_log_Phi_diff(real a, real b) {
  if (b <= a) return negative_infinity();
  if (a >= 0) return swald_log_diff_exp(rdm_log_Phi(-a), rdm_log_Phi(-b));
  if (b <= 0) return swald_log_diff_exp(rdm_log_Phi(b), rdm_log_Phi(a));
  // straddling: both O(1), exact directly; the midpoint fallback keeps b - a >= 1e-4
  return log(Phi(b) - Phi(a));
}

// log(y Phi(y) + phi(y)) given log Phi(y): direct above -1, log space to -12, Mills-ratio series below
real rdm_log_g(real y, real log_Phi_y) {
  if (y < -12) {
    real y2 = square(y);
    real term = 1;
    real total = 0;
    for (j in 1:12) {
      term *= (2.0 * j - 1) / y2;
      total += (j % 2 == 1) ? term : -term;
    }
    return -0.5 * y2 - 0.9189385332046727 + log(total);
  }
  real log_phi_y = -0.5 * square(y) - 0.9189385332046727;
  if (y < -1) return swald_log_diff_exp(log_phi_y, log(-y) + log_Phi_y);
  return log(y * exp(log_Phi_y) + exp(log_phi_y));
}

// log density at decision time t > 0 (Tillman et al., 2020, Eq. 5); midpoint plain Wald for small A
real rdm_log_pdf(real t, real drift, real gap, real A, real s) {
  real st = sqrt(t);
  if (A < 1e-4 * s * st) return swald_lpdf(t | drift, gap + 0.5 * A, 0, s);

  real denom = s * st;
  real alpha = (gap - drift * t) / denom;
  real beta = (gap + A - drift * t) / denom;
  real l1 = log(drift) + rdm_log_Phi_diff(alpha, beta);
  real lpa = -0.5 * square(alpha) - 0.9189385332046727;
  real lpb = -0.5 * square(beta) - 0.9189385332046727;
  real lnum;
  if (abs(alpha) <= abs(beta)) {
    lnum = log_sum_exp(l1, swald_log_diff_exp(lpa, lpb) + log(s / st));
  } else {
    lnum = swald_log_diff_exp(l1, swald_log_diff_exp(lpb, lpa) + log(s / st));
  }
  return lnum - log(A);
}

// log survivor at decision time t > 0 (Tillman et al., 2020, Appendix A); midpoint plain Wald for small A
real rdm_log_surv(real t, real drift, real gap, real A, real s) {
  real st = sqrt(t);
  if (A < 1e-4 * s * st) return swald_lccdf(t | drift, gap + 0.5 * A, 0, s);

  real denom = s * st;
  real b = gap + A;
  real alpha = (gap - drift * t) / denom;
  real beta = (b - drift * t) / denom;
  real lPa = rdm_log_Phi(alpha);
  real lPb = rdm_log_Phi(beta);
  real lD1 = log(denom)
             + swald_log_diff_exp(rdm_log_g(beta, lPb), rdm_log_g(alpha, lPa));

  real s2 = square(s);
  real lEb = 2 * drift * b / s2 + rdm_log_Phi(-(b + drift * t) / denom);
  real lEk = 2 * drift * gap / s2 + rdm_log_Phi(-(gap + drift * t) / denom);
  // the lower-half difference reuses the two log Phi already in hand
  real lP = (beta <= 0) ? swald_log_diff_exp(lPb, lPa) : rdm_log_Phi_diff(alpha, beta);
  real lD2;
  if (lEb >= lEk) {
    lD2 = log_sum_exp(lP, swald_log_diff_exp(lEb, lEk));
  } else {
    lD2 = swald_log_diff_exp(log_sum_exp(lP, lEb), lEk);
  }
  lD2 += 2 * log(s) - log(2 * drift);

  return fmin(swald_log_diff_exp(lD1, lD2) - log(A), 0);
}
