// Normal drift truncated at d > 0 (rtdists posdrift); derivation by .lba_normal_log_M() in R/distributions.R

real lba_normal_g(real z) {
  return z * Phi(z) + exp(std_normal_lpdf(z |));
}

real lba_normal_h(real z) {
  return z * Phi(-z) - exp(std_normal_lpdf(z |));
}

// log M, assembled from signed terms in log space
real lba_normal_log_M(real t, real v, real b, real A, real s) {
  real z_lo = ((b - A) / t - v) / s;
  real z_hi = (b / t - v) / s;
  real delta = A / (t * s);
  if (delta < 1e-4) {
    // midpoint rule; u_m = v + s z_m written without its cancellation
    real z_m = 0.5 * (z_lo + z_hi);
    real u_m = (b - 0.5 * A) / t;
    real bracket = u_m + square(delta) / 24 * (u_m * (square(z_m) - 1) - 2 * s * z_m);
    if (bracket <= 0) return lba_log_floor();
    return log(delta) + std_normal_lpdf(z_m |) + log(bracket);
  }
  {
    real log_dPhi = lba_log_Phi_diff(z_lo, z_hi);
    real sum_z = z_lo + z_hi;
    // x = 0 guarded: log1m_exp(0) would put +Inf on the tape
    real x = 0.5 * delta * abs(sum_z);
    real log_dphi = x > 0
      ? log(s) + std_normal_lpdf((sum_z >= 0 ? z_lo : z_hi) |) + log1m_exp(-x)
      : negative_infinity();
    if (v > 0) {
      real l1 = log(v) + log_dPhi;
      if (sum_z >= 0) return log_sum_exp(l1, log_dphi);
      return l1 > log_dphi ? log_diff_exp(l1, log_dphi) : lba_log_floor();
    }
    if (v < 0) {
      // z_lo > 0 here, so the phi term is the larger one
      real l1 = log(-v) + log_dPhi;
      return log_dphi > l1 ? log_diff_exp(log_dphi, l1) : lba_log_floor();
    }
    return log_dphi;
  }
}

// Single-accumulator log-PDF for LBA with normal drift
real lba_normal_single_lpdf(real t, real v, real b, real A, real s) {
  return lba_normal_log_M(t, v, b, A, s) - log(A) - std_normal_lcdf(v / s |);
}

// Single-accumulator log-survival for LBA with normal drift, in probability space via exact-gradient Phi()
real lba_normal_single_lccdf(real t, real v, real b, real A, real s) {
  real z_lo = ((b - A) / t - v) / s;
  real z_hi = (b / t - v) / s;
  real delta = A / (t * s);
  real log_denom = std_normal_lcdf(v / s |);
  real num;
  if (delta < 1e-4) {
    // (1/delta) int Phi = Phi(z_m) - delta^2 / 24 z_m phi(z_m) + O(delta^4)
    real z_m = 0.5 * (z_lo + z_hi);
    real corr = square(delta) / 24 * z_m * exp(std_normal_lpdf(z_m |));
    num = (v >= 0 ? Phi(z_m) + expm1(log_denom) : exp(log_denom) - Phi(-z_m)) - corr;
  } else if (v >= 0) {
    num = (lba_normal_g(z_hi) - lba_normal_g(z_lo)) / delta + expm1(log_denom);
  } else {
    num = exp(log_denom) - (lba_normal_h(z_hi) - lba_normal_h(z_lo)) / delta;
  }
  return lba_log_clip(num) - log_denom;
}
