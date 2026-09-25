// Normal drift d ~ N(v, s^2) conditional on d > 0 (the posdrift convention of
// rtdists). In z-units, z = (u - v) / s, with delta = z_hi - z_lo = A / (t s):
//   M = int_{z_lo}^{z_hi} (v + s z) phi(z) dz = v [Phi(z_hi) - Phi(z_lo)] + s [phi(z_lo) - phi(z_hi)]
//   S - q = (1/delta) int_{z_lo}^{z_hi} [Phi(z) - Phi(-v/s)] dz,  q = Phi(-v/s)
// and the truncated density and survivor are M / (A Phi(v/s)) and
// (S - q) / Phi(v/s). Both Phi differences vanish with delta and both phi
// differences vanish with delta (z_lo + z_hi); the sum for M can also be a
// difference of two terms of the same size (v < 0, or v > 0 with the interval
// left of v). Every piece is therefore assembled in log space from a signed
// term, or by the midpoint rule below delta = 1e-4, where the direct form has
// lost four digits and the expansion has ~1e-12 left.

real lba_normal_g(real z) {
  return z * Phi(z) + exp(std_normal_lpdf(z |));
}

real lba_normal_h(real z) {
  return z * Phi(-z) - exp(std_normal_lpdf(z |));
}

// log M, see above
real lba_normal_log_M(real t, real v, real b, real A, real s) {
  real z_lo = ((b - A) / t - v) / s;
  real z_hi = (b / t - v) / s;
  real delta = A / (t * s);
  if (delta < 1e-4) {
    // midpoint rule with its second-order term; u_m = v + s z_m is the
    // midpoint of (lo, hi), written without the cancellation of v + s z_m
    real z_m = 0.5 * (z_lo + z_hi);
    real u_m = (b - 0.5 * A) / t;
    real bracket = u_m + square(delta) / 24 * (u_m * (square(z_m) - 1) - 2 * s * z_m);
    if (bracket <= 0) return lba_log_floor();
    return log(delta) + std_normal_lpdf(z_m |) + log(bracket);
  }
  {
    real log_dPhi = lba_log_Phi_diff(z_lo, z_hi);
    real sum_z = z_lo + z_hi;
    // phi(near) - phi(far) = phi(near) (1 - exp(-x)); the term is exactly
    // zero when z_lo = -z_hi, and log1m_exp(0) would put +Inf on the tape
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
      // z_lo > 0 here, so the phi term is positive and the larger one
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

// Single-accumulator log-survival for LBA with normal drift. S - q is formed
// in probability space from the exact-gradient Phi(): from the lower tails
// through g(z) = z Phi(z) + phi(z) when v >= 0 (q < 1/2 is the small term),
// from the upper tails through h(z) = z Phi(-z) - phi(z) when v < 0 (1 - q
// and F are the small terms). The remaining cancellation is intrinsic to the
// closed form and bounded by (S - q) / S, which stays above 1e-5 for any
// decision time under 30 s with gap >= 1e-3.
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
