// LBA helpers shared by all drift distributions; derivation by the R twins in R/distributions.R
// No branch may return -Inf with an infinite partial: Stan keeps the adjoint of a discarded branch

// Last-resort floor for a numerator that has rounded to zero or below
real lba_log_clip(real x) {
  return log(fmax(x, 1e-300));
}

real lba_log_floor() {
  return log(1e-300);
}

// log(Phi(hi) - Phi(lo)), hi > lo, via std_normal_lcdf(-z) in the right half, never lccdf
real lba_log_Phi_diff(real lo, real hi) {
  if (lo + hi >= 0) {
    return log_diff_exp(std_normal_lcdf(-lo |), std_normal_lcdf(-hi |));
  }
  return log_diff_exp(std_normal_lcdf(hi |), std_normal_lcdf(lo |));
}

// log int phi(z) dz over a narrow interval: midpoint rule with its second-order term
real lba_log_phi_int_narrow(real z_lo, real z_hi) {
  real dz = z_hi - z_lo;
  real z_m = 0.5 * (z_lo + z_hi);
  return log(dz) + std_normal_lpdf(z_m |) + log1p(square(dz) / 24 * (square(z_m) - 1));
}
