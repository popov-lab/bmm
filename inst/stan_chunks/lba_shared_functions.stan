// LBA helpers shared by all drift distributions; derivation by the R twins in R/distributions.R
// No branch may return -Inf with an infinite partial: Stan keeps the adjoint of a discarded branch

// Last-resort floor for a numerator that has rounded to zero or below
real lba_log_clip(real x) {
  return log(fmax(x, 1e-300));
}

real lba_log_floor() {
  return log(1e-300);
}

// log Phi(z) with the exact gradient: erfc, not the approximate-derivative std_normal_lcdf; the upper
// tail via log1m so that expm1() of the result keeps 1 - Phi(z) at full relative precision
// (erfc(26.5) ~ 1e-306: beyond z = -37 fall back to the lcdf's asymptotic branch)
real lba_log_Phi(real z) {
  if (z >= 0) return log1m(0.5 * erfc(0.7071067811865476 * z));
  if (z > -37) return log(0.5 * erfc(-0.7071067811865476 * z));
  return std_normal_lcdf(z |);
}

// log(Phi(hi) - Phi(lo)), hi > lo, from whichever tail keeps both probabilities away from 1
real lba_log_Phi_diff(real lo, real hi) {
  if (lo + hi >= 0) return log_diff_exp(lba_log_Phi(-lo), lba_log_Phi(-hi));
  return log_diff_exp(lba_log_Phi(hi), lba_log_Phi(lo));
}

// log int phi(z) dz over a narrow interval: midpoint rule with its second-order term
real lba_log_phi_int_narrow(real z_lo, real z_hi) {
  real dz = z_hi - z_lo;
  real z_m = 0.5 * (z_lo + z_hi);
  return log(dz) + std_normal_lpdf(z_m |) + log1p(square(dz) / 24 * (square(z_m) - 1));
}
