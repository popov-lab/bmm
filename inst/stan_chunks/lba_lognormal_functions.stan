// Lognormal drift d ~ LN(v, s), so z = (log u - v) / s and F(u) = Phi(z). The
// width of the interval in z-units is dz = log(b / (b - A)) / s, taken from
// log1m(A / b) rather than as a difference of two logs so that a tiny A keeps
// its digits. M = exp(v + s^2 / 2) [Phi(z_hi - s) - Phi(z_lo - s)].

real lba_lognormal_log_dPhi(real z_lo, real z_hi, real dz) {
  if (dz < 1e-4) return lba_log_phi_int_narrow(z_lo, z_hi);
  return lba_log_Phi_diff(z_lo, z_hi);
}

// Single-accumulator log-PDF for LBA with lognormal drift
real lba_lognormal_single_lpdf(real t, real v, real b, real A, real s) {
  real z_hi = (log(b / t) - v) / s;
  real dz = -log1m(A / b) / s;
  return v + 0.5 * square(s) + lba_lognormal_log_dPhi(z_hi - dz - s, z_hi - s, dz) - log(A);
}

// Single-accumulator log-survival for LBA with lognormal drift
real lba_lognormal_single_lccdf(real t, real v, real b, real A, real s) {
  real z_hi = (log(b / t) - v) / s;
  real dz = -log1m(A / b) / s;
  real z_lo = z_hi - dz;
  real log_tM = log(t) + v + 0.5 * square(s) + lba_lognormal_log_dPhi(z_lo - s, z_hi - s, dz);
  real log_u = log_sum_exp(log(A) + std_normal_lcdf(z_lo |), log(b) + lba_lognormal_log_dPhi(z_lo, z_hi, dz));
  return (log_u > log_tM ? log_diff_exp(log_u, log_tM) : lba_log_floor()) - log(A);
}
