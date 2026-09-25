// Helpers shared by the four LBA drift distributions. One accumulator with
// threshold b, start point k ~ U(0, A) and drift d finishes at (b - k) / d, so
// at decision time t the drifts that finish exactly then are u = (b - k) / t
// in (lo, hi) = ((b - A) / t, b / t). The density is M / A with
// M = int_{lo}^{hi} u f(u) du, and the survivor is (u_num - t M) / A with
// u_num = b F(hi) - (b - A) F(lo) = A F(lo) + b (F(hi) - F(lo)), the second
// form being a sum of positive terms and the one used here.
//
// The differences of CDFs below are taken from the tail where both terms are
// small, and where an interval is too narrow for a difference to carry
// digits the integral is replaced by the midpoint rule with its second-order
// term. No branch evaluates an expression that can return -Inf with an
// infinite partial: Stan keeps the adjoint of a discarded branch, so
// log_diff_exp(x, x) or log(0) poison the gradient of the whole trial even
// when the value is not used.
//
// cogmod 0.3.2 (Makowski) exposed the Phi saturation and the -690 floor of the
// previous normal-drift kernel; the layouts here were derived independently.

// Last-resort floor for a numerator that has rounded to zero or below. On the
// verification grid (t - ndt down to 1 ms, |z| to 40, A from 1e-8 to 2, gap
// from 1e-3 to 1.5) no branch reaches it; it is left for parameter values
// beyond that range (v / s below -37, or a decision time so long that the
// truncated survivor falls under 1e-300).
real lba_log_clip(real x) {
  return log(fmax(x, 1e-300));
}

real lba_log_floor() {
  return log(1e-300);
}

// log(Phi(hi) - Phi(lo)) for hi > lo, from the upper tails when the pair sits
// in the right half. Never -Inf for finite arguments: std_normal_lcdf stays
// finite where log(Phi()) and std_normal_lccdf (Stan Math 5.3) do not.
real lba_log_Phi_diff(real lo, real hi) {
  if (lo + hi >= 0) {
    return log_diff_exp(std_normal_lcdf(-lo |), std_normal_lcdf(-hi |));
  }
  return log_diff_exp(std_normal_lcdf(hi |), std_normal_lcdf(lo |));
}

// int_{z_lo}^{z_hi} phi(z) dz for a narrow interval, by the midpoint rule with
// its second-order term (relative error O((dz max(|z|, 1))^4 / 1920))
real lba_log_phi_int_narrow(real z_lo, real z_hi) {
  real dz = z_hi - z_lo;
  real z_m = 0.5 * (z_lo + z_hi);
  return log(dz) + std_normal_lpdf(z_m |) + log1p(square(dz) / 24 * (square(z_m) - 1));
}
