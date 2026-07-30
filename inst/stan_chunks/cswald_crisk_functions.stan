// log-PDF of competing risks shifted Wald model. With sndt > 0 the two
// accumulators receive independent non-decision-time draws (a race between
// total finishing times); choice probabilities are unaffected and the density
// differs from the shared-NDT convolution by less than ~2% even at sndt = 0.3.
real cswald_crisk_lpdf(real rt, real mu, real drift, real bound, real ndt, real zr, real s, real sndt, int response) {
  // compute bounds for upper and lower response
  real bound_upper = bound - zr*bound;
  real bound_lower = zr*bound;

  // compute lpdf dependent on response type
  if (response == 1) {
    return swald_sndt_lpdf(rt | drift, bound_upper, ndt, sndt, s) + swald_sndt_lccdf(rt | -drift, bound_lower, ndt, sndt, s);
  } else {
    return swald_sndt_lpdf(rt | -drift, bound_lower, ndt, sndt, s) + swald_sndt_lccdf(rt | drift, bound_upper, ndt, sndt, s);
  }
}

// vectorized overload used as a loop = FALSE family when sndt is fixed at 0:
// returns the summed log-likelihood over all observations. The winning
// accumulator's density is a fully vectorized closed-form shifted Wald
// (no CDF calls); the losing accumulator's survivor vectorizes the z-score
// preparation and keeps a scalar loop only for the std_normal_lcdf pairs
// (Stan has no elementwise log-CDF). With sndt > 0 the convolution is
// lccdf-bound and the scalar per-observation path is faster, so the general
// path delegates to the scalar helpers (and configure_model keeps
// loop = TRUE in that case)
real cswald_crisk_lpdf(vector rt, vector mu, vector drift, vector bound,
                       vector ndt, vector zr, vector s, vector sndt,
                       array[] int dec) {
  int N = rows(rt);
  real lp = 0;
  if (min(sndt) < 0) return negative_infinity();
  if (max(sndt) >= 1e-8) {
    for (n in 1:N) {
      real bu = bound[n] - zr[n] * bound[n];
      real bl = zr[n] * bound[n];
      if (dec[n] == 1) {
        lp += swald_sndt_lpdf(rt[n] | drift[n], bu, ndt[n], sndt[n], s[n])
              + swald_sndt_lccdf(rt[n] | -drift[n], bl, ndt[n], sndt[n], s[n]);
      } else {
        lp += swald_sndt_lpdf(rt[n] | -drift[n], bl, ndt[n], sndt[n], s[n])
              + swald_sndt_lccdf(rt[n] | drift[n], bu, ndt[n], sndt[n], s[n]);
      }
    }
    return lp;
  }
  // both accumulators share rt - ndt, so a single rt <= ndt makes the
  // winner's density (and thus the summed target) -Inf
  vector[N] t = rt - ndt;
  if (min(t) <= 0) return negative_infinity();
  // winner = the accumulator matching the decision (drift toward its bound),
  // loser = the opposite accumulator with mirrored drift; selecting via the
  // 0/1 data vector w keeps everything vectorized
  vector[N] w = to_vector(dec);
  vector[N] bu = bound - zr .* bound;
  vector[N] bl = zr .* bound;
  vector[N] dw = (2 * w - 1) .* drift;
  vector[N] bw = w .* bu + (1 - w) .* bl;
  vector[N] bL = bound - bw;
  vector[N] sg_sq = square(s);
  vector[N] resid = bw - dw .* t;
  lp += sum(log(bw))
        - 0.5 * (N * log(2 * pi()) + 2 * sum(log(s)) + 3 * sum(log(t)))
        - 0.5 * sum(square(resid) ./ (sg_sq .* t));
  vector[N] denom = s .* sqrt(t);
  vector[N] dxt = -(dw .* t);
  vector[N] z1 = (dxt - bL) ./ denom;
  vector[N] z2 = -(dxt + bL) ./ denom;
  vector[N] log_c = -2 * (bL .* dw) ./ sg_sq;
  for (n in 1:N) {
    lp += swald_log_diff_exp(std_normal_lcdf(-z1[n] | ),
                             log_c[n] + std_normal_lcdf(z2[n] | ));
  }
  return lp;
}
