// log-PDF of competing risks shifted Wald model
real cswald_crisk_lpdf(real rt, real mu, real drift, real bound, real ndt, real zr, real s, int response) {
  // compute bounds for upper and lower response
  real bound_upper = bound - zr*bound;
  real bound_lower = zr*bound;

  // compute lpdf dependent on response type
  if (response == 1) {
    return swald_lpdf(rt | drift, bound_upper, ndt, s) + swald_lccdf(rt | -drift, bound_lower, ndt, s);
  } else {
    return swald_lpdf(rt | -drift, bound_lower, ndt, s) + swald_lccdf(rt | drift, bound_upper, ndt, s);
  }
}

// vectorized overload used by the loop = FALSE family: returns the summed
// log-likelihood, winner's density plus loser's survivor per observation
real cswald_crisk_lpdf(vector rt, vector mu, vector drift, vector bound,
                       vector ndt, vector zr, vector s, array[] int dec) {
  int N = rows(rt);
  // both accumulators share rt - ndt, so a single rt <= ndt makes the winner's
  // density (and thus the summed target) -inf
  vector[N] t = rt - ndt;
  if (min(t) <= 0) return negative_infinity();

  // winner = the accumulator matching the decision (drift toward its bound),
  // loser = the opposite accumulator with mirrored drift; selecting via the
  // 0/1 data vector w keeps everything vectorized
  vector[N] w = to_vector(dec);
  vector[N] bound_upper = bound - zr .* bound;
  vector[N] bound_lower = zr .* bound;
  vector[N] drift_win = (2 * w - 1) .* drift;
  vector[N] bound_win = w .* bound_upper + (1 - w) .* bound_lower;
  vector[N] bound_lose = bound - bound_win;

  return sum(swald_log_dens_vec(t, drift_win, bound_win, s))
         + sum(swald_log_surv_vec(t, -drift_win, bound_lose, s));
}
