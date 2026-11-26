// log density of one response time under the competing-risks shifted Wald model with variable non-decision time
real cswald_crisk_lpdf(real rt, real mu, real drift, real bound, real ndt,
                       real zr, real s, real sndt, int response) {
  real bound_upper = bound - zr*bound;
  real bound_lower = zr*bound;

  if (response == 1) {
    return swald_sndt_lpdf(rt | drift, bound_upper, ndt, sndt, s)
           + swald_sndt_lccdf(rt | -drift, bound_lower, ndt, sndt, s);
  } else {
    return swald_sndt_lpdf(rt | -drift, bound_lower, ndt, sndt, s)
           + swald_sndt_lccdf(rt | drift, bound_upper, ndt, sndt, s);
  }
}

// the same log density summed over all observations, valid only while sndt is fixed at 0
real cswald_crisk_lpdf(vector rt, vector mu, vector drift, vector bound,
                       vector ndt, vector zr, vector s, vector sndt,
                       array[] int dec) {
  int N = rows(rt);
  vector[N] t = rt - ndt;
  if (min(t) <= 0) return negative_infinity();

  vector[N] w = to_vector(dec);
  vector[N] bound_upper = bound - zr .* bound;
  vector[N] bound_lower = zr .* bound;
  vector[N] drift_win = (2 * w - 1) .* drift;
  vector[N] bound_win = w .* bound_upper + (1 - w) .* bound_lower;
  vector[N] bound_lose = bound - bound_win;

  return sum(swald_log_dens_vec(t, drift_win, bound_win, s))
         + sum(swald_log_surv_vec(t, -drift_win, bound_lose, s));
}
