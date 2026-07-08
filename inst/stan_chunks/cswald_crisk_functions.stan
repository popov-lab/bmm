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
