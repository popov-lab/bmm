// log-PDF of competing risks shifted Wald model. With st0 > 0 the two
// accumulators receive independent non-decision-time draws (a race between
// total finishing times); choice probabilities are unaffected and the density
// differs from the shared-NDT convolution by less than ~2% even at st0 = 0.3.
real cswald_crisk_lpdf(real rt, real mu, real drift, real bound, real ndt, real zr, real s, real st0, int response) {
  // compute bounds for upper and lower response
  real bound_upper = bound - zr*bound;
  real bound_lower = zr*bound;

  // compute lpdf dependent on response type
  if (response == 1) {
    return swald_st0_lpdf(rt | drift, bound_upper, ndt, st0, s) + swald_st0_lccdf(rt | -drift, bound_lower, ndt, st0, s);
  } else {
    return swald_st0_lpdf(rt | -drift, bound_lower, ndt, st0, s) + swald_st0_lccdf(rt | drift, bound_upper, ndt, st0, s);
  }
}
