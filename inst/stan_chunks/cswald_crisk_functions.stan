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
