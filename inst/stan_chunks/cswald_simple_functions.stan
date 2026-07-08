// log-PDF of the censored shifted Wald model
real cswald_lpdf(real rt, real mu, real drift, real bound, real ndt, real s, real st0, int response) {
  if (response == 1) {
    return swald_st0_lpdf(rt | drift, bound, ndt, st0, s);
  } else {
    return swald_st0_lccdf(rt | drift, bound, ndt, st0, s);
  }
}
