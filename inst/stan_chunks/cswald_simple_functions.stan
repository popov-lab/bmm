// log-PDF of the censored shifted Wald model
real cswald_lpdf(real rt, real mu, real drift, real bound, real ndt, real s, int response) {
  if (response == 1) {
    return swald_lpdf(rt | drift, bound, ndt, s);
  } else {
    return swald_lccdf(rt | drift, bound, ndt, s);
  }
}
