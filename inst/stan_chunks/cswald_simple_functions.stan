// log-PDF of the censored shifted Wald model
real cswald_lpdf(real rt, real mu, real drift, real bound, real ndt, real s, real sndt, int response) {
  if (response == 1) {
    return swald_sndt_lpdf(rt | drift, bound, ndt, sndt, s);
  } else {
    return swald_sndt_lccdf(rt | drift, bound, ndt, sndt, s);
  }
}
