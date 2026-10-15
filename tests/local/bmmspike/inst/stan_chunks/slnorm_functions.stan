// shifted log-normal density, identical to the one brms uses for
// shifted_lognormal(): the shift ndt is bounded above by min(Y) in tdata
real slnorm_lpdf(real y, real mu, real sigma, real ndt) {
  return lognormal_lpdf(y - ndt | mu, sigma);
}
