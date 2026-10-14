  // x is the latent stick of a simplex group and inv_link(x) ~ Beta(a, b);
  // Beta(1, K - k) on stick k of K gives the group a symmetric Dirichlet(1)
  real logitbeta_lpdf(real x, real a, real b) {
    return a * log_inv_logit(x) + b * log1m_inv_logit(x) - lbeta(a, b);
  }
  real logitbeta_lpdf(vector x, real a, real b) {
    return a * sum(log_inv_logit(x)) + b * sum(log1m_inv_logit(x))
           - rows(x) * lbeta(a, b);
  }
  real logitbeta_rng(real a, real b) {
    return logit(beta_rng(a, b));
  }
  // lcdf and lccdf stay finite where log(Phi(x)) and log1m(Phi(x)) underflow
  real probitbeta_lpdf(real x, real a, real b) {
    return (a - 1) * std_normal_lcdf(x) + (b - 1) * std_normal_lccdf(x)
           + std_normal_lpdf(x) - lbeta(a, b);
  }
  real probitbeta_lpdf(vector x, real a, real b) {
    real lp = std_normal_lpdf(x) - rows(x) * lbeta(a, b);
    for (n in 1:rows(x)) {
      lp += (a - 1) * std_normal_lcdf(x[n]) + (b - 1) * std_normal_lccdf(x[n]);
    }
    return lp;
  }
  real probitbeta_rng(real a, real b) {
    return inv_Phi(beta_rng(a, b));
  }
