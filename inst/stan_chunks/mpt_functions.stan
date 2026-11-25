  // log density of x when inv_logit(x) follows a Beta(a, b) distribution
  real logitbeta_lpdf(real x, real a, real b) {
    return a * log_inv_logit(x) + b * log1m_inv_logit(x) - lbeta(a, b);
  }
  // the same log density summed over a vector
  real logitbeta_lpdf(vector x, real a, real b) {
    return a * sum(log_inv_logit(x)) + b * sum(log1m_inv_logit(x))
           - rows(x) * lbeta(a, b);
  }
  // random draw of x such that inv_logit(x) follows a Beta(a, b) distribution
  real logitbeta_rng(real a, real b) {
    return logit(beta_rng(a, b));
  }
  // log density of x when Phi(x) follows a Beta(a, b) distribution
  real probitbeta_lpdf(real x, real a, real b) {
    return (a - 1) * std_normal_lcdf(x) + (b - 1) * std_normal_lccdf(x)
           + std_normal_lpdf(x) - lbeta(a, b);
  }
  // the same log density summed over a vector
  real probitbeta_lpdf(vector x, real a, real b) {
    real lp = std_normal_lpdf(x) - rows(x) * lbeta(a, b);
    for (n in 1:rows(x)) {
      lp += (a - 1) * std_normal_lcdf(x[n]) + (b - 1) * std_normal_lccdf(x[n]);
    }
    return lp;
  }
  // random draw of x such that Phi(x) follows a Beta(a, b) distribution
  real probitbeta_rng(real a, real b) {
    return inv_Phi(beta_rng(a, b));
  }
