// m-AFC SDT likelihood (Green & Swets, 1966; DeCarlo, 2012): P(correct) is P(signal is the largest of m)
// returns logit P(correct), because P(correct) rounds to 1 at large d' and the likelihood goes flat

real mafc_logit_pc(real d, int m, int dist_type,
                   data vector gh_nodes, data vector gh_weights,
                   data vector gl_nodes, data vector gl_weights) {
  if (dist_type == 3)                                      // gumbel_max: softmax
    return d - log(m - 1);
  if (dist_type == 2) {                                    // gumbel_min: gamma ratio
    // telescoped to prod(k / (k + e)), so log1p resolves the tiny e where the lgammas cancel
    real e = exp(-d);
    real log_pc = 0;
    for (k in 1:(m - 1))
      log_pc -= log1p(e / k);
    return log_pc - log(-expm1(log_pc));
  }

  if (dist_type == 1) {                                    // normal: Gauss-Hermite
    if (m == 2)
      return sdt_log_cumprob(d / sqrt(2.0), dist_type) -
             sdt_log_one_minus_cumprob(d / sqrt(2.0), dist_type);
    vector[40] log_terms;
    real q = 0;
    for (i in 1:40) {
      real log_cdf = (m - 1) * sdt_log_cumprob(gh_nodes[i] + d, dist_type);
      log_terms[i] = log(gh_weights[i]) + log_cdf;
      q += gh_weights[i] * (-expm1(log_cdf));
    }
    // the weights sum to 1 only to rounding, so the complement can land above it
    return log_sum_exp(log_terms) - log(fmin(q, 1));
  }

  // logistic (dist_type == 4): Gauss-Legendre on [0, 1]
  vector[64] log_terms;
  real q = 0;
  for (i in 1:64) {
    real log_cdf = (m - 1) *
      sdt_log_cumprob(sdt_quantile(gl_nodes[i], dist_type) + d, dist_type);
    log_terms[i] = log(gl_weights[i]) + log_cdf;
    q += gl_weights[i] * (-expm1(log_cdf));
  }
  return log_sum_exp(log_terms) - log(fmin(q, 1));
}

real sdt_mafc_lpmf(int y, real mu, real d, int m, int dist_type, int trials,
                   data vector gh_nodes, data vector gh_weights,
                   data vector gl_nodes, data vector gl_weights) {
  real logit_pc = mafc_logit_pc(d, m, dist_type,
                                gh_nodes, gh_weights, gl_nodes, gl_weights);
  // binomial_logit_lpmf errors on an infinite logit, where every trial is correct or none is
  if (is_inf(logit_pc))
    return (y == (logit_pc > 0 ? trials : 0)) ? 0.0 : negative_infinity();
  return binomial_logit_lpmf(y | trials, logit_pc);
}
