// logit of P(correct) in m-alternative forced choice (Green & Swets, 1966; DeCarlo, 2012)
real mafc_logit_pc(real d, int m, int dist_type,
                   data vector gh_nodes, data vector gh_weights,
                   data vector gl_nodes, data vector gl_weights) {
  if (dist_type == 3)
    return d - log(m - 1);
  if (dist_type == 2) {
    real e = exp(-d);
    real log_pc = 0;
    for (k in 1:(m - 1))
      log_pc -= log1p(e / k);
    return log_pc - log(-expm1(log_pc));
  }

  if (dist_type == 1) {
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
    return log_sum_exp(log_terms) - log(fmin(q, 1));
  }

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

// binomial log likelihood of the number of correct responses, summed over all observations
real sdt_mafc_lpmf(array[] int y, vector mu, vector d, array[] int m, int dist_type,
                   array[] int trials, data vector gh_nodes, data vector gh_weights,
                   data vector gl_nodes, data vector gl_weights) {
  int N = size(y);
  vector[N] logit_pc;
  array[N] int finite;
  int n_finite = 0;
  for (n in 1:N) {
    logit_pc[n] = mafc_logit_pc(d[n], m[n], dist_type,
                                gh_nodes, gh_weights, gl_nodes, gl_weights);
    if (is_inf(logit_pc[n])) {
      if (y[n] != (logit_pc[n] > 0 ? trials[n] : 0))
        return negative_infinity();
    } else {
      n_finite += 1;
      finite[n_finite] = n;
    }
  }
  return binomial_logit_lpmf(y[finite[1:n_finite]] | trials[finite[1:n_finite]],
                             logit_pc[finite[1:n_finite]]);
}
