// m-AFC SDT likelihood: P(correct | m, d') = integral f(x - d') * F(x)^(m-1) dx
// over the m latent alternatives (one signal shifted by d', m-1 noise), where
// the observer picks the maximum. F is the noise CDF (see sdt_dist_funs.stan).
//
//   normal     : 40-point Gauss-Hermite quadrature (closed form Phi(d'/sqrt(2)) at m = 2)
//   gumbel_max : exact softmax 1 / (1 + (m-1) * exp(-d'))
//   gumbel_min : exact gamma ratio Gamma(1 + e^-d') * Gamma(m) / Gamma(m + e^-d')
//   logistic   : 64-point Gauss-Legendre on [0, 1] of F(Q(u) + d')^(m-1)
//
// The quadrature tables are passed in from transformed data
// (sdt_mafc_tdata.stan); the family loops over rows, so building them here
// would rebuild 208 doubles per row per gradient evaluation.
//
// Reference: DeCarlo (2012); Green & Swets (1966)

real mafc_pc(real d, int m, int dist_type,
             data vector gh_nodes, data vector gh_weights,
             data vector gl_nodes, data vector gl_weights) {
  if (dist_type == 3)                                      // gumbel_max: softmax
    return 1.0 / (1 + (m - 1) * exp(-d));
  if (dist_type == 2)                                      // gumbel_min: gamma ratio
    return exp(lgamma(1 + exp(-d)) + lgamma(m) - lgamma(m + exp(-d)));

  if (dist_type == 1) {                                    // normal: Gauss-Hermite
    if (m == 2) return Phi(d / sqrt(2.0));
    vector[40] log_terms;
    for (i in 1:40)
      log_terms[i] = log(gh_weights[i]) + (m - 1) * std_normal_lcdf(gh_nodes[i] + d);
    return exp(log_sum_exp(log_terms));
  }

  // logistic (dist_type == 4): Gauss-Legendre on [0, 1]
  real p = 0;
  for (i in 1:64)
    p += gl_weights[i] *
         pow(sdt_cumprob(sdt_quantile(gl_nodes[i], dist_type) + d, dist_type), m - 1);
  return p;
}

real sdt_mafc_lpmf(int y, real mu, real d, int m, int dist_type, int trials,
                   data vector gh_nodes, data vector gh_weights,
                   data vector gl_nodes, data vector gl_weights) {
  return binomial_lpmf(y | trials,
                       mafc_pc(d, m, dist_type,
                               gh_nodes, gh_weights, gl_nodes, gl_weights));
}
