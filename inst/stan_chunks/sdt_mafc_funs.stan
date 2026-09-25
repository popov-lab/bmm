// m-AFC SDT likelihood: P(correct | m, d') = integral f(x - d') * F(x)^(m-1) dx
// over the m latent alternatives (one signal shifted by d', m-1 noise), where
// the observer picks the maximum. F is the noise CDF (see sdt_dist_funs.stan).
//
//   normal     : 40-point Gauss-Hermite quadrature (closed form at m = 2)
//   gumbel_max : exact softmax 1 / (1 + (m-1) * exp(-d'))
//   gumbel_min : exact gamma ratio Gamma(1 + e^-d') * Gamma(m) / Gamma(m + e^-d')
//   logistic   : 64-point Gauss-Legendre on [0, 1] of F(Q(u) + d')^(m-1)
//
// The binomial is taken on the logit scale, so this returns logit P(correct)
// rather than P(correct). On the probability scale P(correct) rounds to 1 as
// soon as its complement falls below the double epsilon, and the likelihood
// stops responding to d': at y = 90 of 100 trials and m = 4 the normal branch
// was constant at -312.04 from d' = 12 and gumbel_min, gumbel_max and logistic
// returned -inf from d' = 37, 38 and 40.
//
// gumbel_max's logit is analytic and exact for every d'. The two quadrature
// branches accumulate the COMPLEMENT 1 - P(correct): every term of that sum is
// positive, so the small quantity keeps its relative precision instead of
// being read off a number that has rounded to 1. gumbel_min is bounded by the
// rounding of 1 + e^-d' inside lgamma rather than by the parameterization, so
// its range is unchanged.
//
// The quadrature tables are passed in from transformed data
// (sdt_mafc_tdata.stan); the family loops over rows, so building them here
// would rebuild 208 doubles per row per gradient evaluation.
//
// Reference: DeCarlo (2012); Green & Swets (1966)

real mafc_logit_pc(real d, int m, int dist_type,
                   data vector gh_nodes, data vector gh_weights,
                   data vector gl_nodes, data vector gl_weights) {
  if (dist_type == 3)                                      // gumbel_max: softmax
    return d - log(m - 1);
  if (dist_type == 2) {                                    // gumbel_min: gamma ratio
    real log_pc = fmin(lgamma(1 + exp(-d)) + lgamma(m) - lgamma(m + exp(-d)), 0);
    return log_pc - log1m_exp(log_pc);
  }

  real q = 0;                                              // 1 - P(correct)
  if (dist_type == 1) {                                    // normal: Gauss-Hermite
    if (m == 2)
      return sdt_log_cumprob(d / sqrt(2.0), dist_type) -
             sdt_log_one_minus_cumprob(d / sqrt(2.0), dist_type);
    for (i in 1:40)
      q += gh_weights[i] *
           (-expm1((m - 1) * sdt_log_cumprob(gh_nodes[i] + d, dist_type)));
  } else {                                                 // logistic
    for (i in 1:64)
      q += gl_weights[i] *
           (-expm1((m - 1) *
                   sdt_log_cumprob(sdt_quantile(gl_nodes[i], dist_type) + d,
                                   dist_type)));
  }
  // the weights sum to 1 only to rounding, so q can land just above it
  q = fmin(q, 1);
  return log1m(q) - log(q);
}

real sdt_mafc_lpmf(int y, real mu, real d, int m, int dist_type, int trials,
                   data vector gh_nodes, data vector gh_weights,
                   data vector gl_nodes, data vector gl_weights) {
  real logit_pc = mafc_logit_pc(d, m, dist_type,
                                gh_nodes, gh_weights, gl_nodes, gl_weights);
  // the quadrature runs out of resolution before the logit scale does, and
  // binomial_logit_lpmf raises a domain error on a non-finite argument rather
  // than saturating; past that point every trial is correct, or none is
  if (is_inf(logit_pc))
    return (y == (logit_pc > 0 ? trials : 0)) ? 0.0 : negative_infinity();
  return binomial_logit_lpmf(y | trials, logit_pc);
}
