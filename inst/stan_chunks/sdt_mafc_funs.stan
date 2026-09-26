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
// rather than P(correct): on the probability scale P(correct) rounds to 1 as
// soon as its complement falls under the double epsilon, and the likelihood
// stops responding to d' well inside the range a log link on d can reach.
// Each branch therefore reads log P(correct) and log(1 - P(correct)) off
// whichever side still resolves it and subtracts. The two quadrature branches
// get both from one sweep of the nodes: log_sum_exp of the weighted log-CDFs
// keeps its precision as P(correct) -> 0, and the sum of the complements --
// every term positive, so nothing cancels -- as P(correct) -> 1. Reading
// either side off the other loses the opposite tail.
//
// The range that buys is not the range it responds over: against a 1024-point
// Gauss-Hermite reference at m = 4, the normal branch is accurate to about
// d' = 20 (relative error 7.0e-07; 8.0e-03 at 25, 4.5e-02 from 30 on) while it
// keeps decreasing out to d' = 48. Beyond 20 it is therefore a slightly wrong
// gradient rather than a right one -- still the better trade, because what it
// replaces is a flat plateau from d' = 12 with no gradient at all, which a
// sampler random-walks through instead of rejecting.
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
    // the gamma ratio telescopes to prod(k / (k + e)), and log1p still resolves
    // the tiny e = exp(-d') at which the difference of lgammas has cancelled to
    // zero -- off the telescoped form this branch dies where the probability
    // scale it is meant to rescue does
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
  // binomial_logit_lpmf saturates on a large finite logit but raises a domain
  // error on an infinite one; past the point where a branch still resolves the
  // complement, every trial is correct, or none is
  if (is_inf(logit_pc))
    return (y == (logit_pc > 0 ? trials : 0)) ? 0.0 : negative_infinity();
  return binomial_logit_lpmf(y | trials, logit_pc);
}
