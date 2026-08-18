// Ranking SDT category log-probabilities for the native multinomial family.
// Each rank position is one multinomial category; sdt_ranking_logmu returns
// log p(rank) so brms' softmax recovers the rank distribution. Ranks beyond a
// row's set size (cat > max_rank) are switched off with a finite -100 logit,
// the same idiom the mixture models use for absent set-size components. The set
// size max_rank arrives as a real covariate (brms passes data covariates into
// non-linear formulas as reals), so the gamma terms use lgamma rather than the
// integer-only choose().

// Gumbel-min ranking: closed form via gamma-function ratios.
// Meyer-Grant et al. (2025), based on extreme-value (min) order statistics.
//   cat:      rank position (1 = most likely target, max_rank = least)
//   max_rank: number of ranked items (m) on this row
//   d:   sensitivity as d_a; for the Gumbel branch it is the equal-variance
//        g' of Meyer-Grant et al., since that model carries no sdratio
real sdt_ranking_logp(int cat, real max_rank, real d) {
  real g = d;
  real e_neg_g = exp(-g);
  return -g + lgamma(max_rank) + lgamma(cat - 1 + e_neg_g)
         - lgamma(cat) - lgamma(max_rank + e_neg_g);
}

// Gaussian UV-SDT ranking: Gauss-Hermite quadrature over the target
// distribution. Retains smooth gradients without adaptive integration.
//
// The node count is NOT fixed here: .ranking_fill_quadrature() substitutes the
// three doubled-brace tokens below with the rule chosen by .ranking_gh_n(),
// which sizes the quadrature by the largest set size and by whether sdratio is
// free. R and Stan therefore share one source of truth (.gh_rule()) rather than
// two hand-copied constant tables. Do not repeat the tokens in comments.
//   sdratio: log ratio of signal to noise SD (exp(sdratio) = sigma_s / sigma_n)
real sdt_ranking_uv_logp(int cat, real max_rank, real d, real sdratio) {
  real sigma = exp(sdratio);
  int N_GH = {{N_GH}};
  vector[N_GH] gh_nodes = to_vector({{GH_NODES}});
  vector[N_GH] gh_weights = to_vector({{GH_WEIGHTS}});
  real log_choose = lgamma(max_rank) - lgamma(cat) - lgamma(max_rank - cat + 1);
  real p = 0;

  for (i in 1:N_GH) {
    // `d` is d_a; sdt_rms_scale() converts it to noise-SD units, and is 1
    // when sigma is 1, so the equal-variance case is unchanged
    real eta = d * sdt_rms_scale(sigma) + sigma * gh_nodes[i];
    // Probability-space formulation: avoids log-CDF underflow (→ -Inf) and the
    // 0 * Inf = NaN that arises when multiplying a zero coefficient by a -Inf
    // log-CDF. Boundary guards use 1.0 so pow(1, 0) = 1.
    real cdf  = (max_rank > cat) ? Phi(eta)         : 1.0;
    real ccdf = (cat > 1)        ? (1.0 - Phi(eta)) : 1.0;
    p += gh_weights[i] * pow(cdf, max_rank - cat) * pow(ccdf, cat - 1);
  }

  return log_choose + log(p);
}

// Multinomial-logit value for rank category `cat`. Returns the finite -100
// sentinel when the rank exceeds this row's set size; otherwise dispatches on
// the noise distribution (id 2 = gumbel_min, id 1 = normal; see .sdt_dists).
real sdt_ranking_logmu(int cat, real max_rank, real d, real sdratio,
                       int dist_type) {
  if (cat > max_rank) return -100;
  if (dist_type == 2) return sdt_ranking_logp(cat, max_rank, d);
  return sdt_ranking_uv_logp(cat, max_rank, d, sdratio);
}
