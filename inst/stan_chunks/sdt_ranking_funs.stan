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
//   dprime:   sensitivity (g' in the Gumbel parameterization)
real sdt_ranking_logp(int cat, real max_rank, real dprime) {
  real g = dprime;
  real e_neg_g = exp(-g);
  return -g + lgamma(max_rank) + lgamma(cat - 1 + e_neg_g)
         - lgamma(cat) - lgamma(max_rank + e_neg_g);
}

// Gaussian UV-SDT ranking: 20-point Gauss-Hermite quadrature over the target
// distribution. Retains smooth gradients without adaptive integration.
//   sdratio: log ratio of signal to noise SD (exp(sdratio) = sigma_s / sigma_n)
real sdt_ranking_uv_logp(int cat, real max_rank, real dprime, real sdratio) {
  real sigma = exp(sdratio);
  int N_GH = 20;
  vector[N_GH] gh_nodes = to_vector({
    -7.6190485416797546e+00, -6.5105901570136488e+00,
    -5.5787388058932059e+00, -4.7345813340460463e+00,
    -3.9439673506573110e+00, -3.1890148165533843e+00,
    -2.4586636111723603e+00, -1.7452473208141255e+00,
    -1.0429453488027509e+00, -3.4696415708135458e-01,
     3.4696415708135830e-01,  1.0429453488027574e+00,
     1.7452473208141317e+00,  2.4586636111723683e+00,
     3.1890148165533900e+00,  3.9439673506573163e+00,
     4.7345813340460552e+00,  5.5787388058932033e+00,
     6.5105901570136551e+00,  7.6190485416797591e+00
  });
  vector[N_GH] gh_weights = to_vector({
    1.2578006724378954e-13, 2.4820623623151972e-10,
    6.1274902599825256e-08, 4.4021210902309806e-06,
    1.2882627996193093e-04, 1.8301031310804826e-03,
    1.3997837447100857e-02, 6.1506372063977507e-02,
    1.6173933398399959e-01, 2.6079306344955683e-01,
    2.6079306344955305e-01, 1.6173933398399776e-01,
    6.1506372063977438e-02, 1.3997837447101162e-02,
    1.8301031310805052e-03, 1.2882627996193072e-04,
    4.4021210902309052e-06, 6.1274902599829068e-08,
    2.4820623623151936e-10, 1.2578006724379269e-13
  });
  real log_choose = lgamma(max_rank) - lgamma(cat) - lgamma(max_rank - cat + 1);
  real p = 0;

  for (i in 1:N_GH) {
    real eta = dprime + sigma * gh_nodes[i];
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
// the noise distribution (id 2 = gumbel_min, id 1 = normal; see .SDT_DISTS).
real sdt_ranking_logmu(int cat, real max_rank, real dprime, real sdratio,
                       int dist_type) {
  if (cat > max_rank) return -100;
  if (dist_type == 2) return sdt_ranking_logp(cat, max_rank, dprime);
  return sdt_ranking_uv_logp(cat, max_rank, dprime, sdratio);
}
