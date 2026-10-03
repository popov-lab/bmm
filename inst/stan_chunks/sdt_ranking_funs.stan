// ranking SDT log-probabilities for the native multinomial family, one category per rank position
// max_rank is real because brms passes data covariates into non-linear formulas as reals

// Gumbel-min ranking in closed form (Meyer-Grant et al., 2026); d is their equal-variance g'
// the lgamma difference loses accuracy below d = -20, which the default prior never reaches
real sdt_ranking_logp(int cat, real max_rank, real d) {
  real g = d;
  real e_neg_g = exp(-g);
  return -g + lgamma(max_rank) + lgamma(cat - 1 + e_neg_g)
         - lgamma(cat) - lgamma(max_rank + e_neg_g);
}

// Gaussian unequal-variance ranking by Gauss-Hermite quadrature; sdratio is the log SD ratio
// .ranking_fill_quadrature() fills in the node count and tables, so never repeat its tokens in a comment
real sdt_ranking_uv_logp(int cat, real max_rank, real d, real sdratio) {
  real sigma = exp(sdratio);
  int N_GH = {{N_GH}};
  vector[N_GH] gh_nodes = to_vector({{GH_NODES}});
  vector[N_GH] gh_weights = to_vector({{GH_WEIGHTS}});
  real log_choose = lgamma(max_rank) - lgamma(cat) - lgamma(max_rank - cat + 1);
  real p = 0;

  for (i in 1:N_GH) {
    // d is d_a; sdt_rms_scale() converts it to noise-SD units
    real eta = d * sdt_rms_scale(sigma) + sigma * gh_nodes[i];
    // probability space avoids 0 * -Inf = NaN, and Phi(-eta) keeps the upper tail's precision
    real cdf  = (max_rank > cat) ? Phi(eta)  : 1.0;
    real ccdf = (cat > 1)        ? Phi(-eta) : 1.0;
    p += gh_weights[i] * pow(cdf, max_rank - cat) * pow(ccdf, cat - 1);
  }

  return log_choose + log(p);
}

// ranks beyond the row's set size get a finite -100 logit, as absent components do in the mixture models
real sdt_ranking_logmu(int cat, real max_rank, real d, real sdratio,
                       int dist_type) {
  if (cat > max_rank) return -100;
  if (dist_type == 2) return sdt_ranking_logp(cat, max_rank, d);
  return sdt_ranking_uv_logp(cat, max_rank, d, sdratio);
}
