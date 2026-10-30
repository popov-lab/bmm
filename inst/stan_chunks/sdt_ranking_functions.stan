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

// Gaussian unequal-variance ranking by quadrature; sdratio is the log SD ratio
// .ranking_fill_quadrature() fills in the node counts and tables, so never repeat its tokens in a comment
real sdt_ranking_uv_logp(int cat, real max_rank, real d, real sdratio) {
  real sigma = exp(sdratio);
  // d is d_a; sdt_rms_scale() converts it to noise-SD units
  real delta = d * sdt_rms_scale(sigma);
  real log_choose = lgamma(max_rank) - lgamma(cat) - lgamma(max_rank - cat + 1);
  real p = 0;

  if (sigma <= 1) {
    // Gauss-Hermite weighted by the target density, the narrower of the two
    int N_GH = {{N_GH}};
    vector[N_GH] gh_nodes = to_vector({{GH_NODES}});
    vector[N_GH] gh_weights = to_vector({{GH_WEIGHTS}});
    for (i in 1:N_GH) {
      real eta = delta + sigma * gh_nodes[i];
      // probability space avoids 0 * -Inf = NaN, and Phi(-eta) keeps the upper tail's precision
      real cdf  = (max_rank > cat) ? Phi(eta)  : 1.0;
      real ccdf = (cat > 1)        ? Phi(-eta) : 1.0;
      p += gh_weights[i] * pow(cdf, max_rank - cat) * pow(ccdf, cat - 1);
    }
  } else {
    // Gauss-Legendre on the noise axis, where the lure terms vary; the mirrored node's Phi is Phi(-x)
    int N_GL = {{N_GL}};
    real half_width = {{GL_HALF_WIDTH}};
    vector[N_GL] gl_nodes = to_vector({{GL_NODES}});
    vector[N_GL] gl_weights = to_vector({{GL_WEIGHTS}});
    vector[N_GL] gl_cdf = to_vector({{GL_CDF}});
    real inv_sigma = inv(sigma);
    for (i in 1:N_GL) {
      real lure = pow(gl_cdf[i], max_rank - cat) * pow(gl_cdf[N_GL + 1 - i], cat - 1);
      p += gl_weights[i] * lure * exp(-0.5 * square((gl_nodes[i] - delta) * inv_sigma));
    }
    p *= inv_sigma;
    // past the nodes the lure terms are 1 for rank 1 (above) and rank max_rank (below), 0 otherwise
    if (cat == 1) p += Phi((delta - half_width) * inv_sigma);
    if (cat == max_rank) p += Phi((-half_width - delta) * inv_sigma);
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
