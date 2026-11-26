// log probability of target rank cat under Gumbel-min noise (Meyer-Grant et al., 2026)
real sdt_ranking_logp(int cat, real max_rank, real d) {
  real g = d;
  real e_neg_g = exp(-g);
  return -g + lgamma(max_rank) + lgamma(cat - 1 + e_neg_g)
         - lgamma(cat) - lgamma(max_rank + e_neg_g);
}

// log probability of target rank cat under unequal-variance Gaussian noise, by quadrature
real sdt_ranking_uv_logp(int cat, real max_rank, real d, real sdratio) {
  real sigma = exp(sdratio);
  real delta = d * sdt_rms_scale(sigma);
  real log_choose = lgamma(max_rank) - lgamma(cat) - lgamma(max_rank - cat + 1);
  real p = 0;

  if (sigma <= 1) {
    int N_GH = {{N_GH}};
    vector[N_GH] gh_nodes = to_vector({{GH_NODES}});
    vector[N_GH] gh_weights = to_vector({{GH_WEIGHTS}});
    for (i in 1:N_GH) {
      real eta = delta + sigma * gh_nodes[i];
      real cdf  = (max_rank > cat) ? Phi(eta)  : 1.0;
      real ccdf = (cat > 1)        ? Phi(-eta) : 1.0;
      p += gh_weights[i] * pow(cdf, max_rank - cat) * pow(ccdf, cat - 1);
    }
  } else {
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
    if (cat == 1) p += Phi((delta - half_width) * inv_sigma);
    if (cat == max_rank) p += Phi((-half_width - delta) * inv_sigma);
  }

  return log_choose + log(p);
}

// log probability of target rank cat; ranks beyond the set size of a trial get a very low value
real sdt_ranking_logmu(int cat, real max_rank, real d, real sdratio,
                       int dist_type) {
  if (cat > max_rank) return -100;
  if (dist_type == 2) return sdt_ranking_logp(cat, max_rank, d);
  return sdt_ranking_uv_logp(cat, max_rank, d, sdratio);
}
