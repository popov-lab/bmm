// TCC change detection: log probability of a response to a probe at distance dist from the target
real tcc_cd_lpmf(int y, real mu, real d, real criterion, int dist,
                 vector similarity) {
  real eta = d * similarity[dist + 1] - criterion;
  if (y == 1) {
    return sdt_log_cumprob(-eta, 1);
  }
  return sdt_log_one_minus_cumprob(-eta, 1);
}

// the same, with the criterion measured from the midpoint between no similarity and the target
real tcc_cd_bias_lpmf(int y, real mu, real d, real bias, int dist,
                      vector similarity) {
  return tcc_cd_lpmf(y | mu, d, bias + d / 2, dist, similarity);
}
