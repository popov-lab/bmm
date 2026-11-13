// Single-probe change detection under TCC; needs sdt_dist_funs.stan.
// sdt_yn_lpmf with the probe's similarity to the target, a real strength in
// [0, 1], in place of its 0/1 stimulus, sdratio = 1 and normal noise: the
// probe's familiarity is normal(d * similarity, 1), and the observer answers
// "same" when it passes the criterion, measured from the familiarity of a probe
// with no similarity to the target. One response, coded 1 = "change"; dist is
// the probe's distance from the target in steps of the similarity, 0 = the target
real tcc_cd_lpmf(int y, real mu, real d, real criterion, int dist,
                 vector similarity) {
  real eta = d * similarity[dist + 1] - criterion;
  if (y == 1) {
    return sdt_log_cumprob(-eta, 1);
  }
  return sdt_log_one_minus_cumprob(-eta, 1);
}

// the same likelihood with the criterion measured from the midpoint between
// the familiarity of a probe with no similarity (0) and of the target (d)
real tcc_cd_bias_lpmf(int y, real mu, real d, real bias, int dist,
                      vector similarity) {
  return tcc_cd_lpmf(y | mu, d, bias + d / 2, dist, similarity);
}
