// yes/no SDT likelihood for aggregated counts; needs sdt_dist_funs.stan
// d is d_a, sdratio the SD ratio on the natural scale, stimulus 0 = new and 1 = old
real sdt_yn_lpmf(int y, real mu, real d, real criterion, real sdratio,
                 int stimulus, int dist_type, int trials) {
  real scale = stimulus == 1 ? sdratio : 1.0;
  real eta = (d * sdt_rms_scale(sdratio) / 2.0 * (2 * stimulus - 1)
              - criterion) / scale;
  // P(old) is the survivor at -eta, which differs from F(eta) for the asymmetric Gumbels
  real log_p = sdt_log_one_minus_cumprob(-eta, dist_type);
  real log_q = sdt_log_cumprob(-eta, dist_type);
  real out = lchoose(trials, y);

  if (y > 0) {
    out += y * log_p;
  }
  if (trials > y) {
    out += (trials - y) * log_q;
  }
  return out;
}
