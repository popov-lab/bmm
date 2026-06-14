// Binary SDT likelihood with binomial on aggregated data
// Distribution dispatch (sdt_log_cumprob, sdt_log_one_minus_cumprob) lives in
// sdt_dist_funs.stan, which configure_model.sdt_binary loads before this chunk.
// Parameters:
//   y:         count of "old"/"signal" responses
//   dprime:    sensitivity (on link scale, transformed by brms)
//   criterion: response bias (identity link)
//   sdratio:   log ratio of signal to noise SD (0 = equal variance)
//   stimulus:  0 = noise/new, 1 = signal/old
//   dist_type: noise distribution (1-4, see sdt_cumprob)
//   trials:    total number of trials in this cell
real sdt_binary_lpmf(int y, real mu, real dprime, real criterion, real sdratio,
                     int stimulus, int dist_type, int trials) {
  real scale = stimulus == 1 ? exp(sdratio) : 1.0;
  real eta = (dprime / 2.0 * (2 * stimulus - 1) - criterion) / scale;
  real log_p = sdt_log_cumprob(eta, dist_type);
  real log_q = sdt_log_one_minus_cumprob(eta, dist_type);
  real out = lchoose(trials, y);

  if (y > 0) {
    out += y * log_p;
  }
  if (trials > y) {
    out += (trials - y) * log_q;
  }
  return out;
}
