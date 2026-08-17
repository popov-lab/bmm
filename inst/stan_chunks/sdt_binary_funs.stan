// Binary SDT likelihood with binomial on aggregated data
// Distribution dispatch (sdt_log_cumprob, sdt_log_one_minus_cumprob) lives in
// sdt_dist_funs.stan, which configure_model.sdt_binary loads before this chunk.
// Parameters:
//   y:         count of "old"/"signal" responses
//   d:         sensitivity as d_a (on link scale, transformed by brms); the
//              separation in noise units is d * sdt_rms_scale(sdratio)
//   criterion: response bias (identity link), on the noise-standardized axis
//   sdratio:   ratio of signal to noise SD (log link, so 1 = equal variance)
//   stimulus:  0 = noise/new, 1 = signal/old
//   dist_type: noise distribution (1-4, see sdt_dist_funs.stan)
//   trials:    total number of trials in this cell
real sdt_binary_lpmf(int y, real mu, real d, real criterion, real sdratio,
                     int stimulus, int dist_type, int trials) {
  real scale = stimulus == 1 ? sdratio : 1.0;
  real eta = (d * sdt_rms_scale(sdratio) / 2.0 * (2 * stimulus - 1)
              - criterion) / scale;
  // "old" means the evidence exceeded the criterion, so P(old) is the survival
  // function of the evidence distribution, evaluated at -eta. For symmetric
  // distributions S(-eta) == F(eta); for the extreme-value ones it does not,
  // and using F(eta) would fit the mirror distribution.
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
