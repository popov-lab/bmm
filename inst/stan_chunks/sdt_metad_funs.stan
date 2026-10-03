// meta-d' SDT (Maniscalco & Lau, 2012) category log-probability; needs sdt_rating_funs.stan
// thresholds follow metad, each side rescaled so its mass matches type-1 d; metad = d recovers sdt_rating
// the criterion is threshold K/2, so K must be even; sdt_rating() refuses the metad version at odd K
real sdt_metad_logmu_cat(int cat, vector thresholds,
                         real d, real metad, real sdratio, real stimulus,
                         int dist_type) {
  int K_full = num_elements(thresholds) + 1;
  int mid = K_full %/% 2;
  real sigma = exp(sdratio);
  real rms = sdt_rms_scale(sigma);
  real scale = stimulus > 0.5 ? sigma : 1.0;
  real d_shift = d * rms / 2.0 * (2 * stimulus - 1);
  real metad_shift = metad * rms / 2.0 * (2 * stimulus - 1);
  real crit = thresholds[mid];
  real log_norm;

  if (cat <= mid) {
    log_norm = sdt_log_cumprob((crit - d_shift) / scale, dist_type) -
      sdt_log_cumprob((crit - metad_shift) / scale, dist_type);
  } else {
    log_norm = sdt_log_one_minus_cumprob((crit - d_shift) / scale, dist_type) -
      sdt_log_one_minus_cumprob((crit - metad_shift) / scale, dist_type);
  }

  return sdt_rating_logmu_cat(cat, thresholds, metad, sdratio, stimulus,
                              dist_type) + log_norm;
}
