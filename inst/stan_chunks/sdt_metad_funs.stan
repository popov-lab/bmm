// Meta-d' SDT (Maniscalco & Lau, 2012) category log-probability for the
// multinomial family. Requires sdt_rating_logmu_cat() from sdt_rating_funs.stan
// and the noise-CDF dispatch from sdt_dist_funs.stan. Confidence thresholds are
// read off the metacognitive sensitivity metad; each side of the central
// criterion is rescaled so the summed "old"/"new" mass matches type-1 dprime.
// metad = dprime makes log_norm = 0 and recovers the standard rating likelihood.
real sdt_metad_logmu_cat(int cat, vector thresholds,
                         real dprime, real metad, real sdratio, real stimulus,
                         int dist_type) {
  int K_full = num_elements(thresholds) + 1;
  int mid = (K_full - 1) %/% 2 + 1;
  real scale = stimulus > 0.5 ? exp(sdratio) : 1.0;
  real d_shift = dprime / 2.0 * (2 * stimulus - 1);
  real metad_shift = metad / 2.0 * (2 * stimulus - 1);
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
