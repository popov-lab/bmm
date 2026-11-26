// log probability of a rating category under the meta-d' model (Maniscalco & Lau, 2012)
real sdt_metad_logmu_cat(int cat, vector thresholds,
                         real d, real metad, real sdratio, real stimulus,
                         int dist_type) {
  int K_full = num_elements(thresholds) + 1;
  int mid = K_full %/% 2;
  real sigma = exp(sdratio);
  real rms = sdt_rms_scale(sigma);
  real scale = stimulus > 0.5 ? sigma : 1.0;
  real sep_d = d * rms;
  real sep_metad = metad * rms;
  vector[K_full - 1] thr_metad = thresholds - (sep_d - sep_metad) / 2.0;
  real eta_d = (thresholds[mid] - sep_d * stimulus) / scale;
  real eta_metad = (thr_metad[mid] - sep_metad * stimulus) / scale;
  real log_norm;

  if (cat <= mid) {
    log_norm = sdt_log_cumprob(eta_d, dist_type) -
      sdt_log_cumprob(eta_metad, dist_type);
  } else {
    log_norm = sdt_log_one_minus_cumprob(eta_d, dist_type) -
      sdt_log_one_minus_cumprob(eta_metad, dist_type);
  }

  return sdt_rating_logmu_cat(cat, thr_metad, metad, sdratio, stimulus,
                              dist_type) + log_norm;
}
