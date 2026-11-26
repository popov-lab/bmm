// log probability of rating category cat; d is d_a and sdratio the log of the SD ratio
real sdt_rating_logmu_cat(int cat, vector thresholds,
                          real d, real sdratio, real stimulus,
                          int dist_type) {
  real sigma = exp(sdratio);
  real shift = d * sdt_rms_scale(sigma) * stimulus;
  real scale = stimulus > 0.5 ? sigma : 1.0;
  return ordinal_log_prob_cat(cat, thresholds, shift, scale, dist_type);
}
