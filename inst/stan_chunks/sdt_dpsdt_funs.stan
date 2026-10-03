// dual-process SDT (Yonelinas, 1994) category log-probability; needs sdt_rating_funs.stan
// recollection adds mass to the most confident category; d is the d_a of the familiarity process only
real sdt_dpsdt_logmu_cat(int cat, vector thresholds,
                         real d, real sdratio, real stimulus,
                         int dist_type, real Ro, real Rn) {
  int K_full = num_elements(thresholds) + 1;
  real log_base = sdt_rating_logmu_cat(cat, thresholds, d, sdratio,
                                       stimulus, dist_type);

  if (stimulus > 0.5) {
    if (cat == K_full) {
      return log_sum_exp(log1m_inv_logit(Ro) + log_base, log_inv_logit(Ro));
    }
    return log1m_inv_logit(Ro) + log_base;
  }

  if (cat == 1) {
    return log_sum_exp(log1m_inv_logit(Rn) + log_base, log_inv_logit(Rn));
  }
  return log1m_inv_logit(Rn) + log_base;
}
