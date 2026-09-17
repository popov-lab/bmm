// Dual-process SDT (Yonelinas, 1994) category log-probability for the
// multinomial family. Requires sdt_rating_logmu_cat() from sdt_rating_funs.stan.
// Recollection adds mass to the most-confident category: old items recollected
// as old (Ro) load the top category, new items recall-rejected (Rn) the bottom.
// Ro/Rn arrive on the linear scale; inv_logit maps them to probabilities, so
// fixing them near zero recovers the standard rating likelihood.
// d is the familiarity distributions' d_a: recollection is a separate threshold
// process, so the d_a scaling belongs to the familiarity SDT process alone and
// is applied by sdt_rating_logmu_cat() below.
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
