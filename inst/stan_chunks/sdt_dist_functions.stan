// CDF of the noise distribution; dist_type: 1 = normal, 2 = gumbel_min, 3 = gumbel_max, 4 = logistic
real sdt_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return Phi(eta);
  if (dist_type == 2) return 1 - exp(-exp(eta));
  if (dist_type == 3) return exp(-exp(-eta));
  if (dist_type == 4) return inv_logit(eta);
  reject("sdt_cumprob: unknown dist_type: ", dist_type);
}

// log CDF of the noise distribution
real sdt_log_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return std_normal_lcdf(eta);
  if (dist_type == 2) return log1m_exp(-exp(eta));
  if (dist_type == 3) return -exp(-eta);
  if (dist_type == 4) return -log1p_exp(-eta);
  reject("sdt_log_cumprob: unknown dist_type: ", dist_type);
}

// converts d_a to units of the noise SD, given the ratio of the signal SD to the noise SD
real sdt_rms_scale(real sdratio) {
  return sqrt((1 + square(sdratio)) / 2);
}

// log of 1 - CDF of the noise distribution
real sdt_log_one_minus_cumprob(real eta, int dist_type) {
  if (dist_type == 1) return std_normal_lcdf(-eta);
  if (dist_type == 2) return -exp(eta);
  if (dist_type == 3) return log1m_exp(-exp(-eta));
  if (dist_type == 4) return -eta - log1p_exp(-eta);
  reject("sdt_log_one_minus_cumprob: unknown dist_type: ", dist_type);
}

// vectorized version of sdt_rms_scale
vector sdt_rms_scale(vector sdratio) {
  return sqrt((1 + square(sdratio)) / 2);
}

// vectorized version of sdt_log_cumprob
vector sdt_log_cumprob(vector eta, int dist_type) {
  vector[rows(eta)] out;
  for (i in 1:rows(eta)) out[i] = sdt_log_cumprob(eta[i], dist_type);
  return out;
}

// vectorized version of sdt_log_one_minus_cumprob
vector sdt_log_one_minus_cumprob(vector eta, int dist_type) {
  vector[rows(eta)] out;
  for (i in 1:rows(eta)) out[i] = sdt_log_one_minus_cumprob(eta[i], dist_type);
  return out;
}

// sum of the log CDF over a vector
real sdt_log_cumprob_sum(vector eta, int dist_type) {
  if (dist_type == 1) return std_normal_lcdf(eta);
  return sum(sdt_log_cumprob(eta, dist_type));
}

// sum of the log of 1 - CDF over a vector
real sdt_log_one_minus_cumprob_sum(vector eta, int dist_type) {
  if (dist_type == 1) return std_normal_lcdf(-eta);
  return sum(sdt_log_one_minus_cumprob(eta, dist_type));
}

// quantile function (inverse CDF) of the noise distribution
real sdt_quantile(real u, int dist_type) {
  if (dist_type == 1) return inv_Phi(u);
  if (dist_type == 2) return log(-log1m(u));
  if (dist_type == 3) return -log(-log(u));
  if (dist_type == 4) return logit(u);
  reject("sdt_quantile: unknown dist_type: ", dist_type);
}
