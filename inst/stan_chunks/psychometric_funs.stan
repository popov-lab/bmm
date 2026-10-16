// psychometric function psi = guess + (1 - guess) * (1 - lapse) * F(z); needs sdt_dist_funs.stan
// F(z) = 0.5 at midpoint, 0.05 and 0.95 one width apart; log_x == 1 reads x and midpoint on the log scale
real psychometric_log_lik(int y, int trials, real x, real midpoint, real width,
                          real guess, real lapse, int dist_type, int log_x,
                          real zmid, real zspan) {
  real dx = log_x == 1 ? log(x) - log(midpoint) : x - midpoint;
  real z = zmid + zspan * dx / width;
  real log1m_guess = log1m(guess);
  real out = lchoose(trials, y);

  // log(0) of a guess or lapse fixed at zero is -inf, which log_sum_exp absorbs
  if (y > 0) {
    out += y * log_sum_exp(log(guess),
                           log1m_guess + log1m(lapse) + sdt_log_cumprob(z, dist_type));
  }
  if (trials > y) {
    out += (trials - y) * (log1m_guess +
             log_sum_exp(log(lapse), log1m(lapse) + sdt_log_one_minus_cumprob(z, dist_type)));
  }
  return out;
}

// one family per combination of estimated (dpar) and fixed (data) guess and lapse rates
real psychometric_lpmf(int y, real mu, real midpoint, real width, real guess, real lapse,
                       int trials, real x, real guess_fixed, real lapse_fixed,
                       int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_log_lik(y, trials, x, midpoint, width, guess, lapse,
                              dist_type, log_x, zmid, zspan);
}

real psychometric_fixguess_lpmf(int y, real mu, real midpoint, real width, real lapse,
                                int trials, real x, real guess_fixed, real lapse_fixed,
                                int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_log_lik(y, trials, x, midpoint, width, guess_fixed, lapse,
                              dist_type, log_x, zmid, zspan);
}

real psychometric_fixlapse_lpmf(int y, real mu, real midpoint, real width, real guess,
                                int trials, real x, real guess_fixed, real lapse_fixed,
                                int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_log_lik(y, trials, x, midpoint, width, guess, lapse_fixed,
                              dist_type, log_x, zmid, zspan);
}

real psychometric_fixboth_lpmf(int y, real mu, real midpoint, real width,
                               int trials, real x, real guess_fixed, real lapse_fixed,
                               int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_log_lik(y, trials, x, midpoint, width, guess_fixed, lapse_fixed,
                              dist_type, log_x, zmid, zspan);
}
