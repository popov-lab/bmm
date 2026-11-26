// binomial log likelihood with success probability guess + (1 - guess) * (1 - lapse) * F(z)
real psychometric_function_log_lik(int y, int trials, real x, real midpoint, real width,
                                   real guess, real lapse, int dist_type, int log_x,
                                   real zmid, real zspan) {
  real dx = log_x == 1 ? log(x) - log(midpoint) : x - midpoint;
  real z = zmid + zspan * dx / width;
  if (dist_type == 2 || dist_type == 3) z = fmin(fmax(z, -600), 600);
  real log1m_guess = log1m(guess);
  real out = lchoose(trials, y);

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

// psychometric function with estimated guess and lapse rates
real psychometric_function_lpmf(int y, real mu, real midpoint, real width, real guess, real lapse,
                                int trials, real x, real guess_fixed, real lapse_fixed,
                                int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_function_log_lik(y, trials, x, midpoint, width, guess, lapse,
                                       dist_type, log_x, zmid, zspan);
}

// psychometric function with a fixed guess rate and an estimated lapse rate
real psychometric_function_fixguess_lpmf(int y, real mu, real midpoint, real width, real lapse,
                                         int trials, real x, real guess_fixed, real lapse_fixed,
                                         int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_function_log_lik(y, trials, x, midpoint, width, guess_fixed, lapse,
                                       dist_type, log_x, zmid, zspan);
}

// psychometric function with an estimated guess rate and a fixed lapse rate
real psychometric_function_fixlapse_lpmf(int y, real mu, real midpoint, real width, real guess,
                                         int trials, real x, real guess_fixed, real lapse_fixed,
                                         int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_function_log_lik(y, trials, x, midpoint, width, guess, lapse_fixed,
                                       dist_type, log_x, zmid, zspan);
}

// psychometric function with fixed guess and lapse rates
real psychometric_function_fixboth_lpmf(int y, real mu, real midpoint, real width,
                                        int trials, real x, real guess_fixed, real lapse_fixed,
                                        int dist_type, int log_x, real zmid, real zspan) {
  return psychometric_function_log_lik(y, trials, x, midpoint, width, guess_fixed, lapse_fixed,
                                       dist_type, log_x, zmid, zspan);
}
