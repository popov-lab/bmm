// Ranking SDT likelihood functions
// Implements ranking model where participants rank m items by perceived strength.
// The likelihood uses y * log(p) (multinomial kernel, constant absorbed).

// Gumbel-min ranking: closed-form via gamma function ratios
// Meyer-Grant et al. (2025), based on extreme-value ordering statistics.
// Parameters:
//   y:        count of times target received this rank
//   mu:       internal parameter (fixed to 0, required by brms)
//   dprime:   sensitivity parameter (g' in Gumbel parameterization)
//   rank_pos: rank position (1 = most likely target, m = least)
//   max_rank: total number of ranked items (m)
real sdt_ranking_lpmf(int y, real mu, real dprime,
                      int rank_pos, int max_rank) {
  if (y == 0) return 0;
  real g = dprime;
  real e_neg_g = exp(-g);
  real log_p = -g + lgamma(max_rank) + lgamma(rank_pos - 1 + e_neg_g)
               - lgamma(rank_pos) - lgamma(max_rank + e_neg_g);
  return y * log_p;
}

// Gaussian UV-SDT ranking: numerical integration via integrate_1d
// Parameters:
//   x:     integration variable
//   xc:    complement (unused, required by Stan)
//   theta: {dprime, sigma} where sigma = exp(sdratio)
//   x_r:   empty real array (required by integrate_1d)
//   x_i:   {rank_pos, max_rank}
real ranking_integrand(real x, real xc,
                       array[] real theta,
                       array[] real x_r, array[] int x_i) {
  real d = theta[1];
  real sigma = theta[2];
  int r = x_i[1];
  int m = x_i[2];
  return pow(Phi(x), m - r)
         * exp(normal_lpdf(x | d, sigma))
         * pow(1 - Phi(x), r - 1);
}

// Gaussian UV-SDT ranking likelihood
// Parameters:
//   y:        count of times target received this rank
//   mu:       internal parameter (fixed to 0, required by brms)
//   dprime:   sensitivity (d')
//   sdratio:  log ratio of signal to noise SD (exp(sdratio) = sigma_s/sigma_n)
//   rank_pos: rank position (1..m)
//   max_rank: total number of ranked items (m)
//   x_r:      empty real array (required by integrate_1d)
//   x_i:      empty int array (required by integrate_1d)
real sdt_ranking_uv_lpmf(int y, real mu, real dprime, real sdratio,
                          int rank_pos, int max_rank,
                          data array[] real x_r, data array[] int x_i) {
  if (y == 0) return 0;
  real sigma = exp(sdratio);
  real p = choose(max_rank - 1, rank_pos - 1)
           * integrate_1d(ranking_integrand,
                          negative_infinity(), positive_infinity(),
                          {dprime, sigma}, x_r, {rank_pos, max_rank});
  return y * log(p);
}
