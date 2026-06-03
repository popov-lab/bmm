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

// Gaussian UV-SDT ranking: fixed Gauss-Hermite quadrature over the target
// distribution. This replaces the slower adaptive integrate_1d path while
// retaining smooth gradients.
// Parameters:
//   y:        count of times target received this rank
//   mu:       internal parameter (fixed to 0, required by brms)
//   dprime:   sensitivity (d')
//   sdratio:  log ratio of signal to noise SD (exp(sdratio) = sigma_s/sigma_n)
//   rank_pos: rank position (1..m)
//   max_rank: total number of ranked items (m)
real sdt_ranking_uv_lpmf(int y, real mu, real dprime, real sdratio,
                         int rank_pos, int max_rank) {
  if (y == 0) return 0;
  real sigma = exp(sdratio);
  int N_GH = 20;
  vector[N_GH] gh_nodes = to_vector({
    -7.6190485416797546e+00, -6.5105901570136488e+00,
    -5.5787388058932059e+00, -4.7345813340460463e+00,
    -3.9439673506573110e+00, -3.1890148165533843e+00,
    -2.4586636111723603e+00, -1.7452473208141255e+00,
    -1.0429453488027509e+00, -3.4696415708135458e-01,
     3.4696415708135830e-01,  1.0429453488027574e+00,
     1.7452473208141317e+00,  2.4586636111723683e+00,
     3.1890148165533900e+00,  3.9439673506573163e+00,
     4.7345813340460552e+00,  5.5787388058932033e+00,
     6.5105901570136551e+00,  7.6190485416797591e+00
  });
  vector[N_GH] gh_weights = to_vector({
    1.2578006724378954e-13, 2.4820623623151972e-10,
    6.1274902599825256e-08, 4.4021210902309806e-06,
    1.2882627996193093e-04, 1.8301031310804826e-03,
    1.3997837447100857e-02, 6.1506372063977507e-02,
    1.6173933398399959e-01, 2.6079306344955683e-01,
    2.6079306344955305e-01, 1.6173933398399776e-01,
    6.1506372063977438e-02, 1.3997837447101162e-02,
    1.8301031310805052e-03, 1.2882627996193072e-04,
    4.4021210902309052e-06, 6.1274902599829068e-08,
    2.4820623623151936e-10, 1.2578006724379269e-13
  });
  vector[N_GH] log_terms;
  real log_choose = log(choose(max_rank - 1, rank_pos - 1));

  for (i in 1:N_GH) {
    real eta = dprime + sigma * gh_nodes[i];
    log_terms[i] = log(gh_weights[i])
      + (max_rank - rank_pos) * std_normal_lcdf(eta)
      + (rank_pos - 1) * std_normal_lccdf(eta);
  }

  real log_p = log_choose + log_sum_exp(log_terms);
  return y * log_p;
}

// Equal-variance normal ranking: delegates to UV with sdratio = 0 (sigma = 1).
// Separate family avoids passing sdratio as a brms constant dpar, which
// causes degenerate HMC geometry.
real sdt_ranking_ev_lpmf(int y, real mu, real dprime,
                         int rank_pos, int max_rank) {
  return sdt_ranking_uv_lpmf(y | mu, dprime, 0.0, rank_pos, max_rank);
}
