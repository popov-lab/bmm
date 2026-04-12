// m-AFC SDT likelihood using Gauss-Hermite quadrature
// P(correct | m, d') = integral phi(x - d') * Phi(x)^(m-1) dx
//                     = integral phi(u) * Phi(u + d')^(m-1) du
//
// For m = 2: closed form Phi(d' / sqrt(2))
// For m >= 3: 20-point Gauss-Hermite quadrature
//
// Reference: DeCarlo (2012); Green & Swets (1966)

real mafc_pc(real dprime, int m) {
  if (m == 2) return Phi(dprime / sqrt(2.0));

  // 20-point Gauss-Hermite nodes/weights for standard normal integration
  // integral f(x) * dnorm(x) dx ≈ sum_k weights[k] * f(nodes[k])
  int N_QUAD = 20;
  vector[N_QUAD] gh_nodes = to_vector({
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
  vector[N_QUAD] gh_weights = to_vector({
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

  vector[N_QUAD] log_terms;
  for (i in 1:N_QUAD) {
    log_terms[i] = log(gh_weights[i]) +
      (m - 1) * std_normal_lcdf(gh_nodes[i] + dprime);
  }
  return exp(log_sum_exp(log_terms));
}

real sdt_mafc_lpmf(int y, real mu, real dprime, int m, int trials) {
  real pc = mafc_pc(dprime, m);
  return binomial_lpmf(y | trials, pc);
}
