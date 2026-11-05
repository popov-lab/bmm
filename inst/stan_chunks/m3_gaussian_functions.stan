// log Phi with exact gradient; std_normal_lcdf only where 1 / Phi would overflow the reverse pass
real m3_gauss_log_Phi(real z) {
  real p = Phi(z);
  return p >= 1e-300 ? log(p) : std_normal_lcdf(z | );
}
// Gaussian choice rule: log P(k) = log n_k + log E_z[Phi(z)^(n_k - 1) prod_{j != k} Phi(z + A_k - A_j)^n_j] by Gauss-Hermite quadrature
real m3_gauss_logp_vec(int k, vector A, vector n, data vector nodes, data vector log_w, data vector log_Phi_nodes) {
  int K = num_elements(A);
  int Q = num_elements(nodes);
  // a category without options carries n = 0.0001; -100 keeps 0 * value finite
  if (n[k] < 0.5) return -100;
  vector[Q] terms = log_w + (n[k] - 1) * log_Phi_nodes;
  for (j in 1:K) {
    if (j == k || n[j] < 0.5) continue;
    real d = A[k] - A[j];
    for (i in 1:Q) terms[i] += n[j] * m3_gauss_log_Phi(nodes[i] + d);
  }
  return log(n[k]) + log_sum_exp(terms);
}
