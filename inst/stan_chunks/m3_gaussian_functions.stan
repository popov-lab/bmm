// log of the standard normal CDF, stable in the far lower tail
real m3_gauss_log_Phi(real z) {
  real p = Phi(z);
  return p >= 1e-300 ? log(p) : std_normal_lcdf(z | );
}
// log probability of response category k under the Gaussian choice rule, by quadrature
real m3_gauss_logp_vec(int k, vector A, vector n, data vector nodes, data vector log_w, data vector log_Phi_nodes) {
  int K = num_elements(A);
  int Q = num_elements(nodes);
  if (n[k] < 0.5) return -100;
  vector[Q] terms = log_w + (n[k] - 1) * log_Phi_nodes;
  for (j in 1:K) {
    if (j == k || n[j] < 0.5) continue;
    real d = A[k] - A[j];
    for (i in 1:Q) terms[i] += n[j] * m3_gauss_log_Phi(nodes[i] + d);
  }
  return log(n[k]) + log_sum_exp(terms);
}
