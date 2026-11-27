  // TCC with forced-choice option sets: log probability of each offered option, given its similarity f
  vector tcc_set_logp(real d, vector f, int noise, data vector gl_nodes, data vector gl_logw) {
    int m = rows(f);
    vector[m] means = d * f;
    if (noise == 2) return means - log_sum_exp(means);
    real top = max(means);
    real half = (top + 8 - fmin(top - 8, 0.5 * (min(means) + top) - 6)) / 2;
    vector[rows(gl_nodes)] x = top + 8 - half + half * gl_nodes;
    vector[rows(gl_nodes)] logw = gl_logw + log(half) - 0.5 * log(2 * pi());
    vector[m] out = rep_vector(negative_infinity(), m);
    for (q in 1:rows(gl_nodes)) {
      vector[m] z = x[q] - means;
      vector[m] lphi;
      for (k in 1:m) lphi[k] = z[k] > -37 ? log(Phi(z[k])) : std_normal_lcdf(z[k] | );
      out = log_sum_exp(out, logw[q] + sum(lphi) - 0.5 * square(z) - lphi);
    }
    return out;
  }

  // adds nothing: the responses enter the likelihood through tcc_options_cells()
  real tcc_options_lpmf(array[] int y, vector mu, vector d) {
    return 0;
  }

  // log likelihood of the cells whose first row lies in start to end, weighted by their choices
  real tcc_options_cells(vector d, int start, int end, int noise, data array[] int cell_first,
                         data array[] int cell_set, data array[,] int cell_n,
                         data array[] int set_size, data array[,] int set_off,
                         data vector similarity, data vector gl_nodes, data vector gl_logw) {
    real lp = 0;
    for (g in 1:size(cell_first)) {
      int n = cell_first[g];
      if (n >= start && n <= end) {
        int s = cell_set[g];
        int m = set_size[s];
        vector[m] f;
        for (j in 1:m) f[j] = similarity[abs(set_off[s, j]) + 1];
        lp += dot_product(to_vector(cell_n[g, 1:m]),
                          tcc_set_logp(d[n - start + 1], f, noise, gl_nodes, gl_logw));
      }
    }
    return lp;
  }
