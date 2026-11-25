  // log number of options at each distance: one at the target and at the far end, two elsewhere
  vector tcc_log_mult(int D) {
    return append_row(append_row(0, rep_vector(log2(), D - 2)), 0);
  }

  // TCC with Gumbel noise: log likelihood of the responses, signed option offsets from the target
  real tcc_gumbel_lpmf(array[] int y, vector mu, vector d, data vector similarity) {
    int N = size(y);
    vector[rows(similarity)] log_mult = tcc_log_mult(rows(similarity));
    real lp = 0;
    int i = 1;
    while (i <= N) {
      int j = i;
      while (j < N && d[j + 1] == d[i]) j += 1;
      real f_sum = 0;
      for (n in i:j) f_sum += similarity[abs(y[n]) + 1];
      lp += d[i] * f_sum - (j - i + 1) * log_sum_exp(log_mult + d[i] * similarity);
      i = j + 1;
    }
    return lp;
  }
