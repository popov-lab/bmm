  // y: signed option offsets; runs of equal (d, par) share one similarity and normalising sum.
  // Scaled to 1 at the target and 0 opposite it, both shapes are c + expm1(-rate b(delta)) / -expm1(-rate b(pi)),
  // b = delta and rate = 1 / lambda (Laplace) or b = 1 - cos(delta) and rate = kappa (von Mises); the
  // softmax drops c, which leaves one vector expression per run instead of tcc_similarity()'s five
  real tcc_gumbel_parametric_lpmf(array[] int y, vector mu, vector d, vector par, data vector delta,
                                  data int shape) {
    int N = size(y);
    int D = rows(delta);
    vector[D] log_mult = tcc_log_mult(D);
    vector[D] b = shape == 1 ? delta : 2 * square(sin(0.5 * delta));
    real lp = 0;
    int i = 1;
    while (i <= N) {
      int j = i;
      while (j < N && d[j + 1] == d[i] && par[j + 1] == par[i]) j += 1;
      real rate = shape == 1 ? inv(par[i]) : par[i];
      vector[D] h = expm1(-rate * b);
      real a = d[i] / -expm1(-rate * b[D]);
      real h_sum = 0;
      for (n in i:j) h_sum += h[abs(y[n]) + 1];
      lp += a * h_sum - (j - i + 1) * log_sum_exp(log_mult + a * h);
      i = j + 1;
    }
    return lp;
  }
