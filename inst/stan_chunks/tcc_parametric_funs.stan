  // TCC similarity at distances delta (radians), 1 at 0 and 0 at pi: shape 1 Laplace (par = lambda), 2 von Mises (par = kappa)
  vector tcc_similarity(vector delta, real par, int shape) {
    if (shape == 1) return exp(-delta / par) .* expm1(-(pi() - delta) / par) / expm1(-pi() / par);
    vector[rows(delta)] s = square(sin(0.5 * delta));
    return exp(-2 * par * s) .* expm1(-2 * par * (1 - s)) / expm1(-2 * par);
  }

  // y: signed option offsets; runs of equal (d, par) share one read of the table over (d, log par)
  real tcc_parametric_lpmf(array[] int y, vector mu, vector d, vector par, data array[] matrix tab,
                           data real grid_lo, data real grid_step, data real u_lo, data real u_step) {
    int N = size(y);
    real lp = 0;
    int i = 1;
    while (i <= N) {
      int j = i;
      while (j < N && d[j + 1] == d[i] && par[j + 1] == par[i]) j += 1;
      real gd = (d[i] - grid_lo) / grid_step + 2;
      real gu = (log(par[i]) - u_lo) / u_step + 2;
      int fd = tcc_first_row(gd, rows(tab[1]));
      int fu = tcc_first_row(gu, size(tab));
      vector[4] wd = tcc_row_weights(gd - fd - 1);
      vector[4] wu = tcc_row_weights(gu - fu - 1);
      for (b in 1:4) lp += wu[b] * dot_product(wd, tcc_row_sums(y, tab[fu + b - 1], fd, i, j));
      i = j + 1;
    }
    return lp;
  }
