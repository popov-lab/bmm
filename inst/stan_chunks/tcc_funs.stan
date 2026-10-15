  // TCC (Schurgin, Wixted & Brady, 2020): log P of one option at each distance, by Gauss-Legendre in x
  vector tcc_logp_by_distance(real d, vector similarity, vector gl_nodes, vector gl_logw) {
    int D = rows(similarity);
    vector[D] mult = append_row(append_row(1, rep_vector(2, D - 2)), 1);
    vector[D] means = d * similarity;
    real top = max(means);
    // see .tcc_logp(): a far option under a sharp similarity peaks halfway to the largest mean
    real half = (top + 8 - fmin(top - 8, 0.5 * (min(means) + top) - 6)) / 2;
    vector[rows(gl_nodes)] x = top + 8 - half + half * gl_nodes;
    vector[rows(gl_nodes)] logw = gl_logw + log(half) - 0.5 * log(2 * pi());
    vector[D] out = rep_vector(negative_infinity(), D);
    for (q in 1:rows(gl_nodes)) {
      vector[D] z = x[q] - means;
      vector[D] lphi = log(Phi(z));
      out = log_sum_exp(out, logw[q] + dot_product(mult, lphi) - 0.5 * square(z) - lphi);
    }
    return out;
  }

  // first of the four table rows that interpolate grid position u (row 2 holds the lowest grid value)
  int tcc_first_row(real u, int n_rows) {
    int lo = 1;
    int hi = n_rows - 3;
    while (lo < hi) {
      int mid = (lo + hi + 1) %/% 2;
      if (mid + 1 <= u) lo = mid; else hi = mid - 1;
    }
    return lo;
  }

  // Catmull-Rom weights of the four rows at t in [0, 1]; outside the grid, linear in the end slope
  vector tcc_row_weights(real t) {
    if (t < 0) return [-0.5 * t, 1, 0.5 * t, 0]';
    if (t > 1) return [0, -0.5 * (t - 1), 1, 0.5 * (t - 1)]';
    real t2 = square(t);
    real t3 = t2 * t;
    return 0.5 * [-t + 2 * t2 - t3, 2 - 5 * t2 + 3 * t3, t + 4 * t2 - 3 * t3, t3 - t2]';
  }

  // sum over rows i:j of the four table entries at each row's distance; data only, so no autodiff,
  // and scalars, because a vector slice per row cost 2.5 ms per gradient at 15,200 rows
  vector tcc_row_sums(data array[] int y, data matrix tab, int first, int i, int j) {
    real s1 = 0;
    real s2 = 0;
    real s3 = 0;
    real s4 = 0;
    for (n in i:j) {
      int col = abs(y[n]) + 1;
      s1 += tab[first, col];
      s2 += tab[first + 1, col];
      s3 += tab[first + 2, col];
      s4 += tab[first + 3, col];
    }
    return [s1, s2, s3, s4]';
  }

  // y: signed option offsets; runs of equal d share one interpolation (a run of length 1 is one trial)
  real tcc_lpmf(array[] int y, vector mu, vector d, data matrix tab, data real grid_lo, data real grid_step) {
    int N = size(y);
    real lp = 0;
    int i = 1;
    while (i <= N) {
      int j = i;
      while (j < N && d[j + 1] == d[i]) j += 1;
      real u = (d[i] - grid_lo) / grid_step + 2;
      int first = tcc_first_row(u, rows(tab));
      lp += dot_product(tcc_row_weights(u - first - 1), tcc_row_sums(y, tab, first, i, j));
      i = j + 1;
    }
    return lp;
  }
