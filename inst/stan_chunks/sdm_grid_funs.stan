  // sdm() with response options (#558): the response is one of K equally
  // spaced options, y its offset from the target in steps, and its probability
  // the softmax of c* s(x - mu) over the options the trial offered. Every
  // log-likelihood is reported on the density scale, log P - log(2 pi / K),
  // as in the continuous model; sdm_kernel_funs.stan defines the cell

  // the numerator: c* s(y - mu) per trial
  real sdm_grid_lpmf(array[] int y, vector mu, vector c, vector kappa, data int K,
                     data int activation, data int scaling) {
    vector[size(y)] delta = to_vector(y) * (2 * pi() / K) - mu;
    return sum(exp(sdm_log_cstar(c, kappa, activation, scaling) + sdm_log_s(delta, kappa, activation)));
  }

  // The log of the sum over the options of a cell of exp(c* s(x - mu)), as
  // peak + log sum exp(peak expm1(log s)), which loses no digits at a large
  // peak. Set 0 is the full wheel. There the von Mises sum is the continuous
  // normaliser over the step wherever the activation spans 2.25 steps, w = 1 /
  // sqrt(max(c*, 1) kappa) >= 2.25 * 2 pi / K: measured within 5e-12 in log
  // over 2,745 cells, where the grid sum would cost K terms per cell. The
  // trapezoid error falls as exp(-2.2 (w / step)^2), not faster, because the
  // powers of s in exp(c* s) are narrower than s. With mu = 0 the full wheel
  // takes each distance once, with its multiplicity
  real sdm_grid_log_sum(real mu, real c, real kappa, int activation, int scaling, int K, int set,
                        data array[] int set_size, data array[,] int set_off,
                        data vector grid_half, data vector grid_logmult, data vector grid_full,
                        data vector gl_x, data vector gl_w, data vector wide_cosm1,
                        data vector wide_w) {
    real log_peak = sdm_log_cstar(c, kappa, activation, scaling);
    real peak = exp(log_peak);
    if (set > 0) {
      int m = set_size[set];
      vector[m] delta = to_vector(set_off[set, 1:m]) * (2 * pi() / K) - mu;
      return peak + log_sum_exp(peak * expm1(sdm_log_s(delta, rep_vector(kappa, m), activation)));
    }
    if (activation == 1 && fmax(log_peak, 0) + log(kappa) + 2 * log(2.25 * 2 * pi() / K) <= 0) {
      return log(K) - sdm_vm_ldenom_peak(log_peak, kappa, gl_x, gl_w, wide_cosm1, wide_w);
    }
    if (mu == 0) {
      int n = rows(grid_half);
      return peak + log_sum_exp(grid_logmult
                                + peak * expm1(sdm_log_s(grid_half, rep_vector(kappa, n), activation)));
    }
    return peak + log_sum_exp(peak * expm1(sdm_log_s(grid_full - mu, rep_vector(kappa, K), activation)));
  }

  // the normalisers of all cells, each at its first row, weighted by its rows
  real sdm_grid_cells_lse(vector mu, vector c, vector kappa, int activation, int scaling, int K,
                          int G_cells, array[] int cell_first, array[] int cell_count,
                          array[] int cell_set, data array[] int set_size,
                          data array[,] int set_off, data vector grid_half,
                          data vector grid_logmult, data vector grid_full, data vector gl_x,
                          data vector gl_w, data vector wide_cosm1, data vector wide_w) {
    real out = 0;
    for (g in 1:G_cells) {
      int n = cell_first[g];
      out += cell_count[g] * sdm_grid_log_sum(mu[n], c[n], kappa[n], activation, scaling, K,
                                              cell_set[g], set_size, set_off, grid_half,
                                              grid_logmult, grid_full, gl_x, gl_w, wide_cosm1,
                                              wide_w);
    }
    return out;
  }

  // the cells whose first row lies in the slice start:end
  real sdm_grid_cells_lse_slice(vector mu, vector c, vector kappa, int activation, int scaling,
                                int K, int start, int end, int G_cells,
                                array[] int cell_first, array[] int cell_count,
                                array[] int cell_set, data array[] int set_size,
                                data array[,] int set_off, data vector grid_half,
                                data vector grid_logmult, data vector grid_full,
                                data vector gl_x, data vector gl_w, data vector wide_cosm1,
                                data vector wide_w) {
    real out = 0;
    for (g in 1:G_cells) {
      if (cell_first[g] >= start && cell_first[g] <= end) {
        int n = cell_first[g] - start + 1;
        out += cell_count[g] * sdm_grid_log_sum(mu[n], c[n], kappa[n], activation, scaling, K,
                                                cell_set[g], set_size, set_off, grid_half,
                                                grid_logmult, grid_full, gl_x, gl_w,
                                                wide_cosm1, wide_w);
      }
    }
    return out;
  }
