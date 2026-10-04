  // log peak of the sdm activation; c arrives on the log scale, as in sdm_simple_funs.stan
  real sdm_cd_log_peak(real c, real kappa) {
    return c + 0.5 * log(kappa * inv(2 * pi()));
  }

  // half-width of the sdm spike, beyond which the activation has settled onto its floor
  real sdm_cd_spike_width(real c, real kappa) {
    return fmin(pi(), sqrt(2 * fmax(sdm_cd_log_peak(c, kappa), 1) * inv(kappa)));
  }

  // log integral of the unnormalised sdm density over [lo, hi], split at the peak images and spike-width multiples
  real sdm_cd_log_int(real lo, real hi, real mu, real c, real kappa,
                      data vector gl_x, data vector gl_w) {
    int n_nodes = num_elements(gl_x);
    real w = sdm_cd_spike_width(c, kappa);
    real peak = exp(sdm_cd_log_peak(c, kappa));
    array[9] real mult = {-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3};
    vector[29] cuts;
    vector[n_nodes] log_gl_w = log(gl_w);
    real out = negative_infinity();
    int k = 3;

    cuts[1] = lo;
    cuts[2] = hi;
    for (p in -1:1) {
      for (j in 1:9) {
        cuts[k] = fmin(fmax(mu + p * 2 * pi() + mult[j] * w, lo), hi);
        k += 1;
      }
    }
    cuts = sort_asc(cuts);

    for (i in 1:28) {
      if (cuts[i + 1] > cuts[i]) {
        real half = (cuts[i + 1] - cuts[i]) / 2;
        vector[n_nodes] u = (cuts[i] + cuts[i + 1]) / 2 + half * gl_x - mu;
        out = log_sum_exp(out, log_sum_exp(log(half) + log_gl_w + peak * exp(kappa * (cos(u) - 1))));
      }
    }
    return out;
  }

  // half-width of the "same" arc: where the sdm density equals exp(-criterion) / (2 pi)
  real sdm_cd_crit_angle(real c, real kappa, real criterion, real log_z) {
    real log_peak = sdm_cd_log_peak(c, kappa);
    real thresh = log_z - criterion - log(2 * pi());
    real log_thresh;
    if (thresh <= 0) return pi();
    log_thresh = log(thresh);
    if (log_thresh <= log_peak - 2 * kappa) return pi();
    if (log_thresh >= log_peak) return 0;
    return acos(fmin(fmax(1 + (log_thresh - log_peak) * inv(kappa), -1), 1));
  }

  // c, kappa and mu are constant within a run (data sorted by predictors), so log Z is computed once per run
  real sdm_cd_lpmf(array[] int y, vector mu, vector c, vector kappa,
                   vector criterion, data array[] real probe, int G_runs,
                   data array[] int run_start, data array[] int run_count,
                   data vector gl_x, data vector gl_w) {
    real out = 0;
    for (g in 1:G_runs) {
      int s = run_start[g];
      real log_z = sdm_cd_log_int(-pi(), pi(), 0, c[s], kappa[s], gl_x, gl_w);
      for (n in s:(s + run_count[g] - 1)) {
        real hw = sdm_cd_crit_angle(c[n], kappa[n], criterion[n], log_z);
        real p_same = 1;
        if (hw <= 0) {
          p_same = 0;
        } else if (hw < pi()) {
          p_same = exp(sdm_cd_log_int(probe[n] - hw, probe[n] + hw, mu[n], c[n],
                                      kappa[n], gl_x, gl_w) - log_z);
        }
        out += cd_bernoulli_lpmf(y[n] | p_same);
      }
    }
    return out;
  }
