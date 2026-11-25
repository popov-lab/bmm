  #include 'fun_tan_half.stan'

  // log of the numerator of the SDM likelihood, summed over all observations
  real sdm_simple_lpdf(vector y, vector mu, vector c, vector kappa) {
    int N = size(y);
    vector[N] num = exp(fma(kappa,cos(y-mu)-1,c)) ;
    real out = dot_product(num, sqrt(kappa));
    out *= inv(sqrt2()) * inv_sqrt(pi());
    return(out);
  }

  // normalising term of the SDM log likelihood for one combination of c and kappa
  real sdm_simple_ldenom(real c, real kappa, data vector gl_x, data vector gl_w,
                         data vector wide_cosm1, data vector wide_w) {
    int n_nodes = num_elements(gl_x);
    real log_peak = c + 0.5 * log(kappa / (2 * pi()));
    real peak = exp(log_peak);
    if (kappa * (peak + 1) <= 4000) {
      return log(pi()) - peak - log(dot_product(wide_w, exp(peak * expm1(kappa * wide_cosm1))));
    }
    real w = fmin(pi(), sqrt(2 * fmax(log_peak, 1) / kappa));
    real s = exp(-0.5 * (log_peak + log(kappa)));
    vector[9] cuts = [0, pi(), 0.5 * w, w, 2 * w, 5 * w, s, 3 * s, 9 * s]';
    real log_int = negative_infinity();
    for (i in 1:9) {
      cuts[i] = fmin(fmax(cuts[i], 0), pi());
    }
    cuts = sort_asc(cuts);
    for (i in 1:8) {
      if (cuts[i + 1] > cuts[i]) {
        real half = (cuts[i + 1] - cuts[i]) / 2;
        vector[n_nodes] u = (cuts[i] + cuts[i + 1]) / 2 + half * gl_x;
        vector[n_nodes] d = peak * expm1(-2 * kappa * square(sin(u / 2)));
        real d_max = max(d);
        log_int = log_sum_exp(log_int, log(half) + d_max + log(dot_product(gl_w, exp(d - d_max))));
      }
    }
    return log(pi()) - peak - log_int;
  }

  // normalising term summed over runs of observations that share c and kappa
  real sdm_simple_run_ldenom(vector c, vector kappa, data vector gl_x,
                             data vector gl_w, data vector wide_cosm1,
                             data vector wide_w, int G_runs,
                             array[] int run_start, array[] int run_count) {
    real out = 0;
    for (g in 1:G_runs) {
      int n = run_start[g];
      out += run_count[g] * sdm_simple_ldenom(c[n], kappa[n], gl_x, gl_w, wide_cosm1, wide_w);
    }
    return(out);
  }

  // the same sum for the observations from start to end, used with threading
  real sdm_simple_run_ldenom_slice(vector c, vector kappa, data vector gl_x,
                                   data vector gl_w, data vector wide_cosm1,
                                   data vector wide_w, int start, int end,
                                   int G_runs, array[] int run_start,
                                   array[] int run_count) {
    real out = 0;
    for (g in 1:G_runs) {
      if (run_start[g] >= start && run_start[g] <= end) {
        int n = run_start[g] - start + 1;
        out += run_count[g] * sdm_simple_ldenom(c[n], kappa[n], gl_x, gl_w, wide_cosm1, wide_w);
      }
    }
    return(out);
  }
