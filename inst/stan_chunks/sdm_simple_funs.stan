  #include 'fun_tan_half.stan'

  // utility function trick for converting real to integer type
  int bin_search(real x, int min_val, int max_val) {
      int mid_p;
      int L = min_val;
      int R = max_val;
      while(L < R) {
        mid_p = (R-L) %/% 2;
        if (L + mid_p < x) {
          L += mid_p + 1;
        } else if (L + mid_p > x) {
          R = L + mid_p - 1;
        } else {
          return(L + mid_p);
        }
      }
      return(L);
    }

  // utility function for determining optimal number of chebyshev points for the denominator approximation
  int get_m(real c, real kappa) {
    real m = floor(2 * exp(0.4*c) * kappa^(fma(c,0.0145,0.7)) + 0.5)+2;
    int M = bin_search(m, 2, 200);
    return(M);
  }

  // log of the numerator of the sdm likelihood
  real sdm_simple_lpdf(vector y, vector mu, vector c, vector kappa) {
    int N = size(y);
    vector[N] num = exp(fma(kappa,cos(y-mu)-1,c)) ;
    real out = dot_product(num, sqrt(kappa));
    out *= inv(sqrt2()) * inv_sqrt(pi());
    return(out);
  }

  // log of the normalization constant, approximated by chebyshev quadrature
  real sdm_simple_ldenom_chquad_adaptive(real c, real kappa, matrix CN) {
    int m = get_m(c,kappa);
    vector[m] cosn = CN[1:m,m];
    vector[m] fn = exp(fma(kappa,cosn,c)) * sqrt(kappa) * inv(sqrt2()) * inv_sqrt(pi());
    real out = -log_sum_exp(fn)+log(m);
    return(out);
  }

  real sdm_simple_run_ldenom(vector c, vector kappa, matrix CN,
                             int G_runs, array[] int run_start,
                             array[] int run_count) {
    real out = 0;
    for (g in 1:G_runs) {
      int n = run_start[g];
      real z = sdm_simple_ldenom_chquad_adaptive(c[n], kappa[n], CN);
      out += run_count[g] * z;
    }
    return(out);
  }
