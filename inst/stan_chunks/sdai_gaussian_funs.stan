  // Integrand for Gaussian SDAI source confusion
  // theta = {d1, s1, d2, s2}
  real sdai_int_inst_uvg(real x, real xc,
                         array[] real theta,
                         array[] real x_r,
                         array[] int x_i) {
    real d1 = theta[1];
    real s1 = theta[2];
    real d2 = theta[3];
    real s2 = theta[4];
    return normal_cdf(x | d2, s2) * exp(normal_lpdf(x | d1, s1));
  }

  // Gaussian discrimination probability (both sources standard normal)
  real sdai_p_hit_gaussian(real l, real u) {
    return std_normal_cdf(u) * std_normal_cdf(u) -
           std_normal_cdf(l) * std_normal_cdf(l);
  }
