// Specify additional hyperbolic functions
real csch (real x) {
  return 1/sinh(x);
}

real coth (real x) {
  return cosh(x) / sinh(x);
}

// Specify likelihood for ezDM
real ezdm_4par_lpdf (real mrt_upper, real mu, real drift, real bound, real ndt, real zr, real s, real mrt_lower, real vrt_upper, real vrt_lower, int hits, int trials) {
  // compute misses
  int misses = trials - hits;

  // re-scale dm parameters for moments calculation
  real v = drift;
  real a = bound/2;
  real z = zr * bound - a;

  // get standardized boundary and start point
  real k_z = v * a / square(s);
  real k_x = v * z / square(s);

  // get proportion correct
  real pC = 1 - (exp(-2 * k_x) - exp(-2 * k_z)) / (exp(2 * k_z) - exp(-2 * k_z));

  // drift is small use functions not relying on drift rate
  if(abs(drift) < 1e-8) {
    real mdt_upper_drift0 = (4*a^2 - square(a + z))/(3*square(s));
    real vrt_upper_drift0 = (32*a^4 - 2*(a + z)^4)/(45*s^4);
    real mdt_lower_drift0 = (4*a^2 -(a - z)^2)/(3*square(s));
    real vrt_lower_drift0 = (32*a^4 - 2*(a - z)^4)/(45*s^4);

    if (misses >= 2 && hits >= 2) {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_drift0, sqrt(vrt_upper_drift0/hits)) +
             gamma_lpdf(vrt_upper | ((hits-1)/2.0), ((hits-1)/(2*vrt_upper_drift0))) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_drift0, sqrt(vrt_lower_drift0/misses)) +
             gamma_lpdf(vrt_lower | ((misses-1)/2.0), ((misses-1)/(2*vrt_lower_drift0)));
    } else if(misses < 2) {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_drift0, sqrt(vrt_upper_drift0/hits)) +
             gamma_lpdf(vrt_upper | ((hits-1)/2.0), ((hits-1)/(2*vrt_upper_drift0)));
    } else {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_drift0, sqrt(vrt_lower_drift0/misses)) +
             gamma_lpdf(vrt_lower | ((misses-1)/2.0), ((misses-1)/(2*vrt_lower_drift0)));
    }
  } else {
    real mdt_upper_implied = square(s)/(v^2)*(2*k_z*coth(2*k_z) - (k_x+k_z)*coth(k_x+k_z));
    real vrt_upper_implied = s^4/(v^4)*(4*k_z^2*square(csch(2*k_z)) + 2*k_z*coth(2*k_z) - (k_z+k_x)^2*square(csch(k_z+k_x)) - (k_z+k_x)*coth(k_z+k_x));
    real mdt_lower_implied = square(s)/(v^2)*(2*k_z*coth(2*k_z) - (-k_x+k_z)*coth(-k_x+k_z));
    real vrt_lower_implied = s^4/(v^4)*(4*k_z^2*square(csch(2*k_z)) + 2*k_z*coth(2*k_z) - (k_z-k_x)^2*square(csch(k_z-k_x)) - (k_z-k_x)*coth(k_z+-k_x));

    // return sum of sample statistics distributions log-likelihood
    if (misses >= 2 && hits >= 2) {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_implied, sqrt(vrt_upper_implied/hits)) +
             gamma_lpdf(vrt_upper | ((hits-1)/2.0), ((hits-1)/(2*vrt_upper_implied))) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_implied, sqrt(vrt_lower_implied/misses)) +
             gamma_lpdf(vrt_lower | ((misses-1)/2.0), ((misses-1)/(2*vrt_lower_implied)));
    } else if(misses < 2) {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_implied, sqrt(vrt_upper_implied/hits)) +
             gamma_lpdf(vrt_upper | ((hits-1)/2.0), ((hits-1)/(2*vrt_upper_implied)));
    } else {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_implied, sqrt(vrt_lower_implied/misses)) +
             gamma_lpdf(vrt_lower | ((misses-1)/2.0), ((misses-1)/(2*vrt_lower_implied)));
    }
  }
}
