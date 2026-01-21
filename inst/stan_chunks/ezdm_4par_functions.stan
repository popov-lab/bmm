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
  // translate into Srivastava labelling of dm parameters
  real a = drift;             // drift rate
  real z = bound/2;           // boundary
  real x0 = (zr * bound) - z; // starting point

  // get standardized boundary and start point
  real k_z = (a * z) / square(s);
  real k_x = (a * x0) / square(s);

  // get proportion correct
  real pC = 1 - ( (exp(-2 * k_x) - exp(-2 * k_z)) / (exp(2 * k_z) - exp(-2 * k_z)) );

  // drift is small: use functions not relying on drift rate
  if(abs(drift) < 1e-6) {
    real mdt_upper_drift0 = (4*z^2 - (z + x0)^2) / (3*s^2);
    real mdt_lower_drift0 = (4*z^2 - (z - x0)^2) / (3*s^2);

    real vrt_upper_drift0 = (32*z^4 - 2*(z + x0)^4)/(45*s^4);
    real vrt_lower_drift0 = (32*z^4 - 2*(z - x0)^4)/(45*s^4);

    if (misses >= 2 && hits >= 2) {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_drift0, sqrt(vrt_upper_drift0/hits)) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_drift0, sqrt(vrt_lower_drift0/misses)) +
             gamma_lpdf(vrt_upper | ((hits-1)/2.0), ((hits-1)/(2*vrt_upper_drift0))) +
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
    real mdt_upper_implied = (s^2)/(a^2) * (2*k_z*coth(2*k_z) - ( k_x+k_z) * coth( k_x+k_z));
    real mdt_lower_implied = (s^2)/(a^2) * (2*k_z*coth(2*k_z) - (-k_x+k_z) * coth(-k_x+k_z));

    real vrt_upper_implied = (s^4)/(a^4) * (4*k_z^2*(csch(2*k_z))^2 + 2*k_z*coth(2*k_z) - ( k_x+k_z)^2*(csch( k_x+k_z))^2 - ( k_x+k_z)*coth( k_x+k_z));
    real vrt_lower_implied = (s^4)/(a^4) * (4*k_z^2*(csch(2*k_z))^2 + 2*k_z*coth(2*k_z) - (-k_x+k_z)^2*(csch(-k_x+k_z))^2 - (-k_x+k_z)*coth(-k_x+k_z));

    // return sum of sample statistics distributions log-likelihood
    if (misses >= 2 && hits >= 2) {
      return binomial_lpmf(hits | trials, pC) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_implied, sqrt(vrt_upper_implied/hits)) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_implied, sqrt(vrt_lower_implied/misses)) +
             gamma_lpdf(vrt_upper | ((hits-1)/2.0), ((hits-1)/(2*vrt_upper_implied))) +
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
