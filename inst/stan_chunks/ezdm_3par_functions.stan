// Specify likelihood for ezDM - zr fixed to 0.5
real ezdm_3par_lpdf (real mrt, real mu, real drift, real bound, real ndt, real s, real vrt, int hits, int trials) {
  if (abs(drift) < 1e-8) {
    return binomial_lpmf(hits | trials, 0.5) +
    normal_lpdf(mrt | ndt + (square(bound) / (4 * square(s))), sqrt(pow(bound,4) / (24 * pow(s,4))/trials)) +
    gamma_lpdf(vrt | ((trials-1)/2.0), ((trials-1)/(pow(bound,4) / (24 * pow(s,4)))));
  }
  // compute helper variables
  real y   = -(bound * drift) / square(s);
  real expy   = exp(y);

  // proportion correct
  real pC  = 1 / (1 + expy);

  // Mean decision time
  real MDT = (bound / (2 * drift)) * ((1 - expy) / (1 + expy));

  // Variance of decision time
  real VRT = ((bound * square(s)) / (2 * pow(drift, 3))) *
             (2 * y * expy - exp(2*y) + 1) / (square(expy + 1));

  // return sum of sample statistics distributions log-likelihood
  return binomial_lpmf(hits | trials, pC) +
         normal_lpdf(mrt | ndt + MDT, sqrt(VRT/trials)) +
         gamma_lpdf(vrt | ((trials-1)/2.0), ((trials-1)/(2*VRT)));
}
