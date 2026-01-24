// Specify likelihood for ezDM - zr fixed to 0.5
real ezdm_3par_lpdf (real mrt, real mu, real drift, real bound, real ndt, real s, real vrt, int hits, int trials) {
  // Cache common calculations
  real bound_sq = square(bound);
  real s_sq = square(s);

  // drift is small: use zero-drift formulas
  if (abs(drift) < 1e-6) {
    // Cache variance term (used in both normal and gamma)
    real var_rt = (bound_sq * bound_sq) / (24 * s_sq * s_sq);
    real sd_mrt = sqrt(var_rt / trials);
    real mean_dt = bound_sq / (4 * s_sq);

    return binomial_lpmf(hits | trials, 0.5) +
           normal_lpdf(mrt | ndt + mean_dt, sd_mrt) +
           gamma_lpdf(vrt | ((trials - 1) / 2.0), ((trials - 1) / (2 * var_rt)));
  }

  // compute helper variables
  real y = -(bound * drift) / s_sq;
  real expy = exp(y);
  real expy_plus_1 = expy + 1;
  real exp2y = expy * expy;  // exp(2*y) = exp(y)^2, faster than exp(2*y)

  // proportion correct (probability of hitting upper boundary)
  real pC = 1 / expy_plus_1;

  // Mean decision time
  real MDT = (bound / (2 * drift)) * ((1 - expy) / expy_plus_1);

  // Variance of decision time
  real VRT = ((bound * s_sq) / (2 * drift * drift * drift)) *
             (2 * y * expy - exp2y + 1) / square(expy_plus_1);

  // return sum of sample statistics distributions log-likelihood
  return binomial_lpmf(hits | trials, pC) +
         normal_lpdf(mrt | ndt + MDT, sqrt(VRT / trials)) +
         gamma_lpdf(vrt | ((trials - 1) / 2.0), ((trials - 1) / (2 * VRT)));
}
