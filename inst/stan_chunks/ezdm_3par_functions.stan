// Specify likelihood for ezDM - zr fixed to 0.5
real ezdm_3par_lpdf (real mrt, real mu, real drift, real bound, real ndt, real s, real vrt, int hits, int trials) {
  // Cache common calculations
  real bound_sq = square(bound);
  real s_sq = square(s);

  // Declare all variables at the top
  real drift_abs;
  real pC;
  real y;
  real expy;
  real expy_plus_1;
  real exp2y;
  real MDT;
  real VRT;

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


  // proportion correct (probability of hitting upper boundary)
  pC = 1 / (exp(-bound * drift / s_sq) + 1);

  // Use soft absolute value: sqrt(drift^2 + tau^2) with larger tau
  // This avoids extreme curvature while maintaining smooth gradients
  drift_abs = sqrt(drift * drift + 0.0001);  // tau = 0.01, tau^2 = 0.0001
  y = -bound * drift_abs / s_sq;
  expy = exp(y);
  expy_plus_1 = expy + 1;
  exp2y = expy * expy;

  // Mean decision time
  MDT = (bound / (2 * drift_abs)) * ((1 - expy) / expy_plus_1);

  // Variance of decision time
  VRT = ((bound * s_sq) / (2 * drift_abs * drift_abs * drift_abs)) *
             (2 * y * expy - exp2y + 1) / square(expy_plus_1);

  // return sum of sample statistics distributions log-likelihood
  return binomial_lpmf(hits | trials, pC) +
         normal_lpdf(mrt | ndt + MDT, sqrt(VRT / trials)) +
         gamma_lpdf(vrt | ((trials - 1) / 2.0), ((trials - 1) / (2 * VRT)));
}
