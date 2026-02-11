// Numerically stable hyperbolic functions
real csch_stable(real x) {
  // For large |x|, sinh(x) overflows but csch(x) -> 0
  if (abs(x) > 20) {
    return 2 * exp(-abs(x)) * ((x > 0) ? 1 : -1);
  }
  return 1 / sinh(x);
}

real coth_stable(real x) {
  // For large |x|, coth(x) -> sign(x)
  if (abs(x) > 20) {
    return (x > 0) ? 1 : -1;
  }
  return cosh(x) / sinh(x);
}

// Specify likelihood for ezDM
real ezdm_4par_lpdf(real mrt_upper, real mu, real drift, real bound, real ndt, real zr, real s, real mrt_lower, real vrt_upper, real vrt_lower, int hits, int trials) {
  // Declare all variables at the top
  int misses;
  real s_sq;
  real z;
  real x0;
  
  // compute misses
  misses = trials - hits;

  // Cache common calculations
  s_sq = square(s);
  z = bound / 2;           // boundary (half of bound separation)
  x0 = (zr * bound) - z;   // starting point relative to center

  // drift is small: use functions not relying on drift rate
  if (abs(drift) < 1e-6) {
    // Cache intermediate values
    real z_sq = square(z);
    real s_sq_sq = s_sq * s_sq;
    real z_plus_x0 = z + x0;
    real z_minus_x0 = z - x0;

    real mdt_upper_drift0 = (4 * z_sq - square(z_plus_x0)) / (3 * s_sq);
    real mdt_lower_drift0 = (4 * z_sq - square(z_minus_x0)) / (3 * s_sq);

    real z_sq_sq = z_sq * z_sq;
    real vrt_upper_drift0 = (32 * z_sq_sq - 2 * square(square(z_plus_x0))) / (45 * s_sq_sq);
    real vrt_lower_drift0 = (32 * z_sq_sq - 2 * square(square(z_minus_x0))) / (45 * s_sq_sq);

    // When drift is small, pC -> zr (limit as drift approaches zero)
    if (misses >= 2 && hits >= 2) {
      return binomial_lpmf(hits | trials, zr) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_drift0, sqrt(vrt_upper_drift0 / hits)) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_drift0, sqrt(vrt_lower_drift0 / misses)) +
             gamma_lpdf(vrt_upper | ((hits - 1) / 2.0), ((hits - 1) / (2 * vrt_upper_drift0))) +
             gamma_lpdf(vrt_lower | ((misses - 1) / 2.0), ((misses - 1) / (2 * vrt_lower_drift0)));
    } else if (misses < 2) {
      return binomial_lpmf(hits | trials, zr) +
             normal_lpdf(mrt_upper | ndt + mdt_upper_drift0, sqrt(vrt_upper_drift0 / hits)) +
             gamma_lpdf(vrt_upper | ((hits - 1) / 2.0), ((hits - 1) / (2 * vrt_upper_drift0)));
    } else {
      return binomial_lpmf(hits | trials, zr) +
             normal_lpdf(mrt_lower | ndt + mdt_lower_drift0, sqrt(vrt_lower_drift0 / misses)) +
             gamma_lpdf(vrt_lower | ((misses - 1) / 2.0), ((misses - 1) / (2 * vrt_lower_drift0)));
    }
  }

  // Non-zero drift: use full formulas with cached hyperbolic functions
  // Use signed drift for pC calculation
  real k_z_signed = (drift * z) / s_sq;
  real k_x_signed = (drift * x0) / s_sq;

  // Compute pC with cached exponentials for numerical stability
  real exp_2kz = exp(2 * k_z_signed);
  real exp_neg2kz = 1.0 / exp_2kz;  // More stable than exp(-2*k_z)
  real exp_neg2kx = exp(-2 * k_x_signed);
  real pC = 1 - ((exp_neg2kx - exp_neg2kz) / (exp_2kz - exp_neg2kz));
  
  // Use soft absolute value: sqrt(drift^2 + tau^2) with larger tau
  // This provides smooth gradients without extreme curvature  
  real a = sqrt(drift * drift + 0.0001);  // tau = 0.01, tau^2 = 0.0001
  real k_z = (a * z) / s_sq;
  real k_x = (a * x0) / s_sq;

  real a_sq = square(a);

  // Cache scaling factors
  real scale_mdt = s_sq / a_sq;
  real scale_vrt = square(s_sq) / square(a_sq);

  // Cache hyperbolic function values (expensive operations)
  real two_kz = 2 * k_z;
  real coth_2kz = coth_stable(two_kz);
  real csch_2kz = csch_stable(two_kz);
  real csch_2kz_sq = square(csch_2kz);
  real k_z_sq = square(k_z);

  // Upper boundary: k_x + k_z
  real kxz_upper = k_x + k_z;
  real coth_kxz_upper = coth_stable(kxz_upper);
  real csch_kxz_upper = csch_stable(kxz_upper);
  real csch_kxz_upper_sq = square(csch_kxz_upper);
  real kxz_upper_sq = square(kxz_upper);

  // Lower boundary: -k_x + k_z
  real kxz_lower = -k_x + k_z;
  real coth_kxz_lower = coth_stable(kxz_lower);
  real csch_kxz_lower = csch_stable(kxz_lower);
  real csch_kxz_lower_sq = square(csch_kxz_lower);
  real kxz_lower_sq = square(kxz_lower);

  // Compute mean decision times
  real common_mdt = two_kz * coth_2kz;
  real mdt_upper_implied = scale_mdt * (common_mdt - kxz_upper * coth_kxz_upper);
  real mdt_lower_implied = scale_mdt * (common_mdt - kxz_lower * coth_kxz_lower);

  // Compute variance of decision times
  real common_vrt = 4 * k_z_sq * csch_2kz_sq + common_mdt;
  real vrt_upper_implied = scale_vrt * (common_vrt - kxz_upper_sq * csch_kxz_upper_sq - kxz_upper * coth_kxz_upper);
  real vrt_lower_implied = scale_vrt * (common_vrt - kxz_lower_sq * csch_kxz_lower_sq - kxz_lower * coth_kxz_lower);

  // return sum of sample statistics distributions log-likelihood
  if (misses >= 2 && hits >= 2) {
    return binomial_lpmf(hits | trials, pC) +
           normal_lpdf(mrt_upper | ndt + mdt_upper_implied, sqrt(vrt_upper_implied / hits)) +
           normal_lpdf(mrt_lower | ndt + mdt_lower_implied, sqrt(vrt_lower_implied / misses)) +
           gamma_lpdf(vrt_upper | ((hits - 1) / 2.0), ((hits - 1) / (2 * vrt_upper_implied))) +
           gamma_lpdf(vrt_lower | ((misses - 1) / 2.0), ((misses - 1) / (2 * vrt_lower_implied)));
  } else if (misses < 2) {
    return binomial_lpmf(hits | trials, pC) +
           normal_lpdf(mrt_upper | ndt + mdt_upper_implied, sqrt(vrt_upper_implied / hits)) +
           gamma_lpdf(vrt_upper | ((hits - 1) / 2.0), ((hits - 1) / (2 * vrt_upper_implied)));
  } else {
    return binomial_lpmf(hits | trials, pC) +
           normal_lpdf(mrt_lower | ndt + mdt_lower_implied, sqrt(vrt_lower_implied / misses)) +
           gamma_lpdf(vrt_lower | ((misses - 1) / 2.0), ((misses - 1) / (2 * vrt_lower_implied)));
  }
}