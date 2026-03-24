// Continuous Dual-Process (CDP) SDT likelihood functions
// Wixted & Mickes (2010, Psychological Review)
//
// Two independent continuous signal-detection dimensions:
//   Familiarity (F) and Recollection (R)
// Old/new decisions based on F + R (aggregated memory strength)
// Remember/Know(/Guess) decisions based on R (and optionally F) criteria
//
// Category probabilities computed via Gauss-Hermite quadrature over R
// (more robust than integrate_1d for near-zero probabilities)

// PDF dispatcher matching sdt_cumprob (CDF) from sdt_binary_funs.stan
// dist_type: 1=Gaussian, 2=Gumbel_min, 3=Gumbel_max, 4=Logistic
real sdt_pdf(real eta, int dist_type) {
  if (dist_type == 1) return exp(std_normal_lpdf(eta));
  if (dist_type == 2) {
    real neg_eta = -eta;
    return exp(neg_eta - exp(neg_eta));
  }
  if (dist_type == 3) {
    return exp(eta - exp(eta));
  }
  // logistic
  real e = exp(-eta);
  real d = 1 + e;
  return e / (d * d);
}

// Parsimonious thresholds (Selker et al., 2019)
vector cdp_thresholds_parsimonious(real criterion, real spacing,
                                   int K_full) {
  int n_thresh = K_full - 1;
  vector[n_thresh] thresholds;
  for (k in 1:n_thresh) {
    real gk = log(k * 1.0 / (K_full - k));
    thresholds[k] = criterion + exp(spacing) * gk;
  }
  return thresholds;
}

// Equidistant thresholds
vector cdp_thresholds_equidistant(real criterion, real spacing,
                                  int K_full) {
  int n_thresh = K_full - 1;
  real mid = (n_thresh + 1) / 2.0;
  vector[n_thresh] thresholds;
  for (k in 1:n_thresh) {
    thresholds[k] = criterion + (k - mid) * exp(spacing);
  }
  return thresholds;
}

// Evaluate the CDP integrand at a given R value
// Returns: [F(F_upper) - F(F_lower)] where bounds encode confidence
// band and kcrit constraints
real cdp_kernel(real R, real mu_F, real c_lo, real c_hi,
                real f_lo, real f_hi, int dist_type) {
  real F_lower = fmax(c_lo - R, f_lo) - mu_F;
  real F_upper = fmin(c_hi - R, f_hi) - mu_F;
  if (F_upper <= F_lower) return 0.0;
  return sdt_cumprob(F_upper, dist_type)
         - sdt_cumprob(F_lower, dist_type);
}

// Compute CDP category probability via Gauss-Hermite quadrature
// For Gaussian R: exact substitution R = mu_R + sd_R*sqrt(2)*t
//   gives int f(R)*phi_R(R)dR = (1/sqrt(pi)) * sum w_j * f(node_j)
// For non-Gaussian R: use importance sampling correction
//   weight_j *= pdf_R(node_j) / phi(t_j) * sd_R * sqrt(2)
real cdp_prob_full(real mu_F, real mu_R, real sd_R,
                   real c_lo, real c_hi,
                   real f_lo, real f_hi,
                   real R_lo, real R_hi,
                   int dist_type) {
  // 30-point Gauss-Hermite nodes and weights
  int N_GH = 30;
  array[30] real nodes = {
    -6.86334529, -6.13827922, -5.53314715, -4.98891897,
    -4.48305536, -4.00390860, -3.54444387, -3.09997053,
    -2.66713212, -2.24339147, -1.82674114, -1.41552780,
    -1.00833827, -0.60392106, -0.20112858,
     0.20112858,  0.60392106,  1.00833827,  1.41552780,
     1.82674114,  2.24339147,  2.66713212,  3.09997053,
     3.54444387,  4.00390860,  4.48305536,  4.98891897,
     5.53314715,  6.13827922,  6.86334529
  };
  array[30] real log_weights = {
    -39.07158, -29.85787, -23.51483, -18.73484,
    -14.94853, -11.85985,  -9.27833,  -7.07908,
     -5.17427,  -3.49643,  -2.00076,  -0.65389,
      0.56834,   1.68001,   2.69010,
      2.69010,   1.68001,   0.56834,  -0.65389,
     -2.00076,  -3.49643,  -5.17427,  -7.07908,
     -9.27833, -11.85985, -14.94853, -18.73484,
    -23.51483, -29.85787, -39.07158
  };
  real sqrt2 = sqrt(2.0);
  real log_inv_sqrt_pi = -0.5 * log(pi());

  real result = 0.0;
  for (j in 1:N_GH) {
    real R_j = mu_R + sd_R * sqrt2 * nodes[j];

    // Skip if R_j is outside the R integration bounds
    if (R_j < R_lo || R_j > R_hi) continue;

    real kernel = cdp_kernel(R_j, mu_F, c_lo, c_hi,
                             f_lo, f_hi, dist_type);
    if (kernel <= 0) continue;

    if (dist_type == 1) {
      // Gaussian R: GH quadrature is exact
      result += exp(log_weights[j] + log_inv_sqrt_pi)
                * kernel;
    } else {
      // Non-Gaussian R: importance sampling correction
      // weight * pdf_R(R_j) / phi_GH(t_j) * sd_R * sqrt(2)
      real t_j = nodes[j];
      real phi_gh = exp(-t_j * t_j);  // unnormalized GH kernel
      real eta_R = (R_j - mu_R) / sd_R;
      real pdf_R = sdt_pdf(eta_R, dist_type) / sd_R;
      // GH: (1/sqrt(pi)) * w_j * f(t_j) where w_j includes
      // exp(t^2). We need f(R)*pdf_R(R)*dR, so:
      // contribution = w_j/sqrt(pi) * kernel * pdf_R *
      //               sd_R*sqrt(2) / exp(-t^2)
      // = w_j/sqrt(pi) * kernel * pdf_R * sd_R*sqrt(2)*exp(t^2)
      // But log_weights already has log(w_j) where w_j absorbs
      // the exp(t^2), so:
      // = exp(log_w + log_inv_sqrt_pi) * kernel
      //   * pdf_R * sd_R * sqrt2 / phi_GH
      // Since phi_GH = exp(-t^2) and the original GH weight
      // w_j = W_j * exp(t_j^2), log_weights[j] = log(W_j*exp(t^2)):
      result += exp(log_weights[j] + log_inv_sqrt_pi)
                * kernel * pdf_R * sd_R * sqrt2 / phi_gh;
    }
  }

  // Adjust for truncated R range (renormalize if needed)
  // For Gaussian: the GH nodes outside [R_lo, R_hi] are skipped,
  // which naturally handles the truncation
  return fmax(result, 1e-20);
}

// CDP UV likelihood (2-way R/K model, sigmar as dpar)
real sdt_cdp_uv_lpmf(int y, real mu, real dprimef, real dprimer,
                     real criterion, real spacing, real rcrit,
                     real sigmar,
                     int cat_type, int cat_conf, int stim,
                     int n_conf, int has_guess, int dist_type,
                     int thresh_type,
                     data array[] real x_r,
                     data array[] int x_i) {
  if (y == 0) return 0;

  int K_full = 2 * n_conf;
  vector[K_full - 1] thresholds;
  if (thresh_type == 1)
    thresholds = cdp_thresholds_parsimonious(criterion, spacing,
                                             K_full);
  else
    thresholds = cdp_thresholds_equidistant(criterion, spacing,
                                            K_full);

  int global_k;
  if (cat_type == 1)
    global_k = cat_conf;
  else
    global_k = n_conf + cat_conf;

  real c_lo = global_k == 1
              ? negative_infinity() : thresholds[global_k - 1];
  real c_hi = global_k == K_full
              ? positive_infinity() : thresholds[global_k];

  real R_lo;
  real R_hi;
  real f_lo = negative_infinity();
  real f_hi = positive_infinity();

  if (cat_type == 1) {
    // new: all R, all F
    R_lo = negative_infinity();
    R_hi = positive_infinity();
  } else if (cat_type == 4) {
    // remember: R > rcrit
    R_lo = rcrit;
    R_hi = positive_infinity();
  } else {
    // know (cat_type == 3): R < rcrit
    R_lo = negative_infinity();
    R_hi = rcrit;
  }

  real mu_F = stim == 1 ? dprimef : 0.0;
  real mu_R = stim == 1 ? dprimer : 0.0;
  real sd_R = stim == 1 ? exp(sigmar) : 1.0;

  real p = cdp_prob_full(mu_F, mu_R, sd_R, c_lo, c_hi,
                         f_lo, f_hi, R_lo, R_hi, dist_type);

  return y * log(p);
}

// CDP UV likelihood with kcrit as additional dpar (3-way R/K/G)
real sdt_cdp_uv3_lpmf(int y, real mu, real dprimef, real dprimer,
                      real criterion, real spacing, real rcrit,
                      real sigmar, real kcrit,
                      int cat_type, int cat_conf, int stim,
                      int n_conf, int has_guess, int dist_type,
                      int thresh_type,
                      data array[] real x_r,
                      data array[] int x_i) {
  if (y == 0) return 0;

  int K_full = 2 * n_conf;
  vector[K_full - 1] thresholds;
  if (thresh_type == 1)
    thresholds = cdp_thresholds_parsimonious(criterion, spacing,
                                             K_full);
  else
    thresholds = cdp_thresholds_equidistant(criterion, spacing,
                                            K_full);

  int global_k;
  if (cat_type == 1)
    global_k = cat_conf;
  else
    global_k = n_conf + cat_conf;

  real c_lo = global_k == 1
              ? negative_infinity() : thresholds[global_k - 1];
  real c_hi = global_k == K_full
              ? positive_infinity() : thresholds[global_k];

  real R_lo;
  real R_hi;
  real f_lo;
  real f_hi;

  if (cat_type == 1) {
    R_lo = negative_infinity(); R_hi = positive_infinity();
    f_lo = negative_infinity(); f_hi = positive_infinity();
  } else if (cat_type == 4) {
    R_lo = rcrit; R_hi = positive_infinity();
    f_lo = negative_infinity(); f_hi = positive_infinity();
  } else if (cat_type == 3) {
    R_lo = negative_infinity(); R_hi = rcrit;
    f_lo = kcrit; f_hi = positive_infinity();
  } else {
    R_lo = negative_infinity(); R_hi = rcrit;
    f_lo = negative_infinity(); f_hi = kcrit;
  }

  real mu_F = stim == 1 ? dprimef : 0.0;
  real mu_R = stim == 1 ? dprimer : 0.0;
  real sd_R = stim == 1 ? exp(sigmar) : 1.0;

  real p = cdp_prob_full(mu_F, mu_R, sd_R, c_lo, c_hi,
                         f_lo, f_hi, R_lo, R_hi, dist_type);

  return y * log(p);
}
