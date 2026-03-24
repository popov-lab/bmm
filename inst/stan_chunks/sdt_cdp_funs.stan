// Continuous Dual-Process (CDP) SDT likelihood functions
// Wixted & Mickes (2010, Psychological Review)
//
// Two independent continuous signal-detection dimensions:
//   Familiarity (F) and Recollection (R)
// Old/new decisions based on F + R (aggregated memory strength)
// Remember/Know(/Guess) decisions based on R (and optionally F) criteria
//
// Category probabilities computed via integrate_1d over R dimension

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

// Integrand for CDP category probabilities
// Integrates over recollection dimension R:
//   phi_R(R) * [F(F_upper) - F(F_lower)]
//
// theta = {mu_F, mu_R, sd_R, c_lo, c_hi, f_lo, f_hi}
// x_i = {dist_type}
real cdp_integrand(real R, real xc,
                   array[] real theta,
                   array[] real x_r, array[] int x_i) {
  real mu_F = theta[1];
  real mu_R = theta[2];
  real sd_R = theta[3];
  real c_lo = theta[4];
  real c_hi = theta[5];
  real f_lo = theta[6];
  real f_hi = theta[7];
  int dist_type = x_i[1];

  real eta_R = (R - mu_R) / sd_R;
  real phi_R = sdt_pdf(eta_R, dist_type) / sd_R;

  real F_lower = fmax(c_lo - R, f_lo) - mu_F;
  real F_upper = fmin(c_hi - R, f_hi) - mu_F;

  if (F_upper <= F_lower) return 0.0;

  return phi_R * (sdt_cumprob(F_upper, dist_type)
                  - sdt_cumprob(F_lower, dist_type));
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

// Compute CDP category probability via integrate_1d
// Uses cdp_dist_type (data int array from transformed data) for the
// distribution type, and a relaxed tolerance (1e-3) to handle
// near-zero probability categories.
real cdp_prob(real mu_F, real mu_R, real sd_R,
              real c_lo, real c_hi,
              real f_lo, real f_hi,
              real R_lo, real R_hi,
              data array[] real x_r,
              data array[] int cdp_dist_type) {
  real p = integrate_1d(cdp_integrand, R_lo, R_hi,
                        {mu_F, mu_R, sd_R, c_lo, c_hi,
                         f_lo, f_hi},
                        x_r, cdp_dist_type, 1e-3);
  return fmax(p, 1e-20);
}

// CDP UV likelihood (2-way R/K model, sigmar as dpar)
// vint layout: cat_type, cat_conf, stim, n_new, n_old,
//              has_guess, dist_type, thresh_type
real sdt_cdp_uv_lpmf(int y, real mu, real dprimef, real dprimer,
                     real criterion, real spacing, real rcrit,
                     real sigmar,
                     int cat_type, int cat_conf, int stim,
                     int n_new, int n_old, int has_guess,
                     int dist_type, int thresh_type,
                     data array[] real x_r,
                     data array[] int cdp_dist_type) {
  if (y == 0) return 0;

  int K_full = n_new + n_old;
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
    global_k = n_new + cat_conf;

  real c_lo = global_k == 1
              ? negative_infinity() : thresholds[global_k - 1];
  real c_hi = global_k == K_full
              ? positive_infinity() : thresholds[global_k];

  real R_lo;
  real R_hi;
  real f_lo = negative_infinity();
  real f_hi = positive_infinity();

  if (cat_type == 1) {
    R_lo = negative_infinity();
    R_hi = positive_infinity();
  } else if (cat_type == 4) {
    R_lo = rcrit;
    R_hi = positive_infinity();
  } else {
    R_lo = negative_infinity();
    R_hi = rcrit;
  }

  real mu_F = stim == 1 ? dprimef : 0.0;
  real mu_R = stim == 1 ? dprimer : 0.0;
  real sd_R = stim == 1 ? exp(sigmar) : 1.0;

  real p = cdp_prob(mu_F, mu_R, sd_R, c_lo, c_hi,
                    f_lo, f_hi, R_lo, R_hi,
                    x_r, cdp_dist_type);

  return y * log(p);
}

// CDP UV likelihood with kcrit as additional dpar (3-way R/K/G)
real sdt_cdp_uv3_lpmf(int y, real mu, real dprimef, real dprimer,
                      real criterion, real spacing, real rcrit,
                      real sigmar, real kcrit,
                      int cat_type, int cat_conf, int stim,
                      int n_new, int n_old, int has_guess,
                      int dist_type, int thresh_type,
                      data array[] real x_r,
                      data array[] int cdp_dist_type) {
  if (y == 0) return 0;

  int K_full = n_new + n_old;
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
    global_k = n_new + cat_conf;

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

  real p = cdp_prob(mu_F, mu_R, sd_R, c_lo, c_hi,
                    f_lo, f_hi, R_lo, R_hi,
                    x_r, cdp_dist_type);

  return y * log(p);
}
