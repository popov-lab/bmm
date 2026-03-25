// Continuous Dual-Process (CDP) SDT likelihood functions
// Wixted & Mickes (2010, Psychological Review)
//
// Two independent continuous signal-detection dimensions:
//   Familiarity (F) and Recollection (R)
// Old/new decisions based on F + R (aggregated memory strength)
// Remember/Know(/Guess) decisions based on R (and optionally F) criteria
//
// For Gaussian: analytical solution via bivariate normal CDF (exact)
// For non-Gaussian: Gauss-Hermite quadrature with sigmoid gating (approximate)

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

// Bivariate normal CDF: P(X <= a, Y <= b) with correlation rho
// Uses Genz (2004) integral representation with 10-point Gauss-Legendre:
//   Phi2(a,b,rho) = Phi(a)*Phi(b)
//     + 1/(2*pi) * integral_0^asin(rho)
//       exp(-(a^2+b^2-2*a*b*sin(t)) / (2*cos(t)^2)) dt
real bvn_Phi2(real a, real b, real rho) {
  // handle infinite bounds: Phi2(-Inf, b) = 0, Phi2(Inf, b) = Phi(b)
  if (is_inf(a) && a < 0) return 0;
  if (is_inf(b) && b < 0) return 0;
  if (is_inf(a) && a > 0) return Phi(b);
  if (is_inf(b) && b > 0) return Phi(a);
  if (abs(rho) < 1e-12) return Phi(a) * Phi(b);

  int N_GL = 10;
  vector[N_GL] gl_nodes = to_vector({
    -9.7390652851717163e-01, -8.6506336668898431e-01,
    -6.7940956829902421e-01, -4.3339539412924699e-01,
    -1.4887433898163108e-01,  1.4887433898163141e-01,
     4.3339539412924721e-01,  6.7940956829902410e-01,
     8.6506336668898431e-01,  9.7390652851717141e-01
  });
  vector[N_GL] gl_weights = to_vector({
    6.6671344308688416e-02, 1.4945134915058092e-01,
    2.1908636251598199e-01, 2.6926671930999607e-01,
    2.9552422471475326e-01, 2.9552422471475270e-01,
    2.6926671930999513e-01, 2.1908636251598207e-01,
    1.4945134915058053e-01, 6.6671344308687819e-02
  });

  real asr = asin(rho);
  real hs = -(a * a + b * b) / 2;
  real ab2 = a * b;
  real result = 0;
  for (i in 1:N_GL) {
    real t = asr * (gl_nodes[i] + 1) / 2;
    real sn = sin(t);
    real cs2 = cos(t) * cos(t);
    result += gl_weights[i] * exp((ab2 * sn + hs) / (1 - sn * sn));
  }
  return Phi(a) * Phi(b) + result * asr / (4 * pi());
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

// CDP category probability for Gaussian (analytical via bivariate normal)
//
// M = F + R ~ N(mu_M, sigma_M^2) where mu_M = mu_F + mu_R, sigma_M = sqrt(1 + sd_R^2)
// Corr(M, R) = sd_R / sigma_M
//
// P(new, bin k) = Phi(a_hi) - Phi(a_lo)
// P(know, bin k) = Phi2(a_hi, b; rho) - Phi2(a_lo, b; rho)
// P(remember, bin k) = P(bin k) - P(know, bin k)
//
// where a = (c - mu_M) / sigma_M, b = (rcrit - mu_R) / sd_R
real cdp_prob_gaussian(real mu_F, real mu_R, real sd_R,
                       real c_lo, real c_hi,
                       real rcrit,
                       int cat_type) {
  real mu_M = mu_F + mu_R;
  real sigma_M = sqrt(1 + sd_R * sd_R);
  real rho = sd_R / sigma_M;
  real b = (rcrit - mu_R) / sd_R;

  // compute Phi(a) and Phi2(a, b, rho) for each bound, avoiding Inf arithmetic
  int lo_inf = is_inf(c_lo);
  int hi_inf = is_inf(c_hi);
  real phi_lo = lo_inf ? 0.0 : Phi((c_lo - mu_M) / sigma_M);
  real phi_hi = hi_inf ? 1.0 : Phi((c_hi - mu_M) / sigma_M);
  real p_bin = phi_hi - phi_lo;

  if (cat_type == 1) return fmax(p_bin, 1e-20);

  real bvn_lo = lo_inf ? 0.0 : bvn_Phi2((c_lo - mu_M) / sigma_M, b, rho);
  real bvn_hi = hi_inf ? Phi(b) : bvn_Phi2((c_hi - mu_M) / sigma_M, b, rho);
  real p_know = bvn_hi - bvn_lo;

  if (cat_type == 4)
    return fmax(p_bin - p_know, 1e-20);
  return fmax(p_know, 1e-20);
}

// CDP category probability for Gaussian 3-way (R/K/G) model
// Guess: R < rcrit AND F < kcrit
// Know:  R < rcrit AND F >= kcrit
// Remember: R >= rcrit
//
// For guess/know, the F constraint introduces a third dimension.
// We decompose using conditional probability:
//   P(know, bin k) = P(M in bin k, R < rcrit) - P(M in bin k, R < rcrit, F < kcrit)
//   P(guess, bin k) = P(M in bin k, R < rcrit, F < kcrit)
//
// P(M in bin, R < rcrit, F < kcrit) requires a trivariate constraint:
//   c_lo < F + R < c_hi, R < rcrit, F < kcrit
// We approximate this with GH quadrature over R (only for know/guess split)
// using a sigmoid gate at rcrit for smoothness.
real cdp_prob_gaussian_3way(real mu_F, real mu_R, real sd_R,
                            real c_lo, real c_hi,
                            real rcrit, real kcrit,
                            int cat_type) {
  real mu_M = mu_F + mu_R;
  real sigma_M = sqrt(1 + sd_R * sd_R);
  real rho = sd_R / sigma_M;
  real b = (rcrit - mu_R) / sd_R;

  int lo_inf = is_inf(c_lo);
  int hi_inf = is_inf(c_hi);
  real phi_lo = lo_inf ? 0.0 : Phi((c_lo - mu_M) / sigma_M);
  real phi_hi = hi_inf ? 1.0 : Phi((c_hi - mu_M) / sigma_M);
  real p_bin = phi_hi - phi_lo;

  if (cat_type == 1) return fmax(p_bin, 1e-20);

  real bvn_lo = lo_inf ? 0.0 : bvn_Phi2((c_lo - mu_M) / sigma_M, b, rho);
  real bvn_hi = hi_inf ? Phi(b) : bvn_Phi2((c_hi - mu_M) / sigma_M, b, rho);
  real p_know_total = bvn_hi - bvn_lo;

  // remember: R >= rcrit (no F constraint)
  if (cat_type == 4) return fmax(p_bin - p_know_total, 1e-20);

  // For guess/know split: need P(M in bin, R < rcrit, F < kcrit)
  // Use 20-point GH quadrature over R with sigmoid gate at rcrit
  int N_GH = 20;
  vector[N_GH] gh_nodes = to_vector({
    -7.6190485416797546e+00, -6.5105901570136488e+00,
    -5.5787388058932059e+00, -4.7345813340460463e+00,
    -3.9439673506573110e+00, -3.1890148165533843e+00,
    -2.4586636111723603e+00, -1.7452473208141255e+00,
    -1.0429453488027509e+00, -3.4696415708135458e-01,
     3.4696415708135830e-01,  1.0429453488027574e+00,
     1.7452473208141317e+00,  2.4586636111723683e+00,
     3.1890148165533900e+00,  3.9439673506573163e+00,
     4.7345813340460552e+00,  5.5787388058932033e+00,
     6.5105901570136551e+00,  7.6190485416797591e+00
  });
  vector[N_GH] gh_weights = to_vector({
    1.2578006724378954e-13, 2.4820623623151972e-10,
    6.1274902599825256e-08, 4.4021210902309806e-06,
    1.2882627996193093e-04, 1.8301031310804826e-03,
    1.3997837447100857e-02, 6.1506372063977507e-02,
    1.6173933398399959e-01, 2.6079306344955683e-01,
    2.6079306344955305e-01, 1.6173933398399776e-01,
    6.1506372063977438e-02, 1.3997837447101162e-02,
    1.8301031310805052e-03, 1.2882627996193072e-04,
    4.4021210902309052e-06, 6.1274902599829068e-08,
    2.4820623623151936e-10, 1.2578006724379269e-13
  });

  // P(M in bin, R < rcrit, F < kcrit) via GH over R with sigmoid gate
  real p_guess = 0;
  for (i in 1:N_GH) {
    real z = gh_nodes[i];
    real gate = inv_logit(200 * (b - z));  // R < rcrit gate
    if (gate < 1e-10) continue;

    real R = mu_R + sd_R * z;
    real F_upper = fmin(c_hi - R, kcrit) - mu_F;
    real F_lower = c_lo - R - mu_F;
    if (F_upper <= F_lower) continue;

    p_guess += gh_weights[i] * gate * (Phi(F_upper) - Phi(F_lower));
  }
  p_guess = fmax(p_guess, 1e-20);

  if (cat_type == 2) return p_guess;
  // know = total_know - guess
  return fmax(p_know_total - p_guess, 1e-20);
}

// Non-Gaussian fallback: GH quadrature with sigmoid gating
real cdp_prob_nongaussian(real mu_F, real mu_R, real sd_R,
                          real c_lo, real c_hi,
                          real f_lo, real f_hi,
                          real z_lo, real z_hi,
                          int dist_type) {
  int N_QUAD = 20;
  vector[N_QUAD] nodes = to_vector({
    -7.6190485416797546e+00, -6.5105901570136488e+00,
    -5.5787388058932059e+00, -4.7345813340460463e+00,
    -3.9439673506573110e+00, -3.1890148165533843e+00,
    -2.4586636111723603e+00, -1.7452473208141255e+00,
    -1.0429453488027509e+00, -3.4696415708135458e-01,
     3.4696415708135830e-01,  1.0429453488027574e+00,
     1.7452473208141317e+00,  2.4586636111723683e+00,
     3.1890148165533900e+00,  3.9439673506573163e+00,
     4.7345813340460552e+00,  5.5787388058932033e+00,
     6.5105901570136551e+00,  7.6190485416797591e+00
  });
  vector[N_QUAD] weights = to_vector({
    1.2578006724378954e-13, 2.4820623623151972e-10,
    6.1274902599825256e-08, 4.4021210902309806e-06,
    1.2882627996193093e-04, 1.8301031310804826e-03,
    1.3997837447100857e-02, 6.1506372063977507e-02,
    1.6173933398399959e-01, 2.6079306344955683e-01,
    2.6079306344955305e-01, 1.6173933398399776e-01,
    6.1506372063977438e-02, 1.3997837447101162e-02,
    1.8301031310805052e-03, 1.2882627996193072e-04,
    4.4021210902309052e-06, 6.1274902599829068e-08,
    2.4820623623151936e-10, 1.2578006724379269e-13
  });

  real result = 0;
  for (i in 1:N_QUAD) {
    real z = nodes[i];
    real gate = 1.0;
    if (z_lo != negative_infinity())
      gate *= inv_logit(100 * (z - z_lo));
    if (z_hi != positive_infinity())
      gate *= inv_logit(100 * (z_hi - z));
    if (gate < 1e-10) continue;

    real R = mu_R + sd_R * z;
    real F_lower = fmax(c_lo - R, f_lo) - mu_F;
    real F_upper = fmin(c_hi - R, f_hi) - mu_F;
    if (F_upper <= F_lower) continue;

    real w = weights[i] * gate;
    w *= sdt_pdf(z, dist_type) / sdt_pdf(z, 1);

    result += w * (sdt_cumprob(F_upper, dist_type)
                   - sdt_cumprob(F_lower, dist_type));
  }
  return fmax(result, 1e-20);
}

// CDP UV likelihood (2-way R/K model, sigmar as dpar)
// vint layout: cat_type, cat_conf, stim, n_new, n_old,
//              has_guess, dist_type, thresh_type
real sdt_cdp_uv_lpmf(int y, real mu, real dprimef, real dprimer,
                     real criterion, real spacing, real rcrit,
                     real sigmar,
                     int cat_type, int cat_conf, int stim,
                     int n_new, int n_old, int has_guess,
                     int dist_type, int thresh_type) {
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

  real mu_F = stim == 1 ? dprimef : 0.0;
  real mu_R = stim == 1 ? dprimer : 0.0;
  real sd_R = stim == 1 ? exp(sigmar) : 1.0;

  real p;
  if (dist_type == 1)
    p = cdp_prob_gaussian(mu_F, mu_R, sd_R, c_lo, c_hi, rcrit, cat_type);
  else {
    real z_lo = negative_infinity();
    real z_hi = positive_infinity();
    if (cat_type == 4) z_lo = (rcrit - mu_R) / sd_R;
    else if (cat_type != 1) z_hi = (rcrit - mu_R) / sd_R;
    p = cdp_prob_nongaussian(mu_F, mu_R, sd_R, c_lo, c_hi,
                             negative_infinity(), positive_infinity(),
                             z_lo, z_hi, dist_type);
  }

  return y * log(p);
}

// CDP UV likelihood with kcrit as additional dpar (3-way R/K/G)
real sdt_cdp_uv3_lpmf(int y, real mu, real dprimef, real dprimer,
                      real criterion, real spacing, real rcrit,
                      real sigmar, real kcrit,
                      int cat_type, int cat_conf, int stim,
                      int n_new, int n_old, int has_guess,
                      int dist_type, int thresh_type) {
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

  real mu_F = stim == 1 ? dprimef : 0.0;
  real mu_R = stim == 1 ? dprimer : 0.0;
  real sd_R = stim == 1 ? exp(sigmar) : 1.0;

  real p;
  if (dist_type == 1) {
    p = cdp_prob_gaussian_3way(mu_F, mu_R, sd_R, c_lo, c_hi,
                               rcrit, kcrit, cat_type);
  } else {
    real f_lo; real f_hi;
    if (cat_type == 1) { f_lo = negative_infinity(); f_hi = positive_infinity(); }
    else if (cat_type == 4) { f_lo = negative_infinity(); f_hi = positive_infinity(); }
    else if (cat_type == 3) { f_lo = kcrit; f_hi = positive_infinity(); }
    else { f_lo = negative_infinity(); f_hi = kcrit; }

    real z_lo = negative_infinity();
    real z_hi = positive_infinity();
    if (cat_type == 4) z_lo = (rcrit - mu_R) / sd_R;
    else if (cat_type != 1) z_hi = (rcrit - mu_R) / sd_R;

    p = cdp_prob_nongaussian(mu_F, mu_R, sd_R, c_lo, c_hi,
                             f_lo, f_hi, z_lo, z_hi, dist_type);
  }

  return y * log(p);
}
