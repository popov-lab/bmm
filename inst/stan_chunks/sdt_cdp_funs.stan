// Continuous Dual-Process (CDP) SDT likelihood (Wixted & Mickes, 2010).
//
// Native-multinomial formulation: every response category's logit is set to
// log(p_cat) so brms' softmax recovers the CDP category probabilities exactly
// (the category probabilities sum to 1 analytically, and softmax is invariant
// to the shared normalising constant).
//
// Two correlated continuous dimensions per item, Familiarity F and Recollection
// R, with corr(F, R) = tanh(rho):
//   target: F ~ N(dprimef, 1),  R ~ N(dprimer, exp(sigmar))
//   lure:   F ~ N(0, 1),        R ~ N(0, 1)
// Old/new confidence is read off the aggregate strength S = F + R; the
// Remember/Know split is read off R against rcrit; the optional Know/Guess split
// is read off F against kcrit.
//
// Normal noise only: the bivariate-normal CDF needed for the Remember/Know split
// is exact and differentiable via Owen's T (owens_t is a Stan primitive).

// Bivariate standard-normal CDF P(Z1 <= z1, Z2 <= z2) with correlation rho,
// via Owen's T (Owen, 1956). z1, z2 are finite (infinite strength bounds are
// handled by the caller); exact zeros are nudged to avoid the z in the
// denominator of the Owen's-T argument.
real cdp_Phi2(real z1, real z2, real rho) {
  real h = z1;
  real k = z2;
  if (h == 0 && k == 0) return 0.25 + asin(rho) / (2 * pi());
  if (abs(h) < 1e-10) h = h >= 0 ? 1e-10 : -1e-10;
  if (abs(k) < 1e-10) k = k >= 0 ? 1e-10 : -1e-10;
  real denom = sqrt((1 + rho) * (1 - rho));
  real a1 = (k / h - rho) / denom;
  real a2 = (h / k - rho) / denom;
  real prod = h * k;
  real delta = (prod < 0 || (prod == 0 && (h + k) < 0)) ? 1 : 0;
  return 0.5 * (Phi(h) + Phi(k) - delta) - owens_t(h, a1) - owens_t(k, a2);
}

// Mass of the guess region for one strength bin: P(R < rcrit, F < kcrit,
// c_lo < F + R < c_hi). No closed form (a strength band clips the F < kcrit,
// R < rcrit quadrant on the diagonal), so integrate over R with F | R normal.
// 20-point Gauss-Legendre on (-inf, rcrit] via the t/(1-t) map.
real cdp_guess_mass(real mu_F, real mu_R, real sd_R, real corr,
                    real c_lo, real c_hi, real rcrit, real kcrit) {
  int N_GL = 20;
  vector[N_GL] gl_nodes = to_vector({
    -9.9312859918509492e-01, -9.6397192727791379e-01,
    -9.1223442825132591e-01, -8.3911697182221882e-01,
    -7.4633190646015087e-01, -6.3605368072651512e-01,
    -5.1086700195082709e-01, -3.7370608871541956e-01,
    -2.2778585114164508e-01, -7.6526521133497324e-02,
     7.6526521133497338e-02,  2.2778585114164508e-01,
     3.7370608871541956e-01,  5.1086700195082709e-01,
     6.3605368072651512e-01,  7.4633190646015087e-01,
     8.3911697182221882e-01,  9.1223442825132591e-01,
     9.6397192727791379e-01,  9.9312859918509492e-01});
  vector[N_GL] gl_weights = to_vector({
    1.7614007139152118e-02, 4.0601429800386941e-02,
    6.2672048334109064e-02, 8.3276741576704749e-02,
    1.0193011981724044e-01, 1.1819453196151842e-01,
    1.3168863844917664e-01, 1.4209610931838205e-01,
    1.4917298647260360e-01, 1.5275338713072585e-01,
    1.5275338713072585e-01, 1.4917298647260360e-01,
    1.4209610931838205e-01, 1.3168863844917664e-01,
    1.1819453196151842e-01, 1.0193011981724044e-01,
    8.3276741576704749e-02, 6.2672048334109064e-02,
    4.0601429800386941e-02, 1.7614007139152118e-02});

  real cond_sd = sqrt(fmax(1 - corr * corr, 1e-12));
  real result = 0;
  for (i in 1:N_GL) {
    real t = 0.5 * (gl_nodes[i] + 1);
    real y = t / (1 - t);
    real R = rcrit - y;                       // maps t in (0,1) to R in (-inf, rcrit]
    real jac = 1 / square(1 - t);
    real cond_mean = mu_F + corr * (R - mu_R) / sd_R;   // E[F | R], sd_F = 1
    real f_lo = c_lo == negative_infinity() ? negative_infinity() : c_lo - R;
    real f_hi = fmin(c_hi == positive_infinity() ? positive_infinity() : c_hi - R,
                     kcrit);
    real mass_f = 0;
    if (f_hi > f_lo) {
      real u_hi = f_hi == positive_infinity() ? 1.0 : Phi((f_hi - cond_mean) / cond_sd);
      real u_lo = f_lo == negative_infinity() ? 0.0 : Phi((f_lo - cond_mean) / cond_sd);
      mass_f = u_hi - u_lo;
    }
    real pdf_R = exp(normal_lpdf(R | mu_R, sd_R));
    result += gl_weights[i] * jac * pdf_R * mass_f;
  }
  return fmax(0.5 * result, 1e-20);
}

// Confidence thresholds on the strength axis S = F + R, anchored so the old/new
// boundary (between bin n_new and bin n_new + 1) sits at `criterion`, with
// n_new - 1 thresholds below and n_old - 1 above. Reduces to the symmetric
// centred construction when n_new == n_old. thresh_type: 1 = parsimonious
// (Selker et al., 2019), 2 = equidistant, 3 = log_distance (free cumulative
// log-distances; `deltas` has length K_full - 2, anchor index n_new carries none).
vector cdp_make_thresholds(real criterion, real spacing, array[] real deltas,
                           int n_new, int n_old, int thresh_type) {
  int K_full = n_new + n_old;
  int n_thresh = K_full - 1;
  vector[n_thresh] thr;
  real s = exp(spacing);
  if (thresh_type == 2) {
    for (k in 1:n_thresh) thr[k] = criterion + (k - n_new) * s;
  } else if (thresh_type == 3) {
    thr[n_new] = criterion;
    for (k in (n_new + 1):n_thresh)
      thr[k] = thr[k - 1] + exp(deltas[k - 1]);
    for (k in 1:(n_new - 1)) {
      int kk = n_new - k;            // Stan counts up: descend n_new-1 .. 1
      thr[kk] = thr[kk + 1] - exp(deltas[kk]);
    }
  } else {
    real anchor = log(n_new * 1.0 / (K_full - n_new));
    for (k in 1:n_thresh)
      thr[k] = criterion + s * (log(k * 1.0 / (K_full - k)) - anchor);
  }
  return thr;
}

// CDP probability for a single response category, ordered
//   new(1..n_new), [guess(1..n_old)], know(1..n_old), remember(1..n_old).
real cdp_category_prob(int cat, vector thresholds,
                       real dprimef, real dprimer, real sigmar, real rho,
                       real rcrit, real kcrit, real stimulus,
                       int n_new, int n_old, int has_guess) {
  int K_full = n_new + n_old;
  int type;
  int conf;
  if (cat <= n_new) {
    type = 1;
    conf = cat;
  } else if (has_guess == 1) {
    int r = cat - n_new;
    int block = (r - 1) %/% n_old;          // 0 guess, 1 know, 2 remember
    type = block == 0 ? 2 : (block == 1 ? 3 : 4);
    conf = r - block * n_old;
  } else {
    int r = cat - n_new;
    int block = (r - 1) %/% n_old;          // 0 know, 1 remember
    type = block == 0 ? 3 : 4;
    conf = r - block * n_old;
  }
  int global_k = type == 1 ? conf : (n_new + conf);
  real c_lo = global_k == 1 ? negative_infinity() : thresholds[global_k - 1];
  real c_hi = global_k == K_full ? positive_infinity() : thresholds[global_k];

  real mu_F = stimulus > 0.5 ? dprimef : 0.0;
  real mu_R = stimulus > 0.5 ? dprimer : 0.0;
  real sd_R = stimulus > 0.5 ? exp(sigmar) : 1.0;
  real corr = tanh(rho);
  real mu_S = mu_F + mu_R;
  real sigma_S = sqrt(square(sd_R + corr) + (1 - square(corr)));
  real rho_RS = (sd_R + corr) / sigma_S;
  real hcrit = (rcrit - mu_R) / sd_R;

  real z_lo = c_lo == negative_infinity() ? negative_infinity() : (c_lo - mu_S) / sigma_S;
  real z_hi = c_hi == positive_infinity() ? positive_infinity() : (c_hi - mu_S) / sigma_S;
  real p_bin = (z_hi == positive_infinity() ? 1.0 : Phi(z_hi))
             - (z_lo == negative_infinity() ? 0.0 : Phi(z_lo));

  if (type == 1) return fmax(p_bin, 1e-20);

  real k_hi = z_hi == positive_infinity() ? Phi(hcrit) : cdp_Phi2(hcrit, z_hi, rho_RS);
  real k_lo = z_lo == negative_infinity() ? 0.0 : cdp_Phi2(hcrit, z_lo, rho_RS);
  real p_know_total = k_hi - k_lo;

  if (type == 4) return fmax(p_bin - p_know_total, 1e-20);
  if (has_guess == 0) return fmax(p_know_total, 1e-20);

  real p_guess = cdp_guess_mass(mu_F, mu_R, sd_R, corr, c_lo, c_hi, rcrit, kcrit);
  if (type == 2) return fmax(p_guess, 1e-20);
  return fmax(p_know_total - p_guess, 1e-20);
}

// The category logit `sdt_cdp_logmu` is code-generated per model by
// .sdt_cdp_logmu_stan() (R/model_sdt_cdp.R) so the per-distance log_distance
// deltas can be passed with a fixed arity, mirroring sdt_rating. It calls the
// helpers above (cdp_make_thresholds + cdp_category_prob).
