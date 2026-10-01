// Continuous Dual-Process (CDP) SDT likelihood (Wixted & Mickes, 2010).
//
// Native-multinomial formulation: every response category's logit is set to
// log(p_cat) so brms' softmax recovers the CDP category probabilities exactly
// (the category probabilities sum to 1 analytically, and softmax is invariant
// to the shared normalising constant).
//
// Two correlated continuous dimensions per item, Familiarity F and Recollection
// R, with corr(F, R) = tanh(rho):
//   target: F ~ N(dfam, 1),  R ~ N(drec, exp(sigmar))
//   lure:   F ~ N(0, 1),        R ~ N(0, 1)
// Old/new confidence is read off the aggregate strength S = F + R; the
// Remember/Know split is read off R against rcrit; the optional Know/Guess split
// is read off F against kcrit.
//
// Normal noise only: the bivariate-normal CDF needed for the Remember/Know split
// is exact and differentiable via Owen's T (owens_t is a Stan primitive).

// Bivariate standard-normal CDF P(Z1 <= h, Z2 <= k) with correlation r, via
// Owen's T (Owen, 1956). Where h and k have opposite signs the constant term
// 0.5 * (Phi(h) + Phi(k) - 1) is written with the complementary Phi of the
// positive argument, so it does not round away with 1 - Phi. An argument on
// its axis takes the limit form below.
real cdp_Phi2(real h, real k, real r) {
  if (h == negative_infinity() || k == negative_infinity()) return 0;
  if (h == positive_infinity()) return Phi(k);
  if (k == positive_infinity()) return Phi(h);
  real denom = sqrt((1 + r) * (1 - r));
  // on an axis the Owen's T argument k / h is a limit, not a value:
  // P(Z1 <= 0, Z2 <= k) = Phi(k) / 2 + T(k, r / denom), and the first-order
  // term carries the gradient phi(0) Phi(k / denom) that a constant stand-in
  // for h would lose (the value error is O(h^2), below 1e-20)
  if (abs(h) < 1e-10) {
    return 0.5 * Phi(k) + owens_t(k, r / denom)
           + h * 0.3989422804014327 * Phi(k / denom);
  }
  if (abs(k) < 1e-10) {
    return 0.5 * Phi(h) + owens_t(h, r / denom)
           + k * 0.3989422804014327 * Phi(h / denom);
  }
  real x = h;
  real y = k;
  real base;
  if (x < 0 && y > 0) {
    base = 0.5 * (Phi(x) - Phi(-y));
  } else if (y < 0 && x > 0) {
    base = 0.5 * (Phi(y) - Phi(-x));
  } else {
    base = 0.5 * (Phi(x) + Phi(y));
  }
  return base - owens_t(x, (y / x - r) / denom) - owens_t(y, (x / y - r) / denom);
}

// P(lo < Z < hi) for a standard normal Z, taken from the tail the interval
// lies in so an interval far above 0 does not cancel as 1 - 1. Infinite
// bounds never reach Phi(): its zero density times an infinite adjoint would
// make the gradient NaN.
real cdp_Phi_interval(real lo, real hi) {
  if (hi <= lo) return 0;
  if (lo > 0) {
    return is_inf(hi) ? Phi(-lo) : Phi(-lo) - Phi(-hi);
  }
  return (is_inf(hi) ? 1 : Phi(hi)) - (is_inf(lo) ? 0 : Phi(lo));
}

// P(a < X < b, Y < k) (above = 0) or P(a < X < b, Y > k) (above = 1) for a
// standard bivariate normal with correlation r. Y > k is evaluated as
// -Y < -k, and a band above 0 as its mirror image below 0, so both CDF
// values sit on the small side of the distribution.
real cdp_rect(real a, real b, real k, real r, int above) {
  real kk = above == 1 ? -k : k;
  real rr = above == 1 ? -r : r;
  real lo = a;
  real hi = b;
  if (a > 0) {
    lo = -b;
    hi = -a;
    rr = -rr;
  }
  return cdp_Phi2(hi, kk, rr) - cdp_Phi2(lo, kk, rr);
}

// Mass of an old-response region inside the strength bin (c_lo, c_hi),
// integrated over the strength S on the bin itself, with R | S = s normal
// (mean mu_R + beta * (s - mu_S), SD sd_c):
//   region 1, Guess:          s - kcrit < R < rcrit  (R < rcrit and F < kcrit)
//   region 2, Know-not-Guess: R < min(rcrit, s - kcrit)
// Guess needs S < rcrit + kcrit, so its range is finite. The range is clipped
// to where the strength density matters, split where the conditional
// probability steps (its conditional mean crosses rcrit or s - kcrit, and at
// rcrit + kcrit), and each segment is cut into pieces of at most one strength
// SD, each on 20 Gauss-Legendre nodes.
real cdp_region_mass(int region, real c_lo, real c_hi, real mu_S, real sigma_S,
                     real mu_R, real beta, real sd_c, real rcrit, real kcrit) {
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

  real hi = region == 1 ? fmin(c_hi, rcrit + kcrit) : c_hi;
  hi = fmin(hi, fmax(mu_S + 12 * sigma_S, c_lo + 4 * sigma_S));
  real lo = fmax(c_lo, fmin(mu_S - 12 * sigma_S, hi - 4 * sigma_S));
  if (hi <= lo) return 0;

  vector[3] cand;
  cand[1] = beta != 0 ? mu_S + (rcrit - mu_R) / beta : hi;
  cand[2] = beta != 1 ? (mu_R - beta * mu_S + kcrit) / (1 - beta) : hi;
  cand[3] = rcrit + kcrit;
  for (c in 1:3) {
    if (!(cand[c] > lo && cand[c] < hi)) cand[c] = hi;
  }
  vector[5] edges = append_row(append_row(lo, sort_asc(cand)), hi);

  real total = 0;
  for (e in 1:4) {
    real len = edges[e + 1] - edges[e];
    // a zero-length piece between coincident breakpoints adds no mass but
    // still carries the derivatives of its two edges, which move with
    // different parameters; skipping it dropped that boundary term
    int np = 1;
    while (np * sigma_S < len) np += 1;
    real half = 0.5 * len / np;
    for (p in 1:np) {
      real mid = edges[e] + (2 * p - 1) * half;
      for (i in 1:N_GL) {
        real s = mid + half * gl_nodes[i];
        real m = mu_R + beta * (s - mu_S);
        real pc = region == 1
                  ? cdp_Phi_interval((s - kcrit - m) / sd_c, (rcrit - m) / sd_c)
                  : Phi((fmin(rcrit, s - kcrit) - m) / sd_c);
        total += gl_weights[i] * half * exp(normal_lpdf(s | mu_S, sigma_S)) * pc;
      }
    }
  }
  return total;
}

// Confidence thresholds on the strength axis S = F + R, anchored so the old/new
// boundary (between bin n_new and bin n_new + 1) sits at `criterion`, with
// n_new - 1 thresholds below and n_old - 1 above. Reduces to the symmetric
// centred construction when n_new == n_old. thresh_type: 1 = parsimonious
// (Selker et al., 2019), 2 = equidistant, 3 = log_distance (free log widths;
// deltas[j] is the interval between thresholds j and j + 1, as in sdt_rating).
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
                       real dfam, real drec, real sigmar, real rho,
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

  real mu_F = stimulus > 0.5 ? dfam : 0.0;
  real mu_R = stimulus > 0.5 ? drec : 0.0;
  real sd_R = stimulus > 0.5 ? exp(sigmar) : 1.0;
  real corr = tanh(rho);
  real mu_S = mu_F + mu_R;
  real sigma_S = sqrt(square(sd_R + corr) + (1 - square(corr)));
  // the outer bins' infinite bounds stay constants: derived from mu_S and
  // sigma_S they would carry infinite partials into the gradient
  real z_lo = is_inf(c_lo) ? negative_infinity() : (c_lo - mu_S) / sigma_S;
  real z_hi = is_inf(c_hi) ? positive_infinity() : (c_hi - mu_S) / sigma_S;

  real p;
  if (type == 1) {
    p = cdp_Phi_interval(z_lo, z_hi);
  } else if (type == 4 || has_guess == 0) {
    // the smaller of Remember and Know comes from its own rectangle, the
    // larger as the rest of the bin, so neither is a difference of two
    // nearly equal numbers
    real rho_RS = (sd_R + corr) / sigma_S;
    real hcrit = (rcrit - mu_R) / sd_R;
    real rem = cdp_rect(z_lo, z_hi, hcrit, rho_RS, 1);
    real kn = cdp_rect(z_lo, z_hi, hcrit, rho_RS, 0);
    if (rem <= kn) {
      kn = cdp_Phi_interval(z_lo, z_hi) - rem;
    } else {
      rem = cdp_Phi_interval(z_lo, z_hi) - kn;
    }
    p = type == 4 ? rem : kn;
  } else {
    real beta = sd_R * (corr + sd_R) / square(sigma_S);
    real sd_c = sd_R * sqrt(fmax(1 - square(corr), 1e-12)) / sigma_S;
    p = cdp_region_mass(type == 2 ? 1 : 2, c_lo, c_hi, mu_S, sigma_S, mu_R,
                        beta, sd_c, rcrit, kcrit);
  }
  // an empty region (e.g. Guess in a bin above rcrit + kcrit) has mass 0;
  // the floor keeps its logit finite, so a zero count there adds 0, not NaN
  return fmax(p, 1e-300);
}

// The category logit `sdt_cdp_logmu` is code-generated per model by
// .sdt_cdp_logmu_stan() (R/model_sdt_cdp.R) so the per-distance log_distance
// deltas can be passed with a fixed arity, mirroring sdt_rating. It calls the
// helpers above (cdp_make_thresholds + cdp_category_prob).
