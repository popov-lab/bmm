// Racing diffusion model (Tillman, Van Zandt & Logan, 2020): Wald accumulators
// with a start point uniform on [0, A] and threshold b = gap + A, i.e. a
// distance to threshold uniform on [gap, gap + A]. Density and survival are
// assembled in log space from the antiderivative of the plain Wald survivor
// in the threshold, so that the tails stay finite where the raw closed forms
// underflow. Keep in step with .dwald_full() / .pwald_full() in
// R/distributions.R, which are the same derivation in R.
//
// Spelling rules. An upper normal tail is std_normal_lcdf(-z), never
// the lccdf of z: the two are equal on Stan Math 5.4 (CmdStan 2.40) but
// on 5.3 (rstan / StanHeaders 2.39) the lccdf is -inf from z = 8.26 with an
// infinite partial, which reached this model as a NaN trial for a loser with
// drift 5 from t = 4 s. cogmod 0.3.2 (Makowski) exposed that flaw in bmm's
// cswald survivor and spells its own RDM with the reflected lcdf; the
// derivation here is bmm's own. log Phi(z) goes through rdm_log_Phi() below:
// erfc where it is representable, for an exact value and an exact gradient,
// the lcdf below that. A log-difference is log_diff_exp of the two
// log-probabilities from the tail that keeps both away from 1, and every
// call is guarded by swald_log_diff_exp(): log_diff_exp(x, x) has infinite
// partials that poison the gradient of the whole trial even in a branch whose
// value is discarded. Stan Math has no inverse Gaussian yet (stan-dev/math
// #3382 adds inv_gaussian_* in log space); nothing here depends on it.

// log Phi(z). Stan's Phi() evaluates 0.5 (1 + erf(z / sqrt 2)) between -5 and
// 0, which cancels to a relative error of 1e-16 / Phi(z) -- five digits gone at
// z = -4.4 -- and the differences g(beta) - g(alpha) and Phi(beta) - Phi(alpha)
// below amplify that by their own cancellation (8.5e-7 in a far-tail log
// survival, measured). erfc has no such cancellation, its value is exact to
// 1e-16 down to the underflow at z = -37.5, and its derivative is the exact
// normal density; the lcdf takes over below, at the cost of its approximate
// (1e-5 relative) partial there.
real rdm_log_Phi(real z) {
  if (z < -37.5) return std_normal_lcdf(z | );
  return log(0.5 * erfc(-0.7071067811865475 * z));
}

// log(Phi(b) - Phi(a)) for b > a
real rdm_log_Phi_diff(real a, real b) {
  if (b <= a) return negative_infinity();
  if (a >= 0) return swald_log_diff_exp(rdm_log_Phi(-a), rdm_log_Phi(-b));
  if (b <= 0) return swald_log_diff_exp(rdm_log_Phi(b), rdm_log_Phi(a));
  // the straddling case has both probabilities O(1) and is exact directly;
  // closer than 1e-8 the subtraction loses digits and the midpoint density,
  // exact to second order in b - a, is the better form
  if (b - a < 1e-8) return log(b - a) + std_normal_lpdf(0.5 * (a + b));
  return log(Phi(b) - Phi(a));
}

// log(y Phi(y) + phi(y)), the antiderivative of Phi, given log Phi(y). Above
// -1 the sum is direct (g(-1) = 0.083, no cancellation, exact gradient). Below
// it the two terms cancel to order y^-2 and the difference is taken in log
// space; below -12 that difference itself loses digits, and g is written as
// phi(y) (1 - |y| M(|y|)) with M the Mills ratio and the bracket expanded
// asymptotically (Abramowitz & Stegun 26.2.12): sum_j (-1)^(j+1) (2j-1)!! /
// y^(2j). Twelve terms and the log-space form agree to about 1e-12 at -12.
real rdm_log_g(real y, real log_Phi_y) {
  if (y < -12) {
    real y2 = square(y);
    real term = 1;
    real total = 0;
    for (j in 1:12) {
      term *= (2.0 * j - 1) / y2;
      total += (j % 2 == 1) ? term : -term;
    }
    return -0.5 * y2 - 0.9189385332046727 + log(total);
  }
  real log_phi_y = -0.5 * square(y) - 0.9189385332046727;
  if (y < -1) return swald_log_diff_exp(log_phi_y, log(-y) + log_Phi_y);
  return log(y * exp(log_Phi_y) + exp(log_phi_y));
}

// Below this ratio of start-point range to the diffusion spread s sqrt(t)
// the difference quotients over the start point cancel; the plain Wald at
// the midpoint threshold is second-order accurate in A and takes over.
// Tillman et al. (2020), Eq. 5, generalised to a diffusion scale s:
// f = (1 / A) [drift (Phi(beta) - Phi(alpha)) + (s / sqrt(t)) (phi(alpha) -
// phi(beta))], alpha and beta the standardised distances at the two ends of
// the start-point range. The first bracket is positive; the second takes the
// sign of |beta| - |alpha|. t is the decision time, > 0.
real rdm_log_pdf(real t, real drift, real gap, real A, real s) {
  real st = sqrt(t);
  if (A < 1e-4 * s * st) return swald_lpdf(t | drift, gap + 0.5 * A, 0, s);

  real denom = s * st;
  real alpha = (gap - drift * t) / denom;
  real beta = (gap + A - drift * t) / denom;
  real l1 = log(drift) + rdm_log_Phi_diff(alpha, beta);
  real lpa = -0.5 * square(alpha) - 0.9189385332046727;
  real lpb = -0.5 * square(beta) - 0.9189385332046727;
  real lnum;
  if (abs(alpha) <= abs(beta)) {
    lnum = log_sum_exp(l1, swald_log_diff_exp(lpa, lpb) + log(s / st));
  } else {
    lnum = swald_log_diff_exp(l1, swald_log_diff_exp(lpb, lpa) + log(s / st));
  }
  return lnum - log(A);
}

// Tillman et al. (2020), Appendix A, from the antiderivative of the plain
// Wald survivor in the threshold. With g(y) = y Phi(y) + phi(y),
// q = s^2 / (2 drift), E(u) = exp(2 u drift / s^2) Phi(-(u + drift t) /
// (s sqrt(t))) and P = Phi(beta) - Phi(alpha):
//   S * A = D1 - D2,  D1 = s sqrt(t) (g(beta) - g(alpha)),
//                     D2 = q (P + E(gap + A) - E(gap)).
// D1, D2 and P are positive (g increases; d/du of the survivor's
// antiderivative pieces are signed by the Wald reflection identity); the
// difference of the two E terms is signed and resolved explicitly. The E
// terms stay single exponents so that they survive where exp(2 u drift / s^2)
// and Phi(.) separately overflow and underflow. D1 - D2 cancels only in the
// far tail, by a factor of order drift t / gap, as the plain Wald survivor
// does. t is the decision time, > 0.
real rdm_log_surv(real t, real drift, real gap, real A, real s) {
  real st = sqrt(t);
  if (A < 1e-4 * s * st) return swald_lccdf(t | drift, gap + 0.5 * A, 0, s);

  real denom = s * st;
  real b = gap + A;
  real alpha = (gap - drift * t) / denom;
  real beta = (b - drift * t) / denom;
  real lPa = rdm_log_Phi(alpha);
  real lPb = rdm_log_Phi(beta);
  real lD1 = log(denom)
             + swald_log_diff_exp(rdm_log_g(beta, lPb), rdm_log_g(alpha, lPa));

  real s2 = square(s);
  real lEb = 2 * drift * b / s2 + rdm_log_Phi(-(b + drift * t) / denom);
  real lEk = 2 * drift * gap / s2 + rdm_log_Phi(-(gap + drift * t) / denom);
  // the lower-half difference reuses the two log Phi already in hand
  real lP = (beta <= 0) ? swald_log_diff_exp(lPb, lPa) : rdm_log_Phi_diff(alpha, beta);
  real lD2;
  if (lEb >= lEk) {
    lD2 = log_sum_exp(lP, swald_log_diff_exp(lEb, lEk));
  } else {
    lD2 = swald_log_diff_exp(log_sum_exp(lP, lEb), lEk);
  }
  lD2 += 2 * log(s) - log(2 * drift);

  return fmin(swald_log_diff_exp(lD1, lD2) - log(A), 0);
}

real rdm_log_lik_one(real rt, array[] real drift, real gap, real ndt,
                     real s, real sp, int response, array[] int n,
                     int use_start_var) {
  int n_cats = size(drift);
  real b = gap + sp;
  real t = rt - ndt;
  real lp;

  if (t <= 0) return negative_infinity();

  if (use_start_var == 0) {
    lp = log(n[response]) + swald_lpdf(rt | drift[response], b, ndt, s);
  } else {
    lp = log(n[response]) + rdm_log_pdf(t, drift[response], gap, sp, s);
  }

  // single-pass race: the winner contributes n_win - 1 survival copies, each
  // loser n_j copies. The strict reps > 0 guard skips the winner's survival
  // when its category has one accumulator (all correct trials; every trial
  // for K = 2) and keeps zero-count or underflowed terms out of the sum
  // (0 * -inf would otherwise poison the likelihood with NaN).
  for (j in 1:n_cats) {
    int reps = (j == response) ? n[j] - 1 : n[j];
    if (reps > 0) {
      real log_surv;
      if (use_start_var == 0) {
        log_surv = swald_lccdf(rt | drift[j], b, ndt, s);
      } else {
        log_surv = rdm_log_surv(t, drift[j], gap, sp, s);
      }
      lp += reps * log_surv;
    }
  }

  return lp;
}
