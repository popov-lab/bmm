  // Lin & Oberauer (2022) change detection on the circmix retrieval mixture; derivations by the R twins in R/helpers-cd.R

  // log crossing level of the decision density against the uniform (Appendix B); +Inf never "same", -Inf always "same"
  real circmix_cd_offset(real criterion, real p_s) {
    // 0 whatever p_s, written with its slope in the criterion (-1 / p_s) so that a criterion at 0 keeps its gradient;
    // nothing stored takes the limit p_s -> 0
    if (criterion == 0) {
      return p_s > 0 ? -criterion / p_s : 0;
    }
    if (p_s <= 0) {
      return criterion < 0 ? positive_infinity() : negative_infinity();
    }
    // log(e^-criterion - (1 - p_s)) without the cancellation of expm1(-criterion) + p_s, which is total for a held item
    if (p_s >= 1) {
      return -criterion;
    }
    real scaled = log1m(p_s) + criterion;
    if (scaled >= 0) {
      return negative_infinity();
    }
    return -criterion + log1m_exp(scaled) - log(p_s);
  }

  // the scaled tanh bound, given the log of the smallest and largest value of 2 pi times the decision density
  real circmix_cd_bound(real criterion, real p_s, real log_min, real log_max) {
    if (criterion > 0) {
      // LLR_max = -log(1 + p_s (m - 1)): log1p keeps the digits of a small p_s (1 - m), which log_sum_exp loses to 1 - m
      // for a flat memory; near p_s (1 - m) = 1 (m underflows for a sharp held memory) only log_sum_exp is exact. p_s = 1
      // (a held slot) has no 1 - p_s term; log1m(1) would put an infinite partial on the stack
      real shift = p_s * expm1(log_min);
      real hi;
      if (p_s >= 1) {
        hi = -log_min;
      } else if (shift > -0.5) {
        hi = -log1p(shift);
      } else {
        hi = -log_sum_exp(log1m(p_s), log(p_s) + log_min);
      }
      return hi * tanh(criterion / hi);
    }
    real lo = log1p(p_s * expm1(log_max));
    return lo * tanh(criterion / lo);
  }

  // criterion kept inside the trial's range (LLR_min, LLR_max) of the log-likelihood ratio by a scaled tanh on each side of
  // 0 (slope 1 at 0, criterion = 0 unchanged), so that both responses stay possible: the range of the decision density
  // over the precision states (kappa, lw), the marginal one for the limited observer; the rich observer passes one state
  // at a time, so that each arc has its own range and vanishes with its weight
  real circmix_cd_bounded_criterion(real criterion, real p_s, vector kappa, vector lw) {
    if (criterion == 0 || p_s <= 0) {
      return criterion;
    }
    int S = num_elements(kappa);
    vector[S] base;
    for (n in 1:S) {
      base[n] = lw[n] - log_modified_bessel_first_kind(0, kappa[n]);
    }
    return circmix_cd_bound(criterion, p_s, log_sum_exp(base - kappa), log_sum_exp(base + kappa));
  }

  // log crossing level for the bounded criterion
  real circmix_cd_level(real criterion, real p_s, vector kappa, vector lw) {
    return circmix_cd_offset(circmix_cd_bounded_criterion(criterion, p_s, kappa, lw), p_s);
  }

  // half-width of the "same" arc around the probe for a von Mises decision density at kappa
  real circmix_cd_crit_angle(real kappa, real level) {
    // before any arithmetic on kappa: an unused node with an infinite partial gives kappa a NaN adjoint (0 * Inf)
    if (is_inf(level)) {
      return level > 0 ? 0 : pi();
    }
    real cos_d = (log_modified_bessel_first_kind(0, kappa) + level) / kappa;
    if (cos_d >= 1) {
      return 0;
    }
    if (cos_d <= -1) {
      return pi();
    }
    return acos(cos_d);
  }

  // Gauss-Legendre mass of vM(0, kappa) on [a, b]
  real circmix_cd_gl_vm(real a, real b, real kappa, real log_norm,
                        data vector gl_x, data vector gl_w) {
    real half = 0.5 * (b - a);
    return half * dot_product(gl_w, exp(kappa * cos(0.5 * (a + b) + half * gl_x) - log_norm));
  }

  // von Mises(kappa) mass of each component inside the arc of half-width hw around the probe; d = probe - location
  vector circmix_cd_arc_mass(vector d, real hw, real kappa, data vector gl_x,
                             data vector gl_w) {
    int K = num_elements(d);
    real log_norm = log_modified_bessel_first_kind(0, kappa) + log(2 * pi());
    real scale = inv_sqrt(kappa);
    array[7] real mult = {-18.0, -6.0, -2.0, 0.0, 2.0, 6.0, 18.0};
    vector[K] mass;
    for (k in 1:K) {
      real centre = atan2(sin(d[k]), cos(d[k]));
      real a = centre - hw;
      real b = centre + hw;
      if (hw <= 6 * scale) {
        mass[k] = circmix_cd_gl_vm(a, b, kappa, log_norm, gl_x, gl_w);
        continue;
      }
      vector[21] candidates;
      int n_in = 0;
      for (p in 1:3) {
        for (m in 1:7) {
          candidates[7 * (p - 1) + m] = 2 * pi() * (p - 2) + mult[m] * scale;
          n_in += candidates[7 * (p - 1) + m] > a && candidates[7 * (p - 1) + m] < b;
        }
      }
      vector[n_in + 2] cuts;
      int j = 1;
      cuts[1] = a;
      for (c in 1:21) {
        if (candidates[c] > a && candidates[c] < b) {
          j += 1;
          cuts[j] = candidates[c];
        }
      }
      cuts[n_in + 2] = b;
      cuts = sort_asc(cuts);
      mass[k] = 0;
      for (i in 1:(n_in + 1)) {
        mass[k] += circmix_cd_gl_vm(cuts[i], cuts[i + 1], kappa, log_norm, gl_x, gl_w);
      }
    }
    return mass;
  }

  // P("same") for a given arc: retrieval mass inside it, guessing contributes hw / pi
  real circmix_cd_psame_at(vector d, vector logw, real logw_guess, real kappa,
                           real hw, data vector gl_x, data vector gl_w) {
    if (hw <= 0) {
      return 0;
    }
    if (hw >= pi()) {
      return 1;
    }
    return exp(logw_guess) * hw / pi()
           + dot_product(exp(logw), circmix_cd_arc_mass(d, hw, kappa, gl_x, gl_w));
  }

  // P("same") at one precision; p_s is the observer's prior that the probed item is stored
  real circmix_cd_psame(vector d, vector logw, real logw_guess, real kappa,
                        real criterion, real p_s, data vector gl_x,
                        data vector gl_w) {
    real level = circmix_cd_level(criterion, p_s, [kappa]', [0.0]');
    return circmix_cd_psame_at(d, logw, logw_guess, kappa,
                               circmix_cd_crit_angle(kappa, level), gl_x, gl_w);
  }

  // knowledge-limited boundary: root in (0, pi) of log sum_n exp(lw_n) vM(hw | kappa_n) * 2 pi = level, safeguarded Newton
  real circmix_cd_limited_hw(vector kappa, vector lw, real level, real start) {
    int S = num_elements(kappa);
    if (is_inf(level)) {
      return level > 0 ? 0.0 : pi();
    }
    vector[S] base;
    for (n in 1:S) {
      base[n] = lw[n] - log_modified_bessel_first_kind(0, kappa[n]);
    }
    if (log_sum_exp(base + kappa) <= level) {
      return 0;
    }
    if (log_sum_exp(base - kappa) >= level) {
      return pi();
    }
    real lo = 0;
    real hi = pi();
    real h = start;
    if (!(h > lo && h < hi)) {
      h = 0.5 * pi();
    }
    for (iter in 1:60) {
      vector[S] lp = base + kappa * cos(h);
      real f = log_sum_exp(lp) - level;
      if (f > 0) {
        lo = h;
      } else {
        hi = h;
      }
      real h_new = h + f / (sin(h) * dot_product(softmax(lp), kappa));
      if (!(h_new > lo && h_new < hi)) {
        h_new = 0.5 * (lo + hi);
      }
      if (abs(h_new - h) < 1e-12) {
        return h_new;
      }
      h = h_new;
    }
    return h;
  }

  // P("same") over precision states (kappa_n, normalised log weight lw_n, -Inf for padding); rich: a boundary per state
  // from that state alone, limited: one shared boundary from the marginal density
  real circmix_cd_mix_psame(vector d, vector logw, real logw_guess, vector kappa,
                            vector lw, real criterion, real p_s, int rich,
                            data vector gl_x, data vector gl_w) {
    int S = num_elements(kappa);
    vector[S] ps = rep_vector(0, S);
    if (rich || S == 1) {
      for (n in 1:S) {
        if (!is_inf(lw[n])) {
          real level = circmix_cd_level(criterion, p_s, [kappa[n]]', [0.0]');
          ps[n] = circmix_cd_psame_at(d, logw, logw_guess, kappa[n],
                                      circmix_cd_crit_angle(kappa[n], level), gl_x, gl_w);
        }
      }
    } else {
      real level = circmix_cd_level(criterion, p_s, kappa, lw);
      real start = circmix_cd_crit_angle(exp(log_sum_exp(lw + log(kappa))), level);
      real hw = circmix_cd_limited_hw(kappa, lw, level, start);
      for (n in 1:S) {
        ps[n] = circmix_cd_psame_at(d, logw, logw_guess, kappa[n], hw, gl_x, gl_w);
      }
    }
    return dot_product(exp(lw), ps);
  }

  // a held item at a whole K / ss: P("same") of the kappa_lo state alone, written as its first-order expansion in extra
  // (extra = 0, so the value is unchanged) to give K the slope from above without log(extra) on the stack, as
  // circmix_slot_averaging_ld() does. Rich: P is linear in extra, each state having its own arc. Limited: the marginal
  // range and the shared arc move with extra too; the arc's slope is the implicit one, from
  // log((1 - extra) vM(hw | kappa_lo) + extra vM(hw | kappa_hi)) = level
  real circmix_cd_whole_slots_psame(vector d, vector logw, real kappa_lo, real kappa_hi,
                                    real extra, real criterion, int rich,
                                    data vector gl_x, data vector gl_w) {
    real hw_lo;
    real hw_hi;
    if (rich) {
      hw_lo = circmix_cd_crit_angle(kappa_lo, circmix_cd_level(criterion, 1, [kappa_lo]', [0.0]'));
      hw_hi = circmix_cd_crit_angle(kappa_hi, circmix_cd_level(criterion, 1, [kappa_hi]', [0.0]'));
    } else {
      real log_i0_lo = log_modified_bessel_first_kind(0, kappa_lo);
      real log_i0_hi = log_modified_bessel_first_kind(0, kappa_hi);
      real log_min = -kappa_lo - log_i0_lo
                     + log1p(extra * expm1(kappa_lo - kappa_hi + log_i0_lo - log_i0_hi));
      real log_max = kappa_lo - log_i0_lo
                     + log1p(extra * expm1(kappa_hi - kappa_lo + log_i0_lo - log_i0_hi));
      hw_lo = circmix_cd_crit_angle(kappa_lo, circmix_cd_offset(
        circmix_cd_bound(criterion, 1, log_min, log_max), 1));
      if (hw_lo > 0 && hw_lo < pi()) {
        hw_lo += extra * expm1((kappa_hi - kappa_lo) * cos(hw_lo) + log_i0_lo - log_i0_hi)
                 / (kappa_lo * sin(hw_lo));
      }
      hw_hi = hw_lo;
    }
    real p_lo = circmix_cd_psame_at(d, logw, negative_infinity(), kappa_lo, hw_lo, gl_x, gl_w);
    return p_lo + extra * (circmix_cd_psame_at(d, logw, negative_infinity(), kappa_hi, hw_hi,
                                               gl_x, gl_w) - p_lo);
  }

  // slot averaging; branch from circmix_slot_averaging_branches(), [held, kappa_lo, kappa_hi, extra]. With a slot held the
  // observer knows the item is stored (p_s = 1) and the two slot counts are the precision states; with none held the item
  // is guessed with 1 - extra and at kappa_hi with extra, which is also the observer's prior of storage. log(extra) at a
  // whole K / set size is log(0), whose NaN adjoint would reach K: with no slot held (K = 0) extra = 0 leaves out the
  // state it weights, with one held the value is written in a form whose gradient is the slope from above
  real circmix_cd_slot_averaging_psame(vector d, vector logw, vector branch,
                                       real criterion, int rich,
                                       data vector gl_x, data vector gl_w) {
    real extra = branch[4];
    if (branch[1] == 0) {
      if (extra == 0) {
        return circmix_cd_psame(d, rep_vector(negative_infinity(), num_elements(d)), 0,
                                branch[3], criterion, 0, gl_x, gl_w);
      }
      return circmix_cd_psame(d, log(extra) + logw, log1m(extra), branch[3],
                              criterion, extra, gl_x, gl_w);
    }
    if (extra == 0) {
      return circmix_cd_whole_slots_psame(d, logw, branch[2], branch[3], extra, criterion, rich,
                                          gl_x, gl_w);
    }
    return circmix_cd_mix_psame(d, logw, negative_infinity(), branch[2:3],
                                [log1m(extra), log(extra)]', criterion, 1, rich,
                                gl_x, gl_w);
  }
