  // Lin & Oberauer (2022) change detection on the circmix retrieval mixture; derivations by the R twins in R/helpers-cd.R

  // log crossing level of the decision density against the uniform (Appendix B); +Inf never "same", -Inf always "same"
  real circmix_cd_offset(real criterion, real p_s) {
    if (criterion == 0) {
      return 0;
    }
    if (p_s <= 0) {
      return criterion < 0 ? positive_infinity() : negative_infinity();
    }
    real scaled = expm1(-criterion) + p_s;
    if (scaled <= 0) {
      return negative_infinity();
    }
    return log(scaled) - log(p_s);
  }

  // half-width of the "same" arc around the probe for a von Mises decision density at kappa
  real circmix_cd_crit_angle(real kappa, real level) {
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
    return circmix_cd_psame_at(d, logw, logw_guess, kappa,
                               circmix_cd_crit_angle(kappa, circmix_cd_offset(criterion, p_s)),
                               gl_x, gl_w);
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

  // P("same") over precision states (kappa_n, normalised log weight lw_n); rich: a boundary per state, limited: one shared boundary
  real circmix_cd_mix_psame(vector d, vector logw, real logw_guess, vector kappa,
                            vector lw, real criterion, real p_s, int rich,
                            data vector gl_x, data vector gl_w) {
    int S = num_elements(kappa);
    real level = circmix_cd_offset(criterion, p_s);
    vector[S] ps;
    if (rich || S == 1) {
      for (n in 1:S) {
        ps[n] = circmix_cd_psame_at(d, logw, logw_guess, kappa[n],
                                    circmix_cd_crit_angle(kappa[n], level), gl_x, gl_w);
      }
    } else {
      real start = circmix_cd_crit_angle(exp(log_sum_exp(lw + log(kappa))), level);
      real hw = circmix_cd_limited_hw(kappa, lw, level, start);
      for (n in 1:S) {
        ps[n] = circmix_cd_psame_at(d, logw, logw_guess, kappa[n], hw, gl_x, gl_w);
      }
    }
    return dot_product(exp(lw), ps);
  }

  // number of precision states: one under constant precision (shape +Inf), else the J grid
  int circmix_cd_n_states(real shape, int nodes) {
    return is_inf(shape) ? 1 : nodes;
  }

  // precision states [kappa, normalised log weight] of a memory with mean kappa; shape from circmix_vp_shape()
  matrix circmix_cd_vp_states(real kappa, real shape, real tau, int nodes,
                              data vector logk, data vector dlogk,
                              data real logJ_min, data real dlogJ) {
    if (is_inf(shape)) {
      return [[kappa, 0]];
    }
    matrix[nodes, 2] grid = circmix_vp_grid(shape, tau, nodes);
    matrix[nodes, 2] states;
    for (i in 1:nodes) {
      states[i, 1] = circmix_kappa(exp(grid[i, 1]), logk, dlogk, logJ_min, dlogJ);
    }
    states[ : , 2] = log_softmax(col(grid, 2));
    return states;
  }

  // P("same") with precision marginalised over gamma(J(kappa) / tau, scale = tau); tau = 0 is constant precision
  real circmix_cd_vp_psame(vector d, vector logw, real logw_guess, real kappa,
                           real tau, real criterion, real p_s, int rich,
                           int nodes, data vector logk, data vector dlogk,
                           data real logJ_min, data real dlogJ,
                           data vector gl_x, data vector gl_w) {
    real shape = circmix_vp_shape(kappa, tau, nodes);
    matrix[circmix_cd_n_states(shape, nodes), 2] states
      = circmix_cd_vp_states(kappa, shape, tau, nodes, logk, dlogk, logJ_min, dlogJ);
    return circmix_cd_mix_psame(d, logw, logw_guess, col(states, 1), col(states, 2),
                                criterion, p_s, rich, gl_x, gl_w);
  }

  // slot averaging; branch from circmix_slot_averaging_branches(); p_s = 1 when a slot is held, else the probability of one
  real circmix_cd_slot_averaging_psame(vector d, vector logw, vector branch,
                                       real tau, real criterion, int rich,
                                       int nodes, data vector logk,
                                       data vector dlogk, data real logJ_min,
                                       data real dlogJ, data vector gl_x,
                                       data vector gl_w) {
    if (branch[1] == 0) {
      return circmix_cd_vp_psame(d, branch[5] + logw, branch[4], branch[3], tau,
                                 criterion, exp(branch[5]), rich, nodes, logk,
                                 dlogk, logJ_min, dlogJ, gl_x, gl_w);
    }
    real shape_lo = circmix_vp_shape(branch[2], tau, nodes);
    real shape_hi = circmix_vp_shape(branch[3], tau, nodes);
    matrix[circmix_cd_n_states(shape_lo, nodes), 2] lo
      = circmix_cd_vp_states(branch[2], shape_lo, tau, nodes, logk, dlogk, logJ_min, dlogJ);
    matrix[circmix_cd_n_states(shape_hi, nodes), 2] hi
      = circmix_cd_vp_states(branch[3], shape_hi, tau, nodes, logk, dlogk, logJ_min, dlogJ);
    return circmix_cd_mix_psame(d, logw, negative_infinity(),
                                append_row(col(lo, 1), col(hi, 1)),
                                append_row(branch[4] + col(lo, 2), branch[5] + col(hi, 2)),
                                criterion, 1, rich, gl_x, gl_w);
  }
