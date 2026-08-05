  // Retrieval weights of the interference measurement model: the target draws
  // context-dependent activation c and context-independent activation a, each
  // non-target draws c attenuated by the spatial generalization gradient over
  // its distance from the probed location, and the background contributes a
  // constant. These are the same activations the continuous reproduction imm
  // builds through its mixture formula.
  real imm_cd_full_lpmf(array[] int y, vector mu, vector kappa, vector c,
                        vector a, vector s, vector criterion,
                        data array[] real probe, data matrix nt_features,
                        data matrix nt_distances, data matrix lure_idx,
                        data vector gl_x, data vector gl_w, data int free_crit) {
    int N = size(y);
    int n_nt = cols(nt_features);
    real out = 0;

    for (n in 1:N) {
      // c, a and s arrive on the natural scale: brms has already applied the
      // log links, so the activations are used directly here
      real w_target = c[n] + a[n];
      vector[n_nt] w_nt = rep_vector(0, n_nt);
      real z = w_target + 1;
      real hw;
      real p_same;

      for (k in 1:n_nt) {
        if (lure_idx[n, k] > 0.5) {
          w_nt[k] = c[n] * exp(-s[n] * nt_distances[n, k]) + a[n];
          z += w_nt[k];
        }
      }

      // Eq. B.21: the target's mixture weight normalised by the sum of all of
      // them. It cancels out of the boundary at an unbiased criterion.
      hw = cd_crit_angle(kappa[n], criterion[n], w_target / z, free_crit);

      if (hw <= 0) {
        p_same = 0;
      } else if (hw >= pi()) {
        p_same = 1;
      } else {
        p_same = w_target *
          cd_vm_arc_mass(mu[n], probe[n], hw, kappa[n], gl_x, gl_w) +
          hw * inv(pi());
        for (k in 1:n_nt) {
          if (lure_idx[n, k] > 0.5) {
            p_same += w_nt[k] *
              cd_vm_arc_mass(nt_features[n, k], probe[n], hw, kappa[n],
                             gl_x, gl_w);
          }
        }
        p_same /= z;
      }

      out += cd_bernoulli_lpmf(y[n] | p_same);
    }

    return out;
  }
