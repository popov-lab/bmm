  real mixture3p_cd_lpmf(array[] int y, vector mu, vector kappa, vector thetat,
                         vector thetant, vector criterion,
                         data array[] real probe, data matrix nt_features,
                         data matrix lure_idx, data vector gl_x,
                         data vector gl_w, data int free_crit) {
    int N = size(y);
    int n_nt = cols(nt_features);
    real out = 0;

    for (n in 1:N) {
      real n_active = sum(lure_idx[n]);
      real w_target = exp(thetat[n]);
      real w_nt = n_active > 0 ? exp(thetant[n]) : 0;
      real z = w_target + w_nt + 1;
      // p_s is the weight of the target component: everything else looks
      // uniform to the decision process (Lin & Oberauer, 2022, Eq. 5)
      real hw = cd_crit_angle(kappa[n], criterion[n], w_target / z, free_crit);
      real p_same;

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
            p_same += (w_nt / n_active) *
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
