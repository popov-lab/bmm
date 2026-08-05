  real mixture2p_cd_lpmf(array[] int y, vector mu, vector kappa, vector thetat,
                         vector criterion, data array[] real probe,
                         data vector gl_x, data vector gl_w,
                         data int free_crit) {
    int N = size(y);
    real out = 0;

    for (n in 1:N) {
      real hw = cd_crit_angle(kappa[n], criterion[n], thetat[n], free_crit);
      real p_same;

      if (hw <= 0) {
        p_same = 0;
      } else if (hw >= pi()) {
        p_same = 1;
      } else {
        p_same = thetat[n] *
          cd_vm_arc_mass(mu[n], probe[n], hw, kappa[n], gl_x, gl_w) +
          (1 - thetat[n]) * hw * inv(pi());
      }

      out += cd_bernoulli_lpmf(y[n] | p_same);
    }

    return out;
  }
