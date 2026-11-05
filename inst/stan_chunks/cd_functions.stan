  // change-detection response y = 1 "change", 0 "same"; kept on the "same" scale so rare "same" responses keep their precision
  real cd_bernoulli_lpmf(int y, real p_same) {
    real p = fmin(fmax(p_same, machine_precision()), 1 - machine_precision());
    return y == 1 ? log1m(p) : log(p);
  }
