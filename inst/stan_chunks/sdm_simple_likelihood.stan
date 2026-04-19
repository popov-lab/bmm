  real current_c = negative_infinity();
  real current_kappa = negative_infinity();
  real current_logz = 0;

  for (n in 1:N) {
    array[1] real item_angles;
    array[1] real item_weights;
    int L;

    item_angles[1] = 0;
    item_weights[1] = exp(c[n]);

    if (n == 1 || c[n] != current_c || kappa[n] != current_kappa) {
      L = sdm_get_L_general(
        kappa[n],
        1,
        item_weights[1],
        item_weights[1],
        item_weights[1],
        pi()
      );
      current_logz = sdm_spectral_logz(
        mu[n],
        kappa[n],
        item_angles,
        item_weights,
        1,
        L,
        SDM_SPECTRAL_GRID_COS,
        SDM_SPECTRAL_GRID_SIN
      );
      current_c = c[n];
      current_kappa = kappa[n];
    }

    target += -current_logz;
  }
