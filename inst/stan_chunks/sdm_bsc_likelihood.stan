  for (n in 1:N) {
    int K = cols(SDM_NT_FEATURES);
    array[K + 1] real item_angles;
    array[K + 1] real item_weights;
    int J = 1;
    real w_sum;
    real w_max;
    real R;
    real delta_min;
    int L;

    item_angles[1] = 0;
    item_weights[1] = exp(c[n]);

    for (k in 1:K) {
      if (SDM_LURE_IDX[n, k] > 0.5) {
        J += 1;
        item_angles[J] = SDM_NT_FEATURES[n, k];
        item_weights[J] = exp(c[n] - s[n] * SDM_NT_DISTANCES[n, k]);
      }
    }

    w_sum = sdm_spectral_weight_sum(item_weights, J);
    w_max = sdm_spectral_weight_max(item_weights, J);
    R = sdm_spectral_resultant_length(item_angles, item_weights, J);
    delta_min = sdm_spectral_min_separation(item_angles, J);
    L = sdm_get_L_general(kappa[n], J, w_sum, w_max, R, delta_min);

    target += -sdm_spectral_logz(
      mu[n],
      kappa[n],
      item_angles,
      item_weights,
      J,
      L,
      SDM_SPECTRAL_GRID_COS,
      SDM_SPECTRAL_GRID_SIN
    );
  }
