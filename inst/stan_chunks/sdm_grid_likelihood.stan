  target += -sdm_grid_cells_lse(mu, c, kappa, sdm_activation, sdm_scaling, sdm_K, G_sdm_cells,
                                sdm_cell_first, sdm_cell_count, sdm_cell_set, sdm_set_size,
                                sdm_set_off, sdm_grid_half, sdm_grid_logmult, sdm_grid_full,
                                sdm_gl_x, sdm_gl_w, sdm_wide_cosm1, sdm_wide_w);
  target += -N * log(2 * pi() / sdm_K);
