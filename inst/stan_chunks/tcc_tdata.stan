  // log probability of one option per distance (columns) at each value of the d grid (rows)
  matrix[tcc_grid_rows, rows(tcc_similarity)] tcc_tab;
  for (g in 1:tcc_grid_rows) {
    tcc_tab[g] = tcc_logp_by_distance(tcc_grid_lo + (g - 2) * tcc_grid_step, tcc_similarity,
                                      tcc_gl_nodes, tcc_gl_logw)';
  }
  tcc_tab[tcc_grid_rows] = fmin(tcc_tab[tcc_grid_rows], tcc_tab[tcc_grid_rows - 2]);
