  // log P of one option per distance (columns) on the d grid with one padding row at each end (rows)
  matrix[tcc_grid_rows, rows(tcc_similarity)] tcc_tab;
  for (g in 1:tcc_grid_rows) {
    tcc_tab[g] = tcc_logp_by_distance(tcc_grid_lo + (g - 2) * tcc_grid_step, tcc_similarity,
                                      tcc_gl_nodes, tcc_gl_logw)';
  }
  // above the d grid a response's log P may fall but never rise (see .tcc_grid()): cap the end slope at 0
  tcc_tab[tcc_grid_rows] = fmin(tcc_tab[tcc_grid_rows], tcc_tab[tcc_grid_rows - 2]);
