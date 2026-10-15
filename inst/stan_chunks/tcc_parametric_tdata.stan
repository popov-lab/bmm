  // log P of one option per distance on the grid of d (rows) and of u = log par (array), one padding row at each end
  array[tcc_u_rows] matrix[tcc_grid_rows, rows(tcc_delta)] tcc_tab2;
  for (j in 1:tcc_u_rows) {
    vector[rows(tcc_delta)] f = tcc_similarity(tcc_delta, exp(tcc_u_lo + (j - 2) * tcc_u_step), tcc_shape);
    for (g in 1:tcc_grid_rows) {
      tcc_tab2[j][g] = tcc_logp_by_distance(tcc_grid_lo + (g - 2) * tcc_grid_step, f, tcc_gl_nodes, tcc_gl_logw)';
    }
  }
  // above the d grid a response's log P may fall but never rise (see .tcc_grid()): cap the end slope at 0
  for (j in 1:tcc_u_rows) {
    tcc_tab2[j][tcc_grid_rows] = fmin(tcc_tab2[j][tcc_grid_rows], tcc_tab2[j][tcc_grid_rows - 2]);
  }
