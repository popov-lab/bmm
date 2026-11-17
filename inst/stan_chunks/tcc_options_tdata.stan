  // each cell's first row and its choices per offered option, from the responses; the tables from R
  // must match the data Stan received
  array[G_tcc_cells] int tcc_cell_first = rep_array(0, G_tcc_cells);
  array[G_tcc_cells, M_tcc_sets] int tcc_cell_n = rep_array(0, G_tcc_cells, M_tcc_sets);
  for (n in 1:N) {
    int g = tcc_row_cell[n];
    if (g < 1 || g > G_tcc_cells || tcc_cell_set[g] < 1 || tcc_cell_set[g] > S_tcc_sets) {
      reject("bmm error: row ", n, " has no cell in the TCC option tables. Please report this ",
             "at https://github.com/popov-lab/bmm/issues");
    }
    int s = tcc_cell_set[g];
    int j = 0;
    for (k in 1:tcc_set_size[s]) if (tcc_set_off[s, k] == Y[n]) j = k;
    if (j == 0) {
      reject("bmm error: the response in row ", n, " is not among the options of its cell. ",
             "Please report this at https://github.com/popov-lab/bmm/issues");
    }
    if (tcc_cell_first[g] == 0) tcc_cell_first[g] = n;
    tcc_cell_n[g, j] += 1;
  }
  for (g in 1:G_tcc_cells) {
    if (tcc_cell_first[g] == 0) {
      reject("bmm error: cell ", g, " of the TCC option tables has no rows. Please report this ",
             "at https://github.com/popov-lab/bmm/issues");
    }
  }
