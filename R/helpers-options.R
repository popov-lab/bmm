# Discrete response options (#558, D3, D4). A model with `n_options` reads the
# response as a choice among n_options equally spaced positions on the circle,
# offset from the target in whole steps; with `options` (columns) or a number
# of equally spaced `n_alternatives`, each trial offers only some of them. The
# offered set per trial, taken from the data, is the general case: nothing
# here assumes that sets are constant over trials or equally spaced, and the
# full wheel and the n_alternatives shorthand are special cases of it. The
# class `option_set` carries these steps for sdm() and, with TCC-K3, tcc()

# signed offset of the nearest option from the target, in (-n_options / 2,
# n_options / 2]: -pi and pi are one option
.bin_options <- function(x, n_options, units) {
  steps <- if (units == "radians") floor(x * n_options / (2 * pi) + 0.5) else x
  half <- n_options / 2
  as.integer(half - (half - steps) %% n_options)
}

# the offset an OptIdx column holds where a trial offers fewer options than
# there are columns; outside the range of .bin_options()
.option_absent <- function(n_options) {
  as.integer(n_options)
}

.option_bin_column <- function(model) {
  if (inherits(model, "tcc")) "tcc_bin" else "sdm_bin"
}

.option_columns <- function(n) {
  if (n == 0) character(0) else paste0("OptIdx", seq_len(n))
}

# the arguments of a constructor with response options, checked once for
# sdm() and the d/r functions
.check_n_options <- function(n_options) {
  stopif(
    !is.numeric(n_options) || length(n_options) != 1 || !is.finite(n_options) || n_options < 4 ||
      n_options %% 2 != 0,
    "n_options must be a single even whole number of at least 4, the number of \\
    equally spaced positions on the response wheel."
  )
  as.integer(n_options)
}

.check_option_args <- function(n_options, options, n_alternatives, units) {
  stopif(
    is.null(n_options) && (!is.null(options) || !is.null(n_alternatives) || units == "bins"),
    "`options`, `n_alternatives` and units = 'bins' need `n_options`, the number of \\
    equally spaced positions on the response wheel that the options sit on."
  )
  if (is.null(n_options)) {
    return(NULL)
  }
  n_options <- .check_n_options(n_options)
  stopif(
    !is.null(options) && (!is.character(options) || length(options) == 0),
    "`options` must be the column names of the offered options, one column per option."
  )
  if (is.null(n_alternatives)) {
    return(nlist(n_options))
  }
  if (is.character(n_alternatives) && length(n_alternatives) == 1) {
    stopif(is.null(options), "A column name in `n_alternatives` needs the `options` columns.")
    return(nlist(n_options, n_alternatives))
  }
  stopif(
    !is.numeric(n_alternatives) || length(n_alternatives) != 1 || is.na(n_alternatives) ||
      n_alternatives < 2 || n_alternatives %% 1 != 0,
    "`n_alternatives` must be a column name or a single whole number of at least 2."
  )
  stopif(
    !is.null(options),
    "Give the offered options either as columns in `options` or as a number of \\
    equally spaced options in `n_alternatives`, not both."
  )
  stopif(
    n_options %% n_alternatives != 0,
    "n_options = {n_options} is not a multiple of {n_alternatives}, so \\
    {n_alternatives} equally spaced options do not sit on its positions."
  )
  nlist(n_options, n_alternatives = as.integer(n_alternatives))
}

############################################################################# !
# CHECK_DATA                                                             ####
############################################################################# !

#' @export
check_data.option_set <- function(model, data, formula) {
  resp <- model$resp_vars[[1]]
  n_options <- model$other_vars$n_options
  units <- model$other_vars$units
  stopif(!is.numeric(data[[resp]]), "The response variable '{resp}' must be numeric.")
  .check_option_units(data[[resp]], units, glue("The response variable '{resp}'"))
  bin_col <- .option_bin_column(model)
  width <- .option_width(model)
  reserved <- intersect(c(bin_col, .option_columns(width), if (width > 0) "opt_set_id"), colnames(data))
  warnif(
    length(reserved) > 0,
    "Column{if (length(reserved) > 1) 's' else ''} {collapse_comma(reserved)} in your \\
    data {if (length(reserved) > 1) 'are' else 'is'} reserved by bmm and will be overwritten."
  )
  data[[bin_col]] <- .bin_options(data[[resp]], n_options, units)
  offsets <- .option_offsets(model, data)
  if (!is.null(offsets)) {
    .check_option_sets(offsets, data[[bin_col]], n_options)
    absent <- .option_absent(n_options)
    for (j in seq_len(ncol(offsets))) {
      column <- offsets[, j]
      column[is.na(column)] <- absent
      data[[.option_columns(ncol(offsets))[j]]] <- column
    }
    data$opt_set_id <- .option_set_id(offsets)
  }
  NextMethod("check_data")
}

.check_option_units <- function(x, units, what) {
  if (units == "bins") {
    stopif(
      any(x %% 1 != 0, na.rm = TRUE),
      "With units = 'bins', {what} must hold whole numbers of option steps."
    )
  } else {
    stopif(
      isTRUE(max(abs(x), na.rm = TRUE) > 2 * pi),
      "{what} must be in radians, but values beyond 2*pi show they are not. \\
      Convert them with deg2rad(), or set units = 'bins' if they count option steps."
    )
  }
}

# The offered options per row as integer offsets from the target, NA where a
# row offers fewer than the widest; NULL for the full wheel. NA cells may sit
# anywhere in a row and are moved to its end
.option_offsets <- function(model, data) {
  n_options <- model$other_vars$n_options
  options <- model$other_vars$options
  n_alternatives <- model$other_vars$n_alternatives
  if (is.null(options)) {
    if (is.null(n_alternatives)) {
      return(NULL)
    }
    spaced <- .bin_options((seq_len(n_alternatives) - 1) * n_options / n_alternatives, n_options, "bins")
    return(matrix(spaced, nrow(data), n_alternatives, byrow = TRUE))
  }
  missing <- setdiff(options, colnames(data))
  stopif(length(missing) > 0, "The option column{if (length(missing) > 1) 's' else ''} \\
         {collapse_comma(missing)} {if (length(missing) > 1) 'are' else 'is'} not in the data.")
  units <- model$other_vars$units
  values <- as.matrix(data[options])
  stopif(!is.numeric(values), "The option columns must be numeric.")
  .check_option_units(values, units, "The option columns")
  used <- !is.na(data[[model$resp_vars[[1]]]])
  steps <- if (units == "radians") values * n_options / (2 * pi) else values
  off_grid <- abs(steps - round(steps))
  off_grid[!used, ] <- 0
  if (isTRUE(max(off_grid, na.rm = TRUE) > 0.1)) {
    at <- which(off_grid > 0.1, arr.ind = TRUE)[1, ]
    stop2("The option in column '{options[at[2]]}' of row {at[1]} lies \\
          {signif(off_grid[at[1], at[2]], 2)} of a step from the nearest of the \\
          {n_options} positions on the wheel. Options must sit on those positions: \\
          check n_options and the units of the option columns.")
  }
  if (isTRUE(max(off_grid, na.rm = TRUE) > 1e-6)) {
    message2("Options up to {signif(max(off_grid, na.rm = TRUE), 2)} of a step off the \\
             grid were rounded onto the {n_options} positions on the wheel.")
  }
  offsets <- matrix(.bin_options(round(steps), n_options, "bins"), nrow(values))
  n_alternatives <- model$other_vars$n_alternatives
  if (is.character(n_alternatives)) {
    .check_option_count(data[[n_alternatives]], rowSums(!is.na(offsets)), used, n_alternatives)
  }
  # NA padding moves to the end of each row
  rows <- lapply(seq_len(nrow(offsets)), function(i) {
    row <- offsets[i, ]
    c(row[!is.na(row)], row[is.na(row)])
  })
  matrix(unlist(rows), nrow(offsets), byrow = TRUE)
}

.check_option_count <- function(count, offered, used, name) {
  stopif(is.null(count), "The column '{name}' named in `n_alternatives` is not in the data.")
  row <- which(used & is.na(count))
  stopif(length(row) > 0, "'{name}' is missing in row {row[1]}.")
  row <- which(used & count != offered)
  stopif(
    length(row) > 0,
    "'{name}' says {count[row[1]]} options in row {row[1]}, but the option columns \\
    hold {offered[row[1]]}."
  )
}

.check_option_sets <- function(offsets, response, n_options) {
  used <- which(!is.na(response))
  offered <- offsets[used, , drop = FALSE]
  target <- rowSums(offered == 0, na.rm = TRUE) == 0
  stopif(
    any(target),
    "The target (an offset of 0) is not among the options in row \\
    {used[which(target)[1]]}. Each trial must offer the target."
  )
  chosen <- rowSums(offered == response[used], na.rm = TRUE) == 0
  stopif(
    any(chosen),
    "The response is not one of the options in row {used[which(chosen)[1]]}."
  )
  twice <- apply(offered, 1, function(row) anyDuplicated(row[!is.na(row)]) > 0)
  stopif(
    any(twice),
    "The data offers the same option twice in row {used[which(twice)[1]]} \\
    (with n_options = {n_options}, offsets are taken modulo the wheel)."
  )
}

# rows that offer the same options, in any order, share an id, numbered by
# first appearance
.option_set_id <- function(offsets) {
  keys <- apply(offsets, 1, function(row) paste(sort(row[!is.na(row)]), collapse = ","))
  match(keys, unique(keys))
}

# the number of OptIdx columns check_data() builds: one per option column, or
# n_alternatives for the equally spaced sets
.option_width <- function(model) {
  options <- model$other_vars$options
  if (!is.null(options)) {
    return(length(options))
  }
  if (is.numeric(model$other_vars$n_alternatives)) model$other_vars$n_alternatives else 0L
}

#' @exportS3Method
built_data_columns.option_set <- function(model) {
  width <- .option_width(model)
  c(.option_bin_column(model), .option_columns(width), if (width > 0) "opt_set_id",
    NextMethod("built_data_columns"))
}

# the binned response is the brms response; the option offsets ride along as
# vint terms, so that log_lik() and posterior_predict() read each trial's set
# from the data, newdata included. The likelihood itself reads the sets from
# the cell tables configure_model() builds
#' @export
bmf2bf.option_set <- function(model, formula) {
  width <- .option_width(model)
  vint <- if (width > 0) glue(" | vint({paste(.option_columns(width), collapse = ', ')})") else ""
  brms::bf(glue("{.option_bin_column(model)}{vint} ~ 1"))
}

# brms keeps the binned response and the option offsets, not the user's
# columns, which are rebuilt at the option centres
#' @exportS3Method
revert_check_data.option_set <- function(model, data) {
  n_options <- model$other_vars$n_options
  to_units <- function(x) if (model$other_vars$units == "radians") x * 2 * pi / n_options else x
  rebuild <- function(data, name, value) {
    if (!name %in% colnames(data)) {
      data[[name]] <- value
      attr(data, "rebuilt") <- c(attr(data, "rebuilt"), name)
    }
    data
  }
  bin_col <- .option_bin_column(model)
  data <- rebuild(data, model$resp_vars[[1]], to_units(data[[bin_col]]))
  columns <- .option_columns(.option_width(model))
  options <- model$other_vars$options
  if (!is.null(options)) {
    offsets <- as.matrix(data[columns])
    offsets[offsets == .option_absent(n_options)] <- NA
    for (j in seq_along(options)) {
      data <- rebuild(data, options[j], to_units(offsets[, j]))
    }
    n_alternatives <- model$other_vars$n_alternatives
    if (is.character(n_alternatives)) {
      data <- rebuild(data, n_alternatives, as.integer(rowSums(!is.na(offsets))))
    }
  }
  data[intersect(c(bin_col, columns, "opt_set_id"), colnames(data))] <- NULL
  NextMethod("revert_check_data")
}

# The normaliser of a trial depends on its parameters and its option set, so
# the rows that share both share one: a cell is a run of equal parameters
# times an option set. Rows of a cell need not be adjacent, since each cell's
# normaliser is computed once, at its first row, and weighted by its count
.option_cells <- function(run_id, set_id) {
  key <- paste(run_id, set_id)
  first <- which(!duplicated(key))
  list(
    first = first,
    count = tabulate(match(key, key[first]), length(first)),
    set = set_id[first]
  )
}

# The option sets the rows use, numbered by first appearance, as sorted
# offsets padded with 0 to the widest; set_id is the number of each row's set
.option_set_table <- function(offsets, n_options) {
  absent <- .option_absent(n_options)
  keys <- apply(offsets, 1, function(row) paste(sort(row[row != absent]), collapse = ","))
  set_id <- match(keys, unique(keys))
  sets <- lapply(which(!duplicated(keys)), function(i) sort(offsets[i, offsets[i, ] != absent]))
  size <- lengths(sets)
  off <- matrix(0L, length(sets), max(size))
  for (s in seq_along(sets)) {
    off[s, seq_len(size[s])] <- sets[[s]]
  }
  nlist(set_id, size, off)
}
