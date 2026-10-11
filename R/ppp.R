############################################################################# !
# PPP.R                                                                  ####
# Posterior predictive p-values (T1, T2) for the count models.          ####
############################################################################# !

#' Posterior predictive p-values for count models
#'
#' Computes the goodness-of-fit statistics \eqn{T_1} and \eqn{T_2} of Klauer
#' (2010) and their posterior predictive p-values (PPP) for models that fit
#' aggregated response counts: [mpt()], [m3()], [sdt_yn()], [sdt_rating()],
#' [sdt_mafc()], [sdt_ranking()] and [sdt_cdp()]. These are the statistics that
#' TreeBUGS' `PPP()` reports.
#'
#' @details
#' Both statistics are computed once per posterior draw, for the observed counts
#' and for counts simulated from that draw with [brms::posterior_predict()],
#' each against the expected counts of the same draw from
#' [brms::posterior_epred()]. The simulated data are new data from the same
#' participants. The PPP is the proportion of draws in which the simulated
#' counts deviate more from the expectation than the observed counts do. A PPP
#' near 0 flags misfit; a PPP away from 0 means that the model reproduces the
#' aspect of the data the statistic measures.
#'
#' **\eqn{T_1}: mean structure.** The rows of the data are summed into cells,
#' and \eqn{T_1} is Pearson's \eqn{X^2} over all cells and response categories,
#' \deqn{T_1 = \sum_{c,k} \frac{(n_{ck} - \hat{n}_{ck})^2}{\hat{n}_{ck}},}
#' where \eqn{n_{ck}} is the observed and \eqn{\hat{n}_{ck}} the expected count
#' of category \eqn{k} in cell \eqn{c}. Klauer (2010, Eq. 17) and TreeBUGS
#' average the counts over participants instead of summing them. With the same
#' number of participants in every cell, the two give the same PPP; otherwise
#' the mean gives the cells with fewer participants more weight.
#'
#' **\eqn{T_2}: covariance structure.** For every participant, the counts are
#' summed per cell, giving one vector per participant over all cells and
#' categories. \eqn{T_2} compares the covariance matrix \eqn{s} of these vectors
#' across participants with the covariance matrix \eqn{\sigma} the model
#' expects (Klauer, 2010, Eq. 18):
#' \deqn{T_2 = \sum_{i,j} \frac{(s_{ij} - \sigma_{ij})^2}{\sqrt{\sigma_{ii}\sigma_{jj}}}.}
#' \eqn{\sigma} is the covariance of the expected counts across participants
#' plus \eqn{(T - 1) / T} times the average multinomial covariance within a
#' participant, where \eqn{T} is the number of participants. The sum runs over
#' the full matrix, so the covariances between trees or conditions, which
#' carry the correlations between the person parameters, are part of the
#' check. When participants differ in the cells they have, as in a
#' between-participants design, they are grouped by their set of cells,
#' \eqn{T_2} is computed per group, and the results are summed. A participant
#' who shares the set of cells with nobody else has no covariance and is left
#' out, with a message.
#'
#' **Cells.** By default, a cell is a combination of the tree ([mpt()]), the
#' stimulus (the SDT models) or the set size ([sdt_ranking()]), and the
#' categorical population-level predictors of the model formula. Numeric
#' predictors are not part of the default, because every distinct value would
#' form a cell of its own. Rows of a cell are summed, so the trees or stimuli
#' must be part of `by` for the statistic to compare each tree's own
#' categories.
#'
#' **Impossible categories.** Categories with an expected count of about 0 in
#' every draw, such as the impossible categories of an [mpt()] or [m3()] model
#' or the ranks beyond the set size of an [sdt_ranking()] model, are left out
#' of both statistics.
#'
#' @param fit A `bmmfit` object returned by [bmm()].
#' @param statistic Character. `"T1"`, `"T2"`, or both (the default). If the
#'   participants cannot be identified, the default computes `"T1"` only.
#' @param by Character vector naming the columns of the model data that define
#'   the cells. `NULL` (default) uses the cells described in Details;
#'   `character(0)` pools all rows into one cell. For [mpt()] fits, the
#'   tree column can be named even if it is not a predictor.
#' @param ndraws Integer. Number of posterior draws to use. `NULL` (default)
#'   uses all draws.
#' @param level Character. `"group"` (default) for the statistics of the whole
#'   data set; `"participant"` for one \eqn{T_1} and PPP per participant, as
#'   TreeBUGS reports them, to find the participants the model fits poorly.
#' @param participant Character. The column that identifies participants,
#'   needed for \eqn{T_2} and for `level = "participant"`. `NULL` (default)
#'   uses the grouping variable of the random effects; when the model has more
#'   than one (e.g. participants and items), name it here.
#' @param ... Unused.
#'
#' @return A data frame of class `bmm_ppp` with one row per statistic (and per
#'   participant for `level = "participant"`) and the columns `statistic`,
#'   `observed` and `replicated` (the means of the statistic over draws, for the
#'   observed and the simulated counts, as reported by Klauer, 2010, and
#'   TreeBUGS), `ppp` and `ndraws`. The draws of the statistics are stored in
#'   the attribute `"draws"`.
#'
#' @references Klauer, K. C. (2010). Hierarchical multinomial processing tree
#'   models: A latent-trait approach. *Psychometrika, 75*(1), 70–98.
#'   https://doi.org/10.1007/s11336-009-9141-0
#'
#' @seealso [pp_check.bmmfit()] for the graphical check.
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' trees <- list(
#'   mpt_tree("old", list(old = "D + (1 - D) * g", new = "(1 - D) * (1 - g)")),
#'   mpt_tree("new", list(old = "(1 - D) * g", new = "D + (1 - D) * (1 - g)"))
#' )
#' dat <- expand.grid(id = factor(1:20), item_type = c("old", "new"),
#'                    stringsAsFactors = FALSE)
#' dat$old <- rbinom(nrow(dat), 50, ifelse(dat$item_type == "old", 0.8, 0.2))
#' dat$new <- 50 - dat$old
#'
#' fit <- bmm(
#'   bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id)), dat, mpt(trees, "item_type"),
#'   cores = 4, backend = "cmdstanr"
#' )
#' ppp(fit)
#' ppp(fit, level = "participant")
ppp <- function(fit, ...) {
  UseMethod("ppp")
}

#' @rdname ppp
#' @export
ppp.bmmfit <- function(fit, statistic = c("T1", "T2"), by = NULL, ndraws = NULL,
                       level = c("group", "participant"), participant = NULL,
                       ...) {
  model <- fit$bmm$model
  stopif(!inherits(model, .ppp_models),
         "ppp() is available for the count models ({collapse_comma(.ppp_models)}), \\
          not for the {utils::tail(class(model), 1)} model.")
  statistic_given <- !missing(statistic)
  statistic <- match.arg(statistic, several.ok = TRUE)
  level <- match.arg(level)
  stopif(level == "participant" && statistic_given && "T2" %in% statistic,
         "T2 is a group-level statistic; level = 'participant' computes T1 only.")
  stopif(!is.null(participant) &&
           (!is.character(participant) || length(participant) != 1L ||
              !participant %in% names(fit$data)),
         "'participant' must name a column of the model data.")
  stopif(!is.null(ndraws) &&
           (!is.numeric(ndraws) || length(ndraws) != 1L || ndraws < 1 ||
              ndraws > brms::ndraws(fit)),
         "'ndraws' must be a number between 1 and {brms::ndraws(fit)}.")

  # one set of draws for posterior_predict() and posterior_epred(), so that
  # each replicated data set is compared with its own expectation
  draw_ids <- if (!is.null(ndraws)) sample.int(brms::ndraws(fit), ndraws)
  counts <- .ppp_counts(model, fit, draw_ids)
  by <- by %||% .ppp_default_by(model, fit, counts$data)
  stopif(!is.character(by) || !all(by %in% names(counts$data)),
         "'by' must name columns of the model data. Available: \\
          {collapse_comma(setdiff(names(counts$data), 'Y'))}.")

  needs_participant <- level == "participant" || "T2" %in% statistic
  participant <- participant %||% if (needs_participant) .ppp_participant(fit)
  if (needs_participant && length(participant) == 0L) {
    requested <- if (level == "participant") "level = 'participant'" else "T2"
    stopif(level == "participant" || statistic_given,
           "{requested} needs to know which rows belong to which participant, \\
            and the model has no random effects to tell. Name the column with \\
            'participant'.")
    message2("The model has no random effects to identify participants, so \\
              ppp() computes T1 only. Name the participant column with \\
              'participant' to get T2.")
    statistic <- "T1"
  }
  person <- if (length(participant) == 1L) as.factor(counts$data[[participant]])

  cell <- .ppp_index(counts$data[by])
  out <- if (level == "participant") {
    .ppp_by_participant(counts, cell, person)
  } else {
    .ppp_by_group(counts, cell, statistic, person)
  }
  structure(out$table, class = c("bmm_ppp", "data.frame"), draws = out$draws,
            level = level, by = by, n_cells = max(cell))
}

#' @export
print.bmm_ppp <- function(x, digits = 3, ...) {
  by <- attr(x, "by")
  cells <- if (length(by) == 0L) {
    "all rows pooled into one cell"
  } else {
    glue("{collapse_comma(by)} ({attr(x, 'n_cells')} cells)")
  }
  cat("Posterior predictive p-values (Klauer, 2010)\n")
  cat("Cells: ", cells, "\n\n", sep = "")
  print(as.data.frame(unclass(x)), digits = digits, row.names = FALSE)
  invisible(x)
}

.ppp_models <- c("mpt", "m3", "sdt_yn", "sdt_rating", "sdt_mafc", "sdt_ranking",
                 "sdt_cdp")


############################################################################# !
# COUNTS PER MODEL                                                       ####
############################################################################# !

# observed: rows x categories; replicated and expected: draws x rows x
# categories, on the count scale; data: the model data the cells come from
.ppp_counts <- function(model, fit, draw_ids) {
  UseMethod(".ppp_counts")
}

.ppp_counts.default <- function(model, fit, draw_ids) {
  list(
    observed = unclass(fit$data$Y),
    replicated = brms::posterior_predict(fit, draw_ids = draw_ids),
    expected = brms::posterior_epred(fit, draw_ids = draw_ids),
    data = fit$data
  )
}

# the response column holds one category; the other is trials - Y
.ppp_counts.sdt_yn <- function(model, fit, draw_ids) {
  response <- model$resp_vars$response
  trials <- fit$data[[model$other_vars$n_trials]]
  y <- fit$data[[response]]
  list(
    observed = matrix(c(y, trials - y), ncol = 2L,
                      dimnames = list(NULL, c(response, "other"))),
    replicated = .ppp_complement(brms::posterior_predict(fit, draw_ids = draw_ids),
                                 trials),
    expected = .ppp_complement(brms::posterior_epred(fit, draw_ids = draw_ids),
                               trials),
    data = fit$data
  )
}

# draws x rows counts of one category -> draws x rows x 2 with trials - count
.ppp_complement <- function(draws, trials) {
  array(c(draws, rep(trials, each = nrow(draws)) - draws), c(dim(draws), 2L))
}

.ppp_counts.sdt_mafc <- .ppp_counts.sdt_yn

# the tree column enters no formula, so brms leaves it out of the model data;
# the tree indicators that do enter it identify the tree of each row
.ppp_counts.mpt <- function(model, fit, draw_ids) {
  out <- NextMethod()
  indicators <- model$other_vars$indicators$tree
  tree_id <- model$other_vars$tree_id
  if (is.null(indicators) || tree_id %in% names(out$data)) {
    return(out)
  }
  tree <- max.col(as.matrix(out$data[indicators]), ties.method = "first")
  out$data[[tree_id]] <- factor(names(indicators)[tree], levels = names(indicators))
  out
}


############################################################################# !
# CELLS                                                                  ####
############################################################################# !

.ppp_default_by <- function(model, fit, data) {
  UseMethod(".ppp_default_by")
}

.ppp_default_by.default <- function(model, fit, data) {
  predictors <- intersect(unique(unlist(.population_vars(fit$bmm$user_formula))),
                          names(data))
  is_categorical <- vapply(data[predictors], function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  }, logical(1))
  predictors[is_categorical]
}

.ppp_default_by.mpt <- function(model, fit, data) {
  tree <- if (!is.null(model$other_vars$indicators$tree)) model$other_vars$tree_id
  unique(c(tree, NextMethod()))
}

.ppp_default_by.sdt <- function(model, fit, data) {
  unique(c(model$other_vars$stimulus, NextMethod()))
}

.ppp_default_by.sdt_ranking <- function(model, fit, data) {
  unique(c("max_rank", NextMethod()))
}

# Only plain grouping variables can identify a participant: an interaction
# such as id:session groups the rows of a participant further
.ppp_participant <- function(fit) {
  candidates <- grep(":", re_group_vars(fit$bmm$user_formula), fixed = TRUE,
                     invert = TRUE, value = TRUE)
  candidates <- intersect(candidates, names(fit$data))
  stopif(length(candidates) > 1L,
         "The model has several grouping variables ({collapse_comma(candidates)}). \\
          Name the one that identifies participants with 'participant'.")
  candidates
}

.ppp_index <- function(df) {
  if (ncol(df) == 0L) {
    return(rep(1L, nrow(df)))
  }
  as.integer(interaction(df, drop = TRUE, lex.order = TRUE))
}

# sums the rows of a draws x rows x categories array within each group
.ppp_rowsum <- function(x, group) {
  d <- dim(x)
  summed <- rowsum(matrix(aperm(x, c(2L, 1L, 3L)), d[2]), group)
  aperm(array(summed, c(nrow(summed), d[1], d[3])), c(2L, 1L, 3L))
}

# A category is possible in a cell if its expected share is not negligible in
# at least one draw. Impossible categories sit at exp(-100) relative to the
# others, far below any probability a fitted category reaches.
.ppp_possible <- function(observed, expected) {
  share <- expected / as.vector(apply(expected, c(1L, 2L), sum))
  possible <- apply(!is.na(share) & share > 1e-10, c(2L, 3L), any)
  impossible_seen <- observed > 0 & !possible
  stopif(any(impossible_seen),
         "The model gives probability 0 to observed responses in categories \\
          {collapse_comma(unique(colnames(observed)[col(observed)[impossible_seen]]))}.")
  possible
}


############################################################################# !
# STATISTICS                                                             ####
############################################################################# !

.ppp_by_group <- function(counts, cell, statistic, person) {
  draws <- list()
  if ("T1" %in% statistic) {
    draws$T1 <- .ppp_t1(counts, cell)
  }
  if ("T2" %in% statistic) {
    draws$T2 <- .ppp_t2(counts, cell, person)
  }
  table <- do.call(rbind, lapply(names(draws), function(stat) {
    .ppp_summary(draws[[stat]][, "observed", drop = FALSE],
                 draws[[stat]][, "replicated", drop = FALSE], stat)
  }))
  list(table = table, draws = draws)
}

# one T1 per participant over all of the participant's cells, as TreeBUGS'
# individual T1; the draws are draws x participants matrices
.ppp_by_participant <- function(counts, cell, person) {
  person_cell <- .ppp_index(data.frame(person, cell))
  terms <- .ppp_t1_terms(counts, person_cell)
  pc_person <- person[match(seq_len(max(person_cell)), person_cell)]
  per_person <- function(x) t(rowsum(t(apply(x, c(1L, 2L), sum)), pc_person))
  draws <- list(T1 = list(observed = per_person(terms$observed),
                          replicated = per_person(terms$replicated)))
  table <- data.frame(
    participant = factor(colnames(draws$T1$observed), levels(person)),
    .ppp_summary(draws$T1$observed, draws$T1$replicated, "T1")
  )
  list(table = table, draws = draws)
}

# observed and replicated are draws x statistics matrices
.ppp_summary <- function(observed, replicated, statistic) {
  data.frame(
    statistic = statistic,
    observed = colMeans(observed),
    replicated = colMeans(replicated),
    ppp = colMeans(replicated > observed),
    ndraws = nrow(observed),
    row.names = NULL
  )
}

.ppp_t1 <- function(counts, cell) {
  terms <- .ppp_t1_terms(counts, cell)
  ndraws <- dim(terms$observed)[1]
  cbind(observed = rowSums(matrix(terms$observed, ndraws)),
        replicated = rowSums(matrix(terms$replicated, ndraws)))
}

# the X^2 terms per draw x group x category, 0 for impossible categories
.ppp_t1_terms <- function(counts, group) {
  observed <- rowsum(counts$observed, group)
  expected <- .ppp_rowsum(counts$expected, group)
  replicated <- .ppp_rowsum(counts$replicated, group)
  ndraws <- dim(expected)[1]
  impossible <- !rep(.ppp_possible(observed, expected), each = ndraws)
  terms <- function(n) {
    x <- (n - expected)^2 / expected
    x[impossible] <- 0
    x
  }
  list(observed = terms(rep(observed, each = ndraws)),
       replicated = terms(replicated))
}

.ppp_t2 <- function(counts, cell, person) {
  possible <- .ppp_possible(rowsum(counts$observed, cell),
                            .ppp_rowsum(counts$expected, cell))
  cell_sets <- tapply(cell, person, function(x) paste(sort(unique(x)), collapse = "-"))
  patterns <- split(names(cell_sets), cell_sets)
  usable <- lengths(patterns) >= 2L
  if (any(!usable)) {
    message2("T2 leaves out {sum(lengths(patterns[!usable]))} participant(s) \\
              whose set of cells no other participant shares.")
  }
  if (!any(usable)) {
    warning2("No two participants share the same set of cells, so T2 cannot \\
              be computed.")
    missing_draws <- rep(NA_real_, dim(counts$expected)[1])
    return(cbind(observed = missing_draws, replicated = missing_draws))
  }
  Reduce(`+`, lapply(patterns[usable], function(members) {
    rows <- which(person %in% members)
    .ppp_t2_pattern(counts, rows, cell[rows], droplevels(person[rows]), possible)
  }))
}

# T2 for participants who have the same cells: the participants are the rows
# and the possible cell x category pairs the columns of the count matrices
.ppp_t2_pattern <- function(counts, rows, cell, person, possible) {
  n_people <- nlevels(person)
  cells <- sort(unique(cell))
  columns <- which(possible[cells, , drop = FALSE], arr.ind = TRUE)
  col_cell <- columns[, "row"]
  col_cat <- columns[, "col"]
  # every participant has every cell, so participant p's row sum for the j-th
  # cell is group (p - 1) * length(cells) + j
  person_cell <- (as.integer(person) - 1L) * length(cells) + match(cell, cells)
  lookup <- cbind(c(outer(seq_len(n_people) - 1L, col_cell,
                          function(p, j) p * length(cells) + j)),
                  rep(col_cat, each = n_people))
  wide <- function(per_group) matrix(per_group[lookup], n_people)

  observed <- counts$observed[rows, , drop = FALSE]
  observed_cov <- stats::cov(wide(rowsum(observed, person_cell)))
  expected <- .ppp_rowsum(counts$expected[, rows, , drop = FALSE], person_cell)
  replicated <- .ppp_rowsum(counts$replicated[, rows, , drop = FALSE], person_cell)
  weight <- 1 / sqrt(rowSums(observed))
  weight[!is.finite(weight)] <- 0
  row_col <- match(cell, cells)

  t2 <- function(s, sigma) {
    sum((s - sigma)^2 / sqrt(outer(diag(sigma), diag(sigma))))
  }
  t(vapply(seq_len(dim(expected)[1]), function(d) {
    expected_wide <- wide(expected[d, , ])
    # the average multinomial covariance of a participant: diag(n) - n n' / N
    # summed over the participant's rows, which are independent
    within <- diag(colMeans(expected_wide), length(col_cell))
    for (j in seq_along(cells)) {
      in_cell <- col_cell == j
      scaled <- matrix(counts$expected[d, rows[row_col == j], col_cat[in_cell]],
                       ncol = sum(in_cell)) * weight[row_col == j]
      within[in_cell, in_cell] <- within[in_cell, in_cell] - crossprod(scaled) / n_people
    }
    sigma <- stats::cov(expected_wide) + (n_people - 1) / n_people * within
    c(observed = t2(observed_cov, sigma),
      replicated = t2(stats::cov(wide(replicated[d, , ])), sigma))
  }, numeric(2)))
}
