############################################################################# !
# DIAGNOSTIC CHECKS FOR MULTIVARIATE SPECIFICATIONS                      ####
############################################################################# !

# Cross-component correlations exist only for random-effects terms that share
# an |ID| label AND a grouping variable across components, and they are
# identified only by grouping values (e.g. persons) observed in several of
# them. Both failures are silent in brms: without a shared label the
# components are estimated jointly but independently, and with zero overlap
# (typically IDs coded differently per dataset) the posterior of the
# correlations just reproduces the LKJ prior
check_shared_random_effects <- function(spec) {
  re_list <- lapply(spec, function(component) formula_re_frame(component$formula))
  pairs <- unique(do.call(rbind, re_list))
  pairs <- pairs[!is.na(pairs$id), , drop = FALSE]

  members <- lapply(seq_len(nrow(pairs)), function(i) {
    which(vapply(seq_along(spec), function(j) {
      any(re_list[[j]]$group == pairs$group[i] & re_list[[j]]$id == pairs$id[i], na.rm = TRUE) &&
        pairs$group[i] %in% names(spec[[j]]$data)
    }, logical(1)))
  })
  warnif(
    !any(lengths(members) >= 2),
    "No random-effects ID is shared between the components. To estimate the \\
    correlations between subject-level parameters across components, use the \\
    |ID| syntax with the same ID in at least two components, e.g. \\
    (1 | p | id). The model will run, but no cross-component correlations \\
    will be estimated."
  )

  for (i in which(lengths(members) >= 2)) {
    group <- pairs$group[i]
    id_label <- pairs$id[i]
    value_sets <- lapply(spec[members[[i]]], function(component) {
      unique(as.character(component$data[[group]]))
    })
    shared <- Reduce(intersect, value_sets)
    total <- unique(unlist(value_sets))
    warnif(
      length(shared) == 0,
      "The components tied by the random-effects ID '{id_label}' have no \\
      values of the grouping variable '{group}' in common. The correlations \\
      between their parameters cannot be estimated from these data, and \\
      their posterior will just reproduce the prior. Check that the values \\
      of '{group}' are coded identically across the component datasets."
    )
    if (length(shared) > 0 && length(shared) < length(total)) {
      message2(
        "{length(shared)} of {length(total)} unique values of '{group}' \\
        appear in all components tied by the random-effects ID '{id_label}'. \\
        The correlations between the components' parameters are informed by \\
        these {length(shared)} groups."
      )
    }
  }
}

# the random-effects terms of a bmmformula as brms parses them: one row per
# (grouping variable, |ID| label) pair, with NA for terms without a label.
# brms keeps the label with its surrounding whitespace
formula_re_frame <- function(formula) {
  forms <- Filter(is_formula, unclass(formula))
  frames <- lapply(forms, function(f) brms::brmsterms(brms::bf(f))$dpars$mu$re)
  frames <- lapply(Filter(Negate(is.null), frames), function(re) {
    data.frame(group = trimws(re$group), id = trimws(re$id))
  })
  empty <- data.frame(group = character(), id = character())
  unique(do.call(rbind, c(list(empty), frames)))
}

# Factors used as fixed-effect predictors must have identical levels across
# components: the stacked data unions the levels, and levels absent from a
# component's rows produce all-zero design columns and rank deficiency
warn_factor_level_mismatch <- function(spec) {
  fe_vars <- lapply(spec, function(component) formula_fe_vars(component$formula))
  for (var in unique(unlist(fe_vars))) {
    levels_list <- lapply(seq_along(spec), function(i) {
      column <- spec[[i]]$data[[var]]
      if (var %in% fe_vars[[i]] && is.factor(column)) levels(column)
    })
    levels_list <- Filter(Negate(is.null), levels_list)
    mismatch <- length(levels_list) > 1 &&
      !all(vapply(levels_list, identical, logical(1), y = levels_list[[1]]))
    warnif(
      mismatch,
      "The factor '{var}' is used as a predictor in several components but \\
      has different levels across their datasets. Levels that do not occur \\
      in a component's data will produce inestimable coefficients for that \\
      component."
    )
  }
}

# fixed-effect variables of a bmmformula, excluding random-effects terms
formula_fe_vars <- function(formula) {
  forms <- Filter(is_formula, unclass(formula))
  unique(unlist(lapply(forms, function(f) {
    trms <- try(stats::terms(f), silent = TRUE)
    if (is_try_error(trms)) {
      return(character(0))
    }
    labels <- attr(trms, "term.labels")
    labels <- labels[!grepl("|", labels, fixed = TRUE)]
    unlist(lapply(labels, function(label) all.vars(str2lang(label))))
  })))
}
