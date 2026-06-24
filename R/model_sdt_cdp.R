############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_cdp <- function(judgment = NULL, confidence = NULL, stimulus = NULL,
                           count = NULL, dist = "normal",
                           threshold_type = "parsimonious",
                           n_new = NULL, n_old = NULL,
                           links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)
  thresh_type_int <- .sdt_threshold_type_id(threshold_type)
  delta_based <- threshold_type %in% c("log_distance", "log_ratio", "softmax")

  # Threshold parameters are either a single `spacing` (parsimonious/equidistant)
  # or per-distance `deltaN` (log_distance), with the anchor at the old/new
  # boundary n_new. They must be declared here so configure_prior/create_initfun
  # see them, hence n_new/n_old are required up front for delta-based types.
  thr_parts <- .sdt_threshold_parameter_parts(
    n_ratings = if (delta_based) n_new + n_old else NULL,
    threshold_type = threshold_type,
    anchor = if (delta_based) n_new else 0L
  )

  parameters <- c(
    list(
      dprimef = glue("Familiarity sensitivity: target mean on the familiarity axis"),
      dprimer = glue("Recollection sensitivity: target mean on the recollection axis"),
      criterion = glue("Response bias: old/new boundary on the strength (F+R) axis")
    ),
    thr_parts$parameters,
    list(
      rcrit = glue("Remember criterion: threshold on the recollection axis"),
      sigmar = glue(
        "Log SD of the recollection target distribution, so exp(sigmar) is the ",
        "SD and 0 means SD = 1"
      ),
      rho = glue(
        "Familiarity-recollection correlation on an unconstrained scale; ",
        "tanh(rho) is the correlation and 0 means independent processes"
      ),
      kcrit = glue(
        "Know criterion: threshold on the familiarity axis that splits Know from ",
        "Guess. Active only when the data include 'guess' judgments"
      )
    )
  )
  default_priors <- c(
    list(
      dprimef = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
      dprimer = list(main = "normal(1, 1)", effects = "normal(0, 0.5)"),
      criterion = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)")
    ),
    thr_parts$default_priors,
    list(
      rcrit = list(main = "normal(0, 1)", effects = "normal(0, 0.5)"),
      sigmar = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)"),
      rho = list(main = "normal(0, 0.5)", effects = "normal(0, 0.3)"),
      kcrit = list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
    )
  )
  param_links <- c(
    list(dprimef = "identity", dprimer = "identity", criterion = "identity"),
    thr_parts$param_links,
    list(rcrit = "identity", sigmar = "identity", rho = "identity",
         kcrit = "identity")
  )
  # sigmar, rho, and kcrit are fixed off by default (equal recollection variance,
  # independent processes, no Know/Guess split) and freed through the formula.
  fixed_parameters <- list(sigmar = 0, rho = 0, kcrit = -100)

  requirements <- glue(
    "Provide tidy long data with the following columns:", "\n\n",
    "  - judgment: 'new', 'remember', 'know' (and optionally 'guess')", "\n",
    "  - confidence: the rating on the unified old/new scale (1 = most ",
    "confident 'new' ... K = most confident 'old')", "\n",
    "  - stimulus: 0 = new/lure, 1 = old/target", "\n",
    "  - count: response counts (omit for trial-level data)"
  )

  init_ranges <- list(
    dprimef = c(0.3, 1.0), dprimer = c(0.3, 1.0), criterion = c(-0.3, 0.3),
    rcrit = c(-0.5, 0.5), sigmar = c(-0.1, 0.1), rho = c(-0.2, 0.2),
    kcrit = c(-1.0, 0.0)
  )
  for (p in names(thr_parts$parameters)) {
    init_ranges[[p]] <- if (p == "spacing") c(-0.7, -0.2) else c(-0.5, 0.2)
  }

  out <- structure(
    list(
      resp_vars = nlist(judgment, confidence, count),
      other_vars = nlist(stimulus, dist, dist_int, threshold_type,
                         thresh_type_int, delta_based, n_new, n_old),
      domain = "Recognition Memory",
      task = "Old/New Recognition with Remember/Know Judgments",
      name = "Continuous Dual-Process Signal Detection Theory (CDP)",
      citation = glue(
        "Wixted, J. T., & Mickes, L. (2010). A continuous dual-process model ",
        "of remember/know judgments. Psychological Review, 117(4), 1025-1054."
      ),
      version = "cdp",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = fixed_parameters,
      default_priors = default_priors,
      init_ranges = init_ranges,
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdt", "sdt_cdp"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Continuous Dual-Process Signal Detection Theory Model
#' @name sdt_cdp
#' @details `r model_info(.model_sdt_cdp())`
#'
#' The continuous dual-process model (Wixted & Mickes, 2010) assumes two
#' continuous memory signals per item, Familiarity (F) and Recollection (R),
#' with correlation `tanh(rho)`. Old/new confidence is read off the aggregate
#' strength S = F + R; the Remember/Know judgment splits "old" responses on R
#' (against `rcrit`); an optional Know/Guess split uses F (against `kcrit`).
#'
#' **Response format.** Data are supplied in tidy long form: a `judgment` column
#' (`"new"`, `"remember"`, `"know"`, and optionally `"guess"`), a `confidence`
#' column on the unified old/new scale (1 = most confident "new" to K = most
#' confident "old"; "new" judgments occupy the low levels, "old" judgments the
#' high levels), a `stimulus` column (0 = new/lure, 1 = old/target), and an
#' optional `count` column (omit for trial-level data, which is aggregated
#' internally). The number of "new" and "old" confidence levels need not be
#' equal (e.g. the 1-new / 5-old scale of Rotello et al., 2005). Whether the
#' Know/Guess split is modelled is **driven by the data**: include a `"guess"`
#' judgment level to fit the three-way Remember/Know/Guess model.
#'
#' **Variants via the formula.** By default `sigmar`, `rho`, and `kcrit` are
#' fixed (equal recollection variance, independent F and R, no Know/Guess
#' split), recovering the classic independent-process CDP. Free them through the
#' formula:
#' \itemize{
#'   \item `sigmar ~ 1` estimates unequal recollection variance.
#'   \item `rho ~ 1` estimates the within-item F-R correlation; `rho` also
#'     accepts predictors (e.g. `rho ~ condition`), bounded to (-1, 1) via an
#'     internal `tanh`. This structural, within-item correlation is distinct
#'     from a between-subject correlation of random effects, which is available
#'     separately through brms syntax such as `(1 |p| id)` on any parameter.
#'   \item `kcrit ~ 1` estimates the Know/Guess criterion (requires `"guess"`
#'     judgments in the data).
#' }
#'
#' When no Remember/Know split is available (confidence ratings only), use
#' [sdt_rating()] instead.
#'
#' @param judgment The name of the column coding the memory judgment, with
#'   values `"new"`, `"remember"`, `"know"`, and optionally `"guess"`.
#' @param confidence The name of the column coding confidence on the unified
#'   old/new scale (integer, 1 = most confident "new" to K = most confident
#'   "old").
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (new/lure) and 1 (old/target).
#' @param count Optional name of a column of response counts. If `NULL`
#'   (default), each row is treated as a single trial and counts are aggregated
#'   internally.
#' @param dist The noise distribution. Only `"normal"` is currently supported
#'   (the CDP model is inherently Gaussian).
#' @param threshold_type Character. Threshold parameterization on the strength
#'   axis: `"parsimonious"` (default) and `"equidistant"` use a single `spacing`
#'   parameter; `"log_distance"` (Meyer-Grant et al., 2025) estimates the
#'   `n_new + n_old - 2` distances between adjacent thresholds freely, each as a
#'   `deltaN` parameter on the log scale (the distance leading into threshold
#'   `N` from the old/new boundary). `"log_distance"` requires `n_new`/`n_old`.
#' @param n_new,n_old Integer numbers of "new" and "old" confidence levels.
#'   Required for `threshold_type = "log_distance"` (the per-distance parameters
#'   are fixed at construction); for `"parsimonious"`/`"equidistant"` they are
#'   ignored and detected from the data.
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Wixted, J. T., & Mickes, L. (2010). A continuous dual-process model of
#'   remember/know judgments. \emph{Psychological Review}, \emph{117}(4),
#'   1025--1054. \doi{10.1037/a0020874}
#'
#' Rotello, C. M., Macmillan, N. A., Reeder, J. A., & Wong, M. (2005). The
#'   remember response: Subject to bias, graded, and not a process-pure
#'   indicator of recollection. \emph{Psychonomic Bulletin & Review},
#'   \emph{12}(5), 865--873. \doi{10.3758/BF03196778}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # Simulate a Remember/Know data set (3 new + 3 old confidence levels)
#' dat <- rsdt_cdp(n_per_cell = 200, n_subjects = 20,
#'                 dprimef = 0.8, dprimer = 1.0, criterion = 0,
#'                 spacing = -0.4, rcrit = 0.7, n_new = 3, n_old = 3)
#'
#' model <- sdt_cdp(
#'   judgment = "judgment", confidence = "confidence",
#'   count = "count", stimulus = "stimulus"
#' )
#'
#' fit <- bmm(
#'   formula = bmf(dprimef ~ 1, dprimer ~ 1, criterion ~ 1,
#'                 spacing ~ 1, rcrit ~ 1),
#'   data = dat, model = model, backend = "cmdstanr"
#' )
#'
#' # Estimate unequal recollection variance and the F-R correlation
#' fit_uv <- bmm(
#'   formula = bmf(dprimef ~ 1, dprimer ~ 1, criterion ~ 1,
#'                 spacing ~ 1, rcrit ~ 1, sigmar ~ 1, rho ~ 1),
#'   data = dat, model = model, backend = "cmdstanr"
#' )
#' }
sdt_cdp <- function(judgment, confidence, stimulus,
                    count = NULL, dist = "normal",
                    threshold_type = c("parsimonious", "equidistant",
                                       "log_distance"),
                    n_new = NULL, n_old = NULL,
                    links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  threshold_type <- match.arg(threshold_type)
  stopif(!identical(dist, "normal"),
         "sdt_cdp currently supports only dist = 'normal'; other noise \\
         distributions are deferred to a future release")

  if (threshold_type == "log_distance") {
    stopif(is.null(n_new) || is.null(n_old),
           "n_new and n_old are required when threshold_type = 'log_distance' \\
           (the per-distance threshold parameters depend on the scale size)")
    stopif(n_new < 1 || n_old < 1 || n_new + n_old < 3,
           "n_new and n_old must be >= 1 and sum to >= 3")
  }

  .model_sdt_cdp(judgment = judgment, confidence = confidence,
                 stimulus = stimulus, count = count, dist = dist,
                 threshold_type = threshold_type, n_new = n_new, n_old = n_old,
                 links = links, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_cdp <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  judg_var <- model$resp_vars$judgment
  conf_var <- model$resp_vars$confidence
  count_var <- model$resp_vars$count

  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stopif(!all(unique(data[[stim_var]]) %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (new/lure) and 1 (old/target)")
  stopif(!judg_var %in% colnames(data),
         "Judgment variable '{judg_var}' missing in the data")
  stopif(!conf_var %in% colnames(data),
         "Confidence variable '{conf_var}' missing in the data")
  stopif(!is.null(count_var) && !count_var %in% colnames(data),
         "Count variable '{count_var}' missing in the data")

  judg <- as.character(data[[judg_var]])
  allowed <- c("new", "remember", "know", "guess")
  bad <- setdiff(unique(judg), allowed)
  stopif(length(bad) > 0,
         "Judgment values {collapse_comma(bad)} not allowed; use {collapse_comma(allowed)}")
  stopif(!all(c("new", "remember", "know") %in% judg),
         "Judgment column must contain 'new', 'remember', and 'know' values")
  has_guess <- "guess" %in% judg

  conf <- data[[conf_var]]
  warnif(any(conf != round(conf), na.rm = TRUE),
         "Confidence values should be integers")
  conf <- as.integer(round(conf))
  stopif(any(conf < 1, na.rm = TRUE), "Confidence values must be >= 1")

  is_new <- judg == "new"
  n_new <- max(conf[is_new])
  K_full <- max(conf)
  n_old <- K_full - n_new
  stopif(n_old < 1,
         "Data must contain 'old' confidence levels above the 'new' levels")
  # For delta-based thresholds the parameter set is fixed at construction from
  # the declared scale, so the data must match it exactly.
  if (isTRUE(model$other_vars$delta_based)) {
    stopif(n_new != model$other_vars$n_new || n_old != model$other_vars$n_old,
           "Declared n_new = {model$other_vars$n_new}, n_old = \\
           {model$other_vars$n_old} do not match the data \\
           (n_new = {n_new}, n_old = {n_old})")
  }
  stopif(!all(conf[is_new] %in% seq_len(n_new)),
         "'new' confidence values must occupy levels 1..{n_new}")
  stopif(!all(conf[!is_new] %in% (n_new + 1L):K_full),
         "old (remember/know/guess) confidence values must occupy levels {n_new + 1L}..{K_full}")

  counts_in <- if (!is.null(count_var)) {
    vals <- data[[count_var]]
    stopif(any(vals < 0, na.rm = TRUE), "Count values must be non-negative")
    as.integer(vals)
  } else {
    rep(1L, nrow(data))
  }

  # Canonical category column: new(1..n_new), [guess], know, remember.
  block <- integer(nrow(data))
  if (has_guess) {
    block[judg == "know"] <- 1L
    block[judg == "remember"] <- 2L
  } else {
    block[judg == "remember"] <- 1L
  }
  col_idx <- ifelse(is_new, conf, n_new + block * n_old + (conf - n_new))
  K_cat <- n_new + (if (has_guess) 3L else 2L) * n_old

  cell_cols <- setdiff(colnames(data),
                       c(judg_var, conf_var, if (!is.null(count_var)) count_var))
  cell_key <- do.call(paste, c(data[cell_cols], sep = "\r"))
  first_idx <- which(!duplicated(cell_key))
  cell_data <- data[first_idx, cell_cols, drop = FALSE]
  rownames(cell_data) <- NULL
  row_of <- match(cell_key, cell_key[first_idx])

  Y <- matrix(0L, nrow = length(first_idx), ncol = K_cat)
  lin <- (col_idx - 1L) * nrow(Y) + row_of
  summed <- tapply(counts_in, lin, sum)
  Y[as.integer(names(summed))] <- as.integer(summed)
  colnames(Y) <- paste0("cdp", seq_len(K_cat))
  stopif(any(rowSums(Y) <= 0), "Each cell must have positive total counts")

  cell_data$Y <- Y
  cell_data$nTrials <- rowSums(Y)
  attr(cell_data, "n_new") <- n_new
  attr(cell_data, "n_old") <- n_old
  attr(cell_data, "has_guess") <- has_guess
  data <- cell_data

  NextMethod("check_data")
}


############################################################################# !
# MULTINOMIAL FORMULA & FAMILY CONSTRUCTION                              ####
############################################################################# !

# Like sdt_rating, CDP uses brms' native multinomial family: each category's
# logit is set to log(p_cat) so softmax recovers the CDP probabilities. n_new,
# n_old, threshold type, and has_guess travel as integer literals; the
# parameters and the stimulus covariate travel by name. kcrit is always passed
# (a fixed -100 sentinel when no Know/Guess split is active). spacing is the
# parameter for parsimonious/equidistant and the literal 0 for log_distance
# (which estimates per-distance deltas appended at the end instead).
.sdt_cdp_logmu_args <- function(model) {
  ov <- model$other_vars
  has_spacing <- "spacing" %in% names(model$parameters)
  c(ov$n_new, ov$n_old, ov$thresh_type_int, as.integer(ov$has_guess),
    "dprimef", "dprimer", "criterion",
    if (has_spacing) "spacing" else "0", "rcrit", "sigmar", "rho",
    "kcrit", ov$stimulus, .sdt_threshold_delta_names(model))
}

#' @export
bmf2bf.sdt_cdp <- function(model, formula) {
  K_cat <- model$other_vars$n_new +
    (if (model$other_vars$has_guess) 3L else 2L) * model$other_vars$n_old
  resp_cats <- paste0("cdp", seq_len(K_cat))
  args <- paste(.sdt_cdp_logmu_args(model), collapse = ", ")

  bform <- brms::bf(
    glue("Y | trials(nTrials) ~ sdt_cdp_logmu(1, {args})"),
    nl = TRUE
  )
  for (k in seq_len(K_cat)[-1]) {
    bform <- bform + brms::nlf(stats::as.formula(
      glue("mu{resp_cats[k]} ~ sdt_cdp_logmu({k}, {args})")
    ))
  }
  bform
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_cdp <- function(model, data, formula) {
  # n_new/n_old/has_guess are detected from the data in check_data; bridge them
  # onto the model so bmf2bf can emit the right number of categories.
  model$other_vars$n_new <- attr(data, "n_new")
  model$other_vars$n_old <- attr(data, "n_old")
  model$other_vars$has_guess <- attr(data, "has_guess")
  K_cat <- model$other_vars$n_new +
    (if (model$other_vars$has_guess) 3L else 2L) * model$other_vars$n_old
  resp_cats <- paste0("cdp", seq_len(K_cat))

  formula <- bmf2bf(model, formula)
  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- resp_cats
  formula$family$dpars <- paste0("mu", resp_cats)

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- paste(read_lines2(paste0(sc_path, "/sdt_cdp_funs.stan")),
                     .sdt_cdp_logmu_stan(model), sep = "\n")
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}

# Code-generate the Stan `sdt_cdp_logmu` category logit. The per-distance
# log_distance deltas arrive as fixed-arity `real deltaN` args (matching the
# order in .sdt_cdp_logmu_args) and are packed into a contiguous `deltas` array
# for cdp_make_thresholds. Mirrors .sdt_rating_logmu_stan.
.sdt_cdp_logmu_stan <- function(model) {
  delta_names <- .sdt_threshold_delta_names(model)
  nd <- length(delta_names)

  signature <- paste(
    c("int cat", "int n_new", "int n_old", "int thresh_type", "int has_guess",
      "real dprimef", "real dprimer", "real criterion", "real spacing",
      "real rcrit", "real sigmar", "real rho", "real kcrit", "real stimulus",
      if (nd) paste("real", delta_names)),
    collapse = ", "
  )
  delta_decl <- if (nd == 0L) {
    "  array[0] real deltas;\n"
  } else {
    paste0("  array[", nd, "] real deltas;\n",
           paste(sprintf("  deltas[%d] = %s;", seq_len(nd), delta_names),
                 collapse = "\n"), "\n")
  }

  paste0(
    "real sdt_cdp_logmu(", signature, ") {\n",
    delta_decl,
    "  vector[n_new + n_old - 1] thr = cdp_make_thresholds(criterion, spacing, ",
    "deltas, n_new, n_old, thresh_type);\n",
    "  return log(cdp_category_prob(cat, thr, dprimef, dprimer, sigmar, rho, ",
    "rcrit, kcrit, stimulus, n_new, n_old, has_guess));\n",
    "}\n"
  )
}


#' @title CDP category log-probability (multinomial logit)
#' @description R companion to the Stan `sdt_cdp_logmu` function. It returns
#'   `log(p_cat)` for response category `cat`, which the [sdt_cdp()] multinomial
#'   formula uses as the category logit (so `softmax` recovers the CDP category
#'   probabilities). `brms` evaluates the non-linear formula in R for
#'   `posterior_predict()` and `posterior_epred()`, so this function must be on
#'   the search path; it is exported for that reason and is not called directly.
#' @param cat Integer response-category index (canonical order: new, then
#'   guess/know/remember blocks).
#' @param n_new,n_old Integer numbers of "new" and "old" confidence levels.
#' @param thresh Integer threshold-parameterization id (1 parsimonious,
#'   2 equidistant, 3 log_distance).
#' @param has_guess Integer flag (1 if the Know/Guess split is active).
#' @param dprimef,dprimer,criterion,spacing,rcrit,sigmar,rho,kcrit Model
#'   parameters (draws-by-observation matrices supplied by brms). `kcrit` is
#'   ignored when `has_guess` is 0. `spacing` is the fixed literal 0 for
#'   log_distance, which uses the per-distance deltas in `...` instead.
#' @param stimulus Stimulus covariate (0 = new/lure, 1 = old/target).
#' @param ... For `log_distance`, the per-distance threshold parameters
#'   (`delta1`, ...) supplied by brms in contiguous order; otherwise unused.
#' @return `log(p_cat)`, matching the shape of `dprimef`.
#' @keywords internal
#' @export
sdt_cdp_logmu <- function(cat, n_new, n_old, thresh, has_guess,
                          dprimef, dprimer, criterion, spacing, rcrit,
                          sigmar, rho, kcrit, stimulus, ...) {
  thresh_name <- .sdt_threshold_type_name(thresh)

  shape <- dim(dprimef)
  dprimef <- as.vector(dprimef)
  n <- length(dprimef)
  dprimer <- as.vector(dprimer)
  criterion <- as.vector(criterion)
  spacing <- rep_len(as.vector(spacing), n)
  rcrit <- as.vector(rcrit)
  sigmar <- rep_len(as.vector(sigmar), n)
  rho <- rep_len(as.vector(rho), n)
  kcrit <- rep_len(as.vector(kcrit), n)
  stimulus <- rep_len(as.vector(stimulus), n)
  kc_active <- has_guess == 1L

  # log_distance passes its per-distance deltas as the trailing formula args,
  # so they arrive through `...` in the same contiguous order as the parameters.
  n_deltas <- if (thresh_name == "log_distance") n_new + n_old - 2L else 0L
  deltas <- lapply(list(...)[seq_len(n_deltas)],
                   function(d) rep_len(as.vector(d), n))

  out <- vapply(seq_len(n), function(j) {
    deltas_j <- if (n_deltas) vapply(deltas, `[`, numeric(1), j) else NULL
    thr <- .cdp_make_thresholds(criterion[j], spacing[j], n_new, n_old,
                                thresh_name, deltas_j)
    kc <- if (kc_active) kcrit[j] else NULL
    probs <- .sdt_cdp_category_probs(thr, dprimef[j], dprimer[j], sigmar[j],
                                     rcrit[j], kc, stimulus[j], n_new, n_old,
                                     "normal", rho[j])
    log(probs[cat])
  }, numeric(1))

  if (!is.null(shape)) dim(out) <- shape
  out
}
