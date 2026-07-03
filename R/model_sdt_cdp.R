############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# Canonical response-category columns in the order the multinomial likelihood
# expects: new(1..n_new), [guess], know, remember, with the old blocks on the
# unified confidence scale n_new+1..K. Single source of truth shared by
# check_data, aggregate_sdt_cdp_data(), and rsdt_cdp().
.sdt_cdp_response_cols <- function(n_new, n_old, has_guess, response = "") {
  old_levels <- n_new + seq_len(n_old)
  paste0(response, c(paste0("new", seq_len(n_new)),
                     if (has_guess) paste0("guess", old_levels),
                     paste0("know", old_levels),
                     paste0("remember", old_levels)))
}

.model_sdt_cdp <- function(response = NULL, stimulus = NULL,
                           n_new = NULL, n_old = NULL, dist = "normal",
                           threshold_type = "parsimonious",
                           links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)
  thresh_type_int <- .sdt_threshold_type_id(threshold_type)

  # Threshold parameters are either a single `spacing` (parsimonious/equidistant)
  # or per-distance `deltaN` (log_distance), with the anchor at the old/new
  # boundary n_new. They must be declared here so configure_prior/create_initfun
  # see them.
  thr_parts <- .sdt_threshold_parameter_parts(
    n_ratings = n_new + n_old,
    threshold_type = threshold_type,
    anchor = n_new
  )

  parameters <- c(
    list(
      dprimef = "Familiarity sensitivity: target mean on the familiarity axis",
      dprimer = "Recollection sensitivity: target mean on the recollection axis",
      criterion = "Response bias: old/new boundary on the strength (F+R) axis"
    ),
    thr_parts$parameters,
    list(
      rcrit = "Remember criterion: threshold on the recollection axis",
      sigmar = paste0(
        "Log SD of the recollection target distribution, so exp(sigmar) is the ",
        "SD and 0 means SD = 1"
      ),
      rho = paste0(
        "Familiarity-recollection correlation on an unconstrained scale; ",
        "tanh(rho) is the correlation and 0 means independent processes"
      ),
      kcrit = paste0(
        "Know criterion: threshold on the familiarity axis that splits Know from ",
        "Guess. Active only when the data include 'guess' counts"
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

  requirements <- paste0(
    "Provide aggregated data with one row per cell (unique combination of ",
    "predictors and stimulus class) and one integer count column per response ",
    "category:", "\n\n",
    "  - new1 ... new<n_new> for 'new' judgments", "\n",
    "  - know<k> and remember<k> (and optionally guess<k>) for 'old' ",
    "judgments, with k on the unified confidence scale n_new+1 ... ",
    "n_new+n_old", "\n",
    "  - an optional common column prefix is set via `response`", "\n",
    "  - stimulus: 0 = new/lure, 1 = old/target", "\n\n",
    "Use aggregate_sdt_cdp_data() to build these columns from long-format ",
    "(trial-level) data"
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
      resp_vars = nlist(response),
      other_vars = nlist(stimulus, dist, dist_int, threshold_type,
                         thresh_type_int, n_new, n_old),
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
#' **Response format.** The model is fit to aggregated counts, like
#' [sdt_rating()]: one row per cell (e.g. participant x stimulus class x
#' condition) and one integer count column per response category. The columns
#' follow a fixed naming scheme -- `new1 ... new<n_new>` for "new" judgments,
#' and `know<k>` / `remember<k>` (optionally `guess<k>`) for "old" judgments,
#' where `k` runs over the unified confidence scale `n_new + 1 ...
#' n_new + n_old` (1 = most confident "new", K = most confident "old") -- so
#' the constructor only needs the numbers of confidence levels, not a long
#' vector of column names. An optional common prefix is set via `response`
#' (e.g. `response = "cdp_"` for columns `cdp_new1`, `cdp_know2`, ...). Use
#' [aggregate_sdt_cdp_data()] to build these columns from long-format
#' (trial-level or count) data; [rsdt_cdp()] generates them directly. The
#' numbers of "new" and "old" confidence levels need not be equal (e.g. the
#' 1-new / 5-old scale of Rotello et al., 2005). Whether the Know/Guess split
#' is modelled is **driven by the data**: include `guess<k>` columns to fit
#' the three-way Remember/Know/Guess model.
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
#' @param response An optional common prefix for the response count columns
#'   (default `""`, i.e. the bare canonical names `new1`, `know2`, ... that
#'   [aggregate_sdt_cdp_data()] and [rsdt_cdp()] produce).
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (new/lure) and 1 (old/target).
#' @param n_new,n_old Integer numbers of "new" and "old" confidence levels.
#'   Together with `response` they determine the response count columns and
#'   the threshold parameters.
#' @param dist The noise distribution. Only `"normal"` is currently supported
#'   (the CDP model is inherently Gaussian).
#' @param threshold_type Character. Threshold parameterization on the strength
#'   axis: `"parsimonious"` (default) and `"equidistant"` use a single `spacing`
#'   parameter; `"log_distance"` (Meyer-Grant et al., 2025) estimates the
#'   `n_new + n_old - 2` distances between adjacent thresholds freely, each as a
#'   `deltaN` parameter on the log scale (the distance leading into threshold
#'   `N` from the old/new boundary).
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
#' # Simulate a Remember/Know data set (3 new + 3 old confidence levels) for
#' # 20 subjects: rsdt_cdp() returns the count columns sdt_cdp() expects
#' dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
#' thresholds <- c(-1.1, -0.5, 0, 0.6, 1.3)
#' dat <- cbind(dat, rsdt_cdp(nrow(dat), 200, dat$stimulus,
#'                            dprimef = 0.8, dprimer = 1.0,
#'                            thresholds = thresholds, rcrit = 0.7,
#'                            n_new = 3))
#'
#' model <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
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
sdt_cdp <- function(response = "", stimulus, n_new, n_old, dist = "normal",
                    threshold_type = c("parsimonious", "equidistant",
                                       "log_distance"),
                    links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  threshold_type <- match.arg(threshold_type)
  stopif(!identical(dist, "normal"),
         "sdt_cdp currently supports only dist = 'normal'; other noise \\
         distributions are deferred to a future release")
  stopif(n_new < 1 || n_old < 1 || n_new + n_old < 3,
         "n_new and n_old must be >= 1 and sum to >= 3")

  .model_sdt_cdp(response = response, stimulus = stimulus, n_new = n_new,
                 n_old = n_old, dist = dist, threshold_type = threshold_type,
                 links = links, call = call, ...)
}


############################################################################# !
# DATA PREPARATION                                                       ####
############################################################################# !

#' @title Aggregate long-format Remember/Know data for [sdt_cdp()]
#' @description Reshapes long-format Remember/Know recognition data -- one row
#'   per trial, or one row per response category with a count column -- into
#'   the aggregated wide format [sdt_cdp()] is fit to: one row per cell and
#'   one count column per response category, named `new1 ... new<n_new>`,
#'   `know<k>` / `remember<k>` (and `guess<k>` when the data contain "guess"
#'   judgments), with `k` on the unified confidence scale
#'   `n_new + 1 ... n_new + n_old`.
#' @param data A data frame in long format. Cells are defined by all columns
#'   other than `judgment`, `confidence`, and `count` (e.g. participant,
#'   stimulus class, and any condition variables), which are carried over to
#'   the output unchanged.
#' @param judgment The name of the column coding the memory judgment, with
#'   values `"new"`, `"remember"`, `"know"`, and optionally `"guess"`.
#' @param confidence The name of the column coding confidence on the unified
#'   old/new scale (integer, 1 = most confident "new" to K = most confident
#'   "old"; "new" judgments occupy the low levels, "old" judgments the high
#'   levels).
#' @param count Optional name of a column of response counts. If `NULL`
#'   (default), each row is treated as a single trial.
#' @param response Optional common prefix for the generated count columns
#'   (default `""`). Pass the same value to the `response` argument of
#'   [sdt_cdp()].
#' @return A data frame with one row per cell: the cell-defining columns
#'   followed by the response count columns, ready to pass to [bmm()] with an
#'   [sdt_cdp()] model. The numbers of "new" and "old" confidence levels are
#'   inferred from the data.
#' @examples
#' dat <- data.frame(
#'   stimulus = c(1L, 1L, 1L, 0L, 0L, 0L),
#'   judgment = c("remember", "know", "new", "new", "know", "new"),
#'   confidence = c(3L, 2L, 1L, 1L, 3L, 1L),
#'   count = c(40L, 25L, 15L, 55L, 10L, 20L)
#' )
#' aggregate_sdt_cdp_data(dat, "judgment", "confidence", "count")
#' @keywords transform
#' @export
aggregate_sdt_cdp_data <- function(data, judgment, confidence, count = NULL,
                                   response = "") {
  stop_missing_args()
  stopif(!judgment %in% colnames(data),
         "Judgment variable '{judgment}' missing in the data")
  stopif(!confidence %in% colnames(data),
         "Confidence variable '{confidence}' missing in the data")
  stopif(!is.null(count) && !count %in% colnames(data),
         "Count variable '{count}' missing in the data")

  judg <- as.character(data[[judgment]])
  allowed <- c("new", "remember", "know", "guess")
  bad <- setdiff(unique(judg), allowed)
  stopif(length(bad) > 0,
         "Judgment values {collapse_comma(bad)} not allowed; use {collapse_comma(allowed)}")
  stopif(!all(c("new", "remember", "know") %in% judg),
         "Judgment column must contain 'new', 'remember', and 'know' values")
  has_guess <- "guess" %in% judg

  conf <- data[[confidence]]
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
  stopif(!all(conf[is_new] %in% seq_len(n_new)),
         "'new' confidence values must occupy levels 1..{n_new}")
  stopif(!all(conf[!is_new] %in% (n_new + 1L):K_full),
         "old (remember/know/guess) confidence values must occupy levels {n_new + 1L}..{K_full}")

  counts_in <- if (!is.null(count)) {
    vals <- data[[count]]
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

  cell_cols <- setdiff(colnames(data), c(judgment, confidence, count))
  cell_key <- do.call(paste, c(data[cell_cols], sep = "\r"))
  first_idx <- which(!duplicated(cell_key))
  cells <- data[first_idx, cell_cols, drop = FALSE]
  rownames(cells) <- NULL
  row_of <- match(cell_key, cell_key[first_idx])

  Y <- matrix(0L, nrow = length(first_idx), ncol = K_cat)
  lin <- (col_idx - 1L) * nrow(Y) + row_of
  summed <- tapply(counts_in, lin, sum)
  Y[as.integer(names(summed))] <- as.integer(summed)
  colnames(Y) <- .sdt_cdp_response_cols(n_new, n_old, has_guess, response)
  cbind(cells, as.data.frame(Y))
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_cdp <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stopif(!all(unique(data[[stim_var]]) %in% c(0, 1)),
         "Stimulus variable '{stim_var}' must be coded as 0 (new/lure) and 1 (old/target)")

  n_new <- model$other_vars$n_new
  n_old <- model$other_vars$n_old
  prefix <- model$resp_vars$response
  guess_cols <- paste0(prefix, "guess", n_new + seq_len(n_old))
  n_guess_found <- sum(guess_cols %in% colnames(data))
  stopif(n_guess_found > 0 && n_guess_found < n_old,
         "Found only some of the 'guess' count columns \\
         ({collapse_comma(guess_cols)}); provide all of them or none")
  has_guess <- n_guess_found == n_old

  resp_cols <- .sdt_cdp_response_cols(n_new, n_old, has_guess, prefix)
  missing <- setdiff(resp_cols, colnames(data))
  stopif(length(missing) > 0,
         "Response columns {collapse_comma(missing)} missing in the data. \\
         Use aggregate_sdt_cdp_data() to build them from long-format data")

  for (col in resp_cols) {
    vals <- data[[col]]
    stopif(any(vals < 0, na.rm = TRUE),
           "Response column '{col}' must contain non-negative counts")
    warnif(any(vals != round(vals), na.rm = TRUE),
           "Response column '{col}' should contain integer counts")
  }

  Y <- as.matrix(data[resp_cols])
  colnames(Y) <- paste0("cdp", seq_len(ncol(Y)))
  stopif(any(rowSums(Y) <= 0, na.rm = TRUE),
         "Row sums of response columns must be positive (no empty rows)")

  data <- data[!colnames(data) %in% resp_cols]
  data$Y <- Y
  data$nTrials <- rowSums(Y)
  attr(data, "has_guess") <- has_guess

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
  # whether the Know/Guess split is active is detected from the data in
  # check_data; bridge it onto the model so bmf2bf can emit the right number
  # of categories.
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
#'   `log(p_cat)` for response category `cat` (unnormalized, exactly like the
#'   Stan function -- the category probabilities sum to 1 analytically and
#'   `softmax` absorbs the shared constant), which the [sdt_cdp()] multinomial
#'   formula uses as the category logit. `brms` evaluates the non-linear
#'   formula in R for `posterior_predict()` and `posterior_epred()`, so this
#'   function must be on the search path; it is exported for that reason and is
#'   not called directly. Vectorized over the draws-by-observations matrices
#'   brms supplies.
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
  n <- length(dprimef)

  # log_distance passes its per-distance deltas as the trailing formula args,
  # so they arrive through `...` in the same contiguous order as the parameters.
  deltas <- if (thresh_name == "log_distance") {
    do.call(cbind, lapply(list(...)[seq_len(n_new + n_old - 2L)],
                          function(d) rep_len(as.vector(d), n)))
  }
  thr <- rbind(.cdp_make_thresholds(as.vector(criterion),
                                    rep_len(as.vector(spacing), n),
                                    n_new, n_old, thresh_name, deltas))

  out <- log(.sdt_cdp_category_prob(
    cat, thr, as.vector(dprimef), as.vector(dprimer), as.vector(sigmar),
    as.vector(rho), as.vector(rcrit), as.vector(kcrit), as.vector(stimulus),
    n_new, n_old, has_guess == 1L
  ))
  if (!is.null(shape)) dim(out) <- shape
  out
}
