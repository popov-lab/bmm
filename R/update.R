#' @title Update a bmm model
#' @description Update an existing bmm mode. This function calls
#'   [brms::update.brmsfit()], but it applies the necessary bmm postprocessing
#'   to the model object before and after the update.
#' @param object An object of class `bmmfit`
#' @param formula. A [bmmformula()]. If missing, the original formula
#'  is used. Currently you have to specify a full `bmmformula`
#' @param newdata An optional data frame containing the variables in the model.
#'   Rows of `fit$data` are accepted as they are, which is how [brms::kfold()]
#'   and [brms::reloo()] pass them.
#' @param recompile Logical, indicating whether the Stan model should be recompiled. If
#'   NULL (the default), update tries to figure out internally, if recompilation
#'   is necessary. Setting it to FALSE will cause all Stan code changing
#'   arguments to be ignored.
#' @param file Either `NULL` or a character string. If a string, the updated
#'   model is saved via [saveRDS] in a file named after the string, as in
#'   [bmm()]. `update()` never writes to the file the original fit was read
#'   from: pass `file` explicitly to save the updated fit.
#' @param file_compress Logical or a character string, specifying one of the
#'   compression algorithms supported by [saveRDS] when saving the updated
#'   model object.
#' @param ... Further arguments passed to [brms::update.brmsfit()]
#' @return An updated `bmmfit` object refit to the new data and/or formula. If
#'   `file` is given, it names the newly written file. If not, the `file` field
#'   is carried over from the original fit, so the updated object still points
#'   at the file it came from -- but that file is *not* rewritten: it still
#'   holds the fit as it was before this update. Pass `file` to save the
#'   updated fit.
#' @details When updating a brmsfit created with the cmdstanr backend in a
#'   different R session, a recompilation will be triggered because by default,
#'   cmdstanr writes the model executable to a temporary directory. To avoid
#'   that, set option "cmdstanr_write_stan_file_dir" to a nontemporary path of
#'   your choice before creating the original bmmfit.
#'
#'   For more information and examples, see [brms::update.brmsfit()]
#' @export
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' # generate artificial data from the Signal Discrimination Model
#' # generate artificial data from the Signal Discrimination Model
#' dat <- data.frame(y = rsdm(2000))
#'
#' # define formula
#' ff <- bmf(c ~ 1, kappa ~ 1)
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = sdm(resp_error = "y"),
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # update the model
#' fit <- update(fit, newdata = data.frame(y = rsdm(2000, kappa = 5)))
#'
update.bmmfit <- function(object, formula., newdata = NULL, recompile = NULL,
                          file = NULL, file_compress = TRUE, ...) {
  dots <- list(...)
  save_file <- check_rds_file(file)
  save_compress <- file_compress
  # NextMethod() forwards this method's formals, so `file` would still reach
  # brms::brm() even though it is not in `dots`: brms would then return the
  # contents of an existing file instead of updating, and write a plain brmsfit
  # before any bmm postprocessing has run. Naming file = NULL in the
  # NextMethod() call errors with "matched by multiple actual arguments", so the
  # formals are blanked here instead
  file <- NULL
  file_compress <- NULL
  # brms::update.brmsfit falls back to the original fit's threading spec only
  # when `threads` is absent from the call -- an explicit NULL means "no
  # threading" -- and the effective spec, not the new request, must drive the
  # option that configure_model reads. The threading(NULL) tail pins the option
  # for fits saved without a `threads` field, so a stray global brms.threads
  # cannot slice a likelihood that brms will run serially
  effective_threads <- if ("threads" %in% names(dots)) dots$threads else object$threads
  local_brms_threads(list(threads = effective_threads %||% brms::threading(NULL)))
  stopif(
    isTRUE(object$version$bmm < "0.3.0"),
    "Updating bmm models works only with models fitted with version 0.3.0 or higher"
  )
  stopif(
    "data" %in% names(dots),
    "Please use argument 'newdata' to update the data."
  )
  stopif(
    "model" %in% names(dots),
    "You cannot update with a different model. Create a new fit with 'bmm()' instead."
  )

  object <- restructure(object)

  model <- object$bmm$model
  old_user_formula <- object$bmm$user_formula
  olddata <- object$data
  old_file <- object$file
  configure_opts <- object$bmm$configure_opts

  # revert some postprocessing changes to brmsfit from postprocess_brm
  object <- revert_postprocess_brm(model, object)

  # use the new configure_opts if they are provided
  if (any(names(dots) %in% names(configure_opts))) {
    new_opts <- names(dots)[names(dots) %in% names(configure_opts)]
    configure_opts[new_opts] <- dots[new_opts]
  }
  opts <- configure_options(configure_opts)

  # reuse or replace formula and data
  if (missing(formula.)) {
    user_formula <- old_user_formula
  } else {
    user_formula <- formula.
  }

  # bmm() resolves constants in check_model() before check_data(); update() never
  # calls check_model(), so without this any change the new formula makes to a
  # parameter's constant is silently ignored: a freed parameter stays pinned by
  # the old constant(), a newly fixed one keeps its old non-constant prior, and a
  # re-valued one keeps the old constant. Comparing the resolved constants rather
  # than the formula covers all three
  old_fixed <- model$fixed_parameters
  model <- update_model_fixed_parameters(model, user_formula)
  # the stored model was checked with the old formula and keeps what that
  # formula derived (mpt: the links switched off for non-linear parameters and
  # their sub-parameters), so a new formula needs a new check
  if (!missing(formula.)) {
    model <- check_model(model, newdata %||% olddata, user_formula)
  }
  changed_pars <- union(names(old_fixed), names(model$fixed_parameters))
  changed_pars <- changed_pars[!vapply(
    changed_pars,
    function(par) identical(old_fixed[[par]], model$fixed_parameters[[par]]),
    logical(1)
  )]

  if (is.null(newdata) || is_stored_frame(model, olddata, newdata)) {
    data <- check_stored_data(model, newdata %||% olddata, user_formula)
    from_stored_frame <- !missing(formula.)
    attr(data, "data_name") <- attr(olddata, "data_name")
  } else {
    data <- check_data(model, newdata, user_formula)
    from_stored_frame <- FALSE
    attr(data, "data_name") <- substitute_name(newdata)
  }

  # standard bmm checks and transformations
  formula <- check_formula(model, data, user_formula)
  config_args <- configure_model(model, data, formula)
  # checked on the configured formula, because only there does brms say which
  # of its symbols are variables
  if (from_stored_frame) {
    stop_unstored_variables(config_args, fit_frame_args(object, dots)$data2)
  }

  # bmm() stores the prior the user gave it, and configure_prior() rebuilds the
  # defaults around it for the new formula and data as bmm() would (#538). A fit
  # saved before that has only its whole prior, in which brms marks bmm's
  # defaults for the old formula as user rows; configure_prior() treats each of
  # those rows as a user prior, so for a parameter whose constant changed each
  # would override the freshly configured one -- the stale constant() of a freed
  # or re-valued parameter and the stale free prior of a newly fixed one alike.
  # brms stores the main dpar without a `dpar` label, so its rows are the ones
  # carrying neither label
  stored_prior <- object$bmm$user_prior
  old_prior <- stored_prior %||% carried_prior(object$prior)
  if (is.null(stored_prior) && length(changed_pars) > 0) {
    main_dpar <- names(brms::brmsterms(config_args$formula)$dpars)[1]
    is_main_dpar_row <- !nzchar(old_prior$dpar) & !nzchar(old_prior$nlpar)
    stale <- old_prior$dpar %in% changed_pars |
      old_prior$nlpar %in% changed_pars |
      (isTRUE(main_dpar %in% changed_pars) &
         old_prior$class %in% c("Intercept", "b") & is_main_dpar_row)
    old_prior <- old_prior[!stale, ]
  }
  # the prior of this call joins the carried rows before configure_prior() runs,
  # not after: the constraints a model sets in its own configure_prior method
  # must win over it, and the message naming the displaced rows is printed
  # there, exactly as for the prior of bmm() (#549). Combined afterwards, a
  # prior on a constrained coefficient would replace the constraint in silence
  user_prior <- combine_prior(old_prior, dots$prior)
  prior <- brms::do_call(
    configure_prior,
    c(list(model, data, config_args$formula, user_prior), fit_frame_args(object, dots))
  )
  dots$prior <- NULL
  new_fit_args <- combine_args(nlist(config_args, dots, prior))

  # configure_model() always returns the complete brmsformula, so brms has
  # nothing to merge. Handed over as `formula.`, it would be rebuilt by
  # update.brmsformula() in another element order than bmm() stores, and
  # brms::combine_models() would then reject the updated fit as having a
  # different formula (#464). As the stored formula, it goes through the same
  # validation in brms as the one bmm() passes to brm(). `formula.` is set to
  # NULL because NextMethod() forwards the frame's value, which would otherwise
  # be the user's bmmformula
  object$formula <- new_fit_args$formula
  formula. <- NULL
  if (!identical(new_fit_args$data, olddata)) {
    newdata <- new_fit_args$data
  }

  # the fit's stored init closure captured the Stan data of the original fit,
  # so a new formula or data needs a new one, built from the prior as brms will
  # read it: brms::update.brmsfit() tags the prior so that rows of the old fit
  # that the new formula or data leave without a parameter are dropped instead
  # of rejected. Built before the NextMethod() call, because an error inside a
  # lazy argument of that call surfaces as "promise already under evaluation"
  # rather than as itself
  attr(prior, "allow_invalid_prior") <- TRUE
  # brms::brm() samples with the deprecated `inits` whenever it is given
  init <- dots$inits %||% if ("init" %in% names(dots)) dots$init else brms::do_call(
    create_initfun, c(list(model, data, config_args$formula, prior), fit_frame_args(object, dots))
  )

  # brms::update.brmsfit() resolves abbreviations such as algorithm = "samp"
  # with match.arg() before it decides whether the run is the fit's, so the
  # control carried over and stored here must be decided on the resolved value,
  # and brms is handed that value. The choices are those of brms's internal
  # algorithm_choices() and backend_choices(), which it does not export
  dots$backend <- match.arg(dots$backend %||% object$backend %||% "rstan",
                            c("rstan", "cmdstanr", "mock"))
  dots$algorithm <- match.arg(
    dots$algorithm %||% object$algorithm %||% "sampling",
    c("sampling", "meanfield", "fullrank", "pathfinder", "laplace", "fixed_param")
  )
  control <- sampler_control(object, dots, configure_control(
    carried_control(object, dots), dots$backend, dots$algorithm
  ))

  # brms::update.brmsfit() adds back every row of the old prior marked "user"
  # that the new prior does not restate, and brms marks bmm's defaults so: the
  # old default of a coefficient the new formula keeps would replace the new
  # one, and the stored L row of a correlation would duplicate the new cor row.
  # The user's own rows are in the stored prior already (#538)
  if (!is.null(stored_prior)) {
    object$prior$source[] <- "default"
  }

  # pass back to brms::update.brmsfit; stanvars must be the freshly configured
  # ones — brms otherwise reuses object$stanvars, whose data values (e.g. the
  # sdm run metadata) were computed for the original data and formula. The
  # named `control` replaces the one in the dots and adds the starting step size
  object <- NextMethod("update", object,
    newdata = newdata,
    prior = prior, recompile = recompile,
    stanvars = new_fit_args$stanvars, init = init, control = control,
    backend = dots$backend, algorithm = dots$algorithm, ...
  )
  # without a recompile brms samples with `init` and `control` but keeps the old
  # fit's stan_args, and a chained update() would carry the control from before
  # this one. Single brackets, because `init = NULL` must keep its slot
  object$stan_args[c("init", "control")] <- list(init, control)

  # bmm postprocessing. A fit saved before bmm stored the user's prior stays
  # without one: its user rows cannot be told from the defaults, and storing the
  # prior of this call alone would drop them at the next update()
  object <- postprocess_brm(model, object,
    fit_args = new_fit_args, user_formula = user_formula,
    user_prior = if (!is.null(stored_prior)) user_prior,
    configure_opts = configure_opts
  )

  # saving here, rather than letting brms write the file from inside
  # brms::update.brmsfit(), is the same order bmm() uses: the object on disk is
  # the postprocessed bmmfit that try_read_bmmfit() can load back. Without
  # `file`, the field is carried over so the fit keeps naming the file it came
  # from -- that file still holds the fit before the update
  if (is.null(save_file)) {
    object$file <- old_file
    return(object)
  }
  try_save_bmmfit(object, save_file, compress = save_compress)
}

# `object$data` is brms's model frame, not the data the user passed to bmm():
# check_data() has already run on it once, and brms keeps only the variables the
# fitted formula references. Running check_data() on it again therefore fails
# for any model whose check_data() consumes or generates columns, so the frame
# is first turned back into something check_data() accepts. A column that can
# only be rebuilt, not recovered, is named in the "rebuilt" attribute and
# dropped again afterwards: it stands in for the user's column while the checks
# run and must not reach the model frame, where a new formula could pick it up
# as a predictor of the wrong type. Dropping it also leaves the column NULL for
# the configure_prior methods of the non-target models, which read it to decide
# whether set size 1 needs a constant prior. That is inert only because a column
# gets rebuilt exactly when the fitted formula does not name it, and the
# constraint helpers then return NULL; a new formula. that names it is refused
# by stop_unstored_variables(), because the column is no longer in the data
check_stored_data <- function(model, data, formula) {
  stored <- revert_check_data(model, data)
  rebuilt <- attr(stored, "rebuilt")
  attr(stored, "rebuilt") <- NULL
  data <- check_data(model, stored, formula)
  for (var in rebuilt) {
    data[[var]] <- NULL
  }
  data
}

# A new formula can only use the columns the fit's model frame kept, which are
# those its own formula used, and none of the ones check_stored_data() rebuilt
# and dropped again. brms refuses the others only after bmm has configured the
# prior, with a message that does not say where the data come from
stop_unstored_variables <- function(config_args, data2) {
  # the variables brms::validate_data() looks up, by its own rule: allvars
  # leaves out the correlation ID p of (1 | p | id), cov matrices and non-linear
  # parameters, and Intercept is a column brms adds itself
  unstored <- setdiff(
    all.vars(brms::brmsterms(config_args$formula)$allvars),
    c(colnames(config_args$data), names(data2), "Intercept")
  )
  stopif(
    length(unstored) > 0,
    "The model frame of this fit does not hold {collapse_comma(unstored)}. It keeps \\
    only the columns the fitted formula used, and none that bmm rebuilds for its \\
    checks, such as the tree column of `mpt()`. Pass the data you gave `bmm()` as \\
    `newdata` to use {if (length(unstored) > 1) 'them' else 'it'} in the new formula."
  )
}

# brms::kfold() and reloo() hand update() slices of the fit's own model frame,
# not data in the layout bmm() was given. Data are such a slice when
# revert_check_data() has something to change in the stored frame, they carry
# every column of it (the generated ones included) and lack the ones it
# restores. A user's own column Y is not enough: the other generated columns
# are missing. A fit that revert_check_data() refuses to rebuild is updated
# with the user's data (its error says so), so a refusal means no model frame
is_stored_frame <- function(model, stored, data) {
  reverted <- tryCatch(revert_check_data(model, stored), error = function(e) NULL)
  if (is.null(reverted) || setequal(names(reverted), names(stored))) {
    return(FALSE)
  }
  restored <- setdiff(names(reverted), names(stored))
  all(names(stored) %in% names(data)) && !any(restored %in% names(data))
}

# Dispatch runs general to specific, so every method chains with NextMethod():
# without it a method on a domain class silently shadows one on a model class
# below it, which is how the sdt_yn symptom would come back the moment the SDT
# stack adds a shared method on `sdt`. Because methods chain, a method that
# rebuilds a column appends its name to the "rebuilt" attribute rather than
# assigning it -- an assignment would erase a name an earlier method in the
# chain recorded, and that column would then reach the model frame
revert_check_data <- function(model, data) {
  UseMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.default <- function(model, data) {
  data
}

#' @exportS3Method
revert_check_data.m3 <- function(model, data) {
  resp_cats <- model$resp_vars$resp_cats
  num_options <- m3_num_options(model)
  data[resp_cats] <- as.data.frame(data$Y[, resp_cats, drop = FALSE])
  if (is.numeric(num_options)) {
    for (var in names(num_options)) {
      data[[var]] <- NULL
    }
  } else {
    # check_data() turned a zero option count into 0.0001, so the Idx_ columns
    # are the only record left of which category a row offered no options for.
    # A fit from before #457 computed them from another category's column, and
    # zeroing by them would rebuild a fit that is wrong in a new way
    no_options <- as.matrix(data[paste0("Idx_", resp_cats)]) == 0
    stopif(
      any(no_options & as.matrix(data[num_options]) != 0.0001),
      "The stored data of this fit pairs the option columns with the wrong response \\
      categories (fitted before the fix for #457 with `num_options` named after the \\
      categories in another order), so `update()` cannot rebuild it. Refit with `bmm()` \\
      on your original data, or pass `newdata`."
    )
    data[num_options][no_options] <- 0
  }
  # check_data() refuses these as user columns, but a category named Y or
  # nTrials was just restored and must stay
  data[c(setdiff(c("Y", "nTrials"), resp_cats), paste0("Idx_", resp_cats))] <- NULL
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.mpt <- function(model, data) {
  resp_cats <- model$resp_vars$resp_cats
  data[resp_cats] <- as.data.frame(unclass(data$Y)[, resp_cats, drop = FALSE])
  data$Y <- NULL
  data$nTrials <- NULL
  tree_id <- model$other_vars$tree_id
  # brms keeps the tree column only when a formula names it; otherwise the
  # one-hot tree indicators are the record of each row's tree
  if (!is.null(tree_id) && !tree_id %in% colnames(data)) {
    idx_vars <- model$other_vars$indicators$tree
    data[[tree_id]] <- names(idx_vars)[
      max.col(as.matrix(data[unname(idx_vars)]), ties.method = "first")
    ]
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), tree_id)
  }
  # check_data() refuses generated indicator columns it finds in the data, and
  # rebuilds them from the tree column
  data[unname(c(model$other_vars$indicators$tree, model$other_vars$indicators$possible))] <- NULL
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.non_targets <- function(model, data) {
  set_size <- model$other_vars$set_size
  # brms keeps the set_size column only when a formula predicts something with
  # it. A circmix family reads the set size from vint(ss_numeric), so the frame
  # holds it. Fits made on brms mixtures hold only LureIdx1..n, the step function
  # check_data() built from it, whose row sums give the set size back. Those
  # names come from nt_features rather than from a pattern match, because a user
  # column called LureIdx9 that a formula kept would be summed in too and the set
  # sizes it shifts are legal integers that nothing downstream rejects
  if (is.character(set_size) && not_in(set_size, colnames(data))) {
    data[[set_size]] <- data[["ss_numeric"]] %||% (1 + rowSums(
      data[paste0("LureIdx", seq_along(model$other_vars$nt_features))]
    ))
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), set_size)
  }
  # brms may have dropped the rows of the largest set size from the frame (#459)
  data <- carry_fit_bounds(model, data)
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.mixture2p <- function(model, data) {
  set_size <- model$other_vars$set_size
  # brms keeps the set_size column only when a formula predicts something with
  # it; ss_numeric is check_data()'s numeric copy, kept for the vint() term
  if (is.character(set_size) && not_in(set_size, colnames(data)) &&
    "ss_numeric" %in% colnames(data)) {
    data[[set_size]] <- data$ss_numeric
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), set_size)
  }
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.sdt_yn <- function(model, data) {
  data$dist_type <- NULL
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.sdt_mafc <- function(model, data) {
  m <- model$other_vars$m
  # brms keeps a set-size column named by `m` only when a formula predicts
  # something with it; m_afc is check_data()'s integer copy of that column
  if (is.character(m) && not_in(m, colnames(data))) {
    data[[m]] <- data$m_afc
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), m)
  }
  data$m_afc <- NULL
  data$dist_type <- NULL
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.sdt_ranking <- function(model, data) {
  resp_cols <- model$resp_vars$response
  m <- model$other_vars$m
  data[resp_cols] <- as.data.frame(unclass(data$Y)[, resp_cols, drop = FALSE])
  # max_rank is check_data()'s numeric copy of the set-size column, which brms
  # keeps only when a formula predicts something with it
  if (is.character(m) && not_in(m, colnames(data))) {
    data[[m]] <- data$max_rank
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), m)
  }
  data$Y <- NULL
  data$nTrials <- NULL
  data$max_rank <- NULL
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.sdt_rating <- function(model, data) {
  resp_cols <- model$resp_vars$response
  data[resp_cols] <- as.data.frame(unclass(data$Y)[, resp_cols, drop = FALSE])
  data$Y <- NULL
  data$nTrials <- NULL
  NextMethod("revert_check_data")
}

#' @exportS3Method
revert_check_data.sdt_cdp <- function(model, data) {
  n_new <- model$other_vars$n_new
  n_old <- model$other_vars$n_old
  # check_data() renamed the count columns cdp1 ... cdpK, and whether guess
  # columns were present is recorded only in how many there are
  Y <- unclass(data$Y)
  has_guess <- ncol(Y) == n_new + 3L * n_old
  resp_cols <- .sdt_cdp_response_cols(n_new, n_old, has_guess,
                                      model$resp_vars$response)
  data[resp_cols] <- as.data.frame(Y, col.names = resp_cols)
  data$Y <- NULL
  data$nTrials <- NULL
  NextMethod("revert_check_data")
}

# Rows brms filled in itself are left out, as brms::update.brmsfit() leaves them
# out: a flat default for the class of an intercept-only parameter would
# otherwise override bmm's prior on the predictors a new formula adds. A row of
# unknown source (NA, or no `source` column) cannot be told from a default and
# is left out too, as brms does; an NA in a logical index yields an all-NA row
carried_prior <- function(prior) {
  prior[which(prior$source != "default"), ]
}

# On the same backend and algorithm brms::update.brmsfit() carries a fit's
# stored control only when the call names none, and update.bmmfit() always names
# one, so the stored control is merged here, key by key under the call's. It is
# the control the fit was sampled with, also after an update() that did not
# recompile, because update.bmmfit() writes it back.
# step_size and stepsize are one argument under two spellings, so a call naming
# either replaces whichever the fit stored; merged by name they would survive as
# two keys and configure_control() would then resolve the duplicate in the fit's
# favour
carried_control <- function(object, dots) {
  if (!is_same_run(object, dots)) {
    return(dots$control)
  }
  stored <- object$stan_args$control %||% list()
  if (any(names(dots$control) %in% c("step_size", "stepsize"))) {
    stored <- stored[not_in(names(stored), c("step_size", "stepsize"))]
  }
  utils::modifyList(stored, dots$control %||% list())
}

# On the same run brms::update.brmsfit() fills the keys the control leaves out
# from the old rstan sampler's own arguments (cmdstanr fits keep none there),
# after configure_control() has set the step size. Filled in here instead, brms
# has nothing left to add, so the control update.bmmfit() stores is the one it
# hands the sampler, as brms stores it when it recompiles. update.bmmfit()
# resolves abbreviated backends and algorithms first, as brms compares the
# resolved values
sampler_control <- function(object, dots, control) {
  if (!is_same_run(object, dots)) {
    return(control)
  }
  sampler <- attr(object$fit@sim$samples[[1]], "args")$control
  c(control, sampler[not_in(names(sampler), names(control))])
}

# A fit without a backend or algorithm field counts as changed, because brms
# resolves the missing field to its first choice and then finds it different
is_same_run <- function(object, dots) {
  !is.null(object$backend) && !is.null(object$algorithm) &&
    identical(dots$backend %||% object$backend, object$backend) &&
    identical(dots$algorithm %||% object$algorithm, object$algorithm)
}
