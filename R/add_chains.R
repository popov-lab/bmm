#' @title Sample more chains for a fitted model
#' @description `add_chains()` samples more chains of a model fitted with
#'   [bmm()] and returns one fit with the chains of `fit` followed by the new
#'   ones. Use it when a fit needs more draws, for example because its
#'   effective sample sizes are too small, and the chains it has should be kept.
#' @param fit A `bmmfit` object sampled with MCMC (`algorithm = "sampling"`) by
#'   the rstan or cmdstanr backend.
#' @param chains The number of chains to add, a positive whole number.
#' @param seed The seed of the new chains. If `NULL` (the default), a seed is
#'   drawn from R's random number generator. bmm draws the default initial
#'   values from R's random number generator too, so call [set.seed()] before
#'   `add_chains()` to make the new chains reproducible, with or without
#'   `seed`. A seed that a chain of `fit` was sampled with is refused, because
#'   the new chains could repeat that chain's draws.
#' @param ... Further arguments passed to [update.bmmfit()], such as `cores`,
#'   `refresh` or `control`. `threads` works only for a fit that was sampled
#'   with threading, because a threaded Stan program differs from an
#'   unthreaded one. Arguments that would change the model, the data, the
#'   draws per chain or what is sampled (`formula.`, `newdata`, `data`,
#'   `data2`, `model`, `prior`, `sample_prior`, `iter`, `warmup`, `thin`,
#'   `algorithm`, `backend` and `normalize`) are refused: use [update()] for
#'   those and sample all chains anew. `file` and `file_compress` are refused
#'   too; save the result with [saveRDS()].
#' @return A `bmmfit` object with the chains of `fit` followed by the new ones,
#'   combined with [brms::combine_models()]. Everything else is taken from
#'   `fit`, except criteria added with [brms::add_criterion()]: they were
#'   computed from the chains of `fit` alone and are removed, so add them again.
#'   The result has no `file`, so [brms::add_criterion()] does not write it
#'   over the file of `fit`. `version` and the other fields
#'   that describe the fit are those of `fit`; the new chains ran the same Stan
#'   program. [fit_info()] reports the time, step size and inverse metric of
#'   every chain.
#' @details The new chains are sampled with [update()], which recompiles the
#'   Stan model only when it needs to, for example for a cmdstanr fit loaded in
#'   a new R session. If the Stan program that bmm now writes for the model
#'   differs from the one `fit` was sampled with, as it can after an update of
#'   bmm or brms, the new chains would not sample the same posterior, and
#'   `add_chains()` stops after sampling them. Refit the model with [bmm()] in
#'   that case.
#' @seealso [fit_info()], [update.bmmfit()]
#' @export
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' fit <- bmm(
#'   formula = bmmformula(c ~ 1, kappa ~ 1),
#'   data = data.frame(y = rsdm(1000)),
#'   model = sdm(resp_error = "y"),
#'   chains = 2
#' )
#'
#' # two more chains, four in all
#' fit <- add_chains(fit, chains = 2)
#' fit_info(fit, "time")
add_chains <- function(fit, chains = 1, seed = NULL, ...) {
  dots <- list(...)
  stopif(!inherits(fit, "bmmfit"), "'fit' must be a fit returned by bmm().")
  stopif(
    !identical(fit$algorithm, "sampling") ||
      !isTRUE(fit$backend %in% c("rstan", "cmdstanr")) ||
      !inherits(fit$fit, "stanfit") || !length(fit$fit@sim$samples),
    "add_chains() needs a fit sampled with MCMC (algorithm = 'sampling') by the \\
    rstan or cmdstanr backend, with at least one chain."
  )
  stopif(
    !is.numeric(chains) || length(chains) != 1 || !isTRUE(chains >= 1 && chains %% 1 == 0),
    "'chains' must be a positive whole number."
  )
  stopif(
    !is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || !isTRUE(seed %% 1 == 0)),
    "'seed' must be a whole number or NULL."
  )
  # an unnamed argument would reach update() as its `formula.`
  stopif(
    length(dots) > 0 && (is.null(names(dots)) || !all(nzchar(names(dots)))),
    "Arguments passed on through '...' must be named."
  )
  # update() and brms::brm() match partial names, so `new` is `newdata` there
  refused <- c(
    "formula.", "newdata", "data", "data2", "model", "prior", "sample_prior",
    "iter", "warmup", "thin", "algorithm", "backend", "normalize", "file",
    "file_compress"
  )
  refused <- names(dots)[vapply(
    names(dots), function(name) any(startsWith(refused, name)), logical(1)
  )]
  stopif(
    length(refused) > 0,
    "add_chains() does not take {collapse_comma(refused)}: the new chains must \\
    sample the same posterior as 'fit', with its model, data, prior, backend and \\
    draws per chain. Use update() to change those and sample all chains anew, and \\
    saveRDS() to save the result."
  )

  # a chain's random numbers follow from its seed and its chain id, and the
  # chains of update() count their ids from 1 again: with a seed of `fit`, a new
  # chain repeats the draws of the old chain with its id whenever their initial
  # values agree as well, as they do for a fixed `init`. A fit without a stored
  # seed (NA) passes, because there is no seed to compare with
  old_seeds <- .chain_seeds(fit)
  # a script that calls set.seed() with the same value before bmm() and before
  # add_chains() draws the seed of `fit` again
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
    while (seed %in% old_seeds) seed <- sample.int(.Machine$integer.max, 1)
  }
  repeated <- match(seed, old_seeds)
  stopif(
    !is.na(repeated),
    "Chain {repeated} of 'fit' was sampled with seed \\
    {format(seed, scientific = FALSE)}. New chains with that seed can repeat its \\
    draws; choose another seed or leave 'seed' as NULL."
  )

  new <- stats::update(fit, chains = chains, seed = seed, ...)
  # the seed a chain ran with is not always the one requested: cmdstanr runs
  # some seeds, such as 100000, as seed 1
  new_seeds <- .chain_seeds(new)
  clash <- match(new_seeds, old_seeds, incomparables = NA)
  first <- which(!is.na(clash))[1]
  stopif(
    !is.na(first),
    "New chain {first} ran with seed {format(new_seeds[first], scientific = FALSE)}, \\
    the seed of chain {clash[first]} of 'fit', and can repeat its draws. Choose another seed."
  )
  # brms recompiles when it has to, and a program written by another version of
  # bmm or brms may then differ from the one `fit` ran. The first line of the
  # code names those versions and is all that may differ
  header <- "^[^\n]*\n"
  stopif(
    !identical(sub(header, "", new$model), sub(header, "", fit$model)),
    "The new chains were sampled from a different Stan program than 'fit', which \\
    bmm {fit$version$bmm} wrote. Refit the model with bmm() and as many chains as \\
    you need."
  )
  out <- brms::combine_models(fit, new)
  out$criteria <- list()
  out$file <- NULL
  out
}

# NA for a chain without a stored seed; numeric, because as.character(100000)
# is "1e+05" and would not match a stored 100000L
.chain_seeds <- function(fit) {
  vapply(
    fit$fit@stan_args,
    function(args) suppressWarnings(as.numeric(args$seed %||% NA)),
    numeric(1)
  )
}
