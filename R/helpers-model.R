#' @title Generic S3 method for configuring the model to be fit by brms
#' @description Called by bmm() to automatically construct the model
#'   formula, family objects and default priors for the model specified by the
#'   user. It will call the appropriate configure_model.* functions based on the
#'   list of classes defined in the .model_* functions. Currently, we have a
#'   method only for the last class listed in the .model_* functions. This is to
#'   keep model configuration as simple as possible. In the future we may add
#'   shared methods for classes of models that share the same configuration.
#' @param model A model list object returned from check_model()
#' @param data The user supplied data.frame containing the data to be checked
#' @param formula The user supplied formula
#' @return A named list containing at minimum the following elements:
#'
#'  - formula: An object of class `brmsformula`. The constructed model formula
#'  - data: the user supplied data.frame, preprocessed by check_data
#'  - family: the brms family object
#'  - prior: the brms prior object
#'  - stanvars: (optional) An object of class `stanvars` (for custom families).
#'   See [brms::custom_family()] for more details.
#'
#' @details A bare bones configure_model.* method should look like this:
#'
#'  ``` r
#'  configure_model.newmodel <- function(model, data, formula) {
#'
#'     # preprocessing - e.g. extract arguments from data check, construct new variables
#'     <preprocessing code>
#'
#'     # construct the formula
#'     formula <- bmf2bf(formula, model)
#'
#'     # construct the family
#'     family <- <code for new family>
#'
#'     # construct the default prior
#'     prior <- <code for new prior>
#'
#'     # return the list
#'     nlist(formula, data, family, prior)
#'  }
#'  ```
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' configure_model.mixture3p <- function(model, data, formula) {
#'   # retrieve arguments from the data check
#'   max_set_size <- attr(data, "max_set_size")
#'   lure_idx <- attr(data, "lure_idx_vars")
#'   nt_features <- model$other_vars$nt_features
#'   set_size_var <- model$other_vars$set_size
#'
#'   # construct initial brms formula
#'   formula <- bmf2bf(model, formula) +
#'     brms::lf(kappa2 ~ 1) +
#'     brms::lf(mu2 ~ 1) +
#'     brms::nlf(theta1 ~ thetat) +
#'     brms::nlf(kappa1 ~ kappa)
#'
#'   # additional internal terms for the mixture model formula
#'   kappa_nts <- paste0("kappa", 3:(max_set_size + 1))
#'   theta_nts <- paste0("theta", 3:(max_set_size + 1))
#'   mu_nts <- paste0("mu", 3:(max_set_size + 1))
#'
#'   for (i in 1:(max_set_size - 1)) {
#'     formula <- formula +
#'       glue_nlf("{kappa_nts[i]} ~ kappa") +
#'       glue_nlf(
#'         "{theta_nts[i]} ~ {lure_idx[i]} * (thetant + log(inv_ss)) + ",
#'         "(1 - {lure_idx[i]}) * (-100)"
#'       ) +
#'       glue_nlf("{mu_nts[i]} ~ {nt_features[i]}")
#'   }
#'
#'   # define mixture family
#'   vm_list <- lapply(1:(max_set_size + 1), function(x) brms::von_mises(link = "identity"))
#'   vm_list$order <- "none"
#'   formula$family <- brms::do_call(brms::mixture, vm_list)
#'
#'   nlist(formula, data)
#' }
#'
#' @export
#' @keywords internal developer
configure_model <- function(model, data, formula) {
  UseMethod("configure_model")
}

############################################################################# !
# CHECK_MODEL methods                                                    ####
############################################################################# !

#' Generic S3 method for checking if the model is supported and model preprocessing
#'
#' In addition for validating the model, specific methods might add information
#' to the model object based on the provided data and formula
#'
#' @param model the model argument supplied by the user
#' @param data the data argument supplied by the user
#' @param formula the formula argument supplied by the user
#'
#' @return An object of type 'bmmodel'
#' @keywords internal developer
check_model <- function(model, data = NULL, formula = NULL) {
  UseMethod("check_model")
}

#' @export
check_model.default <- function(model, data = NULL, formula = NULL) {
  bmm_models <- supported_models(print_call = FALSE)
  if (is.function(model)) {
    fun_name <- as.character(substitute(model))
    stopif(
      fun_name %in% bmm_models,
      "Did you forget to provide the required arguments to the model function?
      See ?{fun_name} for details on properly specifying the model argument"
    )
  }

  stopif(
    !is_supported_bmmodel(model),
    "You provided an object of class `{class(model)}` to the model argument.
    The model argument should be a `bmmodel` function.
    You can see the list of supported models by running `supported_models()`

    {supported_models()}"
  )
  model
}

#' @export
check_model.bmmodel <- function(model, data = NULL, formula = NULL) {
  model <- replace_regex_variables(model, data)
  model <- update_model_fixed_parameters(model, formula)
  NextMethod("check_model")
}

# check if the user has provided a regular expression for any model variables and
# replace the regular expression with the actual variables
replace_regex_variables <- function(model, data) {
  regex <- isTRUE(attr(model, "regex"))
  regex_vars <- attr(model, "regex_vars")

  # check if the regex transformation has already been applied (e.g., if
  # updating a previously fit model)
  regex_applied <- isTRUE(attr(model, "regex_applied"))
  if (regex_applied || !regex || length(regex_vars) == 0) {
    return(model)
  }

  data_cols <- names(data)
  # save original user-provided variables
  user_vars <- c(model$resp_vars, model$other_vars)
  attr(model, "user_vars") <- user_vars

  for (var in regex_vars) {
    var_type <- if (var %in% names(model$other_vars)) "other_vars" else "resp_vars"
    model[[var_type]][[var]] <- get_variables(model[[var_type]][[var]], data_cols, regex)
  }

  attr(model, "regex_applied") <- regex
  model
}

# if the user has provided a constant in the bmmformula, add that info to the
# model object; if they have predicted a parameter that is constant by default,
# remove it from the model object
update_model_fixed_parameters <- function(model, formula) {
  constants <- names(formula)[is_constant(formula)]
  free <- names(formula)[!is_constant(formula)]
  # add new constants to the model object
  if (length(constants) > 0) {
    model$fixed_parameters[constants] <- strip_attributes(formula[constants],
      protect = "names",
      recursive = TRUE
    )
  }
  overwrite <- intersect(names(model$fixed_parameters), free)
  if (length(overwrite) > 0) {
    model$fixed_parameters[overwrite] <- NULL
  }
  model
}

############################################################################# !
# HELPER FUNCTIONS                                                       ####
############################################################################# !

#' Measurement models available in `bmm`
#'
#' @param print_call Logical; If TRUE (default), the function will print
#'   information about how each model function should be called and its required
#'   arguments. If FALSE, the function will return a character vector with the
#'   names of the available models
#' @return A character vector of measurement models available in `bmm`
#' @export
#'
#' @examples
#' supported_models()
supported_models <- function(print_call = TRUE) {
  supported_models <- lsp("bmm", pattern = "^\\.model_")
  supported_models <- sub("^\\.model_", "", supported_models)
  if (!print_call) {
    return(supported_models)
  }

  out <- "The following models are supported:\n\n"
  for (model in supported_models) {
    args <- methods::formalArgs(get(model))
    args <- args[!args %in% c("...")]
    args <- collapse_comma(args)
    args <- gsub("'", "", args)
    out <- glue("{out}- `{model}({args})`\n\n")
  }
  out <- glue("{out}\nType `?modelname` to get information about a specific model, e.g. `?imm`\n")
  out <- gsub("`", " ", out)
  class(out) <- "message"
  out
}


#' @title Generate a markdown list of the measurement models available in `bmm`
#' @description Used internally to automatically populate information in the
#'   README file
#' @return Markdown code for printing the list of measurement models available
#'   in `bmm`
#' @export
#'
#' @examples
#' print_pretty_models_md()
#'
#' @keywords internal
print_pretty_models_md <- function() {
  ok_models <- supported_models(print_call = FALSE)
  domains <- c()
  models <- c()
  for (model in ok_models) {
    m <- get_model(model)()
    domains <- c(domains, m$domain)
    models <- c(models, m$name)
  }
  unique_domains <- unique(domains)
  for (dom in unique_domains) {
    cat("**", dom, "**\n\n", sep = "")
    dom_models <- unique(models[domains == dom])
    for (model in dom_models) {
      cat("*", model, "\n")
    }
    cat("\n")
  }
}

# used to extract well formatted information from the model object to print
# in the @details section for the documentation of each model
model_info <- function(model, components = "all") {
  UseMethod("model_info")
}


#' @export
model_info.bmmodel <- function(model, components = "all") {
  pars <- model$parameters
  par_info <- ""
  if (length(pars) > 0) {
    for (par in names(pars)) {
      par_info <- paste0(par_info, "   - `", par, "`: ", pars[[par]], "\n")
    }
  }

  fixed_pars <- model$fixed_parameters
  fixed_par_info <- ""
  if (length(fixed_pars) > 0) {
    for (fixed_par in names(fixed_pars)) {
      fixed_par_info <- paste0(
        fixed_par_info, "   - `", fixed_par,
        "` = ", fixed_pars[[fixed_par]], "\n"
      )
    }
  }

  links <- model$links
  links_info <- summarise_links(links)

  priors <- model$default_priors
  priors_info <- summarise_default_prior(priors)

  info_all <- list(
    domain = paste0("* **Domain:** ", model$domain, "\n\n"),
    task = paste0("* **Task:** ", model$task, "\n\n"),
    name = paste0("* **Name:** ", model$name, "\n\n"),
    citation = paste0("* **Citation:** \n\n   - ", model$citation, "\n\n"),
    version = paste0("* **Version:** ", model$version, "\n\n"),
    requirements = paste0("* **Requirements:** \n\n  ", model$requirements, "\n\n"),
    parameters = paste0("* **Parameters:** \n\n  ", par_info, "\n"),
    fixed_parameters = paste0("* **Fixed parameters:** \n\n  ", fixed_par_info, "\n"),
    links = paste0("* **Default parameter links:** \n\n     - ", links_info, "\n\n"),
    prior = paste0("* **Default priors:** \n\n", priors_info, "\n")
  )

  if (length(components) == 1 && components == "all") {
    components <- names(info_all)
  }

  if (model$version == "NA" || model$version == "") {
    components <- components[components != "version"]
  }

  # return only the specified components
  collapse(info_all[components])
}



#' @param model A string with the name of the model supplied by the user
#' @return A function of type .model_*
#' @details the returned object is a function. To get the model object, call the
#'   returned function, e.g. `get_model("mixture2p")()`
#' @noRd
get_model <- function(model) {
  get(paste0(".model_", model), mode = "function")
}

# same as get_model2, but with the new model structure for the user facing alias
get_model2 <- function(model) {
  get(model, mode = "function")
}

#' Create a file with a template for adding a new model (for developers)
#'
#' @param model_name A string with the name of the model. The file will be named
#'  `model_model_name.R` and all necessary functions will be created with
#'  the appropriate names and structure. The file will be saved in the `R/`
#'  directory
#' @param testing Logical; If TRUE, the function will return the file content but
#'  will not save the file. If FALSE (default), the function will save the file
#' @param custom_family Logical; Do you plan to define a brms::custom_family()?
#'  If TRUE the function will add a section for the custom family, placeholders
#'  for the stan_vars and corresponding empty .stan files in
#'  `inst/stan_chunks/`, that you can fill For an example, see the sdm
#'  model in `/R/model_sdm.R`. If FALSE (default) the function will
#'  not add the custom family section nor stan files.
#' @param stanvar_blocks A character vector with the names of the blocks that
#'  will be added to the custom family section. See [brms::stanvar()] for more
#'  details. The default lists all the possible blocks, but it is unlikely that
#'  you will need all of them. You can specify a vector of only those that you
#'  need. The function will add a section for each block in the list
#' @param open_files Logical; If TRUE (default), the function will open the
#'  template files that were created in RStudio
#'
#' @return If `testing` is TRUE, the function will return the file content as a
#'  string. If `testing` is FALSE, the function will return NULL
#'
#' @details If you get a warning during check() about non-ASCII characters, this
#'  is often due to the citation field. You can find what the problem is by
#'  running
#'  ```r
#'  remotes::install_github("eddelbuettel/dang")
#'  dang::checkPackageAsciiCode(dir = ".")
#'  ```
#'  usually rewriting the numbers (issue, page numbers) manually fixes it
#' @keywords internal developer
#' @export
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' library(usethis)
#'
#'
#' # create a new model file with a brms::custom_family, three .stan files in
#' # inst/stan_chunks/ and open the files
#' use_model_template("abc",
#'   custom_family = TRUE,
#'   stanvar_blocks = c("functions", "likelihood", "tdata")
#' )
#'
use_model_template <- function(model_name,
                               custom_family = FALSE,
                               stanvar_blocks = c(
                                 "data", "tdata", "parameters",
                                 "tparameters", "model", "likelihood",
                                 "genquant", "functions"
                               ),
                               open_files = TRUE,
                               testing = FALSE) {
  file_name <- paste0("model_", model_name, ".R")

  # check if model exists
  if (model_name %in% supported_models(print_call = FALSE)) {
    stop(paste0("Model ", model_name, " already exists"))
  }
  if (file.exists(paste0("R/", file_name))) {
    stop(paste0("File ", file_name, " already exists"))
  }

  model_header <- glue(
    "#############################################################################!
     # MODELS                                                                 ####
     #############################################################################!
     # see file 'R/model_mixture3p.R' for an example\n\n\n"
  )


  check_data_header <- glue(
    "#############################################################################!
     # CHECK_DATA S3 methods                                                  ####
     #############################################################################!
     # A check_data.* function should be defined for each class of the model.
     # If a model shares methods with other models, the shared methods should be
     # defined in helpers-data.R. Put here only the methods that are specific to
     # the model. See ?check_data for details.
     # (YOU CAN DELETE THIS SECTION IF YOU DO NOT REQUIRE ADDITIONAL DATA CHECKS)\n\n\n"
  )

  bmf2bf_header <- glue(
    "#############################################################################!
     # Convert bmmformula to brmsformla methods                               ####
     #############################################################################!
     # A bmf2bf.* function should be defined if the default method for constructing
     # the brmsformula from the bmmformula does not apply (e.g if aterms are required).
     # The shared method for all `bmmodels` is defined in bmmformula.R.
     # See ?bmf2bf for details.
     # (YOU CAN DELETE THIS SECTION IF YOUR MODEL USES A STANDARD FORMULA WITH 1 RESPONSE VARIABLE)\n\n\n"
  )

  configure_model_header <- glue(
    "#############################################################################!
     # CONFIGURE_MODEL S3 METHODS                                             ####
     #############################################################################!
     # Each model should have a corresponding configure_model.* function. See
     # ?configure_model for more information.\n\n\n"
  )

  postprocess_brm_header <- glue(
    "#############################################################################!
     # POSTPROCESS METHODS                                                    ####
     #############################################################################!
     # A postprocess_brm.* function should be defined for the model class. See
     # ?postprocess_brm for details\n\n\n"
  )


  model_object <- glue('
    .model_<<model_name>> <- function(resp_var1 = NULL, required_arg1 = NULL, required_arg2 = NULL, links = NULL, version = NULL, call = NULL, ...) {
      out <- structure(
        list(
          resp_vars = nlist(resp_var1),
          other_vars = nlist(required_arg1, required_arg2),
          domain = "",
          task = "",
          name = "",
          citation = "",
          version = version,
          requirements = "",
          parameters = list(),
          links = list(),
          fixed_parameters = list(),
          default_priors = list(par1 = list(), par2 = list()),
          void_mu = FALSE
        ),
        class = c("bmmodel", "<<model_name>>"),
        call = call
      )
      if(!is.null(version)) class(out) <- c(class(out), paste0("<<model_name>>_",version))
      out$links[names(links)] <- links
      out
    }\n\n',
    .open = "<<", .close = ">>"
  )

  user_facing_alias <- glue("
    # user facing alias
    # information in the title and details sections will be filled in
    # automatically based on the information in the .model_<<model_name>>()$info\n
    #\' @title `r .model_<<model_name>>()$name`
    #\' @name Model Name,
    #\' @details `r model_info(.model_<<model_name>>())`
    #\' @param resp_var1 A description of the response variable
    #\' @param required_arg1 A description of the required argument
    #\' @param required_arg2 A description of the required argument
    #\' @param links A list of links for the parameters.
    #\' @param version A character label for the version of the model. Can be empty or NULL if there is only one version.
    #\' @param ... used internally for testing, ignore it
    #\' @return An object of class `bmmodel`
    #\' @export
    #\' @examples
    #\' \\dontrun{
    #\' # put a full example here (see 'R/model_mixture3p.R' for an example)
    #\' }
    <<model_name>> <- function(resp_var1, required_arg1, required_arg2, links = NULL, version = NULL, ...) {
       call <- match.call()
       stop_missing_args()
       .model_<<model_name>>(resp_var1 = resp_var1, required_arg1 = required_arg1, required_arg2 = required_arg2,
                    links = links, version = version,call = call, ...)
    }\n\n\n",
    .open = "<<", .close = ">>"
  )

  check_data_method <- glue(
    "#' @export
    check_data.<<model_name>> <- function(model, data, formula) {
       # retrieve required arguments
       required_arg1 <- model$other_vars$required_arg1
       required_arg2 <- model$other_vars$required_arg2\n
       # check the data (required)\n
       # compute any necessary transformations (optional)\n
       # save some variables as attributes of the data for later use (optional)\n
       NextMethod('check_data')
    }\n\n\n",
    .open = "<<", .close = ">>"
  )

  # add bmf2bf method if necessary
  bmf2bf_method <- glue("#' @export
    bmf2bf.<<model_name>> <- function(model, formula) {
       # retrieve required response arguments
       resp_var1 <- model$resp_vars$resp_var1
       resp_var2 <- model$resp_vars$resp_arg2\n
       # set the base brmsformula based
       brms_formula <- brms::bf(paste0(resp_var1, \" | \", vreal(resp_var2), \" ~ 1\"))\n
       # return the brms_formula to add the remaining bmmformulas to it.
       brms_formula
    }\n\n\n",
    .open = "<<", .close = ">>"
  )


  # add custom family section if custom_family is TRUE
  # PS: do not try to replace with glue - already wasted enough time, it doesn't work well
  if (custom_family) {
    family_template <- paste0(
      "   <<model_name>>_family <- brms::custom_family(\n",
      "     '<<model_name>>',\n",
      "     dpars = c(),\n",
      "     links = c(),\n",
      "     lb = c(), # upper bounds for parameters\n",
      "     ub = c(), # lower bounds for parameters\n",
      "     type = '', # real for continous dv, int for discrete dv\n",
      "     loop = TRUE, # is the likelihood vectorized\n",
      "   )\n   formula$family <- <<model_name>>_family\n\n"
    )

    stan_vars_template <- paste0(
      "   # prepare initial stanvars to pass to brms, model formula and priors\n",
      "   sc_path <- system.file('stan_chunks', package='bmm')\n"
    )

    for (stanvar_block in stanvar_blocks) {
      stan_vars_file <- glue("inst/stan_chunks/{model_name}_{stanvar_block}.stan")
      if (!testing) {
        file.create(stan_vars_file)
        if (open_files) {
          usethis::edit_file(stan_vars_file)
        }
      }
      # PS: do not try to replace with glue - already wasted enough time, it doesn't work well
      stan_vars_template <- paste0(
        stan_vars_template,
        "   stan_", stanvar_block, " <- read_lines2(paste0(sc_path, '/", model_name, "_", stanvar_block, ".stan'))\n"
      )
    }
    stan_vars_template <- paste0(stan_vars_template, "\n   stanvars <- ")
    i <- 1
    for (stanvar_block in stanvar_blocks) {
      if (i < length(stanvar_blocks)) {
        stan_vars_template <- paste0(stan_vars_template, "stanvar(scode = stan_", stanvar_block, ", block = '", stanvar_block, "') +\n      ")
        i <- i + 1
      } else {
        stan_vars_template <- paste0(stan_vars_template, "stanvar(scode = stan_", stanvar_block, ", block = '", stanvar_block, "')\n\n")
      }
    }
    out_template <- "   nlist(formula, data, stanvars)\n"
  } else {
    stan_vars_template <- ""
    family_template <- "   formula$family <- NULL\n\n"
    out_template <- "   nlist(formula, data)\n"
  }


  family_comment <- ifelse(custom_family,
    "   # construct the family & add to formula object\n",
    "   # add family to formula object\n"
  )

  # PS: do not try to replace with glue - already wasted enough time, it doesn't work well
  configure_model_method <- glue::glue("#' @export\n",
    "configure_model.<<model_name>> <- function(model, data, formula) {\n",
    "   # retrieve required arguments\n",
    "   required_arg1 <- model$other_vars$required_arg1\n",
    "   required_arg2 <- model$other_vars$required_arg2\n\n",
    "   # retrieve arguments from the data check\n",
    "   my_precomputed_var <- attr(data, 'my_precomputed_var')\n\n",
    "   # construct brms formula from the bmm formula\n",
    "   formula <- bmf2bf(model, formula)\n\n",
    family_comment,
    family_template,
    stan_vars_template,
    "   # return the list\n",
    out_template,
    "}\n\n",
    .open = "<<", .close = ">>"
  )

  postprocess_brm_method <- glue(
    "#' @export
    postprocess_brm.<<model_name>> <- function(model, fit) {
       # any required postprocessing (if none, delete this section)
       fit
    }\n",
    .open = "<<", .close = ">>"
  )

  file_content <- paste0(
    model_header,
    model_object,
    user_facing_alias,
    check_data_header,
    check_data_method,
    bmf2bf_header,
    bmf2bf_method,
    configure_model_header,
    configure_model_method,
    postprocess_brm_header,
    postprocess_brm_method
  )

  if (!testing) {
    writeLines(file_content, paste0("R/", file_name))
    if (open_files) {
      usethis::edit_file(paste0("R/", file_name))
    }
  } else {
    cat(file_content)
  }
}


#' @title Generate Stan code for bmm models
#' @description Given the `model`, the `data` and the `formula` for the model,
#'   this function will return the combined stan code generated by `bmm` and
#'   `brms`
#'
#' @inheritParams bmm
#' @aliases stancode
#' @param object A `bmmformula` object
#' @param ... Further arguments passed to [brms::stancode()]. See the
#'   description of [brms::stancode()] for more details
#'
#' @return A character string containing the fully commented Stan code to fit a
#'   bmm model.
#'
#' @seealso [supported_models()], [brms::stancode()]
#' @keywords extract_info
#' @examples
#' scode1 <- stancode(bmf(c ~ 1, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#' cat(scode1)
#' @importFrom brms stancode
#' @export
stancode.bmmformula <- function(object, data, model, prior = NULL, ...) {
  withr::local_options(bmm.sort_data = FALSE)

  # check model, formula and data, and transform data if necessary
  formula <- object
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)

  # generate the model specification to pass to brms later
  config_args <- configure_model(model, data, formula)

  # configure the default prior and combine with user-specified prior
  prior <- configure_prior(model, data, config_args$formula, prior)

  # extract stan code
  dots <- list(...)
  fit_args <- combine_args(nlist(config_args, dots, prior))
  fit_args$object <- fit_args$formula
  fit_args$formula <- NULL
  code <- brms::do_call(brms::stancode, fit_args)
  add_bmm_version_to_stancode(code)
}

add_bmm_version_to_stancode <- function(stancode) {
  version <- packageVersion("bmm")
  text <- paste0("and bmm ", version)
  brms_comp <- regexpr("brms.*(?=\\n)", stancode, perl = T)
  insert_loc <- brms_comp + attr(brms_comp, "match.length") - 1
  new_stancode <- paste0(
    substr(stancode, 1, insert_loc),
    " ", text,
    substr(stancode, insert_loc + 1, nchar(stancode))
  )
  class(new_stancode) <- class(stancode)
  new_stancode
}

############################################################################# !
# StanCode Helper FUNCTIONS                                              ####
############################################################################# !

# Return integer positions of all occurrences of a fixed pattern in x
# simple wrapper around gregexpr that return a 0-length vector if no matches
# instead of -1
which_positions <- function(x, pat) {
  m <- gregexpr(pat, x, fixed = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L) integer(0) else as.integer(m)
}


#' Find the matching closing brace for an opening brace position
#'
#' Given a string containing brace-like delimiters, return the index of the
#' closing brace that matches the opening brace at \code{open_pos}, accounting
#' for nested braces.
#'
#' The function operates on indices of all occurrences of \code{open_brace} and
#' \code{close_brace} in \code{x} and uses a cumulative-sum depth counter to
#' identify the first position at which nesting depth returns to zero.
#'
#' @param x A length-1 character string.
#' @param open_pos A 1-based integer index into \code{x} indicating the position
#'   of an opening brace character.
#' @param open_brace A length-1 character string giving the opening delimiter
#'   (default \code{"\{"}).
#' @param close_brace A length-1 character string giving the closing delimiter
#'   (default \code{"\}"}).
#'
#' @returns A single integer: the 1-based index of the matching closing brace in
#'   \code{x}.
#' @noRd
#' @examples
#' find_matching_brace("{}", 1L)
#'
#' x <- "{{}{}}"
#' find_matching_brace(x, 1L) # outer brace -> 6
#' find_matching_brace(x, 2L) # inner brace -> 3
#' find_matching_brace(x, 4L) # inner brace -> 5
#'
#' y <- "abc{def{ghi}jkl}mno"
#' find_matching_brace(y, 4L) # -> 16
#' find_matching_brace(y, 8L) # -> 12
find_matching_brace <- function(x, open_pos,
                                open_brace = "{", close_brace = "}") {
  stopifnot(is.character(x), length(x) == 1L)
  n <- nchar(x)
  stopifnot(length(open_pos) == 1L, open_pos >= 1L, open_pos <= n)
  stopif(
    substr(x, open_pos, open_pos) != open_brace,
    "Character at open_pos {open_pos} is not open_brace."
  )

  opens <- which_positions(x, open_brace)
  closes <- which_positions(x, close_brace)

  # combine + sort
  pos <- c(opens, closes)
  stopif(!length(pos), "No braces found.")
  delta <- c(rep.int(1L, length(opens)), rep.int(-1L, length(closes)))

  delta <- delta[order(pos)]
  pos <- pos[order(pos)]

  # locate the opening brace occurrence in the brace stream
  k0 <- match(as.integer(open_pos), pos)
  stopif(is.na(k0), "open_pos was not found among opening brace indices (unexpected).")
  stopif(delta[k0] != 1L, "open_pos does not correspond to an opening brace in the stream.")

  # cumulative sum from that opening brace
  cs <- cumsum(delta[k0:length(delta)])

  # first return to 0 gives the matching close
  j <- which(cs == 0L)[1]
  stopif(is.na(j), "No matching closing brace found (unbalanced braces?).")

  pos[k0 + j - 1L]
}

#' @title Extract code from different STAN program blocks
#'
#' @description
#'   This function extracts the code from the different program blocks of a STAN
#'   program. This can be used in combination with the `stancode` function to
#'   access information about the STAN code generated by `brms` and `bmm`.
#'
#' @param stan_code The STAN code for which the elements should be extracted
#' @param blocks A character vector specifying for which program blocks
#'   the code should be extracted. The default extracks all standard blocks:
#'   "functions", "data", "transformed data", "parameters", "transformed parameters",
#'   "model", and "generated quantities"
#'
#' @return A named list with each element containing the code of one of the STAN
#'   program blocks. If a block
#'
#' @keywords extract_info
#'
#' @examples
#' # generate simple stan code from brms
#' stan_code <- stancode(brms::bf(x ~ 1), data = data.frame(x = rnorm(100)))
#'
#' extracted_program_blocks <- extract_stan_blocks(stan_code)
#'
#' @export
extract_stan_blocks <- function(stan_code, blocks = c("functions", "data", "transformed data", "parameters", "transformed parameters", "model", "generated quantities")) {
  stopifnot(is.character(stan_code), length(stan_code) == 1L)
  blocks <- match.arg(blocks, several.ok = TRUE)
  out <- lapply(blocks, function(b) .extract_stan_block(stan_code, b))
  names(out) <- blocks
  out
}

.extract_stan_block <- function(stan_code, block,
                                include_braces = FALSE,
                                trim = TRUE) {
  stopifnot(is.character(block), length(block) == 1L)

  # Anchor to start-of-line (multiline mode), allow leading spaces/tabs.
  # Match the exact block name, then optional whitespace, then '{'.
  header_pat <- paste0("(?m)^[[:space:]]*", block, "[[:space:]]*\\{")
  m <- regexpr(header_pat, stan_code, perl = TRUE)
  if (m[1] == -1L) {
    return(NULL)
  }

  open_pos <- m[1] + attr(m, "match.length") - 1L # points at '{'
  close_pos <- find_matching_brace(stan_code, open_pos)

  out <- if (include_braces) {
    substr(stan_code, open_pos, close_pos)
  } else {
    substr(stan_code, open_pos + 1L, close_pos - 1L)
  }

  if (trim) out <- sub("^\\s+|\\s+$", "", out)
  out
}


#' @title Extract dimension from parameters in STAN parameter block
#'
#' @description
#'  This functions extracts the names, dimensions, and types from a compiled STAN
#'  parameters blocks generated by bmm or brms. This function is used to specify
#'  initial values for bmm models.
#'
#' @param parameters_block The parameters block extracted via `extract_stan_blocks`
#'
#' @return A list of all parameters, their types, and dimensions as as specified in
#'   the STAN data generated by bmm and brms
#'
#' @keywords extract_info
#'
#' @examples
#' # generate simple stan code from brms
#' stan_code <- stancode(brms::bf(x ~ 1 + cond + (1 + cond | ID)),
#'   data = data.frame(
#'     x = rnorm(100),
#'     ID = rep(1:50, each = 2),
#'     cond = rep(1:2, times = 50)
#'   )
#' )
#'
#' extracted_program_blocks <- extract_stan_blocks(stan_code)
#'
#' par_dims <- extract_parameter_dimensions(extracted_program_blocks$parameters) #'
#'
#' @export
extract_parameter_dimensions <- function(parameters_block) {
  lines <- unlist(strsplit(parameters_block, "\n"))
  lines <- trimws(lines)
  lines <- gsub("//.*", "", lines)
  lines <- lines[nzchar(lines)]

  res <- lapply(lines, parse_parameters_line)
  names(res) <- vapply(res, `[[`, character(1), "name")
  res
}

parse_parameters_line <- function(x) {
  # strip trailing comments and semicolon; normalize whitespace
  x <- sub("//.*$", "", x)
  x <- trimws(x)
  x <- sub(";$", "", x)
  if (!nzchar(x)) stop("Empty or comment-only line.", call. = FALSE)

  # 1) optional leading array[...] prefix
  array_dims <- character(0)
  if (grepl("^array\\s*\\[", x, perl = TRUE)) {
    m <- regexpr("^array\\s*\\[([^\\]]*)\\]\\s*", x, perl = TRUE)
    if (m > 0) {
      dims_str <- sub("^array\\s*\\[([^\\]]*)\\]\\s*.*$", "\\1", regmatches(x, m), perl = TRUE)
      array_dims <- trimws(unlist(strsplit(dims_str, ",")))
      array_dims <- array_dims[nzchar(array_dims)]
      x <- sub("^array\\s*\\[[^\\]]*\\]\\s*", "", x, perl = TRUE)
    }
  }
  is_array <- length(array_dims) > 0

  # 2) detect base type
  base_types <- c(
    "cholesky_factor_corr", "cholesky_factor_cov",
    "corr_matrix", "cov_matrix",
    "row_vector", "unit_vector", "positive_ordered", "simplex", "ordered",
    "vector", "matrix", "real", "int"
  )
  base_type <- NULL
  for (bt in base_types) {
    if (grepl(paste0("^", bt, "\\b"), x, perl = TRUE)) {
      base_type <- bt
      break
    }
  }
  if (is.null(base_type)) stop("Unknown or unsupported base type in: ", x, call. = FALSE)

  # consume the base type token
  x_after_bt <- sub(paste0("^", base_type), "", x, perl = TRUE)

  # 3) constraints <...>
  bounds <- parse_bounds(x_after_bt)
  if (!is.null(bounds)) {
    x_after_bt <- sub("^\\s*<[^>]*>\\s*", "", x_after_bt, perl = TRUE)
  }

  # 4) base-type dims [ ... ] (needed for most non-scalars)
  base_dims <- character(0)
  if (grepl("^\\s*\\[", x_after_bt, perl = TRUE)) {
    m <- regexpr("(?<=\\[)[^\\]]+(?=\\])", x_after_bt, perl = TRUE)
    if (m[1] == -1) stop("Could not parse base dimensions.", call. = FALSE)
    dims_str <- regmatches(x_after_bt, m)
    base_dims <- trimws(strsplit(dims_str, ",", fixed = TRUE)[[1]])
    base_dims <- base_dims[nzchar(base_dims)]
    x_after_bt <- sub("^\\s*\\[[^\\]]*\\]\\s*", "", x_after_bt, perl = TRUE)
  } else {
    if (base_type %in% c(
      "vector", "row_vector", "matrix",
      "simplex", "unit_vector", "ordered", "positive_ordered",
      "corr_matrix", "cov_matrix", "cholesky_factor_corr", "cholesky_factor_cov"
    )) {
      stop("Missing dimensions for type ", base_type, ".", call. = FALSE)
    }
  }

  # 5) name
  name <- trimws(x_after_bt)
  if (!nzchar(name)) stop("Missing parameter name.", call. = FALSE)

  # 6) normalize dims by base type
  dims_by_type <- switch(base_type,
    real                 = 1,
    int                  = 1,
    vector               = base_dims[1],
    row_vector           = base_dims[1],
    simplex              = base_dims[1],
    unit_vector          = base_dims[1],
    ordered              = base_dims[1],
    positive_ordered     = base_dims[1],
    matrix               = base_dims[1:2],
    corr_matrix          = base_dims[1],
    cov_matrix           = base_dims[1],
    cholesky_factor_corr = base_dims[1],
    cholesky_factor_cov  = base_dims[1],
    character(0)
  )
  dims <- c(array_dims, dims_by_type)

  # NEW: dual typing
  types <- if (is_array) c("array", base_type) else base_type

  list(
    name    = name,
    type    = base_type, # backward-compat: base type only
    types   = types, # NEW: includes "array" when applicable
    dims    = dims,
    bounds  = bounds
  )
}

# helper: parse <...> constraints into a named list
parse_bounds <- function(s) {
  if (!grepl("<[^>]*>", s, perl = TRUE)) {
    return(NULL)
  }
  inside <- sub(".*?<([^>]*)>.*", "\\1", s, perl = TRUE)
  parts <- trimws(unlist(strsplit(inside, ",")))
  kvs <- lapply(parts, function(p) {
    if (!grepl("=", p, fixed = TRUE)) {
      return(NULL)
    }
    sp <- strsplit(p, "=", fixed = TRUE)[[1]]
    setNames(list(trimws(sp[2])), trimws(sp[1]))
  })
  # merge into a single named list
  kvs <- Filter(Negate(is.null), kvs)
  if (!length(kvs)) {
    return(list())
  }
  Reduce(function(a, b) c(a, b), kvs)
}

