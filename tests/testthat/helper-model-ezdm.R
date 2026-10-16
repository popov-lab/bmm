# the cases the shared tests run for ezdm; see helper-model-cases.R
model_cases_ezdm <- function() {
  three <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")
  four <- ezdm(
    mean_rt = c("mean_rt_upper", "mean_rt_lower"), var_rt = c("var_rt_upper", "var_rt_lower"),
    n_upper = "n_upper", n_trials = "n_trials", version = "4par"
  )
  four_data <- withr::with_seed(1, rezdm(5, n_trials = 40, drift = 0.3, bound = 1.2,
                                         ndt = 0.3, version = "4par"))
  re_formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1 + (1 | id), ndt ~ 1 + (1 | id))
  list(
    stored_frame = list(
      ezdm = list(model = three, formula = bmf(drift ~ 1, bound ~ 1, ndt ~ 1), data = toy_ezdm_data())
    ),
    sd_default = list(
      ezdm = list(
        model = three, formula = re_formula, data = toy_ezdm_data(),
        sd = c(drift = "exponential(1)", bound = "exponential(2)", ndt = "exponential(2)")
      )
    ),
    init = list(
      ezdm = list(model = three, formula = re_formula, data = toy_ezdm_data()),
      ezdm_4par = list(model = four, formula = bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1),
                       data = four_data)
    ),
    links = list(
      ezdm = list(model = three, settable = c("drift", "bound", "ndt", "s")),
      ezdm_4par = list(model = four, settable = c("drift", "bound", "ndt", "zr", "s"))
    ),
    fixture = list(
      ezdm = list(file = "bmmfit_ezdm3_ppcheck.rds"),
      ezdm_4par = list(file = "bmmfit_ezdm4_ppcheck.rds")
    ),
    observables = list(
      ezdm = list(
        model = ezdm(mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
                     n_trials = "n_trials"),
        data = withr::with_seed(1, rezdm(5, n_trials = 40, drift = 0.3, bound = 1.2, ndt = 0.3)),
        formula = bmf(drift ~ 1)
      ),
      ezdm_4par = list(model = four, data = four_data, formula = bmf(drift ~ 1))
    )
  )
}
