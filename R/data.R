#' Data from Experiment 2 reported by Zhang & Luck (2008)
#'
#' Raw data of 8 subjects for the response error in a continuous reproduction task
#' with set size 1, 2, 3, and 6 reported by Zhang & Luck (2008).
#'
#' @format ## `zhang_luck_2008`
#' A data frame with 4,000 rows and 9 columns:
#' \describe{
#'   \item{subID}{Integer uniquely identifying different subjects}
#'   \item{trial}{Trial identifyier}
#'   \item{setsize}{The set_size of the data in this row}
#'   \item{response_error}{The response error, that is the difference between the response
#'   given and the target color in radians.}
#'   \item{col_lure1, col_Lure2, col_Lure3, col_Lure4, col_Lure5}{Color value of the lure items coded relative to the target color.}
#'
#' }
#' @keywords dataset
#' @source <https://www.nature.com/articles/nature06860>
"zhang_luck_2008"


#' Data from Experiment 1 reported by Oberauer & Lin (2017)
#'
#' Raw data of 19 subjects that completed a continuous reproduction task
#' with set size 1 to 8 reported by Oberauer & Lin (2017).
#'
#' @format ## `oberauer_lin_2017`
#' A data frame with 15,200 rows and 19 columns:
#' \describe{
#'   \item{ID}{Integer uniquely identifying different subjects}
#'   \item{session}{Session number}
#'   \item{trial}{Trial number within each session}
#'   \item{set_size}{The set_size of the data in this row}
#'   \item{dev_rad}{The response error, that is the difference between the response
#'   given and the target color in radians.}
#'   \item{col_nt1, col_nt2, col_nt3, col_nt4, col_nt5, col_nt6, col_nt7}{The non-target items' color value relative to the target.}
#'   \item{dist_nt1, dist_nt2, dist_nt3, dist_nt4, dist_nt5, dist_nt6, dist_nt7, dist_nt8}{The spatial distance between all non-target items and the target item in radians.}
#'
#' }
#'
#' @keywords dataset
#' @source <https://osf.io/m4shu>
"oberauer_lin_2017"

#' Data from Experiment 1 reported by Oberauer & Lewandowsky (2019)
#'
#' Raw data of 40 subjects that completed a verbal memory recall task in three different
#' conditions using different types of distractor words.
#'
#' @format ## `oberauer_lewandowsky_2019_e1`
#' A data frame with 120 rows and 10 columns:
#' \describe{
#'   \item{ID}{Integer uniquely identifying each subject}
#'   \item{cond}{Factor sperating the three experimental conditions: `new distractors` refers to
#'   new words being used as distractors, `old reordered` refers to the to be remembered words
#'   being the distractors, but reordered relative to the serial position, `old same` refers
#'   to the to be remebered words being the distractors, and appearing in the same order as
#'   the to be remembered words.}
#'   \item{corr}{The frequency a subject recalled the correct item}
#'   \item{other}{The frequency a subject recalled one of the other to be remebered words}
#'   \item{dist}{The frequency a subject recalled one of the distractors}
#'   \item{npl}{The frequency a subject recalled a not-presented lure (NPL), that is a word
#'   that was not presented during a trial}
#'   \item{n_corr, n_other, n_dist, n_npl}{The number of candidataes in each of the response categories}
#' }
#' @keywords dataset
"oberauer_lewandowsky_2019_e1"


#' Example Data from a Color Judgement Task
#'
#' Raw data from 50 subjects that completed a color judgement task. Participants
#' had to detect if there were more blue or red dots in a cloud of colored dots.
#' The task included a difficulty manipulation coded by the `condition` variable.
#' In the easy condition 60 percent of the dots were either blue or red, and in
#' the difficult condition 52 percent were either blue or red. In addition the
#' fixation cross preceding the onset of the dot matrix had a variable presentation
#' time coded in the `fix_duration` variable
#'
#' @format ## `data_color_judgement_task`
#' A data frame with 9941 rows and 10 columns:
#' \describe{
#'   \item{ID}{Character String uniquely identifying each subject}
#'   \item{trial}{The trial number of the experimental trials}
#'   \item{condition}{A character variable coding the difficulty condition}
#'   \item{fix_duration}{The duration in ms the fixation cross was shown before onset of the dot matrix.}
#'   \item{target_response}{The color / response that was more frequent in the dot matrix}
#'   \item{response}{The response coded in red / blue given by the participant}
#'   \item{response_correct}{Boolean variable coding if the response was correct or not.}
#'   \item{rt}{The reaction time in seconds.}
#'   \item{pressed_key}{The actual key pressed by the participants}
#'   \item{correct_key}{The correct key that should have been pressed.}
#' }
#' @keywords dataset
"data_color_judgement_task"


#' Recognition ROC data from Broeder & Schuetz (2009, Experiment 3)
#'
#' Binary old/new recognition data from 40 subjects, aggregated to response
#' counts. Each subject was tested under five base-rate conditions: the
#' proportion of old items shifts the decision criterion from conservative
#' (`br1`) to liberal (`br5`) while leaving sensitivity unchanged, tracing a
#' five-point binary ROC per subject. This criterion variation is what makes the
#' unequal-variance ratio (`sdratio`) of [sdt_yn()] identifiable: a single
#' condition yields only one hit/false-alarm pair and cannot separate a wider
#' signal distribution from a larger d'. Counts were digitised from the
#' frequencies reported in the original article.
#'
#' @format ## `broeder_schuetz_2009_e3`
#' A data frame with 400 rows (40 subjects x 5 conditions x 2 stimulus types)
#' and 5 columns:
#' \describe{
#'   \item{id}{Integer uniquely identifying each subject}
#'   \item{condition}{Factor with five base-rate conditions, ordered from the
#'   most conservative (`br1`) to the most liberal (`br5`) induced criterion}
#'   \item{stimulus}{Integer stimulus type: 0 = new/lure, 1 = old/target}
#'   \item{n_old}{Integer count of "old" responses in that cell: hits for old
#'   items (`stimulus == 1`) and false alarms for new items (`stimulus == 0`)}
#'   \item{n_trials}{Integer number of items presented in that cell}
#' }
#' @keywords dataset
#' @source Broeder, A., & Schuetz, J. (2009). Recognition ROCs are curvilinear---or
#'   are they? On premature arguments against the two-high-threshold model of
#'   recognition. \emph{Journal of Experimental Psychology: Learning, Memory, and
#'   Cognition}, 35(3), 587--606. \doi{10.1037/a0015279}
#' @examples
#' \dontrun{
#' # Unequal-variance yes/no SDT: the criterion varies across base-rate
#' # conditions, while sensitivity (d) and the signal/noise SD ratio (sdratio)
#' # are held constant across conditions.
#' model <- sdt_yn(
#'   response = "n_old", stimulus = "stimulus", n_trials = "n_trials"
#' )
#' fit <- bmm(
#'   formula = bmf(
#'     d ~ 1 + (1 | id),
#'     criterion ~ 0 + condition + (1 | id),
#'     sdratio ~ 1
#'   ),
#'   data = broeder_schuetz_2009_e3,
#'   model = model,
#'   backend = "cmdstanr"
#' )
#' }
"broeder_schuetz_2009_e3"
