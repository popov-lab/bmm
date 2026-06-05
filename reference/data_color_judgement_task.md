# Example Data from a Color Judgement Task

Raw data from 50 subjects that completed a color judgement task.
Participants had to detect if there were more blue or red dots in a
cloud of colored dots. The task included a difficulty manipulation coded
by the `condition` variable. In the easy condition 60 percent of the
dots were either blue or red, and in the difficult condition 52 percent
were either blue or red. In addition the fixation cross preceding the
onset of the dot matrix had a variable presentation time coded in the
`fix_duration` variable

## Usage

``` r
data_color_judgement_task
```

## Format

### `data_color_judgement_task`

A data frame with 9941 rows and 10 columns:

- ID:

  Character String uniquely identifying each subject

- trial:

  The trial number of the experimental trials

- condition:

  A character variable coding the difficulty condition

- fix_duration:

  The duration in ms the fixation cross was shown before onset of the
  dot matrix.

- target_response:

  The color / response that was more frequent in the dot matrix

- response:

  The response coded in red / blue given by the participant

- response_correct:

  Boolean variable coding if the response was correct or not.

- rt:

  The reaction time in seconds.

- pressed_key:

  The actual key pressed by the participants

- correct_key:

  The correct key that should have been pressed.
