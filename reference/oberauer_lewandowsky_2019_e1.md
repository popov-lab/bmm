# Data from Experiment 1 reported by Oberauer & Lewandowsky (2019)

Raw data of 40 subjects that completed a verbal memory recall task in
three different conditions using different types of distractor words.

## Usage

``` r
oberauer_lewandowsky_2019_e1
```

## Format

### `oberauer_lewandowsky_2019_e1`

A data frame with 120 rows and 10 columns:

- ID:

  Integer uniquely identifying each subject

- cond:

  Factor sperating the three experimental conditions: `new distractors`
  refers to new words being used as distractors, `old reordered` refers
  to the to be remembered words being the distractors, but reordered
  relative to the serial position, `old same` refers to the to be
  remebered words being the distractors, and appearing in the same order
  as the to be remembered words.

- corr:

  The frequency a subject recalled the correct item

- other:

  The frequency a subject recalled one of the other to be remebered
  words

- dist:

  The frequency a subject recalled one of the distractors

- npl:

  The frequency a subject recalled a not-presented lure (NPL), that is a
  word that was not presented during a trial

- n_corr, n_other, n_dist, n_npl:

  The number of candidataes in each of the response categories
