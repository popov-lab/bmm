# Continuous reproduction tasks (CRT)

## 1 Task description

In research on visual working memory participants are often asked to
remember and reproduce continuous features of visual objects such as
their color or orientation. In the continuous reproduction task
(sometimes also called delayed estimation task), participants encode a
set of visual objects into visual working memory and are then asked to
reproduce a specific feature of one cued object on a continuous scale at
test (see Figure [1.1](#fig:vwmcrt) for an illustration).

Most often the features used in these tasks are colors sampled from a
color wheel (Wilken and Ma 2004) or continuous orientations of a bar or
a triangle (Bays et al. 2011). The set of to-be-remembered objects
typically consists of one up to eight objects spatially distributed over
the screen. Thus, participants must associate the to-be-remembered
features (e.g. color or orientation) with the spatial locations they are
presented at. The precision of the representation of an object’s feature
in visual working memory is measured as the angular deviation from the
true feature presented at encoding.

![A typical continuous reproduction task](assets/vwm-crt.png)

Figure 1.1: A typical continuous reproduction task

## 2 The role of measurement models

In these continuous reproduction tasks, the simplest measure of
performance is the average angle deviation of the response from the true
feature value. In many studies, this average recall error has been the
main dependent variable for evaluating the effect of experimental
manipulations. Yet, the average recall error confounds different
properties of memory representations and does not sufficiently represent
the theoretical processes assumed by current models of visual working
memory. Therefore, different measurement models have been proposed to
formalize distinct aspects of visual working memory models and how they
translate into observed behavior.

A measurement model is a statistical model that describes the
relationship between latent cognitive processes and observed behavior.
For continuous reproduction tasks, measurement models provide a more
refined representation of memory processes because they decompose the
average recall error into several theoretically meaningful parameters.

At the core of these models is the assumption that responses in
continuous reproduction tasks can stem from different distributions
depending on the continuous activation of different memory
representation or the cognitive state a person is in at recall.

## 3 CRT Models in the `bmm` package

The `bmm` package implements several measurement models for analyzing
continuous reproduction data:

##### The two-parameter mixture model (Zhang and Luck 2008)

- see [`?mixture2p`](https://venpopov.com/bmm/reference/mixture2p.md)
  and [the mixture models
  article](https://venpopov.com/bmm/articles/bmm_mixture_models.html)

##### The three-parameter mixture model (Bays et al. 2009)

- see [`?mixture3p`](https://venpopov.com/bmm/reference/mixture3p.md)
  and [the mixture models
  article](https://venpopov.com/bmm/articles/bmm_mixture_models.html)

##### The Interference Measurement Model (Oberauer and Lin 2017)

- see [`?imm`](https://venpopov.com/bmm/reference/imm.md) and [the IMM
  article](https://venpopov.com/bmm/articles/bmm_imm.html)

##### The Signal Discrimination Model (SDM) by (Oberauer 2023)

- see [`?sdm`](https://venpopov.com/bmm/reference/sdm.md) and [the SDM
  article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html)

## 4 Preparing data from half-circular stimulus spaces

As already mentioned, some task require subjects to remember
orientations (e.g. of bars or Gabor patches) without a direction. With
such stimulus material the response error can range only from -90 to 90
degrees (or -pi/2 to pi/2). When using data from such a task, you have
to multiply the `response_error` by 2 when pre-processing the data, so
that the response error ranges from -180 to 180 degrees (or -pi to pi).
The same applies to the `nt_features` relative to the target
orientation.

## References

Bays, Paul M., Raquel F. G. Catalao, and Masud Husain. 2009. “The
Precision of Visual Working Memory Is Set by Allocation of a Shared
Resource.” *Journal of Vision* 9 (10): 7–7.
<https://doi.org/10.1167/9.10.7>.

Bays, Paul M., Nikos Gorgoraptis, Natalie Wee, Louise Marshall, and
Masud Husain. 2011. “Temporal Dynamics of Encoding, Storage, and
Reallocation of Visual Working Memory.” *Journal of Vision* 11 (10):
6–6.

Oberauer, Klaus. 2023. “Measurement Models for Visual Working Memory—a
Factorial Model Comparison.” *Psychological Review* (US) 130 (3):
841–52. <https://doi.org/10.1037/rev0000328>.

Oberauer, Klaus, and Hsuan-Yu Lin. 2017. “An Interference Model of
Visual Working Memory.” *Psychological Review* 124 (1): 21–59.
<https://doi.org/10.1037/rev0000044>.

Wilken, Patrick, and Wei Ji Ma. 2004. “A Detection Theory Account of
Change Detection.” *Journal of Vision* 4 (12): 11–11.

Zhang, Weiwei, and Steven J. Luck. 2008. “Discrete Fixed-Resolution
Representations in Visual Working Memory.” *Nature* 453 (7192): 233–35.
<https://doi.org/10.1038/nature06860>.
