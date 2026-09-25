# Distribution functions for the EZ-Diffusion Model (ezdm)

Density and random generation functions for the EZ-Diffusion Model. The
model operates on aggregated data: mean reaction time, variance of
reaction time, and number of responses to the upper boundary.

## Usage

``` r
dezdm(
  mean_rt,
  var_rt,
  n_upper,
  n_trials,
  drift,
  bound,
  ndt,
  zr = 0.5,
  s = 1,
  version = c("3par", "4par"),
  log = TRUE
)

rezdm(
  n,
  n_trials,
  drift,
  bound,
  ndt,
  zr = 0.5,
  s = 1,
  version = c("3par", "4par")
)
```

## Arguments

- mean_rt:

  Observed mean reaction time(s) in seconds. For version "3par", a
  numeric vector or single value. For version "4par", either a vector of
  length 2 (c(mean_rt_upper, mean_rt_lower)) for single observation, or
  a matrix with 2 columns for multiple observations.

- var_rt:

  Observed variance of reaction times in seconds^2. For version "3par",
  a numeric vector or single value. For version "4par", either a vector
  of length 2 (c(var_rt_upper, var_rt_lower)) for single observation, or
  a matrix with 2 columns for multiple observations.

- n_upper:

  Number of responses to the upper boundary

- n_trials:

  Total number of trials

- drift:

  Drift rate (evidence accumulation rate; can be positive or negative
  for below-chance performance).

- bound:

  Boundary separation (distance between decision thresholds).

- ndt:

  Non-decision time (seconds).

- zr:

  Relative starting point (0 to 1). Only used for version "4par".

- s:

  Diffusion constant (standard deviation of noise), default = 1.

- version:

  Character; either "3par" (default) or "4par"

- log:

  Logical; if `TRUE`, values are returned on the log scale.

- n:

  Number of samples to generate

## Value

`dezdm` gives the log-density of the observed summary statistics under
the EZDM, and `rezdm` generates random summary statistics from the
implied sampling distributions.

## Details

The number of upper-boundary responses is binomial. The two RT summaries
follow the joint sampling distribution of the mean and the variance of
`n` independent decision times, matched to the exact first four
cumulants of the first-passage time. Writing \\\mathrm{MDT}\\ and
\\\mathrm{VRT}\\ for its mean and variance and \\\kappa_3\\,
\\\kappa_4\\ for its third and fourth cumulants, and \\W = \kappa_4 /
n + 2\\\mathrm{VRT}^2 / (n - 1)\\ for the exact variance of the sample
variance, \$\$\mathrm{var\\rt} \sim \mathrm{Gamma}(\mathrm{VRT}^2 / W,
\mathrm{VRT} / W),\$\$ \$\$\mathrm{mean\\rt} \mid \mathrm{var\\rt} \sim
N\left(\mathrm{ndt} + \mathrm{MDT} + \frac{\kappa_3 /
n}{W}(\mathrm{var\\rt} - \mathrm{VRT}),\\ \sqrt{\mathrm{VRT} / n -
(\kappa_3 / n)^2 / W}\right),\$\$ so that
\\\mathrm{Var}(\mathrm{mean\\rt}) = \mathrm{VRT} / n\\ and
\\\mathrm{Cov}(\mathrm{mean\\rt}, \mathrm{var\\rt}) = \kappa_3 / n\\ are
exact. Decision times are right-skewed (with a symmetric start point
their kurtosis is 8.8 at zero drift and falls towards 3 as drift grows),
so the older form that assumes normal reaction times — independent
normal and scaled chi-square \\\mathrm{Gamma}((n - 1)/2, (n -
1)/(2\\\mathrm{VRT}))\\ terms — understates the sampling variance of
`var_rt` by a factor of \\1 + (\mathrm{kurtosis} - 3)(n - 1)/(2n)\\:
about 3.8 with 100 trials per cell and 3.5 to 3.6 with 10, less at
accuracies above .95. With an asymmetric start point in version `"4par"`
it is more at the boundary nearer the start point (up to 5.7 for `zr`
between .3 and .7) and less at the other. It also ignores a correlation
of about 0.7 between the two statistics, making posteriors too narrow.
The terms above reduce to it when \\\kappa_3 = \kappa_4 = 0\\.

For version `"3par"` the start point is symmetric, so the decision time
is independent of which boundary is hit and all `n_trials` responses
inform one set of cumulants. For version `"4par"` the two boundaries
have different decision-time distributions, so each is given its own
summaries and its own response count. A boundary reached fewer than
twice has no sample variance; in `dezdm()` it contributes only through
the binomial term, but
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) currently drops
such cells, because `rezdm()` and
[`ezdm_summary_stats()`](https://venpopov.com/bmm/dev/reference/ezdm_summary_stats.md)
code their missing summaries as `NA` and `brms` excludes rows with
missing values. The per-boundary formulas condition on the realised
counts, which are themselves random.

The two additional cumulants cost sampling time. In two simulated
designs (30 subjects with 200 or 250 trials, 3 seeds each, one machine)
a gradient of the whole model took 1.3 times as long as with the older
form for version `"3par"` and 1.9 times for version `"4par"`. The
likelihood alone took 1.5 to 2.0 and 1.7 to 2.5 times as long, most
where accuracy is near chance, so the ratio grows with the number of
cells.

Simulated `mean_rt` is not truncated at `ndt`. When a summary rests on
few responses (a handful of trials in version `"3par"`, or a rarely
reached boundary in version `"4par"`, whatever `n_trials` is), `rezdm()`
can return `mean_rt` below `ndt` or even `mean_rt <= 0`, which
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) rejects; drop
those rows before fitting.

## References

Wagenmakers, E.-J., Van Der Maas, H. L. J., & Grasman, R. P. P. P.
(2007). An EZ-diffusion model for response time and accuracy.
Psychonomic Bulletin & Review, 14(1), 3-22.

Chávez De la Peña, A. F., & Vandekerckhove, J. (2025). An EZ Bayesian
hierarchical drift diffusion model for response time and accuracy.
Psychonomic Bulletin & Review.

## Examples

``` r
# 3-parameter version (single observation)
dezdm(
  mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
  drift = 2, bound = 1.5, ndt = 0.3
)
#> [1] -30.99022

# 3-parameter version (vectorized)
dezdm(
  mean_rt = c(0.5, 0.55), var_rt = c(0.02, 0.025),
  n_upper = c(80, 75), n_trials = c(100, 100),
  drift = 2, bound = 1.5, ndt = 0.3
)
#> [1] -30.99022 -27.30733

# 4-parameter version (single observation)
dezdm(
  mean_rt = c(0.45, 0.55), var_rt = c(0.018, 0.025),
  n_upper = 80, n_trials = 100,
  drift = 2, bound = 1.5, ndt = 0.3, zr = 0.55, version = "4par"
)
#> [1] -37.10626

# generate random summary statistics
rezdm(n = 100, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3)
#>       mean_rt     var_rt n_upper n_trials
#> 1   0.6291430 0.03027850      96      100
#> 2   0.6492095 0.06582494      93      100
#> 3   0.6374134 0.04844420      97      100
#> 4   0.6246423 0.06195106      95      100
#> 5   0.6442479 0.06416348      97      100
#> 6   0.6582447 0.07978936      86      100
#> 7   0.6146450 0.05987511      95      100
#> 8   0.6319291 0.06053525      93      100
#> 9   0.6427734 0.06728796      95      100
#> 10  0.6409447 0.05690502      93      100
#> 11  0.6214507 0.04894997      94      100
#> 12  0.6625078 0.06973032      94      100
#> 13  0.6415679 0.05052144      97      100
#> 14  0.6098283 0.06584346      97      100
#> 15  0.6540606 0.09346122      97      100
#> 16  0.5956537 0.04299261      96      100
#> 17  0.6588027 0.08463187      97      100
#> 18  0.6548355 0.04276585      93      100
#> 19  0.6603150 0.07401777      94      100
#> 20  0.6219985 0.06426829      93      100
#> 21  0.6525230 0.06463210      94      100
#> 22  0.6524002 0.07795316      94      100
#> 23  0.6383603 0.05601125      90      100
#> 24  0.6232808 0.04159296      97      100
#> 25  0.6124396 0.05957785      94      100
#> 26  0.7168649 0.11293057      91      100
#> 27  0.6400327 0.04395243      97      100
#> 28  0.6115228 0.03963022      98      100
#> 29  0.6070992 0.05355732      97      100
#> 30  0.5984335 0.04713857      96      100
#> 31  0.6611604 0.06462787      91      100
#> 32  0.6370724 0.04166656      90      100
#> 33  0.6296500 0.06143637      96      100
#> 34  0.6482779 0.05982352      94      100
#> 35  0.6613657 0.07306913      93      100
#> 36  0.6657784 0.05625606      93      100
#> 37  0.6461891 0.05701411      94      100
#> 38  0.7051881 0.09011906      95      100
#> 39  0.6128421 0.02929434      93      100
#> 40  0.6622412 0.06875625      97      100
#> 41  0.6673827 0.06939368      96      100
#> 42  0.6393786 0.06073983      97      100
#> 43  0.6713710 0.05403128      98      100
#> 44  0.6520346 0.05920155      95      100
#> 45  0.6547611 0.05273972      94      100
#> 46  0.6543898 0.06830856     100      100
#> 47  0.6834951 0.08005194      94      100
#> 48  0.6456146 0.06550383      97      100
#> 49  0.6337608 0.04058992      94      100
#> 50  0.6269983 0.03760525      94      100
#> 51  0.6627175 0.06084106      91      100
#> 52  0.7007489 0.11105881      96      100
#> 53  0.6544035 0.07092160      95      100
#> 54  0.6142211 0.04421934      96      100
#> 55  0.6272045 0.06404083      98      100
#> 56  0.6404205 0.08109415      93      100
#> 57  0.6440556 0.04912187      94      100
#> 58  0.5999217 0.04408336      97      100
#> 59  0.6580514 0.06680245      96      100
#> 60  0.6221916 0.05288977      94      100
#> 61  0.6550374 0.05700139      93      100
#> 62  0.6605298 0.06514998      98      100
#> 63  0.6019853 0.03445298      97      100
#> 64  0.6923796 0.11841626      98      100
#> 65  0.6058672 0.04289969      97      100
#> 66  0.6198149 0.07331693      98      100
#> 67  0.6394459 0.04227146      98      100
#> 68  0.6584875 0.07289812      95      100
#> 69  0.6447479 0.03979900      98      100
#> 70  0.6341242 0.05185994      96      100
#> 71  0.6886301 0.06827487      98      100
#> 72  0.6201042 0.05201076      92      100
#> 73  0.6127449 0.03901742      99      100
#> 74  0.5902882 0.03053039      96      100
#> 75  0.6442957 0.04189008      94      100
#> 76  0.6367534 0.05546564      92      100
#> 77  0.6470561 0.04079083      94      100
#> 78  0.6466342 0.07105353      91      100
#> 79  0.6728277 0.07139048      97      100
#> 80  0.6848235 0.07617792      96      100
#> 81  0.5999806 0.03823880      96      100
#> 82  0.6110980 0.04934902      95      100
#> 83  0.6511748 0.06412352      94      100
#> 84  0.6220781 0.04592337      97      100
#> 85  0.6459097 0.06048572      94      100
#> 86  0.6403477 0.05700306      94      100
#> 87  0.6351961 0.06126420      95      100
#> 88  0.6610006 0.10102193      98      100
#> 89  0.6453243 0.08172159      94      100
#> 90  0.6120853 0.05954973      97      100
#> 91  0.6267192 0.04742482      95      100
#> 92  0.6254491 0.03562686      95      100
#> 93  0.6595509 0.05883631      97      100
#> 94  0.6099907 0.04270317      95      100
#> 95  0.6531484 0.06075766      95      100
#> 96  0.6000379 0.04489141      95      100
#> 97  0.6138791 0.04395106      93      100
#> 98  0.6321043 0.04573348      95      100
#> 99  0.6464785 0.09181266      97      100
#> 100 0.6348542 0.04545770      96      100
rezdm(
  n = 100, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3,
  zr = 0.55, version = "4par"
)
#>     mean_rt_upper mean_rt_lower var_rt_upper var_rt_lower n_upper n_trials
#> 1       0.6105410     0.6153878   0.04172400 8.964842e-02      95      100
#> 2       0.6064770     0.5654622   0.04749583 1.640271e-05      97      100
#> 3       0.6479446     0.7451780   0.07492969 4.276932e-02      96      100
#> 4       0.5861251     0.6318085   0.03434061 2.899165e-02      93      100
#> 5       0.6195691     0.6494892   0.06955375 1.417763e-03      96      100
#> 6       0.6589714     0.8701973   0.06642816 1.847820e-01      97      100
#> 7       0.5606422     0.4920157   0.04241431 1.034777e-03      97      100
#> 8       0.6201360            NA   0.05966949           NA      99      100
#> 9       0.6370978     0.6525236   0.07589530 8.916151e-02      93      100
#> 10      0.6201676     0.7282987   0.05888262 2.669084e-07      98      100
#> 11      0.6145232     0.6405244   0.05144151 4.632519e-02      97      100
#> 12      0.6025327            NA   0.05363346           NA      99      100
#> 13      0.6119176     0.6720856   0.03899302 1.751221e-01      98      100
#> 14      0.5732375     0.6120217   0.06123352 5.232511e-03      98      100
#> 15      0.5643589     0.4234687   0.02574615 6.908943e-04      96      100
#> 16      0.6122113            NA   0.06037782           NA      99      100
#> 17      0.6082959     0.4691547   0.05198223 5.138635e-03      98      100
#> 18      0.6041515            NA   0.04406634           NA     100      100
#> 19      0.6341519     0.6427258   0.06904517 4.189798e-03      98      100
#> 20      0.5696163     0.4778739   0.04260332 1.840902e-02      95      100
#> 21      0.5887743            NA   0.03575235           NA     100      100
#> 22      0.6349632     0.6092789   0.07158360 1.939911e-03      95      100
#> 23      0.6325510     0.6366945   0.05599486 1.114887e-02      96      100
#> 24      0.5880626     0.4695017   0.03706553 1.442436e-02      97      100
#> 25      0.6417976     0.6995824   0.06050376 2.116798e-07      97      100
#> 26      0.6333740     0.6113200   0.08342104 1.404047e-01      98      100
#> 27      0.6099129     0.6494813   0.05705158 2.839436e-02      95      100
#> 28      0.6343020     0.5254134   0.07705888 1.272284e-01      98      100
#> 29      0.5768736     0.6009186   0.03267448 9.047154e-04      98      100
#> 30      0.6441462     0.7475775   0.07728427 3.413250e-02      96      100
#> 31      0.5757580     0.6908072   0.04556671 5.926274e-04      97      100
#> 32      0.6119220     0.5704954   0.05824639 3.658027e-02      96      100
#> 33      0.5699428     0.7271286   0.02667627 2.231874e-02      97      100
#> 34      0.6077413     0.8173563   0.06320824 2.121883e-03      96      100
#> 35      0.6131613     0.5996415   0.04023633 4.015468e-03      96      100
#> 36      0.6278366     0.8162665   0.05609730 7.799658e-02      97      100
#> 37      0.6617476     0.8525888   0.07985665 1.788829e-01      98      100
#> 38      0.6082052     0.8254431   0.05895344 2.555622e-01      96      100
#> 39      0.6159618     0.7825679   0.06150029 6.626438e-02      95      100
#> 40      0.5543013     0.3822259   0.02844729 6.145950e-04      96      100
#> 41      0.5882336     0.6548785   0.04176370 6.226703e-02      97      100
#> 42      0.6092594     0.7757917   0.02631283 4.799534e-02      96      100
#> 43      0.6038300     0.6925668   0.05738837 3.508513e-02      93      100
#> 44      0.5972670     0.7440663   0.06340468 7.024203e-02      96      100
#> 45      0.5864188     0.7253065   0.03289055 2.067244e-02      91      100
#> 46      0.6387383     0.6582383   0.06294350 4.510694e-02      93      100
#> 47      0.6176930     0.7519802   0.09003310 1.061546e-01      95      100
#> 48      0.6325580     0.7242858   0.08218069 1.018427e-04      98      100
#> 49      0.6183147     0.5012574   0.05662094 6.273373e-04      98      100
#> 50      0.6392508     0.6691435   0.05128741 7.641986e-03      94      100
#> 51      0.5827214     0.7177507   0.04021351 3.460467e-05      97      100
#> 52      0.6197440     0.6063592   0.04059592 6.597293e-02      97      100
#> 53      0.5963650     0.5122099   0.05029105 8.652198e-02      98      100
#> 54      0.5875334     0.6462004   0.05170266 6.565021e-02      97      100
#> 55      0.5991676     0.6479927   0.06490299 3.878416e-02      94      100
#> 56      0.6551721     0.5516797   0.10342628 1.268339e-02      97      100
#> 57      0.6095444     0.5410523   0.05624043 7.627708e-03      94      100
#> 58      0.6374213     0.6967806   0.08802373 1.086419e-02      97      100
#> 59      0.6294178     0.6224703   0.05770797 1.684506e-01      98      100
#> 60      0.6440596            NA   0.08599735           NA     100      100
#> 61      0.6060809     0.7674364   0.05427837 1.383788e-01      97      100
#> 62      0.6299425     0.6528520   0.08544313 1.265284e-03      97      100
#> 63      0.5909394     0.7953506   0.05924569 2.480738e-02      98      100
#> 64      0.5689687     0.7654835   0.04916194 1.731515e-01      97      100
#> 65      0.6046898     0.4989256   0.04763176 6.425726e-02      95      100
#> 66      0.5986981            NA   0.04501262           NA      99      100
#> 67      0.6274647     0.6832746   0.05045769 1.613579e-02      96      100
#> 68      0.6124290     0.6248754   0.05380970 2.802804e-02      96      100
#> 69      0.6240701     0.8416091   0.05647346 1.158669e-01      94      100
#> 70      0.6033102     0.7263990   0.04919170 3.952066e-02      98      100
#> 71      0.5680285            NA   0.02928384           NA      99      100
#> 72      0.6199282     0.5454155   0.06100048 6.728827e-02      94      100
#> 73      0.6473579     0.6639256   0.08661778 4.797269e-02      96      100
#> 74      0.5762139     0.5863706   0.03735191 8.271690e-02      96      100
#> 75      0.6129981     0.4173216   0.05966032 1.129470e-01      98      100
#> 76      0.6396884     0.8708836   0.05201270 1.706749e-01      98      100
#> 77      0.5830684     0.5323943   0.04922296 9.233990e-02      95      100
#> 78      0.6195502     0.8065830   0.04106235 4.523342e-02      96      100
#> 79      0.5804487     0.6881519   0.03166932 1.122860e-01      94      100
#> 80      0.6298899     0.5869249   0.05084530 3.977110e-02      96      100
#> 81      0.5876688     0.4950612   0.04865863 7.836664e-03      98      100
#> 82      0.5790522     0.4696505   0.02762189 2.117349e-03      95      100
#> 83      0.6255699     0.6509479   0.06612254 1.987467e-03      96      100
#> 84      0.6190222     0.4218801   0.06684458 8.349153e-04      97      100
#> 85      0.6318047     0.6635272   0.06646771 2.863020e-04      98      100
#> 86      0.6077563            NA   0.05685653           NA     100      100
#> 87      0.6549499     0.5107883   0.07950133 1.023438e-04      96      100
#> 88      0.5967797            NA   0.04079575           NA      99      100
#> 89      0.6243934            NA   0.05883978           NA      99      100
#> 90      0.6066617     0.7060872   0.07256192 1.782459e-02      98      100
#> 91      0.6140634     0.6939546   0.07137744 3.907729e-02      92      100
#> 92      0.6091063     0.6811034   0.03736169 5.299203e-02      96      100
#> 93      0.6216234     0.6031775   0.05886406 1.347452e-01      92      100
#> 94      0.6112944     0.4471751   0.05979233 1.339185e-03      97      100
#> 95      0.6207863     0.8295591   0.06849903 4.408690e-02      95      100
#> 96      0.6364422     0.9839856   0.06935133 1.708045e-01      98      100
#> 97      0.5705686     0.6619939   0.02204953 5.844770e-02      97      100
#> 98      0.6209488     1.0281930   0.06560860 3.485706e-01      98      100
#> 99      0.5960874     0.6677002   0.07037318 1.684311e-01      95      100
#> 100     0.6298949     0.6841336   0.04600919 1.407285e-02      96      100
```
