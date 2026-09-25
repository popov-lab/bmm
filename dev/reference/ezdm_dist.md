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
#> 1   0.6111635 0.05418579      95      100
#> 2   0.6086589 0.04968604      91      100
#> 3   0.6043066 0.06557410      92      100
#> 4   0.6337035 0.04187681      98      100
#> 5   0.6628371 0.06905988      98      100
#> 6   0.6306668 0.06789205      92      100
#> 7   0.6683561 0.09017106      97      100
#> 8   0.6390170 0.05048177      92      100
#> 9   0.6594366 0.06443427      97      100
#> 10  0.6388977 0.05186053      96      100
#> 11  0.6466595 0.05898557      99      100
#> 12  0.6441727 0.06067506      96      100
#> 13  0.5960057 0.03501029      95      100
#> 14  0.6566931 0.07312851      97      100
#> 15  0.6480473 0.05061371      99      100
#> 16  0.6338738 0.06588902      97      100
#> 17  0.6749657 0.07643939      93      100
#> 18  0.6315851 0.04292302      94      100
#> 19  0.6299139 0.04702333      95      100
#> 20  0.6415729 0.06164678      97      100
#> 21  0.6569962 0.04665327      94      100
#> 22  0.6475104 0.06588247      96      100
#> 23  0.6607518 0.09228922      95      100
#> 24  0.6486437 0.06127875      96      100
#> 25  0.6397690 0.06979176      92      100
#> 26  0.6414065 0.07433980      96      100
#> 27  0.6021627 0.05596795      95      100
#> 28  0.6385743 0.05330900      95      100
#> 29  0.6524296 0.07905578      94      100
#> 30  0.6554120 0.06261139      92      100
#> 31  0.6469188 0.07973757      99      100
#> 32  0.6369412 0.05716800      96      100
#> 33  0.6245654 0.04790737      93      100
#> 34  0.6647938 0.07865629      97      100
#> 35  0.6356469 0.06488684      96      100
#> 36  0.6579946 0.06010943      94      100
#> 37  0.6288020 0.04890998      99      100
#> 38  0.6955779 0.08011255      95      100
#> 39  0.6633387 0.03779677      95      100
#> 40  0.6497766 0.05726640      95      100
#> 41  0.6020050 0.03912827      95      100
#> 42  0.6512410 0.07179975      94      100
#> 43  0.6187174 0.04709756      94      100
#> 44  0.6831167 0.07639139      98      100
#> 45  0.6561083 0.05187717      94      100
#> 46  0.6456383 0.07429988      93      100
#> 47  0.6413068 0.06202381      94      100
#> 48  0.6439319 0.08259685      94      100
#> 49  0.6545809 0.05948838      91      100
#> 50  0.6902880 0.09919671      93      100
#> 51  0.6490158 0.04714586      92      100
#> 52  0.6391610 0.07818092      97      100
#> 53  0.6669382 0.07987881      95      100
#> 54  0.6470221 0.07570689      99      100
#> 55  0.6697892 0.07322692      96      100
#> 56  0.6629982 0.05985320      96      100
#> 57  0.6085790 0.04388921      97      100
#> 58  0.6733613 0.08550546      97      100
#> 59  0.6025789 0.04238616      93      100
#> 60  0.6782462 0.08774496      90      100
#> 61  0.5853326 0.03055388      97      100
#> 62  0.6808352 0.06220229      96      100
#> 63  0.6804366 0.06035126      96      100
#> 64  0.6473923 0.05182805      96      100
#> 65  0.6464125 0.06720487      93      100
#> 66  0.6808136 0.08828020      93      100
#> 67  0.6593359 0.09232922      96      100
#> 68  0.6183544 0.06604877      94      100
#> 69  0.6550676 0.06929830      94      100
#> 70  0.6021303 0.06651822      94      100
#> 71  0.6212957 0.05692183      91      100
#> 72  0.6120274 0.03727352      89      100
#> 73  0.6525679 0.05112530      90      100
#> 74  0.6318721 0.03595612      95      100
#> 75  0.6268748 0.03405663      92      100
#> 76  0.6122607 0.04591200      91      100
#> 77  0.7305996 0.12676027      91      100
#> 78  0.6501535 0.07278958      95      100
#> 79  0.6840221 0.07520687      97      100
#> 80  0.6258505 0.07007770      98      100
#> 81  0.6213556 0.04796828      96      100
#> 82  0.6270388 0.04742668      99      100
#> 83  0.6663180 0.07532801      90      100
#> 84  0.6386077 0.06691333      94      100
#> 85  0.6302745 0.04252136      93      100
#> 86  0.5907167 0.04764884      97      100
#> 87  0.6222910 0.05049253      93      100
#> 88  0.6280339 0.05476211      95      100
#> 89  0.6579823 0.07378763      93      100
#> 90  0.6320513 0.06339803      91      100
#> 91  0.6010266 0.04239233      94      100
#> 92  0.6229369 0.04027990      95      100
#> 93  0.6783554 0.09407671      91      100
#> 94  0.6453377 0.08696251      96      100
#> 95  0.6168118 0.05827374      97      100
#> 96  0.6156420 0.04072509      99      100
#> 97  0.6278741 0.04883660      96      100
#> 98  0.6728986 0.07497164      97      100
#> 99  0.6355790 0.04843782      96      100
#> 100 0.6425947 0.03357066      97      100
rezdm(
  n = 100, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3,
  zr = 0.55, version = "4par"
)
#>     mean_rt_upper mean_rt_lower var_rt_upper var_rt_lower n_upper n_trials
#> 1       0.6253226     0.5823772   0.07911455 5.282973e-03      98      100
#> 2       0.6171880     0.5839276   0.05528031 5.966834e-02      95      100
#> 3       0.6485196     0.7139917   0.08391141 2.858171e-01      97      100
#> 4       0.5829925     1.0441456   0.04721922 1.909124e-01      96      100
#> 5       0.6273172     0.5092828   0.07401407 2.321978e-02      93      100
#> 6       0.6249412     0.6171310   0.07956924 2.244483e-03      97      100
#> 7       0.6114124     0.8940065   0.05081401 2.970588e-03      98      100
#> 8       0.5880939     0.8036448   0.06025104 1.227541e-01      96      100
#> 9       0.6447928     0.6301377   0.06424441 4.105989e-02      96      100
#> 10      0.6082380     0.6723502   0.03238539 1.181740e-03      97      100
#> 11      0.5904479     0.7559707   0.05169755 2.350302e-02      93      100
#> 12      0.6196726     0.6583712   0.07667082 2.837501e-02      97      100
#> 13      0.5899915            NA   0.05189302           NA      99      100
#> 14      0.6181806     0.5681570   0.06644203 1.374754e-02      96      100
#> 15      0.6175320     0.7707452   0.07191115 4.200276e-02      95      100
#> 16      0.5931011     0.7873282   0.03816982 1.337478e-04      96      100
#> 17      0.5673028     0.5018691   0.04601163 2.907999e-02      98      100
#> 18      0.6278858     0.6996673   0.07084437 4.308046e-02      93      100
#> 19      0.5999634     0.6190288   0.04638203 1.800411e-03      94      100
#> 20      0.5716147            NA   0.04551777           NA      99      100
#> 21      0.5933489     0.6675141   0.03567301 1.414062e-04      95      100
#> 22      0.6270457     0.6874140   0.07178717 4.104850e-02      94      100
#> 23      0.6327607     0.6135700   0.08203474 6.503511e-03      96      100
#> 24      0.6026259     0.6489808   0.04942633 4.013743e-03      95      100
#> 25      0.6261162     0.4786807   0.06545772 1.936777e-03      97      100
#> 26      0.6402308            NA   0.07972958           NA      99      100
#> 27      0.5963508            NA   0.05022043           NA      99      100
#> 28      0.6202917     0.6139906   0.04721892 1.773589e-01      97      100
#> 29      0.6292854     0.8031477   0.07546246 1.045227e-01      98      100
#> 30      0.6201752     0.6355989   0.06832560 6.812979e-03      98      100
#> 31      0.6184427     0.5993869   0.04891331 1.497810e-10      98      100
#> 32      0.6200807     0.3899048   0.04971070 1.130956e-02      98      100
#> 33      0.5821060     0.6999998   0.04120795 1.410072e-01      94      100
#> 34      0.5808450            NA   0.04444016           NA      99      100
#> 35      0.5996229     0.7369012   0.05796532 1.158442e-01      97      100
#> 36      0.5896469     0.8406242   0.05602702 3.661818e-01      98      100
#> 37      0.5749020     0.6998572   0.05193127 8.825138e-02      97      100
#> 38      0.6089534     0.7608306   0.05016780 1.185885e-01      97      100
#> 39      0.6091953     0.6913986   0.05791185 1.166646e-01      98      100
#> 40      0.6566303     0.7872351   0.09057143 1.427873e-01      97      100
#> 41      0.6259888     0.5853079   0.07171781 6.039863e-04      97      100
#> 42      0.6524525            NA   0.09050340           NA     100      100
#> 43      0.5883644            NA   0.04354750           NA     100      100
#> 44      0.5955397     0.5841933   0.05831974 5.377575e-03      96      100
#> 45      0.6346774     0.5960504   0.06707157 9.774911e-03      95      100
#> 46      0.6210293     0.7869111   0.05421985 1.308940e-01      95      100
#> 47      0.5643313     0.7206490   0.04127384 4.529126e-02      95      100
#> 48      0.6290337     0.8056297   0.05998544 3.975849e-04      95      100
#> 49      0.5864586     0.5344239   0.06338933 2.497840e-03      94      100
#> 50      0.6017033     0.5552666   0.04908010 7.124369e-03      95      100
#> 51      0.6382134     0.6836493   0.08997426 6.404295e-06      97      100
#> 52      0.5878631     0.7653400   0.03712318 1.468800e-04      96      100
#> 53      0.6099902     0.9212601   0.03800446 2.197637e-01      95      100
#> 54      0.6127636     0.5738889   0.08001620 2.442924e-03      94      100
#> 55      0.5898540            NA   0.04125372           NA      99      100
#> 56      0.6565393     0.5712908   0.05915474 6.186089e-02      98      100
#> 57      0.6036652     0.7524118   0.07486943 1.142665e-03      95      100
#> 58      0.6205222            NA   0.07666748           NA      99      100
#> 59      0.6186504     0.7601851   0.06530945 4.286781e-03      96      100
#> 60      0.5846995     0.6870436   0.02863498 6.128926e-02      96      100
#> 61      0.6007758     0.7520200   0.05632263 3.525548e-04      95      100
#> 62      0.6146287     0.6095731   0.06517262 9.433624e-04      96      100
#> 63      0.6079881     0.7482211   0.03627917 2.192291e-01      95      100
#> 64      0.5717904     0.6237479   0.02249992 1.621140e-02      98      100
#> 65      0.6201225            NA   0.06193363           NA      99      100
#> 66      0.6113292     0.7473014   0.05148123 1.414600e-01      97      100
#> 67      0.6075502     0.4635556   0.04719018 2.600340e-03      96      100
#> 68      0.6120317     0.5996071   0.05008479 6.928488e-03      96      100
#> 69      0.6371032     0.6414161   0.08146846 6.582110e-02      93      100
#> 70      0.6095704     0.6700533   0.06765039 1.487197e-02      96      100
#> 71      0.5794014     0.5437786   0.02781296 1.159084e-03      98      100
#> 72      0.5805540     0.8679745   0.04403614 2.297599e-01      95      100
#> 73      0.5882489     0.6572810   0.03381068 4.042184e-03      95      100
#> 74      0.6337591     0.6219328   0.08497394 3.914663e-02      95      100
#> 75      0.6519035     0.6074190   0.09859535 2.255194e-03      96      100
#> 76      0.6020498            NA   0.06652119           NA     100      100
#> 77      0.5973022     1.0329858   0.04875611 1.740835e-01      97      100
#> 78      0.5637603     0.7237134   0.03017398 7.604240e-03      93      100
#> 79      0.6162830     0.7907933   0.06911791 5.392265e-02      95      100
#> 80      0.6616871     0.6805740   0.09730298 7.029844e-02      96      100
#> 81      0.6310367     0.5784882   0.07253849 4.318314e-02      95      100
#> 82      0.6056327     0.6992810   0.03734506 6.048835e-04      98      100
#> 83      0.6501273     0.8648440   0.08690550 2.984268e-01      98      100
#> 84      0.5939566            NA   0.03903370           NA      99      100
#> 85      0.6073185     0.9779838   0.03715447 2.959167e-02      96      100
#> 86      0.6047717     0.4530562   0.03138902 8.292242e-03      95      100
#> 87      0.5719962     0.6500034   0.04307318 6.181304e-02      95      100
#> 88      0.6053661     0.7520025   0.03703732 1.544416e-03      96      100
#> 89      0.6325180            NA   0.06526068           NA      99      100
#> 90      0.6139106     0.4433208   0.04435218 1.282726e-02      95      100
#> 91      0.6117540     0.6004376   0.07323417 1.341352e-02      94      100
#> 92      0.6574545     0.6501949   0.06220291 1.239173e-03      97      100
#> 93      0.5754018     0.5888569   0.06561530 2.414883e-03      98      100
#> 94      0.5895091     0.6046400   0.03646934 1.080198e-01      97      100
#> 95      0.6030274     0.6675923   0.04760016 7.117333e-02      93      100
#> 96      0.6019021     0.7496666   0.05104677 5.595125e-02      92      100
#> 97      0.6342686     0.4714468   0.06731189 1.101831e-02      97      100
#> 98      0.6026688            NA   0.04723144           NA     100      100
#> 99      0.6025780     0.8252478   0.06193917 8.912872e-03      98      100
#> 100     0.6258862     0.7871865   0.04819153 9.731447e-02      96      100
```
