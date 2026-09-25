  // Generated, not typed. Do not edit by hand.
  //
  // Series of the two functions whose derivatives give the decision-time
  // cumulants at small drift (see ezdm_cumulants.stan for how they are used):
  //   G(x) = log(sinh(sqrt x) / sqrt x) = sum_j a_j x^j
  //   C(y) = log cosh(sqrt y)           = sum_j a_j (4^j - 1) y^j
  // ezdm_log_sinhc_dn is G^(n) and ezdm_log_cosh_dn is C^(n), n = 1..4, each a
  // Horner polynomial of the 16 terms a_j j! / (j - n)!. The coefficients are
  // formed in exact rational arithmetic and rounded once, because a wrong one
  // makes the truncated series diverge instead of failing loudly. At x <= 0.49
  // the truncation error is below 1e-13. The a_j are .EZDM_LOG_SINHC_COEF in
  // R/distributions.R, and a test in tests/testthat/test-model_ezdm.R
  // recomputes every literal below from them.

  real ezdm_log_sinhc_d1(real x) {
    real acc = -1.2336844022586037e-16;
    acc = fma(acc, x, 1.2175977014591684e-15);
    acc = fma(acc, x, -1.2017207666653852e-14);
    acc = fma(acc, x, 1.1860508700116827e-13);
    acc = fma(acc, x, -1.1705853409912441e-12);
    acc = fma(acc, x, 1.1553216299501312e-11);
    acc = fma(acc, x, -1.1402575602296092e-10);
    acc = fma(acc, x, 1.1253923258404497e-09);
    acc = fma(acc, x, -1.1107304394989839e-08);
    acc = fma(acc, x, 1.0962973925936889e-07);
    acc = fma(acc, x, -1.0822021404031986e-06);
    acc = fma(acc, x, 1.0688899577788467e-05);
    acc = fma(acc, x, -0.00010582010582010582);
    acc = fma(acc, x, 0.0010582010582010583);
    acc = fma(acc, x, -0.011111111111111112);
    acc = fma(acc, x, 0.16666666666666666);
    return acc;
  }

  real ezdm_log_sinhc_d2(real x) {
    real acc = -1.8505266033879055e-15;
    acc = fma(acc, x, 1.704636782042836e-14);
    acc = fma(acc, x, -1.5622369966650008e-13);
    acc = fma(acc, x, 1.4232610440140192e-12);
    acc = fma(acc, x, -1.2876438750903686e-11);
    acc = fma(acc, x, 1.1553216299501312e-10);
    acc = fma(acc, x, -1.0262318042066481e-09);
    acc = fma(acc, x, 9.003138606723598e-09);
    acc = fma(acc, x, -7.775113076492888e-08);
    acc = fma(acc, x, 6.577784355562133e-07);
    acc = fma(acc, x, -5.411010702015993e-06);
    acc = fma(acc, x, 4.275559831115387e-05);
    acc = fma(acc, x, -0.00031746031746031746);
    acc = fma(acc, x, 0.0021164021164021165);
    acc = fma(acc, x, -0.011111111111111112);
    return acc;
  }

  real ezdm_log_sinhc_d3(real x) {
    real acc = -2.5907372447430678e-14;
    acc = fma(acc, x, 2.2160278166556864e-13);
    acc = fma(acc, x, -1.874684395998001e-12);
    acc = fma(acc, x, 1.565587148415421e-11);
    acc = fma(acc, x, -1.2876438750903685e-10);
    acc = fma(acc, x, 1.039789466955118e-09);
    acc = fma(acc, x, -8.209854433653185e-09);
    acc = fma(acc, x, 6.302197024706518e-08);
    acc = fma(acc, x, -4.6650678458957325e-07);
    acc = fma(acc, x, 3.2888921777810667e-06);
    acc = fma(acc, x, -2.1644042808063973e-05);
    acc = fma(acc, x, 0.0001282667949334616);
    acc = fma(acc, x, -0.0006349206349206349);
    acc = fma(acc, x, 0.0021164021164021165);
    return acc;
  }

  real ezdm_log_sinhc_d4(real x) {
    real acc = -3.367958418165988e-13;
    acc = fma(acc, x, 2.659233379986824e-12);
    acc = fma(acc, x, -2.0621528355978012e-11);
    acc = fma(acc, x, 1.5655871484154212e-10);
    acc = fma(acc, x, -1.1588794875813318e-09);
    acc = fma(acc, x, 8.318315735640945e-09);
    acc = fma(acc, x, -5.74689810355723e-08);
    acc = fma(acc, x, 3.781318214823911e-07);
    acc = fma(acc, x, -2.3325339229478664e-06);
    acc = fma(acc, x, 1.3155568711124267e-05);
    acc = fma(acc, x, -6.493212842419192e-05);
    acc = fma(acc, x, 0.0002565335898669232);
    acc = fma(acc, x, -0.0006349206349206349);
    return acc;
  }

  real ezdm_log_cosh_d1(real x) {
    real acc = -5.298634160052327e-07;
    acc = fma(acc, x, 1.3073855756453773e-06);
    acc = fma(acc, x, -3.2258446078277153e-06);
    acc = fma(acc, x, 7.959452534664482e-06);
    acc = fma(acc, x, -1.9639161941658416e-05);
    acc = fma(acc, x, 4.8457689784647255e-05);
    acc = fma(acc, x, -0.00011956455712177624);
    acc = fma(acc, x, 0.000295013720472793);
    acc = fma(acc, x, -0.0007279171935256592);
    acc = fma(acc, x, 0.0017960640182862406);
    acc = fma(acc, x, -0.004431617764951099);
    acc = fma(acc, x, 0.010934744268077601);
    acc = fma(acc, x, -0.026984126984126985);
    acc = fma(acc, x, 0.06666666666666667);
    acc = fma(acc, x, -0.16666666666666666);
    acc = fma(acc, x, 0.5);
    return acc;
  }

  real ezdm_log_cosh_d2(real x) {
    real acc = -7.94795124007849e-06;
    acc = fma(acc, x, 1.8303398059035282e-05);
    acc = fma(acc, x, -4.19359799017603e-05);
    acc = fma(acc, x, 9.551343041597378e-05);
    acc = fma(acc, x, -0.0002160307813582426);
    acc = fma(acc, x, 0.00048457689784647253);
    acc = fma(acc, x, -0.001076081014095986);
    acc = fma(acc, x, 0.002360109763782344);
    acc = fma(acc, x, -0.005095420354679614);
    acc = fma(acc, x, 0.010776384109717443);
    acc = fma(acc, x, -0.022158088824755492);
    acc = fma(acc, x, 0.043738977072310406);
    acc = fma(acc, x, -0.08095238095238096);
    acc = fma(acc, x, 0.13333333333333333);
    acc = fma(acc, x, -0.16666666666666666);
    return acc;
  }

  real ezdm_log_cosh_d3(real x) {
    real acc = -0.00011127131736109887;
    acc = fma(acc, x, 0.00023794417476745867);
    acc = fma(acc, x, -0.0005032317588211236);
    acc = fma(acc, x, 0.0010506477345757117);
    acc = fma(acc, x, -0.002160307813582426);
    acc = fma(acc, x, 0.004361192080618252);
    acc = fma(acc, x, -0.008608648112767889);
    acc = fma(acc, x, 0.016520768346476408);
    acc = fma(acc, x, -0.030572522128077685);
    acc = fma(acc, x, 0.05388192054858722);
    acc = fma(acc, x, -0.08863235529902197);
    acc = fma(acc, x, 0.1312169312169312);
    acc = fma(acc, x, -0.1619047619047619);
    acc = fma(acc, x, 0.13333333333333333);
    return acc;
  }

  real ezdm_log_cosh_d4(real x) {
    real acc = -0.0014465271256942852;
    acc = fma(acc, x, 0.002855330097209504);
    acc = fma(acc, x, -0.005535549347032359);
    acc = fma(acc, x, 0.010506477345757118);
    acc = fma(acc, x, -0.019442770322241832);
    acc = fma(acc, x, 0.03488953664494602);
    acc = fma(acc, x, -0.06026053678937523);
    acc = fma(acc, x, 0.09912461007885845);
    acc = fma(acc, x, -0.15286261064038842);
    acc = fma(acc, x, 0.21552768219434887);
    acc = fma(acc, x, -0.2658970658970659);
    acc = fma(acc, x, 0.2624338624338624);
    acc = fma(acc, x, -0.1619047619047619);
    return acc;
  }
