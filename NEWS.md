# cumulcalib 0.2.0

* Added a `ties` argument to `cumulcalib()` and `cumulcalibITE()`, controlling
  how observations with exactly tied predictor values (`p` for `cumulcalib()`,
  `h` for `cumulcalibITE()`) are handled. Options are `"average"` (the new
  default), `"random"`, and `"ignore"`. Ties are only checked for, and acted
  on, when they are actually present, so there is no overhead for the common
  case of continuous, untied predictors.

* **Behavior change:** previously, when predictor values were tied (as
  routinely happens with tree-based models, e.g. causal forests, which often
  produce many repeated predicted values), results depended on the arbitrary
  order in which tied observations appeared in the input data -- the same
  data could yield a different p-value on every run. With the new default,
  `ties = "average"`, results are fully deterministic regardless of input row
  order: for `cumulcalib()`, the outcome is averaged within each tied group of
  predicted risks before computing the cumulative-sum process; for
  `cumulcalibITE()`, this is supported for both the conditional (`p` supplied)
  and marginal (`p = NULL`) approaches. `ties = "ignore"` reproduces the
  previous (order-dependent) behavior, and `ties = "random"` resolves ties via
  random reordering, for users who want an explicit robustness check.

# cumulcalib 0.1.0

* Added `cumulcalibITE()`, extending the cumulative calibration methodology to
  the assessment of moderate calibration of individualized treatment effect
  (ITE) models using data from a randomized trial (Sadatsafavi et al. 2026,
  <doi:10.1002/sim.70724>). A companion vignette demonstrates its use.

* Added `print()` methods for objects returned by `cumulcalib()` and
  `cumulcalibITE()`, and extended `summary()` to also support
  `cumulcalibITE()` output.

* `print()` and `summary()` now report the direction of the maximum cumulative
  calibration error C* (as a worded tag, e.g. "observed benefit < predicted")
  together with the predicted value at its location. Direction is reported only
  for C*, the descriptive metric; the test statistics (S_n, S*, B*) are referred
  to null distributions of absolute deviations and remain unsigned. The
  direction is derived on the fly and is consistent with `plot()`.

* `summary()` now describes the shape of miscalibration around the maximum
  cumulative calibration error (C*): whether the cumulative error reverses at an
  interior peak (opposite directions on either side of that location) or is
  one-directional (a monotone accumulation), with the direction(s) worded in
  terms of observed vs predicted. Reported only when the standardized maximum
  deviation S* is at least `shape_threshold` (default 1.5), so it stays silent
  for unremarkable processes. Added the `crossover` element to the summary
  object and the `shape_threshold` argument to `summary()`.

* **Behavior change:** `plot()` no longer draws significance-threshold lines by
  default; only the statistic line(s) are drawn. Pass
  `stats_config = list(lines = "vh")` to restore the previous default of
  showing both, or `list(lines = "h")` for thresholds only.

* Fixed a bug where the bridge test's significance line was not drawn parallel
  to the line connecting the start and end of the random walk.

* Fixed an error in `summary()` for the `BM2p` method caused by a stray
  `writeLines()` call.

# cumulcalib 0.0.1

* Initial CRAN submission.
