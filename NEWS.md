# cumulcalib 0.0.1.9002 (development)

* `summary()` now describes the shape of miscalibration around the maximum
  cumulative calibration error (C*): whether the cumulative error reverses at an
  interior peak (opposite directions on either side of that location) or is
  one-directional (a monotone accumulation), with the direction(s) worded in
  terms of observed vs predicted. Reported only when the standardized maximum
  deviation S* is at least `shape_threshold` (default 1.5), so it stays silent
  for unremarkable processes. Added the `crossover` element to the summary
  object and the `shape_threshold` argument to `summary()`.

# cumulcalib 0.0.1.9001 (development)

* `print()` and `summary()` now report the direction of the maximum cumulative
  calibration error C* (as a worded tag, e.g. "observed benefit < predicted")
  together with the predicted value at its location. Direction is reported only
  for C*, the descriptive metric; the test statistics (S_n, S*, B*) are referred
  to null distributions of absolute deviations and remain unsigned. The
  direction is derived on the fly and is consistent with `plot()`.

# cumulcalib 0.0.1

* Initial CRAN submission.
