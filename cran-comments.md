## Resubmission

This is a resubmission of cumulcalib, updating it from the previously
published CRAN version (0.0.1) to 0.1.0.

## Summary of changes since 0.0.1

* Added `cumulcalibITE()`, extending the cumulative calibration methodology
  to the assessment of moderate calibration of individualized treatment
  effect (ITE) models using data from a randomized trial. This is described
  in a new preprint, added to the Description field as
  Sadatsafavi et al. (2025) <doi:10.48550/arXiv.2512.08140>.
* Added `print()` methods for objects returned by `cumulcalib()` and
  `cumulcalibITE()`, and extended the existing `summary()` method to also
  support `cumulcalibITE()` output.
* Added a second vignette demonstrating `cumulcalibITE()` on a case study.
* Added a testthat test suite (previously the package had a placeholder test
  only).
* Minor documentation and bug fixes (see NEWS.md).

## Test environments

* Local: Windows 11, R 4.5.2
* GitHub Actions (R-CMD-check workflow): windows-latest (release),
  macOS-latest (release), ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

Thank you for reviewing the package.

Regards,
Mohsen Sadatsafavi
