## Resubmission

This is a resubmission of cumulcalib, updating it from the currently
published CRAN version (0.1.0) to 0.2.0.

## Summary of changes since 0.1.0

* The ITE methodology paper cited in the Description field was previously an
  arXiv preprint (a CRAN reviewer had noted this DOI would become valid once
  the paper was published). It has now been published, and the citation is
  updated accordingly to Sadatsafavi et al. (2026) <doi:10.1002/sim.70724>.

* Added a `ties` argument to `cumulcalib()` and `cumulcalibITE()`, controlling
  how observations with exactly tied predictor values are handled (`"group"`
  default, `"random"`, `"ignore"`). Ties are only checked for, and acted on,
  when actually present, so there is no overhead for the common case of
  continuous, untied predictors.

  **Behavior change:** previously, tied predictor values (which routinely
  occur with tree-based models, e.g. causal forests) caused results to depend
  on the arbitrary input row order. With the new default, `ties = "group"`,
  results are fully deterministic regardless of row order. See NEWS.md for
  details.

## Test environments

* Local: Windows 11, R 4.6.0
* GitHub Actions (R-CMD-check workflow): windows-latest (release),
  macOS-latest (release), ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

Thank you for reviewing the package.

Regards,
Mohsen Sadatsafavi
