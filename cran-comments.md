## Resubmission

This is a resubmission of cumulcalib, updating it from the currently
published CRAN version (0.1.0) to 0.2.0.

This release primarily fixes a reproducibility defect: when predictor values
were tied, results depended on the arbitrary order in which tied observations
happened to appear in the input data, so the same data could yield a different
p-value on every run. Tied predictions are common in practice (tree-based
models such as causal forests routinely produce many repeated values), so this
affected a substantial class of real use cases.

## Summary of changes since 0.1.0

* **Bug fix (behavior change).** Added a `ties` argument to `cumulcalib()` and
  `cumulcalibITE()`, controlling how observations with exactly tied predictor
  values are handled (`"group"` default, `"random"`, `"ignore"`). Under the
  new default, `ties = "group"`, results are fully deterministic regardless of
  input row order; `ties = "ignore"` reproduces the previous, order-dependent
  behavior. Because the default changes how tied data are handled, results for
  such data will differ from 0.1.0; this is intentional and is documented in
  NEWS.md. Ties are only checked for, and acted on, when actually present, so
  there is no overhead for the common case of continuous, untied predictors.

* The ITE methodology paper cited in the Description field was previously an
  arXiv preprint (a CRAN reviewer had noted this DOI would become valid once
  the paper was published). It has now been published, and the citation is
  updated accordingly to Sadatsafavi et al. (2026) <doi:10.1002/sim.70724>.

## Test environments

* Local: Windows 11, R 4.6.0
* GitHub Actions (R-CMD-check workflow): windows-latest (release),
  macOS-latest (release), ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

Thank you for reviewing the package.

Regards,
Mohsen Sadatsafavi
