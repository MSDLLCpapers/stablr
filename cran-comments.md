## Resubmission

This is a resubmission. In this version I have:

* Added references for the implemented methods to the Description field of
  DESCRIPTION in the form authors (year) <doi:...>, and added the
  corresponding references to the help pages of `fit_sb()`,
  `expr_arrhenius_humidity()` and the `flavier_supp` dataset.
* Replaced `F` with `FALSE` in function defaults, and renamed a vector
  element called `T` in the tests.
* Unwrapped the `simulate_kinetic_data()` examples from `\dontrun{}`; they
  run in under 5 seconds.
* `simulate_kinetic_data()` no longer writes to `.GlobalEnv` (it previously
  restored `.Random.seed` on exit). The seed is now only set when the user
  supplies one via the `sdpt` argument, which defaults to `NULL`.
* Increased the version to 1.0.2.

The previous resubmission (1.0.1) fixed a test failure on Debian (R-devel,
clang/flang). The test asserted that the optimiser stops early, which is
platform-dependent. It now checks only the platform-independent behaviour.

## Test environments

* local macOS (aarch64-apple-darwin), R 4.4.1
* win-builder: R-devel, R-release
* Debian (R-devel, clang/flang)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* "Sestak", "Berggren", "Perez", "Maqueda", "Waterman", "Krishnamoorthy" and
  "Mathew", flagged as possibly misspelled in DESCRIPTION, are author names
  of the cited references and are spelled correctly.
