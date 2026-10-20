## Resubmission

This is a resubmission. In this version I have:

* Fixed a test failure on Debian (R-devel, clang/flang). The test asserted
  that the optimiser stops early, which is platform-dependent. It now checks
  only the platform-independent behaviour.
* Increased the version to 1.0.1.

## Test environments

* local macOS (aarch64-apple-darwin), R 4.4.1
* win-builder: R-devel, R-release
* Debian (R-devel, clang/flang)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* "Sestak" and "Berggren", flagged as possibly misspelled in DESCRIPTION, are
  the authors of the Sestak-Berggren kinetic model and are spelled correctly.
