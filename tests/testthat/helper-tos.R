# predict_tos() returns a time-0 baseline row (segment 0) carrying the fitted
# intercept, ahead of one row per schedule segment, so a k-segment schedule
# yields k + 1 rows per variable. Tests that are about the segments themselves
# drop the baseline row; tests about the baseline keep only it.

tos_segments <- function(p) {
  out <- p[p$segment > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

tos_baseline <- function(p) {
  out <- p[p$segment == 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}
