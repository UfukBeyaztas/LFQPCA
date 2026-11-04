lfqpca_continuation_path <- function(tau, reference_tau = 0.5,
                                     fractions = c(0, 0.5, 0.75, 1)) {
  tau <- as.numeric(tau)
  reference_tau <- as.numeric(reference_tau)
  fractions <- sort(unique(as.numeric(fractions)))
  if (length(tau) != 1 || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("'tau' must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (length(reference_tau) != 1 || !is.finite(reference_tau) ||
      reference_tau <= 0 || reference_tau >= 1) {
    stop("'reference_tau' must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!length(fractions) || any(!is.finite(fractions)) ||
      min(fractions) < 0 || max(fractions) > 1) {
    stop("'fractions' must lie in [0, 1].", call. = FALSE)
  }
  if (!any(abs(fractions) < 1e-12)) fractions <- c(0, fractions)
  if (!any(abs(fractions - 1) < 1e-12)) fractions <- c(fractions, 1)
  path <- reference_tau + fractions * (tau - reference_tau)
  path[!duplicated(round(path, 12))]
}
