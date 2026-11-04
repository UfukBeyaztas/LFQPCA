lfqpca_sim_l2_normalize <- function(f, w, norm = 1) {
  f <- as.numeric(f)
  nf <- sqrt(sum(w * f^2))
  if (!is.finite(nf) || nf <= 0) {
    stop("Cannot normalize a zero function.", call. = FALSE)
  }
  norm * f / nf
}