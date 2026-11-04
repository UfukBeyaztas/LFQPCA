lfqpca_weighted_gram <- function(F, weights) {
  F <- as.matrix(F)
  weights <- as.numeric(weights)
  if (nrow(F) != length(weights)) {
    stop("Length of weights must equal nrow(F).", call. = FALSE)
  }
  crossprod(F, F * weights)
}
