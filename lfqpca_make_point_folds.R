lfqpca_make_point_folds <- function(Y, n_folds = 3) {
  Y <- as.matrix(Y)
  n_folds <- as.integer(n_folds)
  if (n_folds < 2) stop("n_folds must be at least 2.", call. = FALSE)
  
  fold_id <- matrix(NA_integer_, nrow(Y), ncol(Y))
  for (r in seq_len(nrow(Y))) {
    observed <- which(is.finite(Y[r, ]))
    if (length(observed) < n_folds + 5) {
      stop(
        "Curve ", r, " has too few observed points for ", n_folds,
        "-fold pointwise CV.",
        call. = FALSE
      )
    }
    assignments <- sample(rep(seq_len(n_folds), length.out = length(observed)))
    fold_id[r, observed] <- assignments
  }
  fold_id
}