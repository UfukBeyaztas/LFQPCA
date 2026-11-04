lfqpca_trapezoid_weights <- function(grid, normalize = TRUE) {
  grid <- as.numeric(grid)
  if (length(grid) < 2 || any(!is.finite(grid)) ||
      is.unsorted(grid, strictly = TRUE)) {
    stop("'grid' must be a strictly increasing finite numeric vector.",
         call. = FALSE)
  }
  delta <- diff(grid)
  w <- numeric(length(grid))
  w[1] <- delta[1] / 2
  w[length(grid)] <- delta[length(delta)] / 2
  if (length(grid) > 2) {
    w[2:(length(grid) - 1)] <- (delta[-length(delta)] + delta[-1]) / 2
  }
  if (normalize) w <- w / sum(w)
  w
}
