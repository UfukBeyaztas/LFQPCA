lfqpca_sim_new_directions <- function(raw, base_space, norms, w,
                                      tol = 1e-10) {
  raw <- as.matrix(raw)
  out <- matrix(0, nrow(raw), ncol(raw))
  space <- as.matrix(base_space)
  for (j in seq_len(ncol(raw))) {
    g <- raw[, j]
    if (ncol(space)) {
      g <- g - lfqpca_sim_project(matrix(g, ncol = 1), space, w, tol)[, 1]
    }
    ng <- sqrt(sum(w * g^2))
    if (!is.finite(ng) || ng <= tol) {
      stop("A new loading direction lies in the existing space.",
           call. = FALSE)
    }
    out[, j] <- norms[j] * g / ng
    space <- cbind(space, out[, j])
  }
  out
}