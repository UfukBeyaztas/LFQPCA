lfqpca_time_basis <- function(time_std, degree = 2) {
  u <- as.numeric(time_std)
  
  degree <- as.integer(degree)
  
  if (degree < 1) {
    stop("'mean_time_degree' must be at least 1.", call. = FALSE)
  }
  
  H <- vapply(0:degree, function(p) u^p, numeric(length(u)))
  
  H <- matrix(H, nrow = length(u), ncol = degree + 1)
  
  colnames(H) <- paste0("u^", 0:degree)
  
  H
}
