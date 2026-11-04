lfqpcr_response_time_basis <- function(
    time_std,
    degree,
    minimum_degree = 1) {
  
  u <- as.numeric(time_std)
  degree <- as.integer(degree)
  minimum_degree <- as.integer(minimum_degree)
  
  if (degree < minimum_degree) {
    stop(
      "response_time_degree must be at least the LFQPCA mean_time_degree ",
      "so that span(H) is contained in span(G).",
      call. = FALSE
    )
  }
  
  G <- vapply(0:degree, function(p) u^p, numeric(length(u)))
  
  G <- matrix(G, nrow = length(u), ncol = degree + 1)
  
  colnames(G) <- paste0("G_u^", 0:degree)
  G
}