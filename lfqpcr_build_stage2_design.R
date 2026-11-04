lfqpcr_build_stage2_design <- function(
    time_std,
    x_B,
    x_U,
    covariates,
    response_time_degree,
    minimum_time_degree) {
  
  n <- length(time_std)
  
  G <- lfqpcr_response_time_basis(
    time_std = time_std,
    degree = response_time_degree,
    minimum_degree = minimum_time_degree
  )
  
  C <- lfqpcr_prepare_covariates(covariates, n)
  
  x_B <- as.matrix(x_B)
  x_U <- as.matrix(x_U)
  
  if (nrow(x_B) != n || nrow(x_U) != n) {
    stop(
      "Functional-coordinate matrices must have one row per visit.",
      call. = FALSE
    )
  }
  
  if (ncol(x_B)) {
    colnames(x_B) <- paste0("xB_", seq_len(ncol(x_B)))
  }
  
  if (ncol(x_U)) {
    colnames(x_U) <- paste0("xU_", seq_len(ncol(x_U)))
  }
  
  X <- cbind(G, C, x_B, x_U)
  
  list(
    X = X,
    G = G,
    C = C,
    x_B = x_B,
    x_U = x_U,
    block_sizes = c(
      alpha = ncol(G),
      beta = ncol(C),
      B = ncol(x_B),
      U = ncol(x_U)
    )
  )
}