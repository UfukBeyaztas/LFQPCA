lfqpca_build_function_design <- function(
    B,
    H,
    subject_id,
    time_std,
    xi,
    zeta,
    mask,
    objective_weights) {
  
  B <- as.matrix(B)
  H <- as.matrix(H)
  mask <- as.matrix(mask)
  objective_weights <- as.matrix(objective_weights)
  
  n <- nrow(mask)
  D <- ncol(mask)
  q <- ncol(B)
  M <- ncol(H)
  Kx <- ncol(xi)
  Ku <- ncol(zeta)
  
  if (!all(dim(objective_weights) == c(n, D))) {
    stop(
      "objective_weights must have the same dimensions as mask.",
      call. = FALSE
    )
  }
  
  p <- q * (M + 2 * Kx + Ku)
  n_obs <- sum(mask)
  
  X <- matrix(0, n_obs, p)
  index <- matrix(0, n_obs, 2)
  case_weights <- numeric(n_obs)
  
  colnames(index) <- c("curve", "grid")
  
  position <- 1
  
  for (r in seq_len(n)) {
    observed <- which(mask[r, ])
    nr <- length(observed)
    
    if (!nr) next
    
    rows <- position:(position + nr - 1)
    Br <- B[observed, , drop = FALSE]
    
    blocks <- vector("list", M + 2 * Kx + Ku)
    b <- 1
    
    # Population part, eta(s, u) = sum_m H_m(u) eta_m(s).
    for (m in seq_len(M)) {
      blocks[[b]] <- Br * H[r, m]
      b <- b + 1
    }
    
    # Subject intercept loadings.
    if (Kx) {
      for (k in seq_len(Kx)) {
        blocks[[b]] <- Br * xi[subject_id[r], k]
        b <- b + 1
      }
      
      # Subject slope loadings.
      for (k in seq_len(Kx)) {
        blocks[[b]] <- Br * xi[subject_id[r], k] * time_std[r]
        b <- b + 1
      }
    }
    
    # Visit loadings.
    if (Ku) {
      for (l in seq_len(Ku)) {
        blocks[[b]] <- Br * zeta[r, l]
        b <- b + 1
      }
    }
    
    X[rows, ] <- do.call(cbind, blocks)
    index[rows, ] <- cbind(r, observed)
    
    # Objective weights of the same observed rows.
    case_weights[rows] <- objective_weights[r, observed]
    
    position <- position + nr
  }
  
  list(X = X, index = index, case_weights = case_weights)
}