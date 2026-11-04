lfqpca_build_joint_score_design_subject <- function(
    Y,
    rows_i,
    time_std,
    H,
    mean_part,
    phi0,
    phi1,
    psi,
    objective_weights) {
  
  Y <- as.matrix(Y)
  objective_weights <- as.matrix(objective_weights)
  
  if (!all(dim(objective_weights) == dim(Y))) {
    stop("objective_weights must have the same dimensions as Y.", call. = FALSE)
  }
  
  Kx <- ncol(phi0)
  Ku <- ncol(psi)
  J <- length(rows_i)
  
  N_i <- if (Ku) {
    lfqpca_constraint_null_basis(H[rows_i, , drop = FALSE])
  } else {
    matrix(0, J, 0)
  }
  
  q_i <- ncol(N_i)
  
  response_list <- vector("list", J)
  design_list <- vector("list", J)
  weight_list <- vector("list", J)
  
  for (jj in seq_len(J)) {
    
    r <- rows_i[jj]
    observed <- which(is.finite(Y[r, ]))
    
    if (!length(observed)) next
    
    y_j <- Y[r, observed] - mean_part[r, observed]
    
    Xi_block <-
      phi0[observed, , drop = FALSE] +
      time_std[r] * phi1[observed, , drop = FALSE]
    
    if (Ku && q_i) {
      
      visit_blocks <- lapply(
        seq_len(Ku),
        function(l) {
          outer(psi[observed, l], N_i[jj,])
        }
      )
      
      Z_block <- do.call(cbind, visit_blocks)
      
    } else {
      
      Z_block <- matrix(0, length(observed), 0)
    }
    
    response_list[[jj]] <- y_j
    
    design_list[[jj]] <- cbind(Xi_block, Z_block)
    
    weight_list[[jj]] <- objective_weights[r, observed]
  }
  
  keep <- vapply(response_list, function(z) !is.null(z), logical(1))
  
  if (!any(keep)) {
    stop("A subject has no observed functional values.", call. = FALSE)
  }
  
  list(
    y = unlist(response_list[keep], use.names = FALSE),
    X = do.call(rbind, design_list[keep]),
    case_weights = unlist(weight_list[keep], use.names = FALSE),
    N = N_i,
    q = q_i,
    Kx = Kx,
    Ku = Ku
  )
}
