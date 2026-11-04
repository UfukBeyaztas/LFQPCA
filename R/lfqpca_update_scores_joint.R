lfqpca_update_scores_joint <- function(
    Y,
    subject_id,
    time_std,
    H,
    eta_functions,
    phi0,
    phi1,
    psi,
    tau,
    objective_weights,
    rq_method = "br") {
  
  I <- max(subject_id)
  n <- nrow(Y)
  Kx <- ncol(phi0)
  Ku <- ncol(psi)
  
  xi <- matrix(0, I, Kx)
  zeta <- matrix(0, n, Ku)
  
  mean_part <- lfqpca_mean_surface(H, eta_functions)
  
  for (i in seq_len(I)) {
    
    rows_i <- which(subject_id == i)
    
    design <- lfqpca_build_joint_score_design_subject(
      Y = Y,
      rows_i = rows_i,
      time_std = time_std,
      H = H,
      mean_part = mean_part,
      phi0 = phi0,
      phi1 = phi1,
      psi = psi,
      objective_weights = objective_weights
    )
    
    p_i <- ncol(design$X)
    
    if (length(design$y) <= p_i) {
      stop(
        "Subject ", i,
        " has too few observed values for the joint score block: ",
        length(design$y),
        " observations and ",
        p_i,
        " coefficients.",
        call. = FALSE
      )
    }
    
    beta <- lfqpca_quantile_fit(
      X = design$X,
      y = design$y,
      tau = tau,
      case_weights = design$case_weights,
      method = rq_method
    )
    
    if (Kx) {
      xi[i, ] <- beta[seq_len(Kx)]
    }
    
    if (Ku && design$q) {
      
      a_index <- Kx + seq_len(design$q * Ku)
      
      A_i <- matrix(beta[a_index], nrow = design$q, ncol = Ku)
      
      zeta[rows_i,] <- design$N %*% A_i
    }
  }
  
  list(xi = xi, zeta = zeta)
}
