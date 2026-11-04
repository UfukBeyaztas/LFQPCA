lfqpca_fit_mean_only <- function(
    Y,
    B,
    H,
    tau,
    objective_weights,
    rq_method = "fnb") {
  
  Y <- as.matrix(Y)
  objective_weights <- as.matrix(objective_weights)
  
  n <- nrow(Y)
  D <- ncol(Y)
  q <- ncol(B)
  M <- ncol(H)
  
  if (!all(dim(objective_weights) == c(n, D))) {
    stop("objective_weights must have the same dimensions as Y.", call. = FALSE)
  }
  
  mask <- is.finite(Y)
  n_obs <- sum(mask)
  
  X <- matrix(0, n_obs, q * M)
  y <- numeric(n_obs)
  case_weights <- numeric(n_obs)
  
  position <- 1
  
  for (r in seq_len(n)) {
    
    observed <- which(mask[r, ])
    nr <- length(observed)
    
    if (!nr) next
    
    rows <- position:(position + nr - 1)
    Br <- B[observed, , drop = FALSE]
    
    X[rows,] <- do.call(cbind, lapply(seq_len(M), function(m) Br * H[r, m]))
    
    y[rows] <- Y[r, observed]
    
    case_weights[rows] <- objective_weights[r, observed]
    
    position <- position + nr
  }
  
  beta <- lfqpca_quantile_fit(
    X = X,
    y = y,
    tau = tau,
    case_weights = case_weights,
    method = rq_method
  )
  
  theta_eta <- matrix(beta, nrow = q, ncol = M)
  
  list(theta_eta = theta_eta, eta_functions = B %*% theta_eta)
}