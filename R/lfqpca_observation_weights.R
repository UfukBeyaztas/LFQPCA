lfqpca_observation_weights <- function(Y, subject_id, quadrature_weights) {
  Y <- as.matrix(Y)
  subject_id <- as.integer(subject_id)
  quadrature_weights <- as.numeric(quadrature_weights)
  
  n <- nrow(Y)
  D <- ncol(Y)
  
  if (length(subject_id) != n) {
    stop("length(subject_id) must equal nrow(Y).", call. = FALSE)
  }
  
  if (length(quadrature_weights) != D) {
    stop("length(quadrature_weights) must equal ncol(Y).", call. = FALSE)
  }
  
  if (any(!is.finite(quadrature_weights)) || any(quadrature_weights < 0) ||
      sum(quadrature_weights) <= 0) {
    stop(
      "quadrature_weights must be finite, nonnegative, and have positive sum.",
      call. = FALSE
    )
  }
  
  if (anyNA(subject_id) || any(subject_id < 1)) {
    stop("subject_id must contain positive integer subject labels.",
         call. = FALSE)
  }
  
  I <- max(subject_id)
  visit_counts <- tabulate(subject_id, nbins = I)
  
  if (any(visit_counts <= 0)) {
    stop("Every subject_id must correspond to at least one visit.",
         call. = FALSE)
  }
  
  domain_weights <- matrix(0, nrow = n, ncol = D)
  objective_weights <- matrix(0, nrow = n, ncol = D)
  
  for (r in seq_len(n)) {
    observed <- which(is.finite(Y[r, ]))
    
    if (!length(observed)) {
      stop(
        paste0("Curve ", r, " has no observed functional measurements."),
        call. = FALSE
      )
    }
    
    denominator <- sum(quadrature_weights[observed])
    
    if (!is.finite(denominator) || denominator <= 0) {
      stop(
        paste0(
          "Observed quadrature weights have nonpositive total for curve ",
          r, "."
        ),
        call. = FALSE
      )
    }
    
    omega_r <- quadrature_weights[observed] / denominator
    
    i <- subject_id[r]
    J_i <- visit_counts[i]
    
    domain_weights[r, observed] <- omega_r
    
    objective_weights[r, observed] <- (1 / I) * (1 / J_i) * omega_r
  }
  
  list(
    domain_weights = domain_weights,
    objective_weights = objective_weights,
    visit_counts = visit_counts,
    I = I
  )
}
