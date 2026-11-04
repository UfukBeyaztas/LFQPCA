lfqpcr_orthonormal_space_basis <- function(
    generating_functions,
    integration_weights,
    tol = 1e-10) {
  
  F <- as.matrix(generating_functions)
  w <- as.numeric(integration_weights)
  
  if (nrow(F) != length(w)) {
    stop(
      "Length of integration_weights must equal nrow(generating_functions).",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop(
      "integration_weights must be finite and strictly positive.",
      call. = FALSE
    )
  }
  
  if (!ncol(F)) {
    return(list(
      basis = matrix(0, nrow(F), 0),
      rank = 0,
      singular_values = numeric(0)
    ))
  }
  
  if (any(!is.finite(F))) {
    stop("generating_functions must contain only finite values.", call. = FALSE)
  }
  
  A <- sqrt(w) * F
  sv <- svd(A, nu = min(dim(A)), nv = min(dim(A)))
  
  if (!length(sv$d)) {
    return(list(
      basis = matrix(0, nrow(F), 0),
      rank = 0,
      singular_values = numeric(0)
    ))
  }
  
  threshold <- max(dim(A)) * max(sv$d) * tol
  keep <- which(sv$d > threshold)
  
  if (!length(keep)) {
    return(list(
      basis = matrix(0, nrow(F), 0),
      rank = 0,
      singular_values = sv$d
    ))
  }
  
  V_keep <- sv$v[, keep, drop = FALSE]
  d_keep <- sv$d[keep]
  
  chi <- F %*% V_keep %*%
    diag(1 / d_keep, nrow = length(d_keep), ncol = length(d_keep))
  
  gram <- crossprod(chi, chi * w)
  err <- max(abs(gram - diag(ncol(chi))))
  
  if (!is.finite(err) || err > 1e-7) {
    stop(
      "Failed to construct an L2-orthonormal predictor-space basis.",
      call. = FALSE
    )
  }
  
  list(
    basis = unname(chi),
    rank = ncol(chi),
    singular_values = sv$d,
    orthonormality_error = err
  )
}