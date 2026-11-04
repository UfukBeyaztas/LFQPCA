lfqpcr_coordinates_from_curves <- function(
    B_curves,
    U_curves,
    chi_B,
    chi_U,
    integration_weights) {
  
  B_curves <- as.matrix(B_curves)
  U_curves <- as.matrix(U_curves)
  chi_B <- as.matrix(chi_B)
  chi_U <- as.matrix(chi_U)
  w <- as.numeric(integration_weights)
  
  if (!all(dim(B_curves) == dim(U_curves))) {
    stop("B_curves and U_curves must have identical dimensions.", call. = FALSE)
  }
  
  if (ncol(B_curves) != length(w)) {
    stop(
      "Curve dimension and integration_weights are incompatible.",
      call. = FALSE
    )
  }
  
  if (nrow(chi_B) != length(w) || nrow(chi_U) != length(w)) {
    stop(
      "Predictor-space bases and integration_weights are incompatible.",
      call. = FALSE
    )
  }
  
  x_B <- if (ncol(chi_B)) {
    B_curves %*% (chi_B * w)
  } else {
    matrix(0, nrow(B_curves), 0)
  }
  
  x_U <- if (ncol(chi_U)) {
    U_curves %*% (chi_U * w)
  } else {
    matrix(0, nrow(U_curves), 0)
  }
  
  list(x_B = unname(x_B), x_U = unname(x_U))
}
