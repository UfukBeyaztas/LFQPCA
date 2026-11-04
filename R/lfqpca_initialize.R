lfqpca_initialize <- function(
    Y,
    subject_id,
    time_std,
    grid,
    B,
    H,
    Kx,
    Ku,
    tau,
    weights,
    objective_weights,
    rq_method = "fnb",
    random_start = FALSE) {
  
  n <- nrow(Y)
  D <- ncol(Y)
  I <- max(subject_id)
  
  mean_fit <- lfqpca_fit_mean_only(
    Y = Y,
    B = B,
    H = H,
    tau = tau,
    objective_weights = objective_weights,
    rq_method = rq_method
  )
  eta_functions <- mean_fit$eta_functions
  mean_part <- lfqpca_mean_surface(H, eta_functions)
  residual <- Y - mean_part
  
  if (random_start) {
    xi <- matrix(rnorm(I * Kx), I, Kx)
    q <- ncol(B)
    phi0 <- B %*% matrix(rnorm(q * Kx), q, Kx)
    phi1 <- B %*% matrix(rnorm(q * Kx), q, Kx)
    pair <- lfqpca_weighted_orthonormalize(rbind(phi0, phi1), rep(weights, 2))
    phi0 <- pair[seq_len(D), , drop = FALSE]
    phi1 <- pair[D + seq_len(D), , drop = FALSE]
    
    if (Ku) {
      psi <- B %*% matrix(rnorm(q * Ku), q, Ku)
      psi <- lfqpca_weighted_orthonormalize(psi, weights)
      zeta <- matrix(rnorm(n * Ku), n, Ku)
      zeta <- lfqpca_project_visit_scores_to_constraints(zeta, subject_id, H)
    } else {
      psi <- matrix(0, D, 0)
      zeta <- matrix(0, n, 0)
    }
  } else {
    intercept_curves <- matrix(NA_real_, I, D)
    slope_curves <- matrix(NA_real_, I, D)
    
    for (i in seq_len(I)) {
      rows_i <- which(subject_id == i)
      ui <- time_std[rows_i]
      
      for (d in seq_len(D)) {
        good <- is.finite(residual[rows_i, d])
        if (sum(good) >= 2 && sd(ui[good]) > 0) {
          fit <- lm.fit(x = cbind(1, ui[good]), y = residual[rows_i, d][good])
          co <- fit$coefficients
          if (all(is.finite(co))) {
            intercept_curves[i, d] <- co[1]
            slope_curves[i, d] <- co[2]
          }
        }
      }
      
      intercept_curves[i, ] <- lfqpca_fill_curve_na(intercept_curves[i, ], grid)
      slope_curves[i, ] <- lfqpca_fill_curve_na(slope_curves[i, ], grid)
    }
    
    subject_matrix <- cbind(intercept_curves, slope_curves)
    subject_matrix <- scale(subject_matrix, center = TRUE, scale = FALSE)
    sv_x <- svd(subject_matrix, nu = min(Kx, nrow(subject_matrix)), nv = Kx)
    
    xi <- sv_x$u[, seq_len(Kx), drop = FALSE] %*%
      diag(sv_x$d[seq_len(Kx)], Kx, Kx)
    pair <- sv_x$v[, seq_len(Kx), drop = FALSE]
    phi0 <- pair[seq_len(D), , drop = FALSE]
    phi1 <- pair[D + seq_len(D), , drop = FALSE]
    phi0 <- B %*% qr.solve(B, phi0)
    phi1 <- B %*% qr.solve(B, phi1)
    
    rot_x <- lfqpca_absorb_and_rotate(
      eta_functions = eta_functions,
      xi = xi,
      phi0 = phi0,
      phi1 = phi1,
      zeta = matrix(0, n, 0),
      psi = matrix(0, D, 0),
      weights = weights
    )
    eta_functions <- rot_x$eta_functions
    xi <- rot_x$xi
    phi0 <- rot_x$phi0
    phi1 <- rot_x$phi1
    
    if (Ku) {
      residual2 <- Y - lfqpca_mean_surface(H, eta_functions) -
        lfqpca_subject_contribution(subject_id, time_std, xi, phi0, phi1)
      residual2_filled <- t(apply(
        residual2,
        1,
        lfqpca_fill_curve_na,
        grid = grid
      ))
      
      sv_u <- svd(residual2_filled, nu = min(Ku, nrow(residual2_filled)),
                  nv = Ku)
      zeta <- sv_u$u[, seq_len(Ku), drop = FALSE] %*%
        diag(sv_u$d[seq_len(Ku)], Ku, Ku)
      psi <- sv_u$v[, seq_len(Ku), drop = FALSE]
      psi <- B %*% qr.solve(B, psi)
      zeta <- lfqpca_project_visit_scores_to_constraints(zeta, subject_id, H)
      
      rot_u <- lfqpca_rotate_block(
        scores = zeta,
        loadings = psi,
        weights = weights,
        center_scores = FALSE
      )
      zeta <- rot_u$scores
      psi <- rot_u$loadings
    } else {
      psi <- matrix(0, D, 0)
      zeta <- matrix(0, n, 0)
    }
  }
  
  list(
    eta_functions = eta_functions,
    xi = xi,
    phi0 = phi0,
    phi1 = phi1,
    zeta = zeta,
    psi = psi
  )
}
