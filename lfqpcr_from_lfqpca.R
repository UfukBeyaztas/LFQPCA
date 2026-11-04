lfqpcr_from_lfqpca <- function(
    lfqpca_fit,
    V,
    covariates = NULL,
    tau_R = 0.5,
    response_time_degree = NULL,
    rq_method = c("fnb", "br", "fn"),
    space_tol = 1e-10,
    rank_tol = 1e-7) {
  
  rq_method <- match.arg(rq_method)
  lfqpcr_validate_lfqpca_fit(lfqpca_fit)
  
  V <- as.numeric(V)
  n <- length(lfqpca_fit$subject_id)
  
  if (length(V) != n) {
    stop(
      "length(V) must equal the number of visits in lfqpca_fit.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(tau_R) || length(tau_R) != 1 || !is.finite(tau_R) ||
      tau_R <= 0 ||
      tau_R >= 1) {
    stop("'tau_R' must be strictly between 0 and 1.", call. = FALSE)
  }
  
  if (is.null(response_time_degree)) {
    response_time_degree <- lfqpca_fit$mean_time_degree
  }
  
  response_time_degree <- as.integer(response_time_degree)
  
  components <- lfqpcr_reconstruct_components(lfqpca_fit)
  
  spaces <- lfqpcr_build_predictor_spaces(
    lfqpca_fit = lfqpca_fit,
    tol = space_tol
  )
  
  coords <- lfqpcr_coordinates_from_curves(
    B_curves = components$B,
    U_curves = components$U,
    chi_B = spaces$chi_B,
    chi_U = spaces$chi_U,
    integration_weights = spaces$integration_weights
  )
  
  design <- lfqpcr_build_stage2_design(
    time_std = lfqpca_fit$time_std,
    x_B = coords$x_B,
    x_U = coords$x_U,
    covariates = covariates,
    response_time_degree = response_time_degree,
    minimum_time_degree = lfqpca_fit$mean_time_degree
  )
  
  scalar_covariates_finite <- if (ncol(design$C)) {
    apply(design$C, 1, function(z) all(is.finite(z)))
  } else {
    rep(TRUE, n)
  }
  
  eligible <-
    is.finite(V) & scalar_covariates_finite &
    apply(design$G, 1, function(z) all(is.finite(z))) &
    apply(design$x_B, 1, function(z) all(is.finite(z))) &
    apply(design$x_U, 1, function(z) all(is.finite(z)))
  
  weight_info <- lfqpcr_subject_weights(
    subject_id = lfqpca_fit$subject_id,
    eligible = eligible
  )
  
  stage2 <- lfqpcr_quantile_fit(
    X = design$X,
    y = V,
    tau_R = tau_R,
    case_weights = weight_info$weights,
    rq_method = rq_method,
    rank_tol = rank_tol
  )
  
  coef_blocks <- lfqpcr_split_coefficients(
    coefficients = stage2$coefficients,
    block_sizes = design$block_sizes
  )
  
  gamma_B <- if (spaces$r_B) {
    as.numeric(spaces$chi_B %*% coef_blocks$a_B)
  } else {
    numeric(length(lfqpca_fit$grid))
  }
  
  gamma_U <- if (spaces$r_U) {
    as.numeric(spaces$chi_U %*% coef_blocks$a_U)
  } else {
    numeric(length(lfqpca_fit$grid))
  }
  
  fitted <- stage2$fitted
  residuals <- V - fitted
  
  objective <- lfqpcr_weighted_objective(
    residual = residuals,
    tau_R = tau_R,
    weights = weight_info$weights
  )
  
  w <- spaces$integration_weights
  
  A0 <- if (spaces$r_B) {
    crossprod(spaces$chi_B, lfqpca_fit$phi0 * w)
  } else {
    matrix(0, 0, ncol(lfqpca_fit$phi0))
  }
  
  A1 <- if (spaces$r_B) {
    crossprod(spaces$chi_B, lfqpca_fit$phi1 * w)
  } else {
    matrix(0, 0, ncol(lfqpca_fit$phi1))
  }
  
  AU <- if (spaces$r_U) {
    crossprod(spaces$chi_U, lfqpca_fit$psi * w)
  } else {
    matrix(0, 0, ncol(lfqpca_fit$psi))
  }
  
  theta0 <- if (length(coef_blocks$a_B)) {
    as.numeric(t(A0) %*% coef_blocks$a_B)
  } else {
    numeric(ncol(lfqpca_fit$phi0))
  }
  
  theta1 <- if (length(coef_blocks$a_B)) {
    as.numeric(t(A1) %*% coef_blocks$a_B)
  } else {
    numeric(ncol(lfqpca_fit$phi1))
  }
  
  delta <- if (length(coef_blocks$a_U)) {
    as.numeric(t(AU) %*% coef_blocks$a_U)
  } else {
    numeric(ncol(lfqpca_fit$psi))
  }
  
  answer <- list(
    tau_F = lfqpca_fit$tau,
    tau_R = tau_R,
    lfqpca_fit = lfqpca_fit,
    V = V,
    covariates = design$C,
    subject_id = lfqpca_fit$subject_id,
    subject_levels = lfqpca_fit$subject_levels,
    visit_time = lfqpca_fit$visit_time,
    time_std = lfqpca_fit$time_std,
    response_time_degree = response_time_degree,
    G = design$G,
    B_hat = components$B,
    U_hat = components$U,
    chi_B = spaces$chi_B,
    chi_U = spaces$chi_U,
    r_B = spaces$r_B,
    r_U = spaces$r_U,
    integration_weights = spaces$integration_weights,
    x_B = coords$x_B,
    x_U = coords$x_U,
    design = design$X,
    block_sizes = design$block_sizes,
    eligible = eligible,
    response_weights = weight_info$weights,
    I_R = weight_info$I_R,
    J_R = weight_info$J_R,
    coefficients = stage2$coefficients,
    kept_columns = stage2$keep,
    dropped_columns = stage2$dropped,
    design_rank = stage2$rank,
    a_alpha = coef_blocks$a_alpha,
    beta = coef_blocks$beta,
    a_B = coef_blocks$a_B,
    a_U = coef_blocks$a_U,
    gamma_B = gamma_B,
    gamma_U = gamma_U,
    A0 = A0,
    A1 = A1,
    AU = AU,
    theta0 = theta0,
    theta1 = theta1,
    delta = delta,
    fitted = fitted,
    residuals = residuals,
    objective = objective,
    B_orthonormality_error = spaces$B_orthonormality_error,
    U_orthonormality_error = spaces$U_orthonormality_error,
    rq_method = rq_method,
    solver_used = stage2$solver_used,
    rank_tol = rank_tol,
    space_tol = space_tol,
    lfqpca_converged = isTRUE(lfqpca_fit$converged),
    lfqpca_nondegenerate = if (is.null(lfqpca_fit$nondegenerate)) {
      NA
    } else {
      isTRUE(lfqpca_fit$nondegenerate)
    }
  )
  
  class(answer) <- "lfqpcr_fit"
  answer
}