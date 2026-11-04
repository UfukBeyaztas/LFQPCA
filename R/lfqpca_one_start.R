lfqpca_one_start <- function(
    Y,
    subject_id,
    time_std,
    grid,
    B,
    H,
    weights,
    objective_weights,
    Kx,
    Ku,
    tau,
    max_iter,
    tol,
    rq_method_functions,
    rq_method_scores,
    random_start,
    init_override = NULL,
    objective_increase_tolerance,
    verbose,
    rank_tol = 1e-8) {
  
  if (is.null(init_override)) {
    init <- lfqpca_initialize(
      Y = Y,
      subject_id = subject_id,
      time_std = time_std,
      grid = grid,
      B = B,
      H = H,
      Kx = Kx,
      Ku = Ku,
      tau = tau,
      weights = weights,
      objective_weights = objective_weights,
      rq_method = rq_method_functions,
      random_start = random_start
    )
  } else {
    init <- init_override
  }
  
  eta_functions <- init$eta_functions
  xi <- init$xi
  phi0 <- init$phi0
  phi1 <- init$phi1
  zeta <- init$zeta
  psi <- init$psi
  
  initial_rotation <- lfqpca_absorb_and_rotate(
    eta_functions, xi, phi0, phi1, zeta, psi, weights
  )
  eta_functions <- initial_rotation$eta_functions
  xi <- initial_rotation$xi
  phi0 <- initial_rotation$phi0
  phi1 <- initial_rotation$phi1
  zeta <- initial_rotation$zeta
  psi <- initial_rotation$psi
  
  state <- lfqpca_state_objective(
    Y = Y,
    H = H,
    eta_functions = eta_functions,
    subject_id = subject_id,
    time_std = time_std,
    xi = xi,
    phi0 = phi0,
    phi1 = phi1,
    zeta = zeta,
    psi = psi,
    tau = tau,
    objective_weights = objective_weights
  )
  objective <- state$objective
  fitted <- state$fitted
  if (!is.finite(objective)) {
    stop("The LFQPCA initialization produced a non-finite objective.",
         call. = FALSE)
  }
  objective_history <- objective
  function_objective_history <- numeric(0)
  score_objective_history <- numeric(0)
  converged <- FALSE
  termination <- "max_iter"
  
  accept_candidate <- function(candidate_objective, current_objective) {
    is.finite(candidate_objective) &&
      candidate_objective <= current_objective + objective_increase_tolerance
  }
  
  rejection_log <- data.frame(
    iteration = integer(0),
    block = character(0),
    objective_increase = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (iter in seq_len(max_iter)) {
    old <- list(
      eta_functions = eta_functions,
      xi = xi,
      phi0 = phi0,
      phi1 = phi1,
      zeta = zeta,
      psi = psi,
      fitted = fitted,
      objective = objective
    )
    rejected_blocks <- character(0)
    
    # Block 1: smooth functions.
    function_state <- function(update) {
      lfqpca_state_objective(
        Y = Y,
        H = H,
        eta_functions = update$eta_functions,
        subject_id = subject_id,
        time_std = time_std,
        xi = xi,
        phi0 = update$phi0,
        phi1 = update$phi1,
        zeta = zeta,
        psi = update$psi,
        tau = tau,
        objective_weights = objective_weights
      )
    }
    
    function_update <- lfqpca_update_functions(
      Y = Y,
      B = B,
      H = H,
      subject_id = subject_id,
      time_std = time_std,
      xi = xi,
      zeta = zeta,
      tau = tau,
      objective_weights = objective_weights,
      rq_method = rq_method_functions
    )
    candidate_fun <- function_state(function_update)
    
    if (accept_candidate(candidate_fun$objective, objective)) {
      eta_functions <- function_update$eta_functions
      phi0 <- function_update$phi0
      phi1 <- function_update$phi1
      psi <- function_update$psi
      fitted <- candidate_fun$fitted
      objective <- candidate_fun$objective
    } else {
      rejected_blocks <- c(rejected_blocks, "functions")
      rejection_log <- rbind(
        rejection_log,
        data.frame(
          iteration = iter,
          block = "functions",
          objective_increase = candidate_fun$objective - objective,
          stringsAsFactors = FALSE
        )
      )
      warning("Function update rejected: the objective increased.",
              call. = FALSE)
    }
    function_objective_history <- c(function_objective_history, objective)
    
    # Block 2: subject and visit scores under H_i' Z_i = 0.
    score_state <- function(update) {
      lfqpca_state_objective(
        Y = Y,
        H = H,
        eta_functions = eta_functions,
        subject_id = subject_id,
        time_std = time_std,
        xi = update$xi,
        phi0 = phi0,
        phi1 = phi1,
        zeta = update$zeta,
        psi = psi,
        tau = tau,
        objective_weights = objective_weights
      )
    }
    
    score_update <- lfqpca_update_scores_joint(
      Y = Y,
      subject_id = subject_id,
      time_std = time_std,
      H = H,
      eta_functions = eta_functions,
      phi0 = phi0,
      phi1 = phi1,
      psi = psi,
      tau = tau,
      objective_weights = objective_weights,
      rq_method = rq_method_scores
    )
    candidate_score <- score_state(score_update)
    
    if (accept_candidate(candidate_score$objective, objective)) {
      xi <- score_update$xi
      zeta <- score_update$zeta
      fitted <- candidate_score$fitted
      objective <- candidate_score$objective
    } else {
      rejected_blocks <- c(rejected_blocks, "scores")
      rejection_log <- rbind(
        rejection_log,
        data.frame(
          iteration = iter,
          block = "scores",
          objective_increase = candidate_score$objective - objective,
          stringsAsFactors = FALSE
        )
      )
      warning("Score update rejected: the objective increased.",
              call. = FALSE)
    }
    score_objective_history <- c(score_objective_history, objective)
    
    # Block 3: identifying rotations.
    rotated <- lfqpca_absorb_and_rotate(
      eta_functions = eta_functions,
      xi = xi,
      phi0 = phi0,
      phi1 = phi1,
      zeta = zeta,
      psi = psi,
      weights = weights
    )
    eta_functions <- rotated$eta_functions
    xi <- rotated$xi
    phi0 <- rotated$phi0
    phi1 <- rotated$phi1
    zeta <- rotated$zeta
    psi <- rotated$psi
    
    rotated_state <- lfqpca_state_objective(
      Y = Y,
      H = H,
      eta_functions = eta_functions,
      subject_id = subject_id,
      time_std = time_std,
      xi = xi,
      phi0 = phi0,
      phi1 = phi1,
      zeta = zeta,
      psi = psi,
      tau = tau,
      objective_weights = objective_weights
    )
    fitted <- rotated_state$fitted
    objective_new <- rotated_state$objective
    
    rotation_error <- abs(objective_new - objective)
    if (!is.finite(rotation_error) ||
        rotation_error > max(1e-8, 10 * objective_increase_tolerance)) {
      warning(
        "The identifying rotation changed the objective more than expected: ",
        format(rotation_error, scientific = TRUE),
        call. = FALSE
      )
    }
    
    objective_history <- c(objective_history, objective_new)
    relative_change <- abs(objective_new - old$objective) /
      (abs(old$objective) + 1e-10)
    
    if (verbose) {
      message(sprintf(
        "LFQPCA iteration %2d | objective = %.8f | relative change = %.3e%s",
        iter, objective_new, relative_change,
        if (length(rejected_blocks)) {
          paste0(" | rejected: ", paste(rejected_blocks, collapse = ", "))
        } else ""
      ))
    }
    
    objective <- objective_new
    
    if (length(rejected_blocks) &&
        (length(rejected_blocks) == 2 || relative_change < tol)) {
      # A rejected block left the state unchanged: stop without convergence.
      termination <- "stalled_rejected_updates"
      break
    }
    
    if (relative_change < tol) {
      converged <- TRUE
      termination <- "converged"
      break
    }
  }
  
  theta_eta <- qr.solve(B, eta_functions)
  theta_phi0 <- qr.solve(B, phi0)
  theta_phi1 <- qr.solve(B, phi1)
  theta_psi <- if (Ku) qr.solve(B, psi) else matrix(0, ncol(B), 0)
  
  lambda <- if (Kx) apply(xi, 2, var) else numeric(0)
  nu <- if (Ku) apply(zeta, 2, var) else numeric(0)
  
  constraint_diagnostics <- if (Ku) {
    lfqpca_constraint_diagnostics(zeta, subject_id, H)
  } else {
    matrix(
      0, max(subject_id), 3,
      dimnames = list(
        NULL,
        c("max_abs_intercept_sum", "max_abs_linear_time_sum",
          "max_abs_all_time_basis")
      )
    )
  }
  
  rank_diagnostics <- lfqpca_rank_diagnostics(
    xi = xi,
    phi0 = phi0,
    phi1 = phi1,
    zeta = zeta,
    psi = psi,
    weights = weights,
    rank_tol = rank_tol
  )
  
  list(
    eta_functions = eta_functions,
    phi0 = phi0,
    phi1 = phi1,
    psi = psi,
    xi = xi,
    zeta = zeta,
    lambda = lambda,
    nu = nu,
    theta_eta = theta_eta,
    theta_phi0 = theta_phi0,
    theta_phi1 = theta_phi1,
    theta_psi = theta_psi,
    fitted = fitted,
    objective = objective,
    objective_history = objective_history,
    function_objective_history = function_objective_history,
    score_objective_history = score_objective_history,
    converged = converged,
    termination = termination,
    iterations = length(objective_history) - 1,
    n_rejected_function_updates = sum(rejection_log$block == "functions"),
    n_rejected_score_updates = sum(rejection_log$block == "scores"),
    rejection_log = rejection_log,
    rank_diagnostics = rank_diagnostics,
    nondegenerate = isTRUE(rank_diagnostics$nondegenerate),
    constraint_diagnostics = constraint_diagnostics
  )
}
