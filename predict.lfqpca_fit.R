predict.lfqpca_fit <- function(object, newY, new_subject, new_visit_time, ...) {
  if (!inherits(object, "lfqpca_fit")) {
    stop("object must inherit from 'lfqpca_fit'.", call. = FALSE)
  }
  
  newY <- as.matrix(newY)
  storage.mode(newY) <- "double"
  if (ncol(newY) != length(object$grid)) {
    stop("newY must have as many grid points as the training data.",
         call. = FALSE)
  }
  if (length(new_subject) != nrow(newY) ||
      length(new_visit_time) != nrow(newY)) {
    stop("new_subject and new_visit_time must have length nrow(newY).",
         call. = FALSE)
  }
  
  sid_factor <- factor(new_subject, levels = unique(new_subject))
  sid <- as.integer(sid_factor)
  u <- (as.numeric(new_visit_time) - object$time_center) / object$time_scale
  H <- lfqpca_time_basis(u, degree = object$mean_time_degree)
  new_weight_info <- lfqpca_observation_weights(
    Y = newY,
    subject_id = sid,
    quadrature_weights = object$weights
  )
  
  new_objective_weights <- new_weight_info$objective_weights
  
  if (object$K_visit > 0) {
    for (i in seq_len(nlevels(sid_factor))) {
      rows_i <- which(sid == i)
      H_i <- H[rows_i, , drop = FALSE]
      if (qr(H_i)$rank < ncol(H_i) || nrow(H_i) - qr(H_i)$rank < 1) {
        stop(
          "Prediction with visit deviations requires at least ",
          ncol(H_i) + 1,
          " distinct visit times per new subject.",
          call. = FALSE
        )
      }
    }
  }
  
  scores <- lfqpca_update_scores_joint(
    Y = newY,
    subject_id = sid,
    time_std = u,
    H = H,
    eta_functions = object$eta_functions,
    phi0 = object$phi0,
    phi1 = object$phi1,
    psi = object$psi,
    tau = object$tau,
    objective_weights = new_objective_weights,
    rq_method = "br"
  )
  
  fitted <- lfqpca_reconstruct(
    time_basis = H,
    eta_functions = object$eta_functions,
    subject_id = sid,
    time_std = u,
    xi = scores$xi,
    phi0 = object$phi0,
    phi1 = object$phi1,
    zeta = scores$zeta,
    psi = object$psi
  )
  
  list(
    fitted = fitted,
    xi = scores$xi,
    zeta = scores$zeta,
    subject_id = sid,
    subject_levels = levels(sid_factor),
    time_std = u,
    objective = lfqpca_weighted_check_loss(
      residual = newY - fitted,
      tau = object$tau,
      objective_weights = new_objective_weights
    ),
    constraint_diagnostics = if (object$K_visit) {
      lfqpca_constraint_diagnostics(scores$zeta, sid, H)
    } else matrix(0, nlevels(sid_factor), 3)
  )
}