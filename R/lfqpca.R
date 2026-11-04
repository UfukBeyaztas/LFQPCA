lfqpca <- function(
    Y,
    subject,
    visit_time,
    grid = NULL,
    tau = 0.5,
    K_subject = 2,
    K_visit = 1,
    spline_df = 8,
    spline_degree = 3,
    periodic = FALSE,
    mean_time_degree = 1,
    max_iter = 30,
    tol = 1e-4,
    n_starts = 2,
    rq_method_functions = c("fnb", "fn", "br"),
    rq_method_scores = c("br", "fnb", "fn"),
    strict_identification = TRUE,
    objective_increase_tolerance = 1e-7,
    verbose = TRUE,
    warm_start_fit = NULL,
    rank_tol = 1e-8) {
  
  rq_method_functions <- match.arg(rq_method_functions)
  rq_method_scores <- match.arg(rq_method_scores)
  
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0 ||
      tau >= 1) {
    stop("'tau' must be a finite scalar strictly between 0 and 1.",
         call. = FALSE)
  }
  
  K_subject <- as.integer(K_subject)
  K_visit <- as.integer(K_visit)
  spline_df <- as.integer(spline_df)
  spline_degree <- as.integer(spline_degree)
  mean_time_degree <- as.integer(mean_time_degree)
  max_iter <- as.integer(max_iter)
  n_starts <- as.integer(n_starts)
  
  if (K_subject < 1) stop("K_subject must be at least 1.", call. = FALSE)
  if (K_visit < 0) stop("K_visit cannot be negative.", call. = FALSE)
  if (mean_time_degree < 1) {
    stop("mean_time_degree must be at least 1.", call. = FALSE)
  }
  if (max_iter < 1 || n_starts < 1 || tol <= 0) {
    stop("Invalid iteration or convergence controls.", call. = FALSE)
  }
  
  dat <- lfqpca_validate_data(
    Y = Y,
    subject = subject,
    visit_time = visit_time,
    grid = grid,
    K_visit = K_visit,
    strict_identification = strict_identification
  )
  
  Y <- dat$Y
  subject_id <- dat$subject_id
  visit_time <- dat$visit_time
  grid <- dat$grid
  
  if (K_subject >= dat$I) {
    stop("K_subject must be smaller than the number of subjects.",
         call. = FALSE)
  }
  if (K_visit >= dat$n && K_visit > 0) {
    stop("K_visit must be smaller than the number of curves.", call. = FALSE)
  }
  
  time_center <- mean(visit_time)
  time_scale <- sd(visit_time)
  if (!is.finite(time_scale) || time_scale <= 0) {
    stop("visit_time must vary.", call. = FALSE)
  }
  time_std <- (visit_time - time_center) / time_scale
  
  B <- lfqpca_function_basis(
    grid = grid,
    df = spline_df,
    degree = spline_degree,
    periodic = periodic
  )
  H <- lfqpca_time_basis(time_std = time_std, degree = mean_time_degree)
  weights <- lfqpca_trapezoid_weights(grid)
  observation_weight_info <- lfqpca_observation_weights(
    Y = Y,
    subject_id = subject_id,
    quadrature_weights = weights
  )
  
  domain_weights <- observation_weight_info$domain_weights
  objective_weights <- observation_weight_info$objective_weights
  
  if (K_subject > min(dat$I - 1, 2 * ncol(B))) {
    stop(
      "K_subject must not exceed the dimension of the centered subject block, ",
      "min(number of subjects - 1, 2 * spline basis dimension).",
      call. = FALSE
    )
  }
  if (K_visit > ncol(B)) {
    stop(
      "K_visit cannot exceed the functional spline-basis dimension.",
      call. = FALSE
    )
  }
  
  total_within_df <- 0
  for (i in seq_len(dat$I)) {
    rows_i <- which(subject_id == i)
    H_i <- H[rows_i, , drop = FALSE]
    rank_H <- qr(H_i)$rank
    
    if (qr(cbind(1, time_std[rows_i]))$rank < 2) {
      stop(
        "Subject ", dat$subject_levels[i],
        " does not have at least two distinct visit times.",
        call. = FALSE
      )
    }
    if (rank_H < ncol(H)) {
      stop(
        "The population time basis is rank deficient for subject ",
        dat$subject_levels[i], ". Reduce mean_time_degree or add visits.",
        call. = FALSE
      )
    }
    within_df_i <- nrow(H_i) - rank_H
    total_within_df <- total_within_df + max(within_df_i, 0)
    if (K_visit > 0 && strict_identification && within_df_i < 1) {
      stop(
        "Subject ", dat$subject_levels[i],
        " has no residual visit-score degree of freedom after removing the ",
        ncol(H), "-dimensional population time basis. At least ",
        ncol(H) + 1, " suitably distinct visits are required when K_visit > 0.",
        call. = FALSE
      )
    }
  }
  if (K_visit > total_within_df) {
    stop(
      "K_visit exceeds the total dimension available after the within-subject ",
      "time-basis constraints.",
      call. = FALSE
    )
  }
  
  warm_init <- NULL
  if (!is.null(warm_start_fit)) {
    if (!inherits(warm_start_fit, "lfqpca_fit")) {
      stop("warm_start_fit must inherit from 'lfqpca_fit'.", call. = FALSE)
    }
    if (n_starts < 2) {
      stop("n_starts must be at least 2 when warm_start_fit is given.",
           call. = FALSE)
    }
    if (warm_start_fit$K_subject != K_subject ||
        warm_start_fit$K_visit != K_visit) {
      stop("warm_start_fit has incompatible K_subject/K_visit.", call. = FALSE)
    }
    if (length(warm_start_fit$grid) != length(grid) ||
        max(abs(as.numeric(warm_start_fit$grid) - grid)) > 1e-10) {
      stop("warm_start_fit must use the same functional grid.",
           call. = FALSE)
    }
    if (length(warm_start_fit$time_std) != length(time_std) ||
        max(abs(as.numeric(warm_start_fit$time_std) - time_std)) > 1e-8) {
      stop("warm_start_fit must use the same subjects and visit times.",
           call. = FALSE)
    }
    if (nrow(warm_start_fit$xi) != dat$I ||
        ncol(warm_start_fit$xi) != K_subject ||
        nrow(warm_start_fit$zeta) != dat$n ||
        ncol(warm_start_fit$zeta) != K_visit ||
        !all(dim(warm_start_fit$phi0) == c(length(grid), K_subject)) ||
        !all(dim(warm_start_fit$phi1) == c(length(grid), K_subject)) ||
        !all(dim(warm_start_fit$psi) == c(length(grid), K_visit)) ||
        !all(dim(warm_start_fit$eta_functions) == c(length(grid), ncol(H)))) {
      stop("warm_start_fit has incompatible function or score dimensions.",
           call. = FALSE)
    }
    warm_init <- list(
      eta_functions = warm_start_fit$eta_functions,
      xi = warm_start_fit$xi,
      phi0 = warm_start_fit$phi0,
      phi1 = warm_start_fit$phi1,
      zeta = warm_start_fit$zeta,
      psi = warm_start_fit$psi
    )
  }
  
  fits <- vector("list", n_starts)
  start_types <- character(n_starts)
  start_errors <- rep(NA_character_, n_starts)
  for (st in seq_len(n_starts)) {
    use_warm <- !is.null(warm_init) && st == 2
    random_start <- if (use_warm) FALSE else st > 1
    start_types[st] <- if (st == 1) {
      "structured_svd"
    } else if (use_warm) {
      "quantile_continuation"
    } else {
      "random"
    }
    
    if (verbose && n_starts > 1) {
      message(
        "Starting LFQPCA fit ", st, " of ", n_starts,
        " [", start_types[st], "]"
      )
    }
    
    fits[[st]] <- tryCatch(
      lfqpca_one_start(
        Y = Y,
        subject_id = subject_id,
        time_std = time_std,
        grid = grid,
        B = B,
        H = H,
        weights = weights,
        objective_weights = objective_weights,
        Kx = K_subject,
        Ku = K_visit,
        tau = tau,
        max_iter = max_iter,
        tol = tol,
        rq_method_functions = rq_method_functions,
        rq_method_scores = rq_method_scores,
        random_start = random_start,
        init_override = if (use_warm) warm_init else NULL,
        objective_increase_tolerance = objective_increase_tolerance,
        verbose = verbose,
        rank_tol = rank_tol
      ),
      error = function(e) {
        start_errors[st] <<- conditionMessage(e)
        NULL
      }
    )
  }
  
  start_ok <- !vapply(fits, is.null, logical(1))
  if (!any(start_ok)) {
    stop(
      "All ", n_starts, " LFQPCA starts failed. ",
      paste(
        sprintf("[start %d, %s] %s", seq_len(n_starts), start_types,
                start_errors),
        collapse = " "
      ),
      call. = FALSE
    )
  }
  
  start_objectives <- vapply(
    fits,
    function(x) if (is.null(x)) NA_real_ else x$objective,
    numeric(1)
  )
  start_converged <- vapply(
    fits,
    function(x) !is.null(x) && isTRUE(x$converged),
    logical(1)
  )
  start_nondegenerate <- vapply(
    fits,
    function(x) !is.null(x) && isTRUE(x$nondegenerate),
    logical(1)
  )
  
  eligible <- which(start_ok & is.finite(start_objectives) &
                      start_nondegenerate)
  if (!length(eligible)) {
    eligible <- which(start_ok & is.finite(start_objectives))
  }
  if (!length(eligible)) {
    stop("No LFQPCA start produced a finite objective.", call. = FALSE)
  }
  best_start <- eligible[which.min(start_objectives[eligible])]
  answer <- fits[[best_start]]
  
  answer$tau <- tau
  answer$K_subject <- K_subject
  answer$K_visit <- K_visit
  answer$grid <- grid
  answer$weights <- weights
  answer$domain_weights <- domain_weights
  answer$objective_weights <- objective_weights
  answer$basis <- B
  answer$time_basis <- H
  answer$time_std <- time_std
  answer$time_center <- time_center
  answer$time_scale <- time_scale
  answer$mean_time_degree <- mean_time_degree
  answer$subject_id <- subject_id
  answer$subject_levels <- dat$subject_levels
  answer$visit_time <- visit_time
  answer$Y <- Y
  answer$mask <- is.finite(Y)
  answer$spline_df <- spline_df
  answer$spline_degree <- spline_degree
  answer$periodic <- periodic
  answer$strict_identification <- strict_identification
  answer$all_start_objectives <- start_objectives
  answer$start_converged <- start_converged
  answer$start_nondegenerate <- start_nondegenerate
  answer$start_errors <- start_errors
  answer$n_failed_starts <- sum(!start_ok)
  answer$rank_tol <- rank_tol
  answer$start_types <- start_types
  answer$best_start <- best_start
  answer$best_start_type <- start_types[best_start]
  answer$warm_start_used <- !is.null(warm_init)
  answer$warm_start_tau <-
    if (!is.null(warm_start_fit)) warm_start_fit$tau else NA_real_
  answer$pve_subject <- if (sum(answer$lambda) > 0) {
    cumsum(answer$lambda) / sum(answer$lambda)
  } else numeric(0)
  answer$pve_visit <- if (sum(answer$nu) > 0) {
    cumsum(answer$nu) / sum(answer$nu)
  } else numeric(0)
  
  class(answer) <- "lfqpca_fit"
  answer
}
