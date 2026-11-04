lfqpcr_cv <- function(
    Y,
    V,
    subject,
    visit_time,
    grid = NULL,
    covariates = NULL,
    tau_F = 0.5,
    tau_R = 0.5,
    K_subject = 2,
    K_visit = 1,
    spline_df = 8,
    spline_degree = 3,
    periodic = FALSE,
    mean_time_degree = 1,
    response_time_degree = NULL,
    n_folds = 5,
    seed = 1,
    lfqpca_max_iter = 30,
    lfqpca_tol = 1e-4,
    lfqpca_n_starts = 2,
    lfqpca_rq_method_functions = c("fnb", "fn", "br"),
    lfqpca_rq_method_scores = c("br", "fnb", "fn"),
    strict_identification = TRUE,
    objective_increase_tolerance = 1e-7,
    lfqpcr_rq_method = c("fnb", "br", "fn"),
    space_tol = 1e-10,
    rank_tol = 1e-7,
    verbose = TRUE) {
  
  Y <- as.matrix(Y)
  V <- as.numeric(V)
  
  if (nrow(Y) != length(V) || length(subject) != nrow(Y) ||
      length(visit_time) != nrow(Y)) {
    stop(
      "Y, V, subject, and visit_time have incompatible lengths.",
      call. = FALSE
    )
  }
  
  C <- lfqpcr_prepare_covariates(covariates, nrow(Y))
  
  subject_levels <- unique(subject)
  I <- length(subject_levels)
  
  n_folds <- as.integer(n_folds)
  
  if (n_folds < 2 || n_folds > I) {
    stop("n_folds must be between 2 and the number of subjects.", call. = FALSE)
  }
  
  restore_rng <- lfqpcr_rng_state_saver()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)
  shuffled <- sample(subject_levels)
  fold_id_subject <- rep(seq_len(n_folds), length.out = I)
  fold_map <- setNames(fold_id_subject, as.character(shuffled))
  fold_id <- unname(fold_map[as.character(subject)])
  
  prediction <- rep(NA_real_, nrow(Y))
  
  for (k in seq_len(n_folds)) {
    if (verbose) {
      message("LFQPCR CV fold ", k, " of ", n_folds)
    }
    
    test <- fold_id == k
    train <- !test
    
    fit_k <- lfqpcr(
      Y = Y[train, , drop = FALSE],
      V = V[train],
      subject = subject[train],
      visit_time = visit_time[train],
      grid = grid,
      covariates = C[train, , drop = FALSE],
      tau_F = tau_F,
      tau_R = tau_R,
      K_subject = K_subject,
      K_visit = K_visit,
      spline_df = spline_df,
      spline_degree = spline_degree,
      periodic = periodic,
      mean_time_degree = mean_time_degree,
      response_time_degree = response_time_degree,
      lfqpca_max_iter = lfqpca_max_iter,
      lfqpca_tol = lfqpca_tol,
      lfqpca_n_starts = lfqpca_n_starts,
      lfqpca_rq_method_functions = lfqpca_rq_method_functions,
      lfqpca_rq_method_scores = lfqpca_rq_method_scores,
      strict_identification = strict_identification,
      objective_increase_tolerance = objective_increase_tolerance,
      lfqpcr_rq_method = lfqpcr_rq_method,
      space_tol = space_tol,
      rank_tol = rank_tol,
      verbose = FALSE
    )
    
    pred_k <- predict(
      fit_k,
      newY = Y[test, , drop = FALSE],
      new_subject = subject[test],
      new_visit_time = visit_time[test],
      new_covariates = C[test, , drop = FALSE]
    )
    
    prediction[test] <- pred_k$fitted
  }
  
  covariates_finite <- if (ncol(C)) {
    apply(C, 1, function(z) all(is.finite(z)))
  } else {
    rep(TRUE, nrow(Y))
  }
  
  eligible <- is.finite(V) & is.finite(prediction) & covariates_finite
  
  sid <- as.integer(factor(subject, levels = unique(subject)))
  
  weight_info <- lfqpcr_subject_weights(subject_id = sid, eligible = eligible)
  
  cv_objective <- lfqpcr_weighted_objective(
    residual = V - prediction,
    tau_R = tau_R,
    weights = weight_info$weights
  )
  
  list(
    tau_F = tau_F,
    tau_R = tau_R,
    n_folds = n_folds,
    fold_id = fold_id,
    prediction = prediction,
    eligible = eligible,
    response_weights = weight_info$weights,
    cv_objective = cv_objective
  )
}