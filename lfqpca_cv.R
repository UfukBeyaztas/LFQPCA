lfqpca_cv <- function(
    Y,
    subject,
    visit_time,
    grid = NULL,
    tuning_grid = expand.grid(
      K_subject = c(1, 2),
      K_visit = c(1, 2),
      spline_df = c(6, 8),
      KEEP.OUT.ATTRS = FALSE
    ),
    tau = 0.5,
    n_folds = 3,
    max_iter = 20,
    n_starts = 1,
    verbose = TRUE,
    ...) {
  
  dat <- lfqpca_validate_data(
    Y = Y,
    subject = subject,
    visit_time = visit_time,
    grid = grid,
    K_visit = max(tuning_grid$K_visit),
    strict_identification = TRUE
  )
  
  required <- c("K_subject", "K_visit", "spline_df")
  if (!all(required %in% names(tuning_grid))) {
    stop("tuning_grid must contain K_subject, K_visit, and spline_df.",
         call. = FALSE)
  }
  
  fold_id <- lfqpca_make_point_folds(dat$Y, n_folds = n_folds)
  rows <- vector("list", nrow(tuning_grid) * n_folds)
  counter <- 1
  
  for (g in seq_len(nrow(tuning_grid))) {
    for (fold in seq_len(n_folds)) {
      if (verbose) {
        message(
          "CV model ", g, "/", nrow(tuning_grid),
          ", fold ", fold, "/", n_folds
        )
      }
      
      held_out <- fold_id == fold
      Y_train <- dat$Y
      Y_train[held_out] <- NA_real_
      
      fit <- lfqpca(
        Y = Y_train,
        subject = dat$subject_levels[dat$subject_id],
        visit_time = dat$visit_time,
        grid = dat$grid,
        tau = tau,
        K_subject = tuning_grid$K_subject[g],
        K_visit = tuning_grid$K_visit[g],
        spline_df = tuning_grid$spline_df[g],
        max_iter = max_iter,
        n_starts = n_starts,
        verbose = FALSE,
        ...
      )
      
      test_objective_weights <- lfqpca_cv_test_weights(
        Y = dat$Y,
        held_out = held_out,
        subject_id = dat$subject_id,
        quadrature_weights = lfqpca_trapezoid_weights(
          dat$grid,
          normalize = TRUE
        )
      )
      
      residual_test <- dat$Y - fit$fitted
      
      loss <- lfqpca_weighted_check_loss(
        residual = residual_test,
        tau = tau,
        objective_weights = test_objective_weights
      )
      rows[[counter]] <- data.frame(
        K_subject = tuning_grid$K_subject[g],
        K_visit = tuning_grid$K_visit[g],
        spline_df = tuning_grid$spline_df[g],
        fold = fold,
        cv_objective = loss
      )
      counter <- counter + 1
    }
  }
  
  fold_results <- do.call(rbind, rows)
  mean_results <- aggregate(
    cv_objective ~ K_subject + K_visit + spline_df,
    data = fold_results,
    FUN = mean
  )
  
  sd_results <- aggregate(
    cv_objective ~ K_subject + K_visit + spline_df,
    data = fold_results,
    FUN = sd
  )
  names(mean_results)[names(mean_results) == "cv_objective"] <-
    "mean_cv_objective"
  
  names(sd_results)[names(sd_results) == "cv_objective"] <- "sd_cv_objective"
  summary_results <- merge(
    mean_results,
    sd_results,
    by = c("K_subject", "K_visit", "spline_df")
  )
  summary_results <-
    summary_results[order(summary_results$mean_cv_objective), , drop = FALSE]
  rownames(summary_results) <- NULL
  
  list(
    fold_results = fold_results,
    summary = summary_results,
    best = summary_results[1, , drop = FALSE],
    fold_id = fold_id
  )
}