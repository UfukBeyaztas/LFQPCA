lfqpcr_quantile_fit <- function(
    X,
    y,
    tau_R,
    case_weights,
    rq_method = c("fnb", "br", "fn"),
    rank_tol = 1e-7) {
  
  rq_method <- match.arg(rq_method)
  X <- as.matrix(X)
  y <- as.numeric(y)
  case_weights <- as.numeric(case_weights)
  
  if (nrow(X) != length(y) || length(y) != length(case_weights)) {
    stop("X, y, and case_weights have incompatible dimensions.", call. = FALSE)
  }
  
  if (!is.numeric(tau_R) || length(tau_R) != 1 || !is.finite(tau_R) ||
      tau_R <= 0 ||
      tau_R >= 1) {
    stop("'tau_R' must be strictly between 0 and 1.", call. = FALSE)
  }
  
  valid <-
    is.finite(y) & is.finite(case_weights) & case_weights > 0 &
    apply(X, 1, function(z) all(is.finite(z)))
  
  Xv <- X[valid, , drop = FALSE]
  yv <- y[valid]
  wv <- case_weights[valid]
  
  if (!nrow(Xv)) {
    stop("No eligible observations remain for LFQPCR stage 2.", call. = FALSE)
  }
  
  wv <- wv / mean(wv)
  Xw_all <- sweep(Xv, MARGIN = 1, STATS = wv, FUN = "*")
  yw <- yv * wv
  
  keep <- lfqpca_safe_independent_columns(Xw_all, tol = rank_tol)
  
  if (!length(keep)) {
    stop("The LFQPCR stage-2 design has rank zero.", call. = FALSE)
  }
  
  Xw <- Xw_all[, keep, drop = FALSE]
  
  if (nrow(Xw) <= ncol(Xw)) {
    stop(
      "The LFQPCR stage-2 regression has ",
      nrow(Xw), " eligible observations but ",
      ncol(Xw), " independent coefficients. ",
      "Reduce the model dimension or add scalar-response visits.",
      call. = FALSE
    )
  }
  
  attempt <- lfqpca_rq_solve(Xw, yw, tau_R, rq_method)
  solver_used <- rq_method
  
  if (!attempt$ok && rq_method != "br") {
    primary_message <- attempt$message
    attempt <- lfqpca_rq_solve(Xw, yw, tau_R, "br")
    solver_used <- "br"
    
    if (!attempt$ok) {
      stop(
        "LFQPCR quantile regression failed. Primary solver (", rq_method, "): ",
        primary_message,
        "; Barrodale-Roberts fallback: ",
        attempt$message,
        call. = FALSE
      )
    }
  } else if (!attempt$ok) {
    stop(
      "LFQPCR quantile regression failed (Barrodale-Roberts): ",
      attempt$message,
      call. = FALSE
    )
  }
  
  beta_reduced <- attempt$coefficients
  
  beta_full <- numeric(ncol(X))
  beta_full[keep] <- beta_reduced
  names(beta_full) <- colnames(X)
  
  fitted_all <- as.numeric(X %*% beta_full)
  
  list(
    coefficients = beta_full,
    coefficients_reduced = beta_reduced,
    keep = keep,
    dropped = setdiff(seq_len(ncol(X)), keep),
    rank = length(keep),
    valid = valid,
    fitted = fitted_all,
    solver_used = solver_used
  )
}