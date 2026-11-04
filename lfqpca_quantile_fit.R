lfqpca_quantile_fit <- function(X, y, tau,
                                case_weights = NULL,
                                method = c("fnb", "br", "fn"),
                                qr_tol = 1e-7,
                                allow_rank_reduction = TRUE) {
  method <- match.arg(method)
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  if (nrow(X) != length(y)) {
    stop("nrow(X) must equal length(y) in quantile regression.", call. = FALSE)
  }
  
  if (is.null(case_weights)) {
    case_weights <- rep(1, length(y))
  } else {
    case_weights <- as.numeric(case_weights)
    
    if (length(case_weights) != length(y)) {
      stop(
        "length(case_weights) must equal length(y) in quantile regression.",
        call. = FALSE
      )
    }
  }
  
  finite_rows <-
    is.finite(y) & is.finite(case_weights) & case_weights > 0 &
    apply(X, 1, function(z) all(is.finite(z)))
  
  X <- X[finite_rows, , drop = FALSE]
  y <- y[finite_rows]
  case_weights <- case_weights[finite_rows]
  
  if (!nrow(X)) {
    stop(
      "No finite observations with positive weight in quantile regression.",
      call. = FALSE
    )
  }
  
  if (!ncol(X)) {
    return(numeric(0))
  }
  
  if (nrow(X) <= ncol(X)) {
    stop(
      "Quantile-regression block has ", nrow(X),
      " observations but ", ncol(X),
      " coefficients. Reduce component counts or increase data density.",
      call. = FALSE
    )
  }
  
  case_weights <- case_weights / mean(case_weights)
  
  Xw_all <- sweep(X, MARGIN = 1, STATS = case_weights, FUN = "*")
  yw <- y * case_weights
  
  keep <- seq_len(ncol(X))
  
  if (allow_rank_reduction) {
    keep <- lfqpca_safe_independent_columns(Xw_all, tol = qr_tol)
  }
  
  if (!length(keep)) {
    return(rep(0, ncol(X)))
  }
  
  Xw <- Xw_all[, keep, drop = FALSE]
  beta_full <- numeric(ncol(X))
  
  attempt <- lfqpca_rq_solve(Xw, yw, tau, method)
  
  if (!attempt$ok && method != "br") {
    primary_message <- attempt$message
    attempt <- lfqpca_rq_solve(Xw, yw, tau, "br")
    
    if (!attempt$ok) {
      stop(
        "Quantile regression failed. Primary solver (", method, "): ",
        primary_message,
        "; Barrodale-Roberts fallback: ",
        attempt$message,
        call. = FALSE
      )
    }
  } else if (!attempt$ok) {
    stop(
      "Quantile regression failed (Barrodale-Roberts): ",
      attempt$message,
      call. = FALSE
    )
  }
  
  beta_full[keep] <- attempt$coefficients
  
  beta_full
}