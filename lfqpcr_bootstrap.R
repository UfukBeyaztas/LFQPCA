lfqpcr_bootstrap <- function(
    Y,
    V,
    subject,
    visit_time,
    grid = NULL,
    covariates = NULL,
    B = 200,
    seed = 1,
    ...,
    verbose = TRUE) {
  
  Y <- as.matrix(Y)
  V <- as.numeric(V)
  C <- lfqpcr_prepare_covariates(covariates, nrow(Y))
  
  subject_levels <- unique(subject)
  I <- length(subject_levels)
  B <- as.integer(B)
  
  if (B < 1) {
    stop("B must be at least 1.", call. = FALSE)
  }
  
  restore_rng <- lfqpcr_rng_state_saver()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)
  fits <- vector("list", B)
  
  for (b in seq_len(B)) {
    if (verbose) {
      message("LFQPCR bootstrap replicate ", b, " of ", B)
    }
    
    sampled_subjects <- sample(subject_levels, size = I, replace = TRUE)
    
    rows_b <- integer(0)
    subject_b <- character(0)
    
    for (m in seq_along(sampled_subjects)) {
      rows_m <- which(subject == sampled_subjects[m])
      rows_b <- c(rows_b, rows_m)
      subject_b <- c(subject_b, rep(paste0("boot_", m), length(rows_m)))
    }
    
    fits[[b]] <- tryCatch(
      lfqpcr(
        Y = Y[rows_b, , drop = FALSE],
        V = V[rows_b],
        subject = subject_b,
        visit_time = visit_time[rows_b],
        grid = grid,
        covariates = C[rows_b, , drop = FALSE],
        ...,
        verbose = FALSE
      ),
      error = function(e) {
        structure(
          list(message = conditionMessage(e)),
          class = "lfqpcr_bootstrap_error"
        )
      }
    )
  }
  
  ok <- vapply(fits, function(x) inherits(x, "lfqpcr_fit"), logical(1))
  
  list(
    fits = fits,
    successful = ok,
    n_successful = sum(ok),
    n_failed = sum(!ok),
    seed = seed
  )
}