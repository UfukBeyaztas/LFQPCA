lfqpca_continuation <- function(Y, subject, visit_time, grid = NULL, tau,
                                reference_tau = 0.5,
                                fractions = c(0, 0.5, 0.75, 1),
                                n_starts = 5, path_starts = 2,
                                reference_fit = NULL, verbose = FALSE, ...) {
  path <- lfqpca_continuation_path(tau, reference_tau, fractions)
  fit_at <- function(tau_now, starts, warm) {
    lfqpca(Y = Y, subject = subject, visit_time = visit_time, grid = grid,
           tau = tau_now, n_starts = starts, warm_start_fit = warm,
           verbose = verbose, ...)
  }
  errors <- character(0)
  steps <- seq_along(path)
  warm <- NULL
  fit <- NULL
  if (!is.null(reference_fit)) {
    if (!inherits(reference_fit, "lfqpca_fit")) {
      stop("'reference_fit' must be an lfqpca_fit object.", call. = FALSE)
    }
    fit <- reference_fit
    warm <- reference_fit
    steps <- steps[-1]
  }
  for (k in steps) {
    last <- k == length(path)
    starts <- if (last || k == 1) n_starts else path_starts
    fit_k <- tryCatch(fit_at(path[k], starts, warm), error = function(e) e)
    if (inherits(fit_k, "error")) {
      if (last) stop(conditionMessage(fit_k), call. = FALSE)
      label <- if (k == 1) " (reference): " else ": "
      errors <- c(errors, paste0("tau=", format(path[k]), label,
                                 conditionMessage(fit_k)))
      next
    }
    fit <- fit_k
    warm <- fit_k
  }
  fit$continuation_errors <- errors
  fit$continuation_tau_path <- path
  fit
}
