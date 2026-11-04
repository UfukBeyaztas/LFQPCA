lfqpcr_validate_lfqpca_fit <- function(lfqpca_fit) {
  if (!inherits(lfqpca_fit, "lfqpca_fit")) {
    stop("lfqpca_fit must inherit from class 'lfqpca_fit'.", call. = FALSE)
  }
  
  required <- c(
    "tau", "grid", "weights",
    "time_std", "time_center", "time_scale",
    "time_basis", "mean_time_degree",
    "subject_id", "subject_levels", "visit_time",
    "xi", "zeta", "phi0", "phi1", "psi"
  )
  
  missing_fields <- setdiff(required, names(lfqpca_fit))
  
  if (length(missing_fields)) {
    stop(
      "The supplied LFQPCA fit is missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  
  n <- length(lfqpca_fit$subject_id)
  D <- length(lfqpca_fit$grid)
  
  if (length(lfqpca_fit$time_std) != n || length(lfqpca_fit$visit_time) != n) {
    stop("LFQPCA visit-level fields have incompatible lengths.", call. = FALSE)
  }
  
  if (nrow(lfqpca_fit$phi0) != D || nrow(lfqpca_fit$phi1) != D ||
      nrow(lfqpca_fit$psi) != D) {
    stop(
      "LFQPCA loading functions are incompatible with the fitted grid.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}